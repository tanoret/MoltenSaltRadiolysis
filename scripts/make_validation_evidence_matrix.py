#!/usr/bin/env python3
"""Generate the publication-style validation hierarchy used in the milestone report.

The figure is intentionally monochrome and table-like. It separates the
question addressed at each evidence level from the strongest claim that can be
made, while keeping operational screening visibly outside the validation
sequence.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO / "manuscript" / "figures"

ROWS = [
    (
        "1",
        "Implementation\nverification",
        "Do units, conservation laws, positivity,\nand numerical tolerances behave as specified?",
        "The implementation solves the stated equations\nfor the tested cases.",
    ),
    (
        "2",
        "Calibration and\nin-sample checks",
        "Can the inferred parameters reproduce the data\nused to calibrate them?",
        "The model is consistent with the calibration data\nunder the stated likelihood and priors.",
    ),
    (
        "3",
        "Conditional\nprediction",
        "Does the kinetic response transfer to withheld\nconditions under an explicit nuisance model?",
        "Predictive transfer is supported only for\nthe withheld conditions and nuisance\ntreatment tested.",
    ),
    (
        "4",
        "External mechanism\nanchors",
        "Are selected rates, channels, or closures\ncompatible with independent\nmeasurements?",
        "Selected components are physically plausible;\nmechanism completeness is not established.",
    ),
    (
        "5",
        "Observation-limited\nbenchmark",
        "Which parameter combination is identifiable\nwhen transport and analytical recovery\nare uncertain?",
        "Only the product $\\eta f_{\\mathrm{rec}}$ is bounded;\nuranium kinetics remain unidentified.",
    ),
    (
        "6",
        "Operational\nscreening",
        "What consequences follow when the upstream\nassumptions and uncertainties are propagated?",
        "Sensitivity, scenario comparison, and\nexperiment planning - not a design-basis\nprediction.",
    ),
]


def add_text(
    ax,
    x: float,
    y: float,
    text: str,
    *,
    size: float = 7.1,
    weight: str = "normal",
    ha: str = "left",
    va: str = "center",
) -> None:
    ax.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        ha=ha,
        va=va,
        color="black",
        linespacing=1.12,
    )


def make_figure(output_dir: Path, png_dpi: int = 300) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = output_dir / "fig_validation_evidence_ladder.pdf"
    out_png = output_dir / "fig_validation_evidence_ladder.png"

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(7.35, 4.15))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    left, right = 0.025, 0.985
    top = 0.955
    header_h = 0.075
    row_h = 0.125
    note_h = 0.09
    x = [left, 0.085, 0.315, 0.640, right]

    ax.add_patch(
        Rectangle(
            (left, top - header_h),
            right - left,
            header_h,
            facecolor="0.90",
            edgecolor="black",
            linewidth=0.9,
        )
    )
    headers = [
        "Level",
        "Evidence category",
        "Question addressed",
        "Strongest admissible claim",
    ]
    for i, label in enumerate(headers):
        add_text(
            ax,
            (x[i] + x[i + 1]) / 2,
            top - header_h / 2,
            label,
            size=7.3,
            weight="bold",
            ha="center",
        )

    y_top = top - header_h
    for idx, (level, category, question, claim) in enumerate(ROWS):
        y_bottom = y_top - row_h
        face = "0.975" if idx % 2 == 0 else "white"
        if idx == 5:
            face = "0.93"
            ax.plot([left, right], [y_top, y_top], color="black", linewidth=1.6)
        ax.add_patch(
            Rectangle(
                (left, y_bottom),
                right - left,
                row_h,
                facecolor=face,
                edgecolor="0.35",
                linewidth=0.55,
            )
        )
        yc = (y_top + y_bottom) / 2
        add_text(ax, (x[0] + x[1]) / 2, yc, level, size=8.0, weight="bold", ha="center")
        add_text(ax, x[1] + 0.010, yc, category, size=7.1, weight="bold")
        add_text(ax, x[2] + 0.010, yc, question, size=6.55)
        add_text(ax, x[3] + 0.010, yc, claim, size=6.55)
        y_top = y_bottom

    table_bottom = top - header_h - len(ROWS) * row_h
    for xv in x[1:-1]:
        ax.plot([xv, xv], [table_bottom, top], color="0.35", linewidth=0.55)

    note_y = table_bottom - note_h / 2
    ax.plot(
        [left, right],
        [table_bottom - 0.018, table_bottom - 0.018],
        color="black",
        linewidth=0.8,
    )
    add_text(
        ax,
        left,
        note_y - 0.008,
        "Boundary: operational screening propagates upstream assumptions and uncertainty; "
        "it is not additional validation evidence.",
        size=6.85,
        weight="bold",
    )

    fig.subplots_adjust(left=0.015, right=0.995, top=0.995, bottom=0.02)
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(out_png, dpi=png_dpi, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return out_pdf, out_png


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--png-dpi",
        type=int,
        default=300,
        help="PNG resolution in dots per inch (default: 300)",
    )
    args = parser.parse_args()
    if args.png_dpi < 72:
        parser.error("--png-dpi must be at least 72")
    out_pdf, out_png = make_figure(args.output_dir, args.png_dpi)
    print(out_pdf)
    print(out_png)


if __name__ == "__main__":
    main()
