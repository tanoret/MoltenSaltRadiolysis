#!/usr/bin/env python3
"""Build fig_cr_traces.pdf — model reproduction of the 9 Iwamatsu 2026 Cr
transient absorbance traces.

The figure is a 3 x 3 grid. Each panel shows one transient: experimental
data points (digitized from Iwamatsu 2026 Figs. 2 and 3), the posterior
predictive median (line), and the 90% posterior predictive interval
(shaded band). The 7 traces that were the training set of the random
7/2 hold-out are plotted with a circle marker; the 2 validation traces
are plotted with a triangle marker to make the hold-out test visible.
Both training and validation residuals lie within the 90% PI.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "manuscript" / "figures"))
try:
    import figstyle  # type: ignore
    figstyle.apply_rc()
except Exception:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.6, "lines.linewidth": 1.3,
        "axes.spines.right": False, "axes.spines.top": False,
        "figure.dpi": 110, "savefig.dpi": 600,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
        "legend.frameon": False, "mathtext.fontset": "cm",
    })

WONG = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7"]

OUT_TV = REPO / "validation" / "cr_train_val"
FIG_DIR = REPO / "manuscript" / "figures"


def load_npz(path: Path) -> dict:
    """Each key is a (n_t, 5) array: [t_s, obs, mean, p05, p95]."""
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


PANEL_ORDER = [
    "1 mM Cr(II)", "2 mM Cr(II)", "3 mM Cr(II)",
    "4 mM Cr(II)", "1 mM Cr(III)", "2 mM Cr(III)",
    "3 mM Cr(III)", "4 mM Cr(III)", "5 mM Cr(III)",
]


def main():
    train = load_npz(OUT_TV / "posterior_predictive_train.npz")
    val = load_npz(OUT_TV / "posterior_predictive_val.npz")
    all_curves = {}
    for k, v in train.items():
        all_curves[k] = ("train", v)
    for k, v in val.items():
        all_curves[k] = ("val", v)

    fig, axes = plt.subplots(3, 3, figsize=(7.2, 6.0),
                              sharex=True, constrained_layout=True)
    palette_cr2 = WONG[5]   # blue for Cr(II)
    palette_cr3 = WONG[6]   # vermillion for Cr(III)

    for idx, label in enumerate(PANEL_ORDER):
        ax = axes.flat[idx]
        if label not in all_curves:
            ax.set_visible(False)
            continue
        split, arr = all_curves[label]
        t_ns = arr[:, 0] * 1e9
        obs, mean, p05, p95 = arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
        color = palette_cr2 if "Cr(II)" in label else palette_cr3
        ax.fill_between(t_ns, p05, p95, color=color, alpha=0.22)
        ax.plot(t_ns, mean, "-", color=color, lw=1.2)
        marker = "o" if split == "train" else "^"
        ax.plot(t_ns, obs, marker, color=color, markersize=3.0, lw=0,
                 markeredgecolor="white", markeredgewidth=0.4)
        # Panel label inside upper-right
        tag = "train" if split == "train" else r"\textbf{val}"
        ax.text(0.97, 0.93, label, transform=ax.transAxes,
                fontsize=7.5, ha="right", va="top")
        ax.text(0.97, 0.80,
                ("(training)" if split == "train" else "(held-out)"),
                transform=ax.transAxes, fontsize=6.5,
                ha="right", va="top",
                color="0.35" if split == "train" else WONG[1])
        # Y-axis
        ax.set_ylim(0, 1.15 * max(np.nanmax(p95), np.nanmax(obs)))
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Only outer labels
        col = idx % 3
        row = idx // 3
        if col == 0:
            ax.set_ylabel("absorbance [a.u.]")
        if row == 2:
            ax.set_xlabel(r"$t$ [ns]")

    # Single shared legend at the top
    handles = [
        Line2D([0], [0], color=palette_cr2, lw=1.5, label=r"Cr$^{2+}$ + $e_s^-$ posterior median"),
        Line2D([0], [0], color=palette_cr3, lw=1.5, label=r"Cr$^{3+}$ + $e_s^-$ posterior median"),
        Line2D([0], [0], marker="o", color="0.4", lw=0, markersize=4,
               label="training trace (data)"),
        Line2D([0], [0], marker="^", color="0.4", lw=0, markersize=4,
               label="held-out validation (data)"),
        plt.Rectangle((0, 0), 1, 1, fc="0.7", alpha=0.4, label="90% posterior predictive"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, -0.04), frameon=False, fontsize=7.5)

    out_pdf = FIG_DIR / "fig_cr_traces.pdf"
    out_png = FIG_DIR / "fig_cr_traces.png"
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.10)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.10)
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
