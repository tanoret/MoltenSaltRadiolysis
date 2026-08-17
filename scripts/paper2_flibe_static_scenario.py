#!/usr/bin/env python3
"""Conditional FLiBe--UF4 static F2 screening scenario.

The display model is

    P_ss(T) proportional to G(F2) / k_rec(T),
    k_rec(T) = A exp[-E_a / (R T)].

The 873.15 K scenario is anchored at 31 Pa with a reported 90% interval
[19, 44] Pa.  The same nominal Arrhenius temperature factor
``E_a = 39 kJ mol^-1`` displays temperature dependence.  Structural
uncertainty from extrapolating the Toth--Felker mechanism from 423 K to reactor
temperature is not represented by the interval and must be interpreted
separately.  The result is a conditional screening scenario, not a reactor
forecast or design margin.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

R_GAS = 8.314462618
REPO = Path(__file__).resolve().parent.parent
DEFAULT_FIGURE_DIR = REPO / "manuscript" / "figures"
DEFAULT_RESULTS_DIR = REPO / "validation" / "reviewer_response"


@dataclass(frozen=True)
class FlibeScenario:
    activation_energy_j_mol: float = 39.0e3
    reference_temperature_k: float = 873.15
    median_pressure_pa: float = 31.0
    lower_pressure_pa: float = 19.0
    upper_pressure_pa: float = 44.0
    screening_pressure_pa: float = 100.0
    lifetime_years: float = 60.0


def pressure_at_temperature(
    temperature_c: np.ndarray | float,
    reference_pressure_pa: float,
    scenario: FlibeScenario = FlibeScenario(),
) -> np.ndarray:
    """Scale a reference pressure with the nominal Arrhenius temperature factor."""
    temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
    factor = np.exp(
        scenario.activation_energy_j_mol
        / R_GAS
        * (1.0 / temperature_k - 1.0 / scenario.reference_temperature_k)
    )
    return float(reference_pressure_pa) * factor


def compute_scenario(
    scenario: FlibeScenario = FlibeScenario(),
) -> dict[str, np.ndarray | FlibeScenario]:
    """Return lifetime and temperature arrays for the static scenario."""
    years = np.linspace(0.0, scenario.lifetime_years, 241)
    temperatures_c = np.linspace(50.0, 700.0, 500)
    return {
        "scenario": scenario,
        "years": years,
        "pressure_median_time": np.full_like(years, scenario.median_pressure_pa),
        "pressure_lower_time": np.full_like(years, scenario.lower_pressure_pa),
        "pressure_upper_time": np.full_like(years, scenario.upper_pressure_pa),
        "temperatures_c": temperatures_c,
        "pressure_median_temperature": pressure_at_temperature(
            temperatures_c, scenario.median_pressure_pa, scenario
        ),
        "pressure_lower_temperature": pressure_at_temperature(
            temperatures_c, scenario.lower_pressure_pa, scenario
        ),
        "pressure_upper_temperature": pressure_at_temperature(
            temperatures_c, scenario.upper_pressure_pa, scenario
        ),
    }


def write_summary_csv(path: Path, scenario: FlibeScenario = FlibeScenario()) -> None:
    """Write the anchor values and selected temperature extrapolations."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = [
        {
            "quantity": "P_F2_median_600C",
            "value": scenario.median_pressure_pa,
            "units_or_note": "Pa; conditional static scenario",
        },
        {
            "quantity": "P_F2_90pct_low_600C",
            "value": scenario.lower_pressure_pa,
            "units_or_note": "Pa",
        },
        {
            "quantity": "P_F2_90pct_high_600C",
            "value": scenario.upper_pressure_pa,
            "units_or_note": "Pa",
        },
        {
            "quantity": "E_rec",
            "value": scenario.activation_energy_j_mol / 1000.0,
            "units_or_note": "kJ/mol; Toth-Felker screening anchor",
        },
    ]
    for temperature_c in (500.0, 600.0, 700.0):
        pressure = float(
            pressure_at_temperature(
                temperature_c, scenario.median_pressure_pa, scenario
            )
        )
        rows.append(
            {
                "quantity": f"P_F2_median_{int(temperature_c)}C",
                "value": pressure,
                "units_or_note": "Pa; conditional Arrhenius extrapolation",
            }
        )

    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["quantity", "value", "units_or_note"]
        )
        writer.writeheader()
        writer.writerows(rows)


def make_figure(
    pdf_path: Path,
    png_path: Path,
    data: dict[str, np.ndarray | FlibeScenario],
) -> None:
    """Create the lifetime and temperature display for the static scenario."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scenario = data["scenario"]
    assert isinstance(scenario, FlibeScenario)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)

    ax = axes[0]
    ax.fill_between(
        data["years"],
        data["pressure_lower_time"],
        data["pressure_upper_time"],
        alpha=0.25,
        label="scenario 90% interval",
    )
    ax.plot(
        data["years"],
        data["pressure_median_time"],
        linewidth=1.8,
        label="scenario median",
    )
    ax.axhline(
        scenario.screening_pressure_pa,
        linestyle=":",
        linewidth=1.2,
        color="0.35",
        label="100 Pa screening anchor",
    )
    ax.set_yscale("log")
    ax.set_xlim(0.0, scenario.lifetime_years)
    ax.set_ylim(1.0, 1.0e3)
    ax.set_xlabel("operating years")
    ax.set_ylabel(r"conditional $P_{\mathrm{F}_2}$ [Pa]")
    ax.set_title(r"(a) Static scenario at 600 $^\circ$C")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    ax.fill_between(
        data["temperatures_c"],
        data["pressure_lower_temperature"],
        data["pressure_upper_temperature"],
        alpha=0.25,
        label="scenario 90% interval",
    )
    ax.plot(
        data["temperatures_c"],
        data["pressure_median_temperature"],
        linewidth=1.8,
        label="scenario median",
    )
    ax.axhline(
        scenario.screening_pressure_pa,
        linestyle=":",
        linewidth=1.2,
        color="0.35",
        label="100 Pa screening anchor",
    )
    ax.axvline(
        150.0,
        linestyle="--",
        linewidth=1.0,
        color="0.45",
        label="423 K experimental anchor",
    )
    ax.set_yscale("log")
    ax.set_xlim(50.0, 700.0)
    ax.set_ylim(10.0, 1.0e6)
    ax.set_xlabel(r"temperature [$^\circ$C]")
    ax.set_ylabel(r"conditional $P_{\mathrm{F}_2}^{\rm ss}$ [Pa]")
    ax.set_title("(b) Arrhenius temperature extrapolation")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=6.7, frameon=False)

    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=DEFAULT_FIGURE_DIR,
        help=f"figure output directory (default: {DEFAULT_FIGURE_DIR})",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"CSV output directory (default: {DEFAULT_RESULTS_DIR})",
    )
    parser.add_argument(
        "--no-figure",
        action="store_true",
        help="write only the CSV audit table",
    )
    args = parser.parse_args()

    scenario = FlibeScenario()
    data = compute_scenario(scenario)
    csv_path = args.results_dir / "flibe_static_scenario_summary.csv"
    pdf_path = args.figure_dir / "fig_flibe_f2_lifetime.pdf"
    png_path = args.figure_dir / "fig_flibe_f2_lifetime.png"

    write_summary_csv(csv_path, scenario)
    if not args.no_figure:
        make_figure(pdf_path, png_path, data)

    print(
        "FLiBe-UF4 conditional F2 scenario at 600 C: "
        f"{scenario.median_pressure_pa:g} Pa "
        f"[{scenario.lower_pressure_pa:g}, {scenario.upper_pressure_pa:g}]"
    )
    for temperature_c in (500.0, 600.0, 700.0):
        pressure = float(
            pressure_at_temperature(
                temperature_c, scenario.median_pressure_pa, scenario
            )
        )
        print(f"T={temperature_c:.0f} C: P_F2={pressure:.4g} Pa")
    print(csv_path)
    if not args.no_figure:
        print(pdf_path)


if __name__ == "__main__":
    main()
