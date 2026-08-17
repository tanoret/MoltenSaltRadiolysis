#!/usr/bin/env python3
"""Recovery-aware reanalysis of the Phillips NaCl--UCl3 non-detection.

This reviewer-driven model replaces the earlier calculation that treated the
nominal approximately 1000 ppm value as a direct upper bound on total Cl2
generation and hard-coded an effective leakage efficiency of 1e-18.  The
Phillips report documents a reactive analytical train: a 100 ppm Cl2/Ar
standard was not recovered by gas mass spectrometry, and the UV signal decayed
in the White cell.  The observation therefore constrains the product

    (gas-release fraction) x (sampling/recovery fraction),

not the gas-release fraction alone and not the total amount of Cl2 generated.

The operational propagation is intentionally a transparent linear
accumulated-dose stress test.  It does not infer the unmeasured U(III)/U(IV)
rate constants or claim a reactor forecast.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

R_GAS = 8.314462618
AVOGADRO = 6.02214076e23
ELECTRON_VOLT_J = 1.602176634e-19
ENERGY_PER_100_EV_J = 100.0 * ELECTRON_VOLT_J

REPO = Path(__file__).resolve().parent.parent
DEFAULT_FIGURE_DIR = REPO / "manuscript" / "figures"
DEFAULT_RESULTS_DIR = REPO / "validation" / "reviewer_response"


@dataclass(frozen=True)
class McfrInputs:
    """Inputs for the recovery-aware laboratory bound and stress test."""

    # Phillips laboratory screening geometry and accumulated dose.
    laboratory_density_kg_m3: float = 2600.0
    laboratory_salt_volume_m3: float = 5.0e-5
    laboratory_cover_volume_m3: float = 5.0e-4
    laboratory_accumulated_dose_gy: float = 31.0e6
    laboratory_temperature_k: float = 873.15
    analytical_pressure_pa: float = 101325.0
    nominal_detection_mole_fraction: float = 1000.0e-6

    # Operational stress-test geometry.
    operational_density_kg_m3: float = 2700.0
    operational_salt_volume_m3: float = 50.0
    operational_cover_volume_m3: float = 5.0
    operational_temperature_k: float = 873.15
    dose_rate_kgy_h: float = 10.0
    lifetime_years: float = 60.0
    henry_constant_mol_m3_pa: float = 2.0e-5
    screening_pressure_pa: float = 100.0

    # Scenario parameters.
    nominal_g_cl: float = 0.2
    g_values: tuple[float, ...] = (0.05, 0.2, 1.0)
    recovery_cases: tuple[float, ...] = (1.0, 0.1, 0.01)


def max_cl2_moles(g_cl: float, deposited_energy_j: np.ndarray | float) -> np.ndarray:
    """Return molecular-Cl2 equivalents for a primary Cl-atom G-value.

    Two primary chlorine atoms are required per molecular Cl2 equivalent.
    """
    return (
        0.5
        * float(g_cl)
        * np.asarray(deposited_energy_j, dtype=float)
        / ENERGY_PER_100_EV_J
        / AVOGADRO
    )


def compute_analysis(inputs: McfrInputs = McfrInputs()) -> dict[str, object]:
    """Compute the recovery bound, operational stress test, and U(III) capacity.

    Returns a dictionary containing scalar audit values and plotting arrays.
    """
    laboratory_mass_kg = (
        inputs.laboratory_density_kg_m3 * inputs.laboratory_salt_volume_m3
    )
    laboratory_energy_j = (
        inputs.laboratory_accumulated_dose_gy * laboratory_mass_kg
    )
    detectable_cl2_moles = (
        inputs.nominal_detection_mole_fraction
        * inputs.analytical_pressure_pa
        * inputs.laboratory_cover_volume_m3
        / (R_GAS * inputs.laboratory_temperature_k)
    )

    eta_times_recovery_bound = {
        g_value: float(
            detectable_cl2_moles / max_cl2_moles(g_value, laboratory_energy_j)
        )
        for g_value in inputs.g_values
    }
    nominal_product_bound = eta_times_recovery_bound[inputs.nominal_g_cl]
    recovery_cases = np.asarray(inputs.recovery_cases, dtype=float)
    eta_bounds = np.minimum(1.0, nominal_product_bound / recovery_cases)

    years = np.concatenate(
        ([0.0], np.geomspace(1.0e-6, inputs.lifetime_years, 500))
    )
    hours = years * 365.25 * 24.0
    accumulated_dose_gy = inputs.dose_rate_kgy_h * 1.0e3 * hours
    operational_mass_kg = (
        inputs.operational_density_kg_m3 * inputs.operational_salt_volume_m3
    )
    operational_energy_j = accumulated_dose_gy * operational_mass_kg
    unbuffered_cl2_moles = max_cl2_moles(
        inputs.nominal_g_cl, operational_energy_j
    )

    denominator_with_henry = (
        inputs.operational_cover_volume_m3
        / (R_GAS * inputs.operational_temperature_k)
        + inputs.henry_constant_mol_m3_pa
        * inputs.operational_salt_volume_m3
    )
    denominator_gas_only = inputs.operational_cover_volume_m3 / (
        R_GAS * inputs.operational_temperature_k
    )
    pressure_with_henry = np.vstack(
        [eta * unbuffered_cl2_moles / denominator_with_henry for eta in eta_bounds]
    )
    pressure_gas_only = np.vstack(
        [eta * unbuffered_cl2_moles / denominator_gas_only for eta in eta_bounds]
    )

    # One-pass U(III) stoichiometric capacity for a 60/40 mol% NaCl/UCl3 salt.
    molar_mass_nacl_kg_mol = 58.44e-3
    molar_mass_ucl3_kg_mol = 344.36e-3
    mean_formula_mass_kg_mol = (
        0.60 * molar_mass_nacl_kg_mol + 0.40 * molar_mass_ucl3_kg_mol
    )
    formula_moles = operational_mass_kg / mean_formula_mass_kg_mol
    u3_moles = 0.40 * formula_moles
    primary_cl_atom_moles_60y = (
        inputs.nominal_g_cl
        * operational_energy_j[-1]
        / ENERGY_PER_100_EV_J
        / AVOGADRO
    )
    capacity_ratio = float(primary_cl_atom_moles_60y / u3_moles)
    one_pass_depletion_years = (
        inputs.lifetime_years / capacity_ratio if capacity_ratio > 0.0 else np.inf
    )

    return {
        "inputs": inputs,
        "years": years,
        "recovery_cases": recovery_cases,
        "eta_bounds": eta_bounds,
        "eta_times_recovery_bound": eta_times_recovery_bound,
        "nominal_product_bound": nominal_product_bound,
        "pressure_with_henry": pressure_with_henry,
        "pressure_gas_only": pressure_gas_only,
        "detectable_cl2_moles": detectable_cl2_moles,
        "u3_moles": u3_moles,
        "primary_cl_atom_moles_60y": primary_cl_atom_moles_60y,
        "capacity_ratio": capacity_ratio,
        "one_pass_depletion_years": one_pass_depletion_years,
    }


def audit_rows(analysis: dict[str, object]) -> list[dict[str, object]]:
    """Convert an analysis result into an auditable long-form table."""
    inputs = analysis["inputs"]
    assert isinstance(inputs, McfrInputs)
    rows: list[dict[str, object]] = []

    bounds = analysis["eta_times_recovery_bound"]
    assert isinstance(bounds, dict)
    for g_value in inputs.g_values:
        rows.append(
            {
                "quantity": f"eta_times_recovery_upper_bound_G={g_value:g}",
                "value": bounds[g_value],
                "units_or_note": (
                    "dimensionless; conditional on geometry and nominal 1000 ppm"
                ),
            }
        )

    recovery_cases = np.asarray(analysis["recovery_cases"], dtype=float)
    eta_bounds = np.asarray(analysis["eta_bounds"], dtype=float)
    p_henry = np.asarray(analysis["pressure_with_henry"], dtype=float)
    p_gas = np.asarray(analysis["pressure_gas_only"], dtype=float)
    for recovery, eta, pressure_henry, pressure_gas in zip(
        recovery_cases, eta_bounds, p_henry[:, -1], p_gas[:, -1]
    ):
        crossing_years = (
            inputs.lifetime_years * inputs.screening_pressure_pa / pressure_henry
            if pressure_henry > 0.0
            else np.inf
        )
        rows.extend(
            [
                {
                    "quantity": f"eta_upper_bound_recovery={recovery:g}",
                    "value": eta,
                    "units_or_note": f"dimensionless; G={inputs.nominal_g_cl:g}",
                },
                {
                    "quantity": f"P60_with_nominal_KH_recovery={recovery:g}",
                    "value": pressure_henry,
                    "units_or_note": "Pa; conditional linear stress test",
                },
                {
                    "quantity": f"P60_gas_only_recovery={recovery:g}",
                    "value": pressure_gas,
                    "units_or_note": "Pa; conditional linear stress test",
                },
                {
                    "quantity": f"time_to_100Pa_recovery={recovery:g}",
                    "value": crossing_years,
                    "units_or_note": "years; analytic linear-stress-test crossing",
                },
            ]
        )

    rows.extend(
        [
            {
                "quantity": "lab_detectable_Cl2_moles",
                "value": analysis["detectable_cl2_moles"],
                "units_or_note": "mol",
            },
            {
                "quantity": "operational_U3_inventory",
                "value": analysis["u3_moles"],
                "units_or_note": "mol",
            },
            {
                "quantity": "primary_Cl_atom_equivalents_60y",
                "value": analysis["primary_cl_atom_moles_60y"],
                "units_or_note": f"mol; G={inputs.nominal_g_cl:g}",
            },
            {
                "quantity": "primary_Cl_over_one_pass_U3_capacity",
                "value": analysis["capacity_ratio"],
                "units_or_note": "dimensionless",
            },
            {
                "quantity": "one_pass_U3_depletion_time",
                "value": analysis["one_pass_depletion_years"],
                "units_or_note": "years; absent regeneration",
            },
        ]
    )
    return rows


def write_audit_csv(path: Path, analysis: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = audit_rows(analysis)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["quantity", "value", "units_or_note"]
        )
        writer.writeheader()
        writer.writerows(rows)


def make_figure(
    pdf_path: Path,
    png_path: Path,
    analysis: dict[str, object],
) -> None:
    """Create the recovery-bound and operational-stress-test figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    inputs = analysis["inputs"]
    assert isinstance(inputs, McfrInputs)
    bounds = analysis["eta_times_recovery_bound"]
    assert isinstance(bounds, dict)

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    ax = axes[0]
    recovery_grid = np.logspace(-4, 0, 300)
    for g_value in inputs.g_values:
        eta = np.minimum(1.0, bounds[g_value] / recovery_grid)
        ax.plot(
            recovery_grid,
            eta,
            label=fr"$G(\mathrm{{Cl}}^\bullet)={g_value:g}$",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("effective Cl$_2$ recovery fraction")
    ax.set_ylabel(r"conditional upper bound on gas-release fraction, $\eta$")
    ax.set_title("(a) What the non-detection constrains")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(True, which="both", alpha=0.2)

    ax = axes[1]
    years = np.asarray(analysis["years"], dtype=float)
    pressure = np.asarray(analysis["pressure_with_henry"], dtype=float)
    for recovery, curve in zip(inputs.recovery_cases, pressure):
        ax.plot(years, curve, label=fr"recovery $f_{{\rm rec}}={recovery:g}$")
    ax.axhline(
        inputs.screening_pressure_pa,
        linestyle=":",
        linewidth=1.2,
        color="0.35",
        label="100 Pa screening anchor",
    )
    ax.set_yscale("log")
    ax.set_xlim(0.0, inputs.lifetime_years)
    ax.set_ylim(1.0e-2, max(1.0e8, 2.0 * float(np.max(pressure[:, -1]))))
    ax.set_xlabel(
        f"operating years at {inputs.dose_rate_kgy_h:g} kGy h$^{{-1}}$ (scenario)"
    )
    ax.set_ylabel(r"conditional $P_{\mathrm{Cl}_2}$ [Pa]")
    ax.set_title("(b) Linear accumulated-dose stress test")
    ax.legend(fontsize=6.6, frameon=False)
    ax.grid(True, which="both", alpha=0.2)

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

    analysis = compute_analysis()
    csv_path = args.results_dir / "mcfr_recovery_aware_summary.csv"
    pdf_path = args.figure_dir / "fig_mcfr_cl2_lifetime_revised.pdf"
    png_path = args.figure_dir / "fig_mcfr_cl2_lifetime_revised.png"

    write_audit_csv(csv_path, analysis)
    if not args.no_figure:
        make_figure(pdf_path, png_path, analysis)

    inputs = analysis["inputs"]
    assert isinstance(inputs, McfrInputs)
    print(
        "Nominal bound on eta*f_rec "
        f"(G={inputs.nominal_g_cl:g}): "
        f"{analysis['nominal_product_bound']:.4e}"
    )
    recovery_cases = np.asarray(analysis["recovery_cases"], dtype=float)
    eta_bounds = np.asarray(analysis["eta_bounds"], dtype=float)
    pressure = np.asarray(analysis["pressure_with_henry"], dtype=float)
    for recovery, eta, p60 in zip(recovery_cases, eta_bounds, pressure[:, -1]):
        crossing_years = (
            inputs.lifetime_years * inputs.screening_pressure_pa / p60
            if p60 > 0.0
            else np.inf
        )
        print(
            f"f_rec={recovery:g}: eta<={eta:.4e}; P60={p60:.4e} Pa; "
            f"100 Pa at {crossing_years:.4g} y"
        )
    print(
        "U(III) one-pass capacity: primary Cl/U(III)="
        f"{analysis['capacity_ratio']:.2f}; depletion "
        f"~{analysis['one_pass_depletion_years']:.2f} y"
    )
    print(csv_path)
    if not args.no_figure:
        print(pdf_path)


if __name__ == "__main__":
    main()
