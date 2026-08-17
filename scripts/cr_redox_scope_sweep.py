#!/usr/bin/env python3
"""Exploratory chromium redox fate balance for the Paper I review response.

This calculation is deliberately not an HBMAE calibration.  It combines the
four measured 400 degC bimolecular rate coefficients from Iwamatsu et al.
(2026, DOI: 10.1039/D4CP04190A) with a stoichiometric fate balance to show
which assumptions are required before a sustained-irradiation redox shift can
be inferred from pulse-radiolysis data.

For one solvated electron captured by chromium, the chromium oxidation state
falls by one unit.  For one dichlorine radical anion captured by Cr(II), the
oxidation state rises by one unit.  Capture by Cr(III) may be oxidative or may
follow a reductive branch producing Cr(II) and Cl2; that branch fraction is not
identified by the transient data.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

K_E_CR2 = 4.1e10
K_E_CR3 = 6.1e10
K_CL2R_CR2 = 7.2e9
K_CL2R_CR3 = 1.4e9
C_CR2 = 4.0e-3  # M, representative mixed-chromium screening composition
C_CR3 = 4.0e-3  # M

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO / "validation" / "reviewer_response" / "cr_redox_scope_sweep.csv"


def capture_probabilities() -> dict[str, float]:
    """Return pseudo-first-order chromium capture probabilities."""
    electron_total = K_E_CR2 * C_CR2 + K_E_CR3 * C_CR3
    radical_total = K_CL2R_CR2 * C_CR2 + K_CL2R_CR3 * C_CR3
    return {
        "P_e_to_CrII": K_E_CR2 * C_CR2 / electron_total,
        "P_e_to_CrIII": K_E_CR3 * C_CR3 / electron_total,
        "P_Cl2rad_to_CrII": K_CL2R_CR2 * C_CR2 / radical_total,
        "P_Cl2rad_to_CrIII": K_CL2R_CR3 * C_CR3 / radical_total,
    }


def sweep(
    source_ratios: Iterable[float] = (0.5, 1.0, 1.5),
    reductive_branches: Iterable[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
) -> list[dict[str, float]]:
    """Return the stoichiometric redox balance over the requested grid.

    ``source_ratio`` is ``G(Cl2 radical anion) / G(solvated electron)``.
    Negative net oxidation units indicate a net-reducing balance.
    """
    probabilities = capture_probabilities()
    p_x_cr2 = probabilities["P_Cl2rad_to_CrII"]
    p_x_cr3 = probabilities["P_Cl2rad_to_CrIII"]

    rows: list[dict[str, float]] = []
    for source_ratio in source_ratios:
        for reductive_branch in reductive_branches:
            if not 0.0 <= reductive_branch <= 1.0:
                raise ValueError("reductive branch fractions must lie in [0, 1]")
            delta = -1.0 + float(source_ratio) * (
                p_x_cr2
                + p_x_cr3 * ((1.0 - reductive_branch) - reductive_branch)
            )
            rows.append(
                {
                    "G_Cl2rad_over_G_electron": float(source_ratio),
                    "CrIII_Cl2rad_reductive_branch": float(reductive_branch),
                    **probabilities,
                    "net_oxidation_units_per_electron": delta,
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"CSV output path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()

    rows = sweep()
    write_csv(args.output, rows)
    probabilities = capture_probabilities()
    print(
        "P(e_s- capture by CrII)="
        f"{probabilities['P_e_to_CrII']:.4f}; "
        f"by CrIII={probabilities['P_e_to_CrIII']:.4f}"
    )
    print(
        "P(Cl2.- capture by CrII)="
        f"{probabilities['P_Cl2rad_to_CrII']:.4f}; "
        f"by CrIII={probabilities['P_Cl2rad_to_CrIII']:.4f}"
    )
    for row in rows:
        if row["G_Cl2rad_over_G_electron"] == 1.0:
            delta = row["net_oxidation_units_per_electron"]
            direction = "reduction" if delta < 0 else ("oxidation" if delta > 0 else "neutral")
            print(
                "equal primary yields, reductive branch="
                f"{row['CrIII_Cl2rad_reductive_branch']:.2f}: "
                f"delta={delta:+.4f} ({direction})"
            )
    print(args.output)


if __name__ == "__main__":
    main()
