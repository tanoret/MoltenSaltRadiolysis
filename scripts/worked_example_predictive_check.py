#!/usr/bin/env python3
"""Posterior-predictive check for Paper I Worked Example I.

The model is the single molten-LiCl--KCl reaction

    e_s^- + Cr(II) -> Cr(I)

using the Iwamatsu et al. (2026) 400 degC rate and the published Arrhenius
pair as the literature prior (DOI: 10.1039/D4CP04190A).  This is an
illustration rather than an independent validation because the scalar datum
and the literature prior come from the same experimental campaign.

The posterior is Gaussian in ``x = [ln(A), E_a]``.  The script evaluates that
posterior analytically, verifies the positive Arrhenius compensation
correlation, and writes exact lognormal posterior-predictive quantiles at
400, 500, and 600 degC.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import numpy as np

R_GAS = 8.314462618  # J mol^-1 K^-1
REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = (
    REPO / "validation" / "reviewer_response" / "worked_example_predictive_check.csv"
)


def posterior_parameters() -> tuple[np.ndarray, np.ndarray]:
    """Return posterior mean and covariance for ``[ln(A), E_a]``.

    ``A`` is in M^-1 s^-1 and ``E_a`` is in J mol^-1.
    """
    prior_mean = np.array([np.log(1.7e13), 33.5e3])
    prior_sd = np.array([0.2 / 1.7, 0.6e3])
    prior_cov = np.diag(prior_sd**2)

    temperature_k = 673.15
    observed_rate = 4.1e10
    sigma_log_rate = 0.05

    design = np.array([[1.0, -1.0 / (R_GAS * temperature_k)]])
    observation = np.array([np.log(observed_rate)])
    observation_cov = np.array([[sigma_log_rate**2]])

    prior_precision = np.linalg.inv(prior_cov)
    observation_precision = design.T @ np.linalg.inv(observation_cov) @ design
    posterior_cov = np.linalg.inv(prior_precision + observation_precision)
    posterior_mean = posterior_cov @ (
        prior_precision @ prior_mean
        + design.T @ np.linalg.inv(observation_cov) @ observation
    )
    return posterior_mean, posterior_cov


def arrhenius_correlation(covariance: np.ndarray) -> float:
    """Return ``corr(ln(A), E_a)`` from a 2x2 covariance matrix."""
    denominator = np.sqrt(covariance[0, 0] * covariance[1, 1])
    return float(covariance[0, 1] / denominator)


def predictive_quantiles(
    temperatures_c: Iterable[float] = (400.0, 500.0, 600.0),
) -> list[dict[str, float]]:
    """Return exact 2.5/50/97.5% posterior-predictive rate quantiles.

    For a Gaussian posterior on ``[ln(A), E_a]``, ``ln(k(T))`` is Gaussian and
    ``k(T)`` is lognormal.  No Monte Carlo sampling is required.
    """
    mean, covariance = posterior_parameters()
    z_975 = 1.959963984540054
    rows: list[dict[str, float]] = []
    for temperature_c in temperatures_c:
        temperature_k = float(temperature_c) + 273.15
        design = np.array([1.0, -1.0 / (R_GAS * temperature_k)])
        mean_log_rate = float(design @ mean)
        variance_log_rate = float(design @ covariance @ design)
        sd_log_rate = np.sqrt(variance_log_rate)
        rows.append(
            {
                "T_C": float(temperature_c),
                "k_2.5pct_M-1_s-1": float(np.exp(mean_log_rate - z_975 * sd_log_rate)),
                "k_median_M-1_s-1": float(np.exp(mean_log_rate)),
                "k_97.5pct_M-1_s-1": float(np.exp(mean_log_rate + z_975 * sd_log_rate)),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    """Write posterior-predictive rows to ``path``."""
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

    mean, covariance = posterior_parameters()
    rows = predictive_quantiles()
    write_csv(args.output, rows)

    print(f"posterior median A = {np.exp(mean[0]):.5e} M^-1 s^-1")
    print(f"posterior median Ea = {mean[1] / 1000.0:.4f} kJ/mol")
    print(f"corr(ln A, Ea) = {arrhenius_correlation(covariance):+.5f}")
    for row in rows:
        print(
            f"T={row['T_C']:.0f} C: k={row['k_median_M-1_s-1']:.5e} "
            f"[{row['k_2.5pct_M-1_s-1']:.5e}, "
            f"{row['k_97.5pct_M-1_s-1']:.5e}]"
        )
    print(args.output)


if __name__ == "__main__":
    main()
