"""Regression tests for the reviewer-response radiolysis calculations."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from cr_redox_scope_sweep import sweep
from paper2_flibe_static_scenario import FlibeScenario, pressure_at_temperature
from paper2_mcfr_recovery_aware import compute_analysis
from worked_example_predictive_check import (
    arrhenius_correlation,
    posterior_parameters,
    predictive_quantiles,
)


def test_worked_example_arrhenius_geometry_and_prediction() -> None:
    mean, covariance = posterior_parameters()
    np.testing.assert_allclose(np.exp(mean[0]), 1.66493e13, rtol=2e-6)
    # The previous manuscript wording incorrectly called this an anti-correlation.
    assert arrhenius_correlation(covariance) > 0.83

    rows = predictive_quantiles((500.0,))
    assert len(rows) == 1
    np.testing.assert_allclose(
        rows[0]["k_median_M-1_s-1"], 8.94557e10, rtol=2e-6
    )


def test_chromium_scope_sweep_requires_unmeasured_branching_assumption() -> None:
    rows = sweep(source_ratios=(1.0,), reductive_branches=(0.0, 1.0))
    np.testing.assert_allclose(rows[0]["net_oxidation_units_per_electron"], 0.0)
    assert rows[1]["net_oxidation_units_per_electron"] < -0.32


def test_mcfr_recovery_bound_and_u3_capacity() -> None:
    analysis = compute_analysis()
    np.testing.assert_allclose(
        analysis["nominal_product_bound"], 1.6708e-4, rtol=5e-5
    )
    p60 = np.asarray(analysis["pressure_with_henry"], dtype=float)[:, -1]
    np.testing.assert_allclose(
        p60, np.array([7.2809e5, 7.2809e6, 7.2809e7]), rtol=5e-5
    )
    np.testing.assert_allclose(analysis["capacity_ratio"], 47.10, rtol=2e-3)
    np.testing.assert_allclose(
        analysis["one_pass_depletion_years"], 1.27, rtol=5e-3
    )


def test_flibe_static_scenario_temperature_scaling() -> None:
    scenario = FlibeScenario()
    pressures = pressure_at_temperature(
        np.array([500.0, 600.0, 700.0]), scenario.median_pressure_pa, scenario
    )
    np.testing.assert_allclose(pressures, [62.10, 31.0, 17.85], rtol=5e-4)
