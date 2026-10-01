from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prospective_msr_operational_studies.py"
SPEC = importlib.util.spec_from_file_location("prospective_msr_operational_studies", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {SCRIPT}")
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class ProspectiveOperationalStudiesTests(unittest.TestCase):
    def test_Phillips_bound_makes_released_source_G_invariant(self) -> None:
        sources = []
        for G in [0.05, 0.2, 1.0]:
            plant = MOD.ChloridePlant(G_cl_atoms_per_100eV=G)
            sources.append(MOD.chloride_stable_cl2_source_mol_h(plant, 0.1))
        self.assertTrue(np.allclose(sources, sources[0], rtol=1.0e-12, atol=0.0))

    def test_no_removal_reproduces_article_II_60_year_pressure(self) -> None:
        plant = MOD.ChloridePlant()
        source = MOD.chloride_stable_cl2_source_mol_h(plant, 1.0)
        pressure = source * 60.0 * MOD.HOURS_PER_YEAR / plant.inventory_capacitance_mol_Pa
        self.assertAlmostEqual(pressure, 728090.9378164454, delta=1.0e-5)

    def test_cleanup_half_life_requirement(self) -> None:
        plant = MOD.ChloridePlant()
        slope = MOD.pressure_source_rate_Pa_h(plant, 1.0)
        half_life = MOD.maximum_half_life_h_for_target(slope, 100.0)
        self.assertAlmostEqual(half_life, 50.0724, delta=0.01)

    def test_U3_capacity_regression(self) -> None:
        plant = MOD.ChloridePlant()
        inv = MOD.U3Inventory()
        depletion = MOD.time_to_U3_depletion_years(plant, inv, 1.0, 0.0)
        self.assertAlmostEqual(depletion, 1.273872556539177, delta=1.0e-10)
        required = MOD.required_regeneration_fraction_for_horizon(plant, inv, 60.0, 1.0)
        self.assertAlmostEqual(required, 0.978768779070462, delta=2.0e-8)

    def test_flibe_static_temperature_regression(self) -> None:
        expected = {500.0: 62.10434733892321, 600.0: 31.0, 700.0: 17.849193993691838}
        for temperature_C, expected_pressure in expected.items():
            pressure = MOD.flibe_pressure_source_Pa_h(1.0) / MOD.k_rec_flibe_h(temperature_C + 273.15)
            self.assertAlmostEqual(pressure, expected_pressure, delta=2.0e-9)

    def test_constant_flibe_state_is_steady(self) -> None:
        t = np.linspace(0.0, 24.0, 241)
        pressure = MOD.simulate_flibe_pressure(
            t,
            dose_fraction=lambda _t: 1.0,
            temperature_K=lambda _t: MOD.T_FLIBE_REF_K,
            P0_Pa=MOD.P_FLIBE_REF,
        )
        self.assertLess(float(np.max(np.abs(pressure - MOD.P_FLIBE_REF))), 1.0e-8)

    def test_linear_inventory_matches_analytic_solution(self) -> None:
        t = np.linspace(0.0, 10.0, 201)
        source = 2.0
        removal = 0.5
        initial = 1.0
        numerical = MOD.simulate_linear_inventory(
            t,
            lambda _t: source,
            lambda _t: removal,
            initial,
        )
        exact = source / removal + (initial - source / removal) * np.exp(-removal * t)
        self.assertLess(float(np.max(np.abs(numerical - exact))), 2.0e-7)


if __name__ == "__main__":
    unittest.main()
