from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "msr_radiolysis_field_studies.py"
SPEC = importlib.util.spec_from_file_location("msr_radiolysis_field_studies", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {MODULE_PATH}")
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class FieldStudyUnitTests(unittest.TestCase):
    def test_G_value_conversion(self) -> None:
        power = 100.0 * MOD.EV_J * MOD.N_A
        self.assertAlmostEqual(float(MOD.g_source(1.0, power)), 1.0, places=13)

    def test_arrhenius_reference_point(self) -> None:
        value = MOD.arrhenius_from_reference(2.2e6, 673.15, 2.6e4, 673.15)
        self.assertAlmostEqual(float(value), 2.2e6, places=7)

    def test_intermediate_balances_close(self) -> None:
        for composition in (MOD.CHLORIDE, MOD.FLUORIDE):
            state = composition.composition_state()
            residuals = MOD.intermediate_balance_residuals(
                composition,
                1000.0,
                40.0e6,
                float(state["halide_mol_m3"]),
                1.0e-4,
                1.0e-4,
            )
            self.assertLess(max(abs(value) for value in residuals.values()), 2.0e-10)

    def test_scavenging_lifetime_controls_intermediate_inventory(self) -> None:
        state = MOD.FLUORIDE.composition_state()
        short = MOD.steady_intermediates(
            MOD.FLUORIDE, 1000.0, 40.0e6, float(state["halide_mol_m3"]), 1.0e-7, 1.0e-7
        )
        long = MOD.steady_intermediates(
            MOD.FLUORIDE, 1000.0, 40.0e6, float(state["halide_mol_m3"]), 1.0e-4, 1.0e-4
        )
        self.assertGreater(float(long["electron_mol_m3"]), float(short["electron_mol_m3"]))
        self.assertGreater(
            float(long["oxidizing_intermediate_mol_m3"]),
            float(short["oxidizing_intermediate_mol_m3"]),
        )

    def test_zero_net_redox_rate_has_zero_potential_shift(self) -> None:
        result = MOD.uranium_redox_state(
            np.array([1000.0]),
            np.array([0.2]),
            np.array([900.0]),
            np.array([0.0]),
            MOD.SECONDS_PER_YEAR,
        )
        self.assertAlmostEqual(float(result["delta_E_V"][0]), 0.0, places=15)
        self.assertAlmostEqual(float(result["u3_fraction_final"][0]), 0.2, places=15)

    def test_burnup_conversion_is_linear(self) -> None:
        t1 = MOD.irradiation_time_s(2.0, 0.003, 500.0)
        t2 = MOD.irradiation_time_s(5.0, 0.003, 500.0)
        self.assertAlmostEqual(t2 / t1, 2.5, places=15)

    def test_composition_interpretation_preserves_requested_compound_fraction(self) -> None:
        chloride = MOD.CHLORIDE.composition_state(2650.0)
        ucl3_moles = float(chloride["total_u_mol_m3"])
        self.assertAlmostEqual(
            ucl3_moles * MOD.M_UCL3 / 2650.0,
            MOD.CHLORIDE.uranium_compound_mass_fraction,
            places=14,
        )
        fluoride = MOD.FLUORIDE.composition_state(4142.0)
        u_moles = float(fluoride["total_u_mol_m3"])
        mean_molar_mass = (
            MOD.FLUORIDE.uf3_mole_fraction * MOD.M_UF3
            + (1.0 - MOD.FLUORIDE.uf3_mole_fraction) * MOD.M_UF4
        )
        self.assertAlmostEqual(
            u_moles * mean_molar_mass / 4142.0,
            MOD.FLUORIDE.uranium_compound_mass_fraction,
            places=14,
        )

    def test_capsule_case_input_table_is_complete(self) -> None:
        self.assertEqual(len(MOD.CAPSULE_CASES), 4)
        keys = {case.key for case in MOD.CAPSULE_CASES}
        self.assertEqual(
            keys,
            {"chloride_I2", "chloride_I3", "fluoride_I2", "fluoride_I3"},
        )


class SuppliedSimCouplingTests(unittest.TestCase):
    def test_one_supplied_case_when_available(self) -> None:
        sim_dir = Path("/mnt/data")
        path = sim_dir / MOD.CAPSULE_CASES[0].file_name
        if not path.exists():
            self.skipTest("supplied STAR-CCM+ integration file is not available")
        data = MOD.extract_capsule_case(MOD.CAPSULE_CASES[0], sim_dir)
        self.assertEqual(data["salt"].n_cells, 2400)
        self.assertTrue(np.allclose(data["liquid_fraction"], 1.0))
        self.assertAlmostEqual(data["salt"].total_volume_m3, 1.3927839083729209e-5, delta=1e-16)
        expected_source = float(
            MOD.g_source(
                MOD.CHLORIDE.stable_g_value,
                MOD.CAPSULE_CASES[0].user_power_W_cc * 1.0e6,
            )
            * data["salt"].total_volume_m3
        )
        self.assertAlmostEqual(data["stable_source_total_mol_s"], expected_source, delta=1e-20)


if __name__ == "__main__":
    unittest.main()
