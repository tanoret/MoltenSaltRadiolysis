from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "starccm_sim_extract.py"
SPEC = importlib.util.spec_from_file_location("starccm_sim_extract", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {MODULE_PATH}")
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class PartitionMapTests(unittest.TestCase):
    def test_ghost_owner_mapping(self) -> None:
        mapping = MOD.PartitionMap(
            n_element=np.array([3, 3]),
            n_interior=np.array([2, 2]),
            ghosts=(((1, 1),), ((0, 0),)),
        )
        self.assertEqual(mapping.global_size, 4)
        self.assertEqual(mapping.global_index(0, 0), 0)
        self.assertEqual(mapping.global_index(1, 1), 3)
        self.assertEqual(mapping.global_index(0, 2), 3)
        self.assertEqual(mapping.global_index(1, 2), 0)

    def test_axisymmetric_rectangle_volume(self) -> None:
        points = np.array([[0.0, 1.0], [2.0, 1.0], [2.0, 3.0], [0.0, 3.0]])
        centroid, area, volume = MOD.StarSimReader._polygon_geometry(points)
        np.testing.assert_allclose(centroid, [1.0, 2.0], rtol=0.0, atol=1.0e-14)
        self.assertAlmostEqual(area, 4.0)
        self.assertAlmostEqual(volume, 16.0 * math.pi)


class SuppliedSimIntegrationTests(unittest.TestCase):
    def test_supplied_salt_mesh_when_available(self) -> None:
        path = Path("/mnt/data/SABRE_Saltfoss_I2_v2.sim")
        if not path.exists():
            self.skipTest("supplied STAR-CCM+ integration file is not available")
        mesh = MOD.StarSimReader(path).extract_salt_mesh("Salt")
        self.assertEqual(mesh.n_cells, 2400)
        self.assertAlmostEqual(mesh.total_volume_m3, 1.2350897003267723e-5, delta=1.0e-16)
        self.assertAlmostEqual(
            mesh.volume_average(mesh.fields["Temperature"]),
            1035.3208649611545,
            delta=1.0e-9,
        )
        self.assertTrue(np.allclose(mesh.fields["Solidity"], 0.0))


if __name__ == "__main__":
    unittest.main()
