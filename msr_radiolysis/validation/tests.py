
import unittest
import numpy as np

from ..builder import build_system
from ..integrator import integrate_system
from ..utils import g_to_source

class TestRadiolysisFramework(unittest.TestCase):

    def test_g_value_conversion(self):
        # Compare against explicit formula in reference_validation.yaml
        Edot = 1.0e6
        G = 0.005
        S = g_to_source(G, Edot)
        expected = 0.005 * 1.0e6 / (100.0 * 1.602176634e-19) / 6.02214076e23
        self.assertAlmostEqual(S, expected, delta=abs(expected)*1e-12)

    def test_pseudo_first_order_electron_capture(self):
        # ODE: d[e]/dt = -k [e] [Zn2+]   -> for constant [Zn2+], e(t)=e0*exp(-k[Zn2+]t)
        cfg = {
            "kernel":"chloride",
            "temperature_K":673.15,
            "liquid_volume_m3":1.0e-3,
            "headspace_volume_m3":0.0,
            "kLa_s^-1":0.0,
            "dose_rate_J_m3_s":0.0,
            "radiation":"gamma",
            "metals":{
                "Zn":{"Zn2+":100.0,"Zn+":0.0}
            },
            "initial_concentrations":{"e_s-":1.0, "Cl-":5000.0},
            "gas_species":[]
        }
        system = build_system(cfg)
        # modify: keep only electron capture reaction active
        system.reactions = [r for r in system.reactions if "e_s- + Zn2+ -> Zn+" in r.name]
        t, C, extra = integrate_system(system, t_final=2.0e-10, n_steps=5)
        e0 = 1.0
        e_end = C[-1, system.species_index["e_s-"]]
        # analytic with k derived from current database.yaml Arrhenius:
        #   e_s- + Zn2+ -> Zn+ uses A and Ea, so k(T) = A * exp(-Ea / R T).
        rxn = [r for r in system.reactions if "e_s- + Zn2+ -> Zn+" in r.name][0]
        k = rxn.k_forward(system.T)
        Zn = 100.0
        ttot = t[-1]
        expected = e0*np.exp(-k*Zn*ttot)
        self.assertAlmostEqual(e_end, expected, delta=abs(expected)*5e-2)

    def test_second_order_disproportionation(self):
        # 2 Cl2•- -> products, C(t) = C0/(1 + k*C0*t)
        cfg = {
            "kernel":"chloride",
            "temperature_K":673.15,
            "liquid_volume_m3":1.0e-3,
            "headspace_volume_m3":0.0,
            "kLa_s^-1":0.0,
            "dose_rate_J_m3_s":0.0,
            "radiation":"gamma",
            "metals":{},
            "initial_concentrations":{"Cl2•-":1.0},
            "gas_species":[]
        }
        system = build_system(cfg)
        # keep only the disproportionation reaction
        system.reactions = [r for r in system.reactions if "2 Cl2•- -> Cl3- + Cl-" in r.name]
        # integrate to reach half-life
        k = 2.2e6
        C0 = 1.0
        t_half = 1.0/(2.0*k*C0)
        t, C, extra = integrate_system(system, t_final=t_half, n_steps=5)
        C_end = C[-1, system.species_index["Cl2•-"]]
        self.assertAlmostEqual(C_end, C0/2.0, delta=C0*5e-2)

    def test_gas_partition_equilibrium(self):
        # Build a minimal system with only gas exchange and no chemistry
        cfg = {
            "kernel":"chloride",
            "temperature_K":673.15,
            "liquid_volume_m3":1.0e-3,
            "headspace_volume_m3":1.0e-3,
            "kLa_s^-1":10.0,
            "dose_rate_J_m3_s":0.0,
            "radiation":"gamma",
            "metals":{},
            "initial_concentrations":{"Cl2_diss":0.0},
            "gas_species":["Cl2"]
        }
        system = build_system(cfg)
        # initialize gas moles to a known value
        idx_g = system.species_index["Cl2_g"]
        y0 = system.initial_concentrations.copy()
        y0[idx_g] = 1.0e-3  # mol
        system.initial_concentrations = y0
        t, C, extra = integrate_system(system, t_final=10.0, n_steps=40)
        C_liq = C[-1, system.species_index["Cl2_diss"]]
        kH = 1.0e-5
        R = system.R_gas
        T = system.T
        Vg = 1.0e-3
        Vl = 1.0e-3
        n_total = 1.0e-3  # mol (initially all in gas)
        expected = (kH * R * T * n_total / Vg) / (1.0 + kH * R * T * Vl / Vg)
        self.assertAlmostEqual(C_liq, expected, delta=abs(expected)*1e-2)

    def test_F2_yield_linear_accumulation(self):
        # With only F2 source and no losses or gas phase
        cfg = {
            "kernel":"fluoride",
            "temperature_K":673.15,
            "liquid_volume_m3":1.0e-3,
            "headspace_volume_m3":0.0,
            "kLa_s^-1":0.0,
            "dose_rate_J_m3_s":1.0e6,
            "radiation":"gamma",
            "G_values_override":{"F2_diss":0.005, "e_s-":0.0, "F•":0.0},
            "metals":{},
            "initial_concentrations":{},
            "gas_species":[]
        }
        system = build_system(cfg)
        # disable fluoride reactions to isolate the source term
        system.reactions = []
        t_final = 10.0
        t, C, extra = integrate_system(system, t_final=t_final, n_steps=10)
        C_end = C[-1, system.species_index["F2_diss"]]
        # expected concentration
        expected = (0.005 * 1.0e6 / (100.0 * 1.602176634e-19) / 6.02214076e23) * t_final
        self.assertAlmostEqual(C_end, expected, delta=abs(expected)*1e-12)

if __name__ == "__main__":
    unittest.main()