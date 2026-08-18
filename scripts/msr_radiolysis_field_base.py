#!/usr/bin/env python3
"""Reactor-scale and STAR-CCM+ field-resolved molten-salt radiolysis studies.

The calculations use the reaction and G-value parameters currently encoded in
``msr_radiolysis/data/database.yaml`` together with the low-order off-gas
closures used in ``prospective_msr_operational_studies.py``.  Two deliberately
separate products are computed:

1. fast local intermediate-species fields from the repository kinetic kernels;
2. slow stable-halogen and uranium-redox envelopes, for which the available
   evidence does not identify all branching, survival, and regeneration terms.

The STAR-CCM+ temperature, density, liquid fraction, heat-source, and plenum
fields are read directly from the supplied STAR-CCM+ 20.06.010 ``.sim`` files by
``starccm_sim_extract.py``.  The CFD state is frozen while chemistry evolves;
the saved files contain solution fields rather than a resolved irradiation
history.

The script writes CSV/NPZ results and publication figures.  It does not promote
placeholder Henry coefficients, scavenging lifetimes, or net-redox branching
parameters to validated design values.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.collections import PolyCollection
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from starccm_sim_types import SaltMesh

R_GAS = 8.31446261815324
FARADAY = 96485.33212
N_A = 6.02214076e23
EV_J = 1.602176634e-19
G_TO_MOL_PER_J = 1.0 / (100.0 * EV_J * N_A)
SECONDS_PER_DAY = 86400.0
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY
HOURS_PER_YEAR = 365.25 * 24.0
GW_DAY_J = 1.0e9 * SECONDS_PER_DAY

# Atomic and molecular masses [kg/mol].
M_U = 238.02891e-3
M_CL = 35.45e-3
M_F = 18.99840316273e-3
M_NA = 22.98976928e-3
M_K = 39.0983e-3
M_UCL3 = M_U + 3.0 * M_CL
M_UF3 = M_U + 3.0 * M_F
M_UF4 = M_U + 4.0 * M_F
M_NACL = M_NA + M_CL
M_KCL = M_K + M_CL
M_NAF = M_NA + M_F
M_KF = M_K + M_F

# Repository G-values [molecules / 100 eV].
G_E_CHLORIDE = 0.30
G_CL_RADICAL = 0.30
G_E_FLUORIDE = 0.20
G_F_RADICAL = 0.20
G_F2_EMPIRICAL = 0.005

# Article-II / Phillips-derived effective stable Cl2 yield.  The product bound
# is eta*f_rec = 1.670782e-4 at G(Cl atom)=0.2.  One Cl2 contains two Cl atoms.
G_CL_ATOM_SCREEN = 0.20
ETA_RECOVERY_PRODUCT_BOUND = 1.670782136734233e-4
G_CL2_EFFECTIVE = 0.5 * G_CL_ATOM_SCREEN * ETA_RECOVERY_PRODUCT_BOUND

# Placeholder Henry coefficients currently encoded in database.yaml.
KH_CL2 = 1.0e-5  # mol m^-3 Pa^-1
KH_F2 = 1.0e-6   # mol m^-3 Pa^-1

# Intermediate closure parameters.  These represent unresolved pseudo-first-
# order capture by uranium couples, impurities, walls, and other scavengers.
TAU_E_REF_S = 1.0e-4
TAU_OX_REF_S = 1.0e-4
BETA_REDOX_CAPSULE = 1.0e-5  # net oxidizing equivalents / primary pair

# Plot style: conventional, restrained, and compatible with journal printing.
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
        "font.size": 8.3,
        "axes.labelsize": 8.3,
        "axes.titlesize": 8.5,
        "legend.fontsize": 7.2,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.25,
        "savefig.dpi": 450,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass(frozen=True)
class SaltComposition:
    key: str
    label: str
    kernel: str
    density_user_kg_m3: float
    uranium_compound_mass_fraction: float
    base_na_mole_fraction: float = 0.5
    uf3_mole_fraction: float = 0.10

    @property
    def stable_g_value(self) -> float:
        return G_CL2_EFFECTIVE if self.kernel == "chloride" else G_F2_EMPIRICAL

    @property
    def radical_g_value(self) -> float:
        return G_CL_RADICAL if self.kernel == "chloride" else G_F_RADICAL

    @property
    def henry(self) -> float:
        return KH_CL2 if self.kernel == "chloride" else KH_F2

    def composition_state(self, density_kg_m3: float | np.ndarray | None = None) -> dict[str, np.ndarray]:
        rho = np.asarray(
            self.density_user_kg_m3 if density_kg_m3 is None else density_kg_m3,
            dtype=float,
        )
        if np.any(rho <= 0.0):
            raise ValueError("density must be positive")
        if self.kernel == "chloride":
            compound_moles = rho * self.uranium_compound_mass_fraction / M_UCL3
            base_mass = rho * (1.0 - self.uranium_compound_mass_fraction)
            base_molar_mass = (
                self.base_na_mole_fraction * M_NACL
                + (1.0 - self.base_na_mole_fraction) * M_KCL
            )
            base_moles = base_mass / base_molar_mass
            total_u = compound_moles
            halide = base_moles + 3.0 * compound_moles
            u_mass = total_u * M_U
            nominal_u3_fraction = np.full_like(rho, 0.99, dtype=float)
        else:
            uranium_salt_molar_mass = (
                self.uf3_mole_fraction * M_UF3
                + (1.0 - self.uf3_mole_fraction) * M_UF4
            )
            compound_moles = (
                rho * self.uranium_compound_mass_fraction / uranium_salt_molar_mass
            )
            base_mass = rho * (1.0 - self.uranium_compound_mass_fraction)
            base_molar_mass = (
                self.base_na_mole_fraction * M_NAF
                + (1.0 - self.base_na_mole_fraction) * M_KF
            )
            base_moles = base_mass / base_molar_mass
            total_u = compound_moles
            halide = base_moles + compound_moles * (
                3.0 * self.uf3_mole_fraction + 4.0 * (1.0 - self.uf3_mole_fraction)
            )
            u_mass = total_u * M_U
            nominal_u3_fraction = np.full_like(rho, self.uf3_mole_fraction, dtype=float)
        return {
            "total_u_mol_m3": total_u,
            "halide_mol_m3": halide,
            "u_mass_kg_m3": u_mass,
            "nominal_u3_fraction": nominal_u3_fraction,
        }


CHLORIDE = SaltComposition(
    key="chloride",
    label=r"NaCl--KCl--UCl$_3$ (8 wt\% UCl$_3$)",
    kernel="chloride",
    density_user_kg_m3=2650.0,
    uranium_compound_mass_fraction=0.08,
)
FLUORIDE = SaltComposition(
    key="fluoride",
    label=r"NaF--KF--UF$_4$/UF$_3$ (10 wt\% uranium fluorides)",
    kernel="fluoride",
    density_user_kg_m3=4142.0,
    uranium_compound_mass_fraction=0.10,
    uf3_mole_fraction=0.10,
)


@dataclass(frozen=True)
class CapsuleCase:
    key: str
    file_name: str
    capsule: str
    composition: SaltComposition
    user_power_W_cc: float
    burnup_3_GWd_tU: float
    gas_3_mol: float
    burnup_4_GWd_tU: float
    gas_4_mol: float


CAPSULE_CASES = (
    CapsuleCase(
        "chloride_I2",
        "SABRE_AFC_Kaeri_I2_v2.sim",
        "I2",
        CHLORIDE,
        33.0,
        3.74,
        1.1e-5,
        4.96,
        1.4e-4,
    ),
    CapsuleCase(
        "chloride_I3",
        "SABRE_AFC_Kaeri_I3_v2.sim",
        "I3",
        CHLORIDE,
        33.8,
        3.77,
        1.1e-5,
        5.00,
        1.4e-4,
    ),
    CapsuleCase(
        "fluoride_I2",
        "SABRE_Saltfoss_I2_v2.sim",
        "I2",
        FLUORIDE,
        39.2,
        2.51,
        8.9e-5,
        3.31,
        1.1e-4,
    ),
    CapsuleCase(
        "fluoride_I3",
        "SABRE_Saltfoss_I3_v2.sim",
        "I3",
        FLUORIDE,
        40.0,
        3.61,
        1.0e-5,
        3.61,
        1.3e-4,
    ),
)


def g_source(G: float, power_density_W_m3: float | np.ndarray) -> np.ndarray:
    """Convert G and deposited power to a molar source [mol m^-3 s^-1]."""
    if G < 0.0:
        raise ValueError("G must be nonnegative")
    q = np.asarray(power_density_W_m3, dtype=float)
    if np.any(q < 0.0):
        raise ValueError("power density must be nonnegative")
    return G * q * G_TO_MOL_PER_J


def arrhenius_from_reference(
    k_ref: float,
    T_ref_K: float,
    Ea_J_mol: float,
    temperature_K: float | np.ndarray,
) -> np.ndarray:
    temperature = np.asarray(temperature_K, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("temperature must be positive")
    return k_ref * np.exp(Ea_J_mol / R_GAS * (1.0 / T_ref_K - 1.0 / temperature))


def _positive_quadratic_root(source: np.ndarray, linear: np.ndarray, quadratic: np.ndarray) -> np.ndarray:
    """Stable positive root of source - linear*x - quadratic*x^2 = 0."""
    source, linear, quadratic = np.broadcast_arrays(
        np.asarray(source, float), np.asarray(linear, float), np.asarray(quadratic, float)
    )
    out = np.empty_like(source)
    linear_only = quadratic <= 0.0
    out[linear_only] = np.divide(
        source[linear_only],
        linear[linear_only],
        out=np.zeros_like(source[linear_only]),
        where=linear[linear_only] > 0.0,
    )
    idx = ~linear_only
    disc = np.sqrt(np.maximum(linear[idx] ** 2 + 4.0 * quadratic[idx] * source[idx], 0.0))
    out[idx] = np.divide(
        2.0 * source[idx],
        linear[idx] + disc,
        out=np.zeros_like(source[idx]),
        where=(linear[idx] + disc) > 0.0,
    )
    return np.maximum(out, 0.0)



__all__ = [name for name in globals() if not name.startswith("__")]
