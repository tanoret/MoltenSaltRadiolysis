
import numpy as np

NA = 6.02214076e23   # 1/mol
EV_TO_J = 1.602176634e-19  # J/eV

def g_to_source(G_molecules_per_100eV: float, dose_rate_J_m3_s: float) -> float:
    """
    Convert a G-value (molecules/100eV) and volumetric dose rate (J/m^3/s)
    into a molar source term S [mol/m^3/s].
    S = G * dose * (1 / (100 eV)) * (1 / N_A)

    Returns S in mol/m^3/s.
    """
    if G_molecules_per_100eV is None:
        return 0.0
    factor = 1.0 / (100.0 * EV_TO_J) / NA  # mol/J per (molecule/100eV)
    return float(G_molecules_per_100eV) * dose_rate_J_m3_s * factor

def M_to_SI_k(k_M_inv_s: float) -> float:
    """Convert a rate constant from (L/mol/s) = M^-1 s^-1 to SI (m^3/mol/s)."""
    return float(k_M_inv_s) * 1.0e-3

def SI_to_M_k(k_m3_mol_s: float) -> float:
    """Convert from SI to L/mol/s (M^-1 s^-1)."""
    return float(k_m3_mol_s) * 1.0e3
