
from dataclasses import dataclass, field
from typing import Dict, List, Callable, Optional, Tuple
import numpy as np

R_GAS = 8.314462618  # J/mol/K

@dataclass
class Species:
    name: str
    phase: str = "liq"  # "liq" or "gas"
    index: int = -1

@dataclass
class Reaction:
    name: str
    reactants: Dict[str, float]  # stoich coeffs (positive), e.g., {"A":1, "B":1}
    products: Dict[str, float]   # stoich coeffs (positive)
    reversible: bool = False
    # kinetic parameters -- either supply (k_ref, T_ref, Ea) or (A, Ea)
    k_ref: Optional[float] = None         # in SI units (m^3/mol/s for bimolecular; s^-1 for unimolecular)
    T_ref: Optional[float] = None         # K
    Ea_J_mol: Optional[float] = None      # J/mol
    A: Optional[float] = None             # pre-exponential in SI units
    order_overrides: Optional[Dict[str, float]] = None  # if kinetic orders differ from stoich
    # For reversible reactions you can also provide reverse parameters (not used in starter)
    notes: str = ""

    def k_forward(self, T: float) -> float:
        if self.A is not None and self.Ea_J_mol is not None:
            return self.A * np.exp(-self.Ea_J_mol/(R_GAS*T))
        if self.k_ref is not None and self.T_ref is not None and self.Ea_J_mol is not None:
            return float(self.k_ref) * np.exp(self.Ea_J_mol/R_GAS * (1.0/self.T_ref - 1.0/T))
        if self.k_ref is not None:
            # T-independent fallback when only k_ref is given (e.g., paper reports k at one T only)
            return float(self.k_ref)
        return 0.0

    def rate_of_progress(self, T: float, C: np.ndarray, species_index: Dict[str,int]) -> float:
        kf = self.k_forward(T)
        # compute mass-action product of reactant concentrations to their kinetic orders
        if self.order_overrides is None:
            # default: orders equal stoichiometric coefficients
            order_dict = self.reactants
        else:
            order_dict = self.order_overrides
        prod = 1.0
        for s, nu in order_dict.items():
            idx = species_index.get(s, None)
            if idx is None:
                # species not tracked -> treat as large constant (not ideal)
                continue
            conc = C[idx]
            if conc < 0.0:
                conc = 0.0
            prod *= conc**nu
        return kf * prod

@dataclass
class GasExchange:
    """Mass transfer between dissolved and gas species for a diatomic X2.

    Form: R = kLa * (C_liq - C_eq), with C_eq = kH * p_gas.
    kLa in s^-1, kH in mol/(m^3·Pa), p_gas from ideal gas.
    """
    name: str                     # e.g., "Cl2"
    liq_species: str             # e.g., "Cl2_diss"
    gas_species: str             # e.g., "Cl2_g"
    kLa_s: float                 # s^-1
    kH_molm3Pa: float            # mol/(m^3·Pa)

@dataclass
class System:
    species: List[Species]
    species_index: Dict[str,int]
    reactions: List[Reaction]
    sources: Dict[str, float]  # constant source terms [mol/m^3/s]
    T: float
    V_liq_m3: float
    V_gas_m3: float
    gas_exchanges: List[GasExchange] = field(default_factory=list)
    R_gas: float = R_GAS

    def rhs(self, t: float, y: np.ndarray) -> Tuple[np.ndarray, Dict[str,float]]:
        """Compute dC/dt for all species + gas amounts.

        Species in 'gas' phase are stored in 'concentration-like' units for convenience:
        we store n_gas (mol) as y[idx] and convert pressures using ideal gas when needed.
        For liquid species, y[idx] is mol/m^3.
        """
        dydt = np.zeros_like(y)
        extra = {}

        # 1) homogeneous reactions (mass action in liquid phase)
        for rxn in self.reactions:
            r = rxn.rate_of_progress(self.T, y, self.species_index)
            # update species using stoichiometric net coefficients
            for s, nu in rxn.products.items():
                idx = self.species_index[s]
                if self.species[idx].phase == "liq":
                    dydt[idx] += nu * r
                else:
                    # gas species store mol, scale by V_liq
                    dydt[idx] += nu * r * self.V_liq_m3
            for s, nu in rxn.reactants.items():
                idx = self.species_index[s]
                if self.species[idx].phase == "liq":
                    dydt[idx] -= nu * r
                else:
                    dydt[idx] -= nu * r * self.V_liq_m3

        # 2) radiolysis sources (0th order, constant)
        for s, Si in self.sources.items():
            if s not in self.species_index:
                continue
            idx = self.species_index[s]
            if self.species[idx].phase == "liq":
                dydt[idx] += Si
            else:
                dydt[idx] += Si * self.V_liq_m3

        # 3) gas-liquid exchange for diatomics
        for gx in self.gas_exchanges:
            idx_liq = self.species_index[gx.liq_species]
            idx_gas = self.species_index[gx.gas_species]
            C_liq = y[idx_liq]  # mol/m^3
            n_gas = y[idx_gas]  # mol
            p_gas = (n_gas * self.R_gas * self.T) / max(self.V_gas_m3, 1e-30) if self.V_gas_m3 > 0 else 0.0
            C_eq = gx.kH_molm3Pa * p_gas
            R = gx.kLa_s * (C_liq - C_eq)
            # liquid loses if above equilibrium (positive R)
            dCdt_liq = -R
            dngdt = +R * self.V_liq_m3
            dydt[idx_liq] += dCdt_liq
            dydt[idx_gas] += dngdt
            extra[f"p_{gx.name}_Pa"] = p_gas
        return dydt, extra
