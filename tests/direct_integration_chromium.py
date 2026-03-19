#!/usr/bin/env python3
"""
Single-file radiolysis model (flattened from your MoltenSaltRadiolysis repo).

- Keeps your config shape (kernel/volumes/kLa/dose_rate/radiation/G overrides/metals/gas species)
- Offers build_system(), integrate_system(), quick_plot(), and a small CLI:
      python msr_radiolysis_onefile.py examples_config.json

Notes:
- Units:
    * Concentrations: mol/L
    * Gas species state variable: mol in headspace
    * k (bimolecular): L/(mol*s); first-order: s^-1
    * Henry: kH in mol/(L*Pa), p in Pa
    * Dose rate: J/(L*s)
    * G-values: molecules / 100 eV
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import json, math, sys, pathlib
import numpy as np
import pandas as pd

# Optional deps
try:
    from scipy.integrate import solve_ivp
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

try:
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False

# ---- constants
NA = 6.02214076e23          # 1/mol
EV_J = 1.602176634e-19      # J
Rgas = 8.314462618          # J/(mol*K)

# ---- data structures
@dataclass
class Species:
    name: str
    phase: str  # "liq" or "gas"
    index: int

@dataclass
class Reaction:
    name: str
    reactants: Dict[str, float]         # species -> stoich coeff (consumed)
    products: Dict[str, float]          # species -> stoich coeff (formed)
    reversible: bool = False            # reverse not parameterized here; add explicit reverse reaction instead
    k_ref: Optional[float] = None
    T_ref: Optional[float] = None
    Ea_J_mol: Optional[float] = None
    A: Optional[float] = None
    order_overrides: Optional[Dict[str, float]] = None
    notes: str = ""

    def k(self, T: float) -> float:
        """Arrhenius: either A/Ea, or k_ref/T_ref/Ea, or just k_ref."""
        if self.A is not None and self.Ea_J_mol is not None:
            return float(self.A) * math.exp(-self.Ea_J_mol/(Rgas*T))
        if self.k_ref is not None and self.T_ref is not None and self.Ea_J_mol is not None:
            # k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref))
            return float(self.k_ref) * math.exp(-self.Ea_J_mol/Rgas * (1.0/T - 1.0/self.T_ref))
        return float(self.k_ref) if self.k_ref is not None else 0.0

@dataclass
class GasExchange:
    """Linear mass transfer with Henry’s law:
       flux [mol/L/s] = kLa * (C_liq - kH * p_gas)
       gas variable is mol in headspace; p_gas = n_gas * R * T / V_gas
    """
    name: str
    liq_species: str
    gas_species: str
    kLa_s: float
    kH_molm3Pa: float

@dataclass
class System:
    species: List[Species]
    species_index: Dict[str, int]
    reactions: List[Reaction]
    sources: Dict[str, float]                 # species -> mol/L/s source into LIQUID pool
    T: float                                  # K
    V_liq_m3: float
    V_gas_m3: float
    gas_exchanges: List[GasExchange] = field(default_factory=list)
    initial_concentrations: Optional[np.ndarray] = None

# ---- sources from G-values
def build_sources_from_G(g_values: Dict[str, float], dose_rate_J_m3_s: float) -> Dict[str, float]:
    """S_i [mol/L/s] = G_i [molecules/100eV] * dose_rate / (100 eV in J) / N_A"""
    if dose_rate_J_m3_s == 0.0 or not g_values:
        return {}
    factor = dose_rate_J_m3_s / (100.0 * EV_J) / NA
    return {sp: float(G) * factor for sp, G in g_values.items()}

# ---- a tiny embedded DB you can extend/replace with your real numbers
# Structure mirrors what your builder expected: kernels, phases, reactions, henry, metals, G_values...
MINI_DB = {
    "kernels": {
        "chloride": {
            "species": ["Cl-", "Cl2•-", "Cl2_diss", "e_s-", "Cl•", "Cl3-"],   # add/adjust as needed
            "phases":  {"Cl-": "liq", "Cl2•-": "liq", "Cl2_diss": "liq", "e_s-": "liq", "Cl•": "liq", "Cl3-": "liq"},
            "reactions": [
                {"name": "e + Cl2_diss -> Cl2•-",
                 "reactants": {"e_s-": 1.0, "Cl2_diss": 1.0},
                 "products":  {"Cl2•-": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 1.0e10}},   # L/(mol*s)
                {"name": "e + Cl2•- -> 2Cl-",
                 "reactants": {"e_s-": 1.0, "Cl2•-": 1.0},
                 "products":  {"Cl-": 2.0},
                 "reversible": False,
                 "params": {"k_ref": 1.0e10}},   # L/(mol*s)
                {"name": "Cl- + Cl• -> Cl2•-",
                 "reactants": {"Cl-": 1.0, "Cl•": 1.0},
                 "products":  {"Cl2•-": 2.0},
                 "reversible": False,
                 "params": {"k_ref": 1.0e10}},   # L/(mol*s)
                 {"name": "Cl2•- + Cl2•- -> Cl- + Cl3-",
                 "reactants": {"Cl2•-": 2.0},
                 "products":  {"Cl-": 1.0, "Cl3-": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 2.2e9}},   # L/(mol*s)
            ],
        },
        "fluoride": {
            "species": ["F-", "F2•-", "F2_diss", "e_s-"],
            "phases":  {"F-": "liq", "F2•-": "liq", "F2_diss": "liq", "e_s-": "liq"},
            "reactions": [
                {"name": "e + F2_diss -> F2•-",
                 "reactants": {"e_s-": 1.0, "F2_diss": 1.0},
                 "products":  {"F2•-": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 1.0e-3}},
            ],
        },
        
    },
    "henry": {  # placeholder Henry constants
        "Cl2_mol_m3_Pa": 2.0e-5,   # mol/(L*Pa)
        "F2_mol_m3_Pa":  1.0e-5,
    },
    # Metal section shows structure only; add your real species/reactions
    "metals": {
        "Zn": {
            "species": ["Zn2+", "Zn+"],
            "chloride_templated_reactions": [],
            "fluoride_templated_reactions": [],
        },
        "U":  {
            "species": ["U4+", "U3+"],
            "chloride_templated_reactions": [],
            "fluoride_templated_reactions": [],
        },
        "chromium": {
            "species": ["Cr2+", "Cr3+", "Cr+", "e_s-"],   # add/adjust as needed
            "phases":  {"Cr2+": "liq", "Cr+": "liq", "Cr+": "liq", "e_s-": "liq"},
            "chloride_templated_reactions": [
                 {"name": "e + Cr3+ -> Cr2+",
                 "reactants": {"e_s-": 1.0, "Cr3+": 1.0},
                 "products":  {"Cr2+": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 6.1e10}},   # L/(mol*s)
                 {"name": "e + Cr2+ -> Cr+",
                 "reactants": {"e_s-": 1.0, "Cr2+": 1.0},
                 "products":  {"Cr+": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 4.1e10}},   # L/(mol*s)
                 {"name": "Cr2+ + Cl2•- -> Cr3+ + 2Cl-",
                 "reactants": {"Cr2+": 1.0, "Cl2•-": 1.0},
                 "products":  {"Cr3+": 1.0, "Cl-": 2.0},
                 "reversible": False,
                 "params": {"k_ref": 7.2e9}},   # L/(mol*s)
                 {"name": "Cr3+ + Cl2•- -> Cr2+ + Cl2_diss",
                 "reactants": {"Cr3+": 1.0, "Cl2•-": 1.0},
                 "products":  {"Cr2+": 1.0, "Cl2_diss": 1.0},
                 "reversible": False,
                 "params": {"k_ref": 1.4e9}},   # L/(mol*s)
                 {"name": "Cr3+ + Cr+ -> 2Cr2+",
                 "reactants": {"Cr3+": 1.0, "Cr+": 1.0},
                 "products":  {"Cr2+": 2.0},
                 "reversible": False,
                 "params": {"k_ref": 1.7e10}},   # L/(mol*s)
            ],
        },
    },
    # Example layout for G-values (put your real ones here or pass via config override)
    "G_values": {
        "gamma": {
            "chloride": {"Cl2_diss": 0.8, "e_s-": 2.7, "Cl2•-": 1.7},   # molecules/100eV
            "fluoride": {"F2_diss": 0.001, "e_s-": 0.001},
        }
    }
}

# ---- builder (flattened)
def _species_list_to_objects(spec_names: List[str], phases: Dict[str, str]) -> List[Species]:
    return [Species(name=s, phase=phases.get(s, "liq"), index=i) for i, s in enumerate(spec_names)]

def build_system(config: Dict, db: Optional[Dict] = None) -> System:
    """
    Build a System from a configuration dictionary.
    Required/recognized keys (matching your README/builder):
      kernel: "chloride" | "fluoride"
      temperature_K, liquid_volume_m3, headspace_volume_m3, kLa_s^-1
      dose_rate_J_m3_s, radiation
      G_values_override: {species -> G (molecules/100eV)}
      initial_concentrations: {species -> mol/L}
      metals: {MetalName -> {species -> mol/L}}
      gas_species: list subset of ["Cl2","F2"] to activate gas exchange
    """
    db = db or MINI_DB
    kernel = config["kernel"]; assert kernel in db["kernels"]
    T_K  = float(config.get("temperature_K", 673.15))
    V_liq = float(config.get("liquid_volume_m3", 1.0e-3))
    V_gas = float(config.get("headspace_volume_m3", 0.0))
    kLa   = float(config.get("kLa_s^-1", 0.0))
    dose_rate = float(config.get("dose_rate_J_m3_s", 0.0))
    radiation = str(config.get("radiation", "gamma"))
    g_override = dict(config.get("G_values_override", {}) or {})
    gas_flags = set(config.get("gas_species", []) or [])

    kdb = db["kernels"][kernel]
    species_list = list(kdb["species"])
    phases = dict(kdb["phases"])
    reactions: List[Reaction] = []

    # base reactions
    for row in kdb.get("reactions", []):
        reactions.append(_reaction_from_row(row))

    # gas phases & exchanges
    gas_exchanges: List[GasExchange] = []
    if "Cl2" in gas_flags and kernel == "chloride":
        if "Cl2_g" not in species_list:
            species_list.append("Cl2_g")
            phases["Cl2_g"] = "gas"
        gas_exchanges.append(GasExchange(
            name="Cl2", liq_species="Cl2_diss", gas_species="Cl2_g",
            kLa_s=kLa, kH_molm3Pa=db["henry"]["Cl2_mol_m3_Pa"]
        ))
    if "F2" in gas_flags and kernel == "fluoride":
        if "F2_g" not in species_list:
            species_list.append("F2_g")
            phases["F2_g"] = "gas"
        gas_exchanges.append(GasExchange(
            name="F2", liq_species="F2_diss", gas_species="F2_g",
            kLa_s=kLa, kH_molm3Pa=db["henry"]["F2_mol_m3_Pa"]
        ))

    # metals (species + template reactions if present)
    metals_cfg = dict(config.get("metals", {}) or {})
    for metal, conc_map in metals_cfg.items():
        if metal not in db["metals"]:
            continue
        md = db["metals"][metal]
        for sp in md.get("species", []):
            if sp not in species_list:
                species_list.append(sp)
                phases[sp] = "liq"
        tmpl_key = f"{kernel}_templated_reactions"
        for row in md.get(tmpl_key, []):
            reactions.append(_reaction_from_row(row))

    # finalize species + index
    species_objs = _species_list_to_objects(species_list, phases)
    idx = {s.name: s.index for s in species_objs}

    # sources from G-values (kernel+radiation, overridden by config)
    g_kernel = db.get("G_values", {}).get(radiation, {}).get(kernel, {}) or {}
    g_all = {**g_kernel, **g_override}
    sources = build_sources_from_G(g_all, dose_rate)

    # system
    sysobj = System(
        species=species_objs, species_index=idx, reactions=reactions, sources=sources,
        T=T_K, V_liq_m3=V_liq, V_gas_m3=V_gas, gas_exchanges=gas_exchanges
    )
    sysobj.initial_concentrations = _assemble_initial(sysobj, config)
    return sysobj

def _reaction_from_row(row: Dict) -> Reaction:
    def _num(v):
        if isinstance(v, str):
            try: return float(v)
            except Exception: return v
        return v
    pars = {k: _num(v) for k, v in row.get("params", {}).items()}
    return Reaction(
        name=row.get("name", "rxn"),
        reactants=dict(row["reactants"]),
        products=dict(row["products"]),
        reversible=bool(row.get("reversible", False)),
        k_ref=pars.get("k_ref"), T_ref=pars.get("T_ref"),
        Ea_J_mol=pars.get("Ea_J_mol"), A=pars.get("A"),
        order_overrides=row.get("orders"), notes=row.get("notes", "")
    )

def _assemble_initial(system: System, config: Dict) -> np.ndarray:
    n = len(system.species)
    C0 = np.zeros(n, dtype=float)
    init = dict(config.get("initial_concentrations", {}) or {})
    for s, v in init.items():
        if s in system.species_index:
            C0[system.species_index[s]] = float(v)
    # metals init
    for _, conc_map in (config.get("metals", {}) or {}).items():
        for s, v in conc_map.items():
            if s in system.species_index:
                C0[system.species_index[s]] = float(v)
    # gas species start at 0 mol (headspace)
    return C0

# ---- RHS assembly
def _rate_law_C(system: System, sp_name: str, y: np.ndarray) -> float:
    """Return concentration for rate law (liquid phases only)."""
    s = system.species[system.species_index[sp_name]]
    if s.phase == "gas":
        return 0.0  # do not involve gas moles in homogeneous kinetics
    return max(0.0, y[s.index])

def _rhs(system: System, t: float, y: np.ndarray) -> np.ndarray:
    dy = np.zeros_like(y)

    # homogeneous reactions (liquid)
    for rxn in system.reactions:
        k = rxn.k(system.T)
        if k == 0.0:
            continue
        # mass action term
        rate = k
        if rxn.reactants:
            for sp, sto in rxn.reactants.items():
                order = rxn.order_overrides.get(sp, sto) if rxn.order_overrides else sto
                Ci = _rate_law_C(system, sp, y)
                if order == 0:
                    continue
                rate *= Ci ** float(order)
        if rate == 0.0 or not math.isfinite(rate):
            continue
        # apply stoichiometry to liquids
        for sp, nu in rxn.reactants.items():
            i = system.species_index.get(sp)
            if i is None: continue
            if system.species[i].phase == "liq":
                dy[i] -= rate * nu
        for sp, nu in rxn.products.items():
            i = system.species_index.get(sp)
            if i is None: continue
            if system.species[i].phase == "liq":
                dy[i] += rate * nu

    # radiolytic sources (liquid only)
    for sp, S in system.sources.items():
        i = system.species_index.get(sp)
        if i is not None and system.species[i].phase == "liq":
            dy[i] += float(S)

    # gas exchange
    if system.gas_exchanges and system.V_gas_m3 > 0.0:
        for gx in system.gas_exchanges:
            il = system.species_index[gx.liq_species]
            ig = system.species_index[gx.gas_species]
            C_liq = max(0.0, y[il])
            n_gas = max(0.0, y[ig])
            p_gas = n_gas * Rgas * system.T / system.V_gas_m3  # Pa
            C_eq = gx.kH_molm3Pa * p_gas
            flux = gx.kLa_s * (C_liq - C_eq)  # mol/L/s (positive = to gas)
            # liquid conc drops by flux
            dy[il] -= flux
            # gas moles increase by flux * liquid volume
            dy[ig] += flux * system.V_liq_m3
    return dy

# ---- integration
def integrate_system(system: System, t_final: float, n_steps: int = 10000,
                     method: str = "auto") -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Return (t, Y, extras) with Y shape [n_times, n_species]."""
    t0 = 0.0
    t = np.linspace(t0, float(t_final), int(n_steps))
    y0 = np.array(system.initial_concentrations, dtype=float)

    if _HAVE_SCIPY and method == "auto":
        sol = solve_ivp(lambda tt, yy: _rhs(system, tt, yy), (t0, t[-1]), y0,
                        t_eval=t, method="Radau", vectorized=False, atol=1e-12, rtol=1e-8)
        Y = sol.y.T
        info = {"scipy": True, "status": int(sol.status), "nfev": getattr(sol, "nfev", None)}
    else:
        # Simple semi-implicit fallback: backward Euler with 1 fixed-point sweep
        Y = np.zeros((t.size, y0.size))
        Y[0] = y0
        for k in range(1, t.size):
            dt = t[k] - t[k-1]
            y_prev = Y[k-1].copy()
            # predictor: explicit Euler
            y = y_prev + dt * _rhs(system, t[k-1], y_prev)
            # one corrector sweep
            y = y_prev + dt * _rhs(system, t[k], y)
            # keep non-negative
            Y[k] = np.maximum(0.0, y)
        info = {"scipy": False, "note": "backward-Euler (1 sweep) fallback"}

    extras = {"species": [s.name for s in system.species],
              "phases": [s.phase for s in system.species]}
    return t, Y, extras


# ---- plotting
def quick_plot(t: np.ndarray, C: np.ndarray, system: System, species_to_plot: Optional[List[str]] = None, absorbance_csv: Optional[str] = None):
    if not _HAVE_MPL:
        print("matplotlib not available; skipping plot.")
        return

    sp_idx = system.species_index
    names = [s.name for s in system.species]
    targets = species_to_plot or [s for s in names if system.species[sp_idx[s]].phase == "liq"]

    fig, ax1 = plt.subplots()

    # Plot concentrations on the left y-axis
    for s in targets:
        i = sp_idx[s]
        y = C[:, i]
        label = f"{s} ({'mol/L' if system.species[i].phase=='liq' else 'mol'})"
        ax1.semilogy(t, y, label=label)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Concentration [mol/L]")
    ax1.legend(loc="upper left")

    # Add absorbance on the right y-axis from CSV
    if absorbance_csv is not None:
        # Load data from CSV
        try:
            absorbance_data = pd.read_csv(absorbance_csv)

            # Convert time from milliseconds to seconds
            absorbance_data['time'] = absorbance_data['time'] * 1e-9

            # Extract time and absorbance columns
            absorbance_time = absorbance_data['time'].to_numpy()  # Time in seconds
            absorbance_values = absorbance_data['absorbance'].to_numpy()  # Absorbance values

            # Plot absorbance data on the right y-axis
            ax2 = ax1.twinx()  # Create a secondary y-axis
            ax2.semilogy(absorbance_time, absorbance_values, color='red', label="Absorbance Data for 3.03 mM Cr2+")
            ax2.set_ylabel("Absorbance")
            ax2.legend(loc="upper right")
        except Exception as e:
            print(f"Error loading or processing absorbance data from CSV: {e}")

    plt.tight_layout()
    plt.savefig('youplot_with_absorbance.png')
    plt.show()


# ---- CLI
def _load_config(path: str) -> Dict:
    with open(path, "r") as f:
        return json.load(f)

def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: python msr_radiolysis_onefile.py <config.json> [t_final] [n_steps]")
        # print("No config provided; running a tiny demo...")
        demo = {
            "kernel": "chloride",
            "metals": {
                    "chromium": {  # Add chromium species with initial concentrations
                        "Cr2+": 3.03e-3,
                        "Cr3+": 0,
                        "Cr+": 0
                        }
                    },
            "temperature_K": 673.15,
            "liquid_volume_m3": 1.0e-3,
            "headspace_volume_m3": 1.0e-3,
            "kLa_s^-1": 0.01,
            "dose_rate_J_m3_s": 1.0e6,
            "radiation": "gamma",
            "G_values_override": {"Cl2_diss": 0.002},  # molecules/100eV
            "initial_concentrations": {
                "Cl-": 20.0, #Very large concentration, essentially constant
                "Cl2•-": 0,
                "Cl2_diss": 0,
                "e_s-": 1e-3, #Small initial concentration modeling pulse
            },
            "gas_species": ["Cl2"]
        }
        system = build_system(demo)
        t, C, extra = integrate_system(system, t_final=20e-9, n_steps=10000)
        print("Species:", extra["species"])
        print("Final state:", C[-1])
        quick_plot(t, C, system, species_to_plot=["e_s-"], absorbance_csv="absorbance3mMCr2.csv")
        return

    cfg = _load_config(argv[0])
    t_final = float(argv[1]) if len(argv) > 1 else 20e-9
    n_steps = int(argv[2]) if len(argv) > 2 else 10000
    system = build_system(cfg)
    t, C, extra = integrate_system(system, t_final=t_final, n_steps=n_steps)

    # write CSV
    out = pathlib.Path("radiolysis_output.csv")
    header = ",".join(["t_s"] + extra["species"])
    data = np.column_stack([t, C])
    np.savetxt(out, data, delimiter=",", header=header, comments="")
    print(f"Wrote {out.resolve()}")

    # plot
    quick_plot(t, C, system)

if __name__ == "__main__":
    main()
