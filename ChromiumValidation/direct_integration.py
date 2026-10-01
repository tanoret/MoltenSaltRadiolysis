#!/usr/bin/env python3
"""
Single-file radiolysis model (flattened from the MoltenSaltRadiolysis repo).

    python direct_integration.py config.json [t_final_s] [n_out]

Units
-----
* Concentrations: mol/L (M)
* Gas species state variable: mol in headspace
* k: first-order s^-1; second-order L mol^-1 s^-1
* G-values: molecules / 100 eV
* Pulse dose: Gy (J/kg); density: g/cm^3

Radiolytic source (revision 2)
------------------------------
Total radiolytic yield of species i per pulse:

    C_i,tot [M] = G_i * D * rho / (100 eV * N_A)

with D in Gy, rho in kg/L (= g/cm^3), 100 eV in J. This yield is either
  * placed in the initial state (pulse_width_s == 0, "instantaneous" pulse), or
  * delivered at the constant rate C_i,tot / tau_p over 0 <= t <= tau_p
    (rectangular pulse of width tau_p), and zero afterwards.
The rectangular pulse is integrated as two segments ([0, tau_p] and
[tau_p, t_final]) so the solver never steps across the source discontinuity.

Background (pseudo-first-order) losses (revision 2)
---------------------------------------------------
Experimental decays at zero solute concentration are not zero. Iwamatsu et al.
(2026) ESI Tables S2-S4 give observed pseudo-first-order rate coefficients vs
[Cr]; the intercepts of k_obs vs [Cr] give the background loss rate of each
transient in each sample series. These are supplied per case through the config
key "background_first_order_s^-1": {"e_s-": k_bg, "Cl2•-": k_bg}. e_s- background
is an open-system sink (unidentified impurity scavenging); Cl2•- background
reduces Cl2•- to 2 Cl- (as in the previous version's 3.3e6 s^-1 term).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import json, math, sys, pathlib
import numpy as np
import pandas as pd

try:
    import scipy
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
    reactants: Dict[str, float]
    products: Dict[str, float]
    reversible: bool = False
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
            return float(self.k_ref) * math.exp(-self.Ea_J_mol/Rgas * (1.0/T - 1.0/self.T_ref))
        return float(self.k_ref) if self.k_ref is not None else 0.0

@dataclass
class GasExchange:
    """Linear mass transfer with Henry's law:
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
    yields: Dict[str, float]                  # total radiolytic yield per pulse [M]
    pulse_width: float                        # s (0 = instantaneous)
    T: float                                  # K
    V_liq_m3: float
    V_gas_m3: float
    gas_exchanges: List[GasExchange] = field(default_factory=list)
    initial_concentrations: Optional[np.ndarray] = None   # pre-pulse state

# ---- radiolytic yields
def radiolytic_yields(g_values: Dict[str, float], dose_Gy: float,
                      density_g_cm3: float) -> Dict[str, float]:
    """C_i,tot [mol/L] = G_i [molecules/100 eV] * D [J/kg] * rho [kg/L] / (100 eV [J] * N_A)."""
    if dose_Gy == 0.0 or not g_values:
        return {}
    factor = dose_Gy * density_g_cm3 / (100.0 * EV_J * NA)
    return {sp: float(G) * factor for sp, G in g_values.items()}

# ---- reaction database
# Rate coefficients at 400 C. R1-R4, R6, R7 are order-of-magnitude estimates (see paper Table 1).
MINI_DB = {
    "kernels": {
        "chloride": {
            "species": ["Cl-", "Cl2•-", "Cl2_diss", "e_s-", "Cl•", "Cl3-"],
            "phases":  {"Cl-": "liq", "Cl2•-": "liq", "Cl2_diss": "liq", "e_s-": "liq", "Cl•": "liq", "Cl3-": "liq"},
            "reactions": [
                {"name": "R1 e_s- + Cl2 -> Cl2•-",
                 "reactants": {"e_s-": 1.0, "Cl2_diss": 1.0}, "products": {"Cl2•-": 1.0},
                 "params": {"k_ref": 1.0e10}},                      # M^-1 s^-1 (estimate)
                {"name": "R2 e_s- + Cl2•- -> 2Cl-",
                 "reactants": {"e_s-": 1.0, "Cl2•-": 1.0}, "products": {"Cl-": 2.0},
                 "params": {"k_ref": 1.0e10}},                      # M^-1 s^-1 (estimate)
                {"name": "R3 Cl- + Cl• -> Cl2•-",
                 "reactants": {"Cl-": 1.0, "Cl•": 1.0}, "products": {"Cl2•-": 1.0},
                 "params": {"k_ref": 1.0e10}},                      # M^-1 s^-1 (estimate)
                {"name": "R4 Cl• + e_s- -> Cl-",
                 "reactants": {"Cl•": 1.0, "e_s-": 1.0}, "products": {"Cl-": 1.0},
                 "params": {"k_ref": 1.0e10}},                      # M^-1 s^-1 (estimate)
                {"name": "R5 Cl2•- + Cl2•- -> Cl- + Cl3-",
                 "reactants": {"Cl2•-": 2.0}, "products": {"Cl-": 1.0, "Cl3-": 1.0},
                 "params": {"k_ref": 2.2e9}},                       # M^-1 s^-1
                {"name": "R6 Cl3- -> Cl- + Cl2",
                 "reactants": {"Cl3-": 1.0}, "products": {"Cl-": 1.0, "Cl2_diss": 1.0},
                 "params": {"k_ref": 1.0e4}},                       # s^-1 (estimate)
                {"name": "R7 Cl- + Cl2 -> Cl3-",
                 "reactants": {"Cl-": 1.0, "Cl2_diss": 1.0}, "products": {"Cl3-": 1.0},
                 "params": {"k_ref": 1.0e3}},                       # M^-1 s^-1 (estimate)
            ],
        },
        "fluoride": {
            "species": ["F-", "F2•-", "F2_diss", "e_s-"],
            "phases":  {"F-": "liq", "F2•-": "liq", "F2_diss": "liq", "e_s-": "liq"},
            "reactions": [
                {"name": "e + F2_diss -> F2•-",
                 "reactants": {"e_s-": 1.0, "F2_diss": 1.0}, "products": {"F2•-": 1.0},
                 "params": {"k_ref": 1.0e-3}},
            ],
        },
    },
    # Background first-order loss channels; k supplied per case via config.
    "background": {
        "e_s-":  {"name": "Rbg,e e_s- -> (impurity sink)", "products": {}},
        "Cl2•-": {"name": "Rbg,Cl2 Cl2•- -> 2Cl- (impurity)", "products": {"Cl-": 2.0}},
    },
    "henry": {
        "Cl2_mol_m3_Pa": 2.0e-5,
        "F2_mol_m3_Pa":  1.0e-5,
    },
    "metals": {
        "Zn": {"species": ["Zn2+", "Zn+"],
               "chloride_templated_reactions": [], "fluoride_templated_reactions": []},
        "U":  {"species": ["U4+", "U3+"],
               "chloride_templated_reactions": [], "fluoride_templated_reactions": []},
        "chromium": {
            "species": ["Cr2+", "Cr3+", "Cr+"],
            "chloride_templated_reactions": [
                {"name": "R8 Cr2+ + Cl2•- -> Cr3+ + 2Cl-",
                 "reactants": {"Cr2+": 1.0, "Cl2•-": 1.0}, "products": {"Cr3+": 1.0, "Cl-": 2.0},
                 "params": {"k_ref": 7.2e9}},
                {"name": "R9 Cr3+ + Cl2•- -> Cr2+ + Cl2",
                 "reactants": {"Cr3+": 1.0, "Cl2•-": 1.0}, "products": {"Cr2+": 1.0, "Cl2_diss": 1.0},
                 "params": {"k_ref": 1.4e9}},
                {"name": "R10 Cr3+ + Cr+ -> 2Cr2+",
                 "reactants": {"Cr3+": 1.0, "Cr+": 1.0}, "products": {"Cr2+": 2.0},
                 "params": {"k_ref": 1.7e10}},
                {"name": "R11 e_s- + Cr2+ -> Cr+",
                 "reactants": {"e_s-": 1.0, "Cr2+": 1.0}, "products": {"Cr+": 1.0},
                 "params": {"k_ref": 4.1e10}},
                {"name": "R12 e_s- + Cr3+ -> Cr2+",
                 "reactants": {"e_s-": 1.0, "Cr3+": 1.0}, "products": {"Cr2+": 1.0},
                 "params": {"k_ref": 6.1e10}},
            ],
        },
    },
    "G_values": {
        "gamma": {
            "chloride": {"e_s-": 2.8, "Cl2•-": 2.8},
            "fluoride": {"F2_diss": 0.001, "e_s-": 0.001},
        }
    }
}

# ---- builder
def _species_list_to_objects(spec_names: List[str], phases: Dict[str, str]) -> List[Species]:
    return [Species(name=s, phase=phases.get(s, "liq"), index=i) for i, s in enumerate(spec_names)]

def _reaction_from_row(row: Dict) -> Reaction:
    def _num(v):
        if isinstance(v, str):
            try: return float(v)
            except Exception: return v
        return v
    pars = {k: _num(v) for k, v in row.get("params", {}).items()}
    return Reaction(
        name=row.get("name", "rxn"),
        reactants=dict(row["reactants"]), products=dict(row["products"]),
        reversible=bool(row.get("reversible", False)),
        k_ref=pars.get("k_ref"), T_ref=pars.get("T_ref"),
        Ea_J_mol=pars.get("Ea_J_mol"), A=pars.get("A"),
        order_overrides=row.get("orders"), notes=row.get("notes", "")
    )

def build_system(config: Dict, db: Optional[Dict] = None) -> System:
    """
    Recognized config keys:
      kernel, temperature_K, liquid_volume_m3, headspace_volume_m3, kLa_s^-1, radiation,
      pulse_dose_Gy, density_g_cm3, pulse_width_s,
      G_values_override {species: G}, initial_concentrations {species: M},
      metals {metal: {species: M}}, gas_species [..],
      background_first_order_s^-1 {species: k},
      rate_overrides {reaction-tag (e.g. "R11"): k},  disabled_reactions [tags]
    """
    db = db or MINI_DB
    kernel = config["kernel"]; assert kernel in db["kernels"]
    T_K   = float(config.get("temperature_K", 673.15))
    V_liq = float(config.get("liquid_volume_m3", 1.0e-3))
    V_gas = float(config.get("headspace_volume_m3", 0.0))
    kLa   = float(config.get("kLa_s^-1", 0.0))
    radiation = str(config.get("radiation", "gamma"))
    g_override = dict(config.get("G_values_override", {}) or {})
    gas_flags = set(config.get("gas_species", []) or [])
    dose_Gy = float(config.get("pulse_dose_Gy", 0.0))
    rho = float(config.get("density_g_cm3", 1.677))
    tau_p = float(config.get("pulse_width_s", 0.0))
    overrides = dict(config.get("rate_overrides", {}) or {})
    disabled = set(config.get("disabled_reactions", []) or [])

    kdb = db["kernels"][kernel]
    species_list = list(kdb["species"])
    phases = dict(kdb["phases"])
    rows = list(kdb.get("reactions", []))

    gas_exchanges: List[GasExchange] = []
    if "Cl2" in gas_flags and kernel == "chloride":
        if "Cl2_g" not in species_list:
            species_list.append("Cl2_g"); phases["Cl2_g"] = "gas"
        gas_exchanges.append(GasExchange("Cl2", "Cl2_diss", "Cl2_g", kLa, db["henry"]["Cl2_mol_m3_Pa"]))
    if "F2" in gas_flags and kernel == "fluoride":
        if "F2_g" not in species_list:
            species_list.append("F2_g"); phases["F2_g"] = "gas"
        gas_exchanges.append(GasExchange("F2", "F2_diss", "F2_g", kLa, db["henry"]["F2_mol_m3_Pa"]))

    for metal in (config.get("metals", {}) or {}):
        if metal not in db["metals"]:
            continue
        md = db["metals"][metal]
        for sp in md.get("species", []):
            if sp not in species_list:
                species_list.append(sp); phases[sp] = "liq"
        rows += list(md.get(f"{kernel}_templated_reactions", []))

    reactions: List[Reaction] = []
    for row in rows:
        rx = _reaction_from_row(row)
        tag = rx.name.split()[0]
        if tag in disabled:
            continue
        if tag in overrides:
            rx.k_ref = float(overrides[tag]); rx.A = None; rx.Ea_J_mol = None
        reactions.append(rx)

    # background first-order losses
    for sp, kbg in (config.get("background_first_order_s^-1", {}) or {}).items():
        if float(kbg) <= 0.0 or sp not in species_list:
            continue
        bg = db["background"][sp]
        reactions.append(Reaction(name=bg["name"], reactants={sp: 1.0},
                                  products=dict(bg["products"]), k_ref=float(kbg)))

    species_objs = _species_list_to_objects(species_list, phases)
    idx = {s.name: s.index for s in species_objs}

    g_kernel = db.get("G_values", {}).get(radiation, {}).get(kernel, {}) or {}
    g_all = {**g_kernel, **g_override}
    yields = radiolytic_yields(g_all, dose_Gy, rho)

    sysobj = System(species=species_objs, species_index=idx, reactions=reactions,
                    yields=yields, pulse_width=tau_p, T=T_K, V_liq_m3=V_liq, V_gas_m3=V_gas,
                    gas_exchanges=gas_exchanges)
    sysobj.initial_concentrations = _assemble_initial(sysobj, config)
    return sysobj

def _assemble_initial(system: System, config: Dict) -> np.ndarray:
    C0 = np.zeros(len(system.species), dtype=float)
    for s, v in (config.get("initial_concentrations", {}) or {}).items():
        if s in system.species_index:
            C0[system.species_index[s]] = float(v)
    for _, conc_map in (config.get("metals", {}) or {}).items():
        for s, v in conc_map.items():
            if s in system.species_index:
                C0[system.species_index[s]] = float(v)
    return C0

# ---- RHS
def _source_vector(system: System) -> np.ndarray:
    """Constant source rate [M/s] during a rectangular pulse."""
    S = np.zeros(len(system.species))
    if system.pulse_width > 0.0:
        for sp, c in system.yields.items():
            i = system.species_index.get(sp)
            if i is not None and system.species[i].phase == "liq":
                S[i] = c / system.pulse_width
    return S

def _rhs(system: System, t: float, y: np.ndarray, S: Optional[np.ndarray] = None) -> np.ndarray:
    dy = np.zeros_like(y)
    T = system.T
    for rxn in system.reactions:
        k = rxn.k(T)
        if k == 0.0:
            continue
        rate = k
        for sp, sto in rxn.reactants.items():
            order = rxn.order_overrides.get(sp, sto) if rxn.order_overrides else sto
            s = system.species[system.species_index[sp]]
            Ci = 0.0 if s.phase == "gas" else max(0.0, y[s.index])
            if order != 0:
                rate *= Ci ** float(order)
        if rate == 0.0 or not math.isfinite(rate):
            continue
        for sp, nu in rxn.reactants.items():
            i = system.species_index[sp]
            if system.species[i].phase == "liq":
                dy[i] -= rate * nu
        for sp, nu in rxn.products.items():
            i = system.species_index[sp]
            if system.species[i].phase == "liq":
                dy[i] += rate * nu

    if S is not None:
        dy += S

    if system.gas_exchanges and system.V_gas_m3 > 0.0:
        for gx in system.gas_exchanges:
            il = system.species_index[gx.liq_species]
            ig = system.species_index[gx.gas_species]
            C_liq = max(0.0, y[il]); n_gas = max(0.0, y[ig])
            p_gas = n_gas * Rgas * system.T / system.V_gas_m3
            flux = gx.kLa_s * (C_liq - gx.kH_molm3Pa * p_gas)
            dy[il] -= flux
            dy[ig] += flux * system.V_liq_m3
    return dy

# ---- integration
SOLVER_DEFAULTS = dict(method="Radau", rtol=1e-8, atol=1e-14)

def integrate_system(system: System, t_final: float, n_out: int = 100000,
                     t_eval: Optional[np.ndarray] = None,
                     rtol: float = SOLVER_DEFAULTS["rtol"],
                     atol: float = SOLVER_DEFAULTS["atol"],
                     method: str = SOLVER_DEFAULTS["method"]) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Integrate from t = 0 (pulse start) to t_final.
    n_out is the number of OUTPUT points (t_eval); the solver's internal steps are adaptive.
    Returns (t, Y [n_times x n_species], info).
    """
    if not _HAVE_SCIPY:
        raise RuntimeError("SciPy is required.")
    t = np.linspace(0.0, float(t_final), int(n_out)) if t_eval is None else np.asarray(t_eval, float)
    y0 = np.array(system.initial_concentrations, dtype=float)
    tau = system.pulse_width
    segs = []
    if tau <= 0.0:
        for sp, c in system.yields.items():           # instantaneous pulse
            i = system.species_index.get(sp)
            if i is not None:
                y0[i] += c
        segs.append((0.0, t[-1], None, t))
    else:
        S = _source_vector(system)
        segs.append((0.0, tau, S, t[t <= tau]))
        segs.append((tau, t[-1], None, t[t > tau]))

    Ts, Ys, info = [], [], {"python": sys.version.split()[0], "scipy": scipy.__version__,
                            "method": method, "rtol": rtol, "atol": atol,
                            "nfev": 0, "njev": 0, "nlu": 0, "success": True, "segments": []}
    y = y0
    for (ta, tb, S, te) in segs:
        if tb <= ta:
            continue
        # always evaluate at tb as well, so the next segment starts from the exact end state
        te_full = np.append(te[te < tb], tb)
        sol = solve_ivp(lambda tt, yy: _rhs(system, tt, yy, S), (ta, tb), y,
                        t_eval=te_full, method=method, rtol=rtol, atol=atol)
        info["success"] &= bool(sol.success)
        for k in ("nfev", "njev", "nlu"):
            info[k] += int(getattr(sol, k, 0) or 0)
        info["segments"].append({"t0": ta, "t1": tb, "status": int(sol.status), "message": sol.message})
        keep = len(te_full) if (len(te) and te[-1] >= tb) else len(te_full) - 1
        Ts.append(sol.t[:keep]); Ys.append(sol.y.T[:keep])
        y = sol.y[:, -1]
    tt = np.concatenate(Ts); Y = np.vstack(Ys)
    info["species"] = [s.name for s in system.species]
    info["phases"] = [s.phase for s in system.species]
    info["min_value"] = float(Y.min())
    return tt, Y, info

# ---- CSV export
def export_species_csv(t: np.ndarray, C: np.ndarray, system: System,
                       species_names: List[str], path: str = "species_concentration.csv") -> None:
    """Export liquid-phase species concentrations vs time (t_s, <species>...)."""
    names, idxs = [], []
    for n in species_names:
        if n not in system.species_index:
            print(f"Warning: '{n}' not found. Available: {list(system.species_index)}"); continue
        i = system.species_index[n]
        if system.species[i].phase != "liq":
            print(f"Warning: '{n}' is gas-phase (mol, not mol/L) — skipping."); continue
        names.append(n); idxs.append(i)
    if not names:
        print("Error: no valid liquid-phase species to export."); return
    safe = lambda s: s.replace("•", "rad").replace(" ", "_")
    header = "t_s," + ",".join(safe(n) for n in names)
    np.savetxt(pathlib.Path(path), np.column_stack([t] + [C[:, i] for i in idxs]),
               delimiter=",", header=header, comments="")
    print(f"Species concentrations written to {pathlib.Path(path).resolve()} ({len(t)} points)")

# ---- plotting
def quick_plot(t: np.ndarray, C: np.ndarray, system: System,
               species_to_plot: Optional[List[str]] = None, path: str = "species_plot.png"):
    """Diagnostic concentration plot (time in ns). Not a model-data comparison."""
    if not _HAVE_MPL:
        print("matplotlib not available; skipping plot."); return
    targets = species_to_plot or [s.name for s in system.species if s.phase == "liq"]
    fig, ax = plt.subplots()
    for s in targets:
        ax.plot(t * 1e9, C[:, system.species_index[s]], label=s)
    ax.set_xlabel("Time after pulse start (ns)"); ax.set_ylabel("Concentration (M)")
    ax.legend(); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)

# ---- CLI
def _load_config(path: str) -> Dict:
    with open(path, "r") as f:
        return json.load(f)

def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: python direct_integration.py <config.json> [t_final_s] [n_out]")
        sys.exit(0)
    cfg = _load_config(argv[0])
    t_final = float(argv[1]) if len(argv) > 1 else 20e-9
    n_out = int(argv[2]) if len(argv) > 2 else 100000
    system = build_system(cfg)
    t, C, info = integrate_system(system, t_final=t_final, n_out=n_out)
    print({k: v for k, v in info.items() if k not in ("species", "phases")})
    export_species_csv(t, C, system, species_names=["Cl2•-", "e_s-"], path="species_concentration.csv")
    header = ",".join(["t_s"] + info["species"])
    np.savetxt("radiolysis_output.csv", np.column_stack([t, C]), delimiter=",", header=header, comments="")
    quick_plot(t, C, system, species_to_plot=["e_s-", "Cl2•-"])

if __name__ == "__main__":
    main()
