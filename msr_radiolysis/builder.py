
import json
from typing import Dict, List, Optional
import copy
import numpy as np

from .core import Species, Reaction, GasExchange, System
from .sources import build_sources
from .utils import M_to_SI_k

def _load_database():
    import yaml, pkgutil
    data = pkgutil.get_data(__name__.split('.')[0], "data/database.yaml")
    return yaml.safe_load(data)

def _load_reference_validation():
    import yaml, pkgutil
    data = pkgutil.get_data(__name__.split('.')[0], "data/reference_validation.yaml")
    return yaml.safe_load(data)

def _species_objects(spec_names: List[str], phases: Dict[str,str]) -> List[Species]:
    species = []
    for i, s in enumerate(spec_names):
        ph = phases.get(s, "liq")
        species.append(Species(name=s, phase=ph, index=i))
    return species

def build_system(config: Dict) -> System:
    """
    Build a System from a configuration dict.
    Required keys:
     - kernel: "chloride" or "fluoride"
     - temperature_K
     - liquid_volume_m3
     - headspace_volume_m3
     - kLa_s^-1
     - dose_rate_J_m3_s
     - radiation: key for g-values (e.g., "gamma")
     - G_values_override: optional dict species->G
     - metals: optional dict metal -> {species: conc}
     - initial_concentrations: optional dict species->conc (mol/m^3) for any tracked species
     - gas_species: list including any of ["Cl2", "F2"] that should have gas phases active
    """
    db = _load_database()
    kernel = config["kernel"]
    assert kernel in ["chloride", "fluoride"]

    temperature_K = float(config.get("temperature_K", 873.15))
    V_liq = float(config.get("liquid_volume_m3", 1.0e-3))
    V_gas = float(config.get("headspace_volume_m3", 0.0))
    kLa = float(config.get("kLa_s^-1", 0.0))
    dose_rate = float(config.get("dose_rate_J_m3_s", 0.0))
    radiation = config.get("radiation","gamma")
    g_override = config.get("G_values_override", {}) or {}
    gas_species_flag = set(config.get("gas_species", []))

    # pick base species & reactions
    kernel_db = db["kernels"][kernel]
    species_list = copy.deepcopy(kernel_db["species"])
    phases = copy.deepcopy(kernel_db["phases"])

    # Add metals (if given) by including their species in the system and enabling templated reactions
    metals_cfg = config.get("metals", {}) or {}
    metals_db = db["metals"]
    reactions: List[Reaction] = []

    # Base reactions (non-templated)
    for rxn_row in kernel_db["reactions"]:
        reactions.append(_reaction_from_row(rxn_row))

    # Gas phases & exchange
    gas_exchanges = []
    if "Cl2" in gas_species_flag and kernel == "chloride":
        # ensure gas species exists
        if "Cl2_g" not in species_list:
            species_list.append("Cl2_g")
            phases["Cl2_g"] = "gas"
        gx = GasExchange(name="Cl2", liq_species="Cl2_diss", gas_species="Cl2_g",
                         kLa_s=kLa, kH_molm3Pa=db["henry"]["Cl2_mol_m3_Pa"])
        gas_exchanges.append(gx)
    if "F2" in gas_species_flag and kernel == "fluoride":
        if "F2_g" not in species_list:
            species_list.append("F2_g")
            phases["F2_g"] = "gas"
        gx = GasExchange(name="F2", liq_species="F2_diss", gas_species="F2_g",
                         kLa_s=kLa, kH_molm3Pa=db["henry"]["F2_mol_m3_Pa"])
        gas_exchanges.append(gx)

    # Add metal species & templated reactions
    for metal, conc_map in metals_cfg.items():
        if metal not in metals_db:
            continue
        md = metals_db[metal]
        for sp in md["species"]:
            if sp not in species_list:
                species_list.append(sp)
                phases[sp] = "liq"
        # templated reactions dependent on kernel
        tmpl_key = f"{kernel}_templated_reactions"
        for rxn_row in md.get(tmpl_key, []):
            reactions.append(_reaction_from_row(rxn_row))

    # Create Species objects & index
    species_objs = _species_objects(species_list, phases)
    index = {s.name: s.index for s in species_objs}

    # Prepare sources (G-values) as per kernel + override
    kernel_g = db["G_values"].get(radiation, {}).get(kernel, {})
    g_values = dict(kernel_g)
    g_values.update(g_override)  # override/add
    sources = build_sources(g_values, dose_rate)

    # System object
    system = System(species=species_objs,
                    species_index=index,
                    reactions=reactions,
                    sources=sources,
                    T=temperature_K,
                    V_liq_m3=V_liq,
                    V_gas_m3=V_gas,
                    gas_exchanges=gas_exchanges)
    # Attach initial conditions
    system.initial_concentrations = _assemble_initial_concs(system, config)
    return system

def _assemble_initial_concs(system: System, config: Dict) -> (np.ndarray):
    C0 = np.zeros(len(system.species))
    init = config.get("initial_concentrations", {}) or {}
    # initialize specified
    for s, val in init.items():
        if s in system.species_index:
            C0[system.species_index[s]] = float(val)
    # if metals provided, set their inits
    metals_cfg = config.get("metals", {}) or {}
    for metal, conc_map in metals_cfg.items():
        for sp, val in conc_map.items():
            if sp in system.species_index:
                C0[system.species_index[sp]] = float(val)
    # gas species default 0 mol (stored as mol)
    return C0

def _reaction_from_row(row: Dict) -> Reaction:
    def _num(v):
        if isinstance(v, str):
            try:
                return float(v)
            except Exception:
                return v
        return v
    # supports: reactants, products, name, reversible, params: k_ref, T_ref, Ea_J_mol, or A+Ea, and optional order_overrides
    pars = row.get("params", {})
    rr = Reaction(
        name=row.get("name","rxn"),
        reactants=row["reactants"],
        products=row["products"],
        reversible=row.get("reversible", False),
        k_ref=_num(pars.get("k_ref", None)),
        T_ref=_num(pars.get("T_ref", None)),
        Ea_J_mol=_num(pars.get("Ea_J_mol", None)),
        A=_num(pars.get("A", None)),
        order_overrides=row.get("orders", None),
        notes=row.get("notes","")
    )
    return rr
