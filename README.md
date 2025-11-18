# MoltenSaltRadiolysis

A compact, extensible framework for **radiolysis modeling in molten salts** (MSRs).

**Features**
- Data-driven species & reactions (YAML database).
- Kinetic network builder for **chloride** and **fluoride** salts with optional metals (Zn, U as examples).
- Radiolysis **G-value** sources per species and dose rate.
- Gas–liquid coupling for **Cl₂** and **F₂** with Henry law + `k_L a` mass transfer.
- Stiff ODE integration via **SciPy** (if available), with a pure-NumPy implicit-Euler fallback.
- Minimal **front end** (CLI) + plotting utilities.
- **Validation tests** that check unit conversions, rate-law behavior, and simple literature-aligned targets.

> This is a *starter* framework with documented placeholders where the community still lacks definitive numbers.
> Swap in your facility data and extend the YAML database to your chemistry.

## Quick start (inside Python)

```python
from msr_radiolysis.builder import build_system
from msr_radiolysis.integrator import integrate_system
from msr_radiolysis.plotting import quick_plot

# Choose a base salt kernel and metals to include
config = {
  "kernel": "chloride",          # or "fluoride"
  "temperature_K": 673.15,       # e.g., 400 °C
  "liquid_volume_m3": 1.0e-3,    # 1 liter control volume
  "headspace_volume_m3": 1.0e-3, # 1 liter headspace (set 0 for none)
  "kLa_s^-1": 0.0,               # mass-transfer coefficient (0 disables)
  "dose_rate_J_m3_s": 1.0e6,     # volumetric energy deposition
  "radiation": "gamma",
  "G_values_override": {         # optionally override or add G-values here (molecules/100eV)
     "F2_diss": 0.005
  },
  "metals": {
     "Zn": {"Zn2+": 1.0, "Zn+": 0.0},  # mol/m^3 initial
     "U":  {"U4+": 0.1, "U3+": 0.0}
  },
  "initial_concentrations": {     # mol/m^3 for tracked species (others default to 0)
     "Cl-": 5000.0
  },
  "gas_species": ["Cl2","F2"]     # which diatomics have gas phases
}

system = build_system(config)
t, C, extra = integrate_system(system, t_final=5.0, n_steps=200)  # seconds
quick_plot(t, C, system, species_to_plot=["e_s-", "Cl2•-", "Cl2_diss"])
```

## CLI

```bash
python -m msr_radiolysis.frontend_cli examples/config_licl_kcl_zinc.json
```

This writes CSVs and PNGs to the working directory.

## Validation

```bash
python -m msr_radiolysis.validation.tests
```

The tests exercise unit conversions (G-values to sources), kinetic half-lives, second-order decays, and gas–liquid partitioning.

---

### Notes on units
- Concentrations: **mol/m³**.
- Bimolecular rate constants: **m³/(mol·s)**.
- First-order constants: **s⁻¹**.
- Henry law (implemented form): **C = kH · p**, with **kH in mol/(m³·Pa)** and pressure in **Pa**.
- Dose rate: **J/(m³·s)**.
- G-values: **molecules per 100 eV**.

A reliable conversion factor is included for `S_i = G_i · ẏ_vol · 1/(100 eV · N_A)`.

---

### Extend
- Add species and reactions in `msr_radiolysis/data/database.yaml`.
- You can specify a reaction via:
  - `k_ref`, `T_ref`, `Ea_J_mol`  (Arrhenius from a known k at T_ref), or
  - `A` and `Ea_J_mol` directly.
- Reactions are assumed **mass-action**. Orders default to reactant stoichiometry unless overridden.
