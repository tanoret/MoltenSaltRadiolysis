# Literature validation cases

Each subdirectory under this folder is one validation case: a digitized
literature dataset plus a `manifest.yaml` that ties the data to a runnable
model configuration. Cases are grouped first by chemical system (e.g.
`cr_licl_kcl/`, `zn_licl_kcl/`, `flibe_msre/`) and then by paper
(`<first_author>_<year>_<journal>/`).

## Layout

```
validation/
├── README.md
├── CANDIDATES.md                       # candidate papers awaiting digitization
├── <system>/
│   └── <author>_<year>_<journal>/
│       ├── manifest.yaml               # model config + traces + metrics
│       └── data/
│           ├── *.csv                   # digitized traces (time, value)
│           └── Source.txt              # citation + digitization notes
```

## Manifest schema

See [msr_radiolysis/validation/harness.py](../msr_radiolysis/validation/harness.py)
for the loader. The `system:` block is passed straight to `build_system()`.

Pulse-radiolysis cases use `dose_rate_J_m3_s: 0` and put the post-pulse
inventory of primary radicals (typically `e_s-`, `Cl•`) in
`initial_concentrations`. Steady-state irradiation cases set a nonzero
`dose_rate_J_m3_s` and let the G-value sources drive the chemistry.

## Running

```python
from msr_radiolysis.validation.harness import run_case
out = run_case("validation/cr_licl_kcl/iwamatsu_2026_pccp/manifest.yaml")
```

`out["results"]` holds the model and experimental arrays for each trace,
ready to overlay-plot. `out["metrics"]` holds residual statistics.

## Conventions

- All concentrations in `manifest.yaml` are **mol/m³** (SI). Literature
  values in mol/L should be multiplied by `1e3`.
- All rate constants in [database.yaml](../msr_radiolysis/data/database.yaml)
  are **m³/(mol·s)**. Literature M⁻¹ s⁻¹ values multiplied by `1e-3`.
- Time in digitized CSV files keeps the original axis units. The manifest
  declares those units (`s`, `ms`, `µs`, `ns`, `ps`) and the harness converts.
- Cite every figure: `digitization.figures` must reference the figure
  number(s) the CSVs came from, and `data/Source.txt` must hold the full
  citation. A reader should be able to go from any CSV back to the
  original figure in the paper.
