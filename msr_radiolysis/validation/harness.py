"""Literature-comparison validation harness.

Loads a validation case from a YAML manifest, runs the model with the
configured system, and compares the resulting trajectories to digitized
experimental traces.

Manifest layout (see validation/cr_licl_kcl/iwamatsu_2026_pccp/manifest.yaml
for a complete example):

    case_id: str
    title: str
    citation: {authors, journal, year, volume, issue, pages, doi}
    digitization: {tool, figures, date_extracted, notes}
    system:                       # passed to build_system()
      kernel: chloride|fluoride
      temperature_K: float
      liquid_volume_m3: float
      headspace_volume_m3: float
      kLa_s^-1: float
      dose_rate_J_m3_s: float
      radiation: gamma
      G_values_override: {...}
      metals: {...}
      initial_concentrations: {...}
      gas_species: [...]
    integration:
      t_final_s: float
      n_steps: int
    traces:
      - file: data/foo.csv                # path relative to manifest
        label: "1 mM Cr(II)"
        observable:
          kind: concentration|absorbance|absorbance_scale_free
          species: e_s-                   # model species the trace corresponds to
          molar_absorptivity_M_cm: float  # optional, only for kind=absorbance
          path_length_cm: float           # optional, only for kind=absorbance
        time_units: s|ms|us|ns|ps
        initial_overrides:                # optional per-trace tweaks
          Cr2+: 1.0
    metrics:
      - kind: rmse_log|rmse|max_rel_error
        species: e_s-

For pulse-radiolysis cases set dose_rate_J_m3_s=0 and put the post-pulse
inventory in initial_concentrations (e.g., e_s- starts at its pulse value).
"""

from __future__ import annotations

import csv
import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml

from ..builder import build_system
from ..integrator import integrate_system

TIME_UNIT_TO_S: Dict[str, float] = {
    "s": 1.0,
    "ms": 1e-3,
    "us": 1e-6,
    "µs": 1e-6,
    "ns": 1e-9,
    "ps": 1e-12,
}


@dataclass
class Trace:
    file: str
    label: str
    observable: Dict[str, Any] = field(default_factory=dict)
    time_units: str = "s"
    initial_overrides: Dict[str, float] = field(default_factory=dict)


@dataclass
class ValidationCase:
    case_id: str
    title: str
    citation: Dict[str, Any]
    digitization: Dict[str, Any]
    system_config: Dict[str, Any]
    integration: Dict[str, Any]
    traces: List[Trace]
    metrics: List[Dict[str, Any]] = field(default_factory=list)
    manifest_path: Optional[Path] = None
    notes: str = ""

    @property
    def root(self) -> Path:
        return self.manifest_path.parent if self.manifest_path else Path(".")


def load_manifest(path) -> ValidationCase:
    p = Path(path)
    with p.open("r") as f:
        m = yaml.safe_load(f)
    traces = [
        Trace(
            file=t["file"],
            label=t.get("label", t["file"]),
            observable=t.get("observable", {}) or {},
            time_units=t.get("time_units", "s"),
            initial_overrides=t.get("initial_overrides", {}) or {},
        )
        for t in m.get("traces", [])
    ]
    return ValidationCase(
        case_id=m["case_id"],
        title=m.get("title", m["case_id"]),
        citation=m.get("citation", {}) or {},
        digitization=m.get("digitization", {}) or {},
        system_config=m["system"],
        integration=m.get("integration", {"t_final_s": 1.0e-6, "n_steps": 1000}),
        traces=traces,
        metrics=m.get("metrics", []) or [],
        manifest_path=p,
        notes=m.get("notes", ""),
    )


def load_trace_csv(path: Path, time_units: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read a 2-column CSV ('time' first column, value second). Convert time to seconds."""
    t_vals, y_vals = [], []
    with Path(path).open("r") as f:
        reader = csv.reader(f)
        rows = list(reader)
    start = 1 if rows and not _looks_numeric(rows[0][0]) else 0
    for row in rows[start:]:
        if len(row) < 2:
            continue
        try:
            t_vals.append(float(row[0]))
            y_vals.append(float(row[1]))
        except ValueError:
            continue
    if time_units not in TIME_UNIT_TO_S:
        raise ValueError(f"Unknown time_units '{time_units}'. Supported: {list(TIME_UNIT_TO_S)}")
    return np.asarray(t_vals) * TIME_UNIT_TO_S[time_units], np.asarray(y_vals)


def _looks_numeric(s: str) -> bool:
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def _apply_overrides(cfg: Dict[str, Any], overrides: Dict[str, float]) -> Dict[str, Any]:
    """Return a deep-copy config with per-trace initial concentration overrides applied.

    Overrides land in initial_concentrations or in metals.<metal>.<species> depending
    on where the key already lives in the config. New keys go to initial_concentrations.
    """
    new_cfg = copy.deepcopy(cfg)
    ic = dict(new_cfg.get("initial_concentrations", {}) or {})
    metals = dict(new_cfg.get("metals", {}) or {})
    for sp, val in overrides.items():
        placed = False
        for metal, conc_map in metals.items():
            if sp in conc_map:
                metals[metal] = {**conc_map, sp: float(val)}
                placed = True
                break
        if not placed:
            ic[sp] = float(val)
    new_cfg["initial_concentrations"] = ic
    new_cfg["metals"] = metals
    return new_cfg


def run_trace(case: ValidationCase, tr: Trace) -> Dict[str, Any]:
    """Build & integrate the system for one trace; return model + experimental arrays."""
    cfg = _apply_overrides(case.system_config, tr.initial_overrides)
    system = build_system(cfg)

    t_final = float(case.integration.get("t_final_s", 1.0e-6))
    n_steps = int(case.integration.get("n_steps", 1000))
    t_model, C_model, _ = integrate_system(system, t_final=t_final, n_steps=n_steps)

    species = tr.observable.get("species")
    if species and species in system.species_index:
        idx = system.species_index[species]
        model_obs = C_model[:, idx]
    else:
        model_obs = None

    csv_path = case.root / tr.file
    t_exp, y_exp = load_trace_csv(csv_path, tr.time_units)

    return {
        "trace": tr,
        "t_model_s": t_model,
        "C_model": C_model,
        "model_observable": model_obs,
        "species_index": system.species_index,
        "t_exp_s": t_exp,
        "y_exp": y_exp,
    }


def model_observable_to_data_space(result: Dict[str, Any]) -> Optional[np.ndarray]:
    """Convert the model's species trace into the observable's data space.

    - kind=concentration: returns mol/m^3 (same as model).
    - kind=absorbance with molar_absorptivity_M_cm and path_length_cm: returns
      Abs = epsilon (M^-1 cm^-1) * (C / 1000) (M) * path (cm).
      The factor of 1000 converts mol/m^3 -> mol/L.
    - kind=absorbance_scale_free: rescales model trace so its max matches the
      experimental max (used when epsilon is unknown).
    """
    tr: Trace = result["trace"]
    model = result["model_observable"]
    if model is None:
        return None
    kind = tr.observable.get("kind", "concentration")
    if kind == "concentration":
        return model
    if kind == "absorbance":
        eps = tr.observable.get("molar_absorptivity_M_cm")
        ell = tr.observable.get("path_length_cm")
        if eps is None or ell is None:
            return None
        return float(eps) * (model / 1000.0) * float(ell)
    if kind == "absorbance_scale_free":
        y_exp = result["y_exp"]
        if model.max() <= 0 or len(y_exp) == 0:
            return None
        return model * (float(np.nanmax(y_exp)) / float(model.max()))
    return None


def compute_metrics(case: ValidationCase, results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute manifest-specified metrics per trace."""
    out: Dict[str, Any] = {"per_trace": [], "metrics": {}}
    for r in results:
        tr: Trace = r["trace"]
        model_in_data_space = model_observable_to_data_space(r)
        if model_in_data_space is None:
            out["per_trace"].append({"label": tr.label, "rmse": None, "note": "no model observable"})
            continue
        t_model = r["t_model_s"]
        t_exp = r["t_exp_s"]
        y_exp = r["y_exp"]
        # interpolate model onto experimental time grid
        mask = (t_exp >= t_model[0]) & (t_exp <= t_model[-1])
        if mask.sum() < 2:
            out["per_trace"].append({"label": tr.label, "rmse": None, "note": "no overlap"})
            continue
        y_model_at_exp = np.interp(t_exp[mask], t_model, model_in_data_space)
        residual = y_model_at_exp - y_exp[mask]
        rmse = float(np.sqrt(np.mean(residual ** 2)))
        out["per_trace"].append({"label": tr.label, "rmse": rmse, "n_points": int(mask.sum())})
    return out


def run_case(manifest_path) -> Dict[str, Any]:
    """End-to-end: load manifest, run every trace, compute metrics."""
    case = load_manifest(manifest_path)
    results = [run_trace(case, tr) for tr in case.traces]
    metrics = compute_metrics(case, results)
    return {"case": case, "results": results, "metrics": metrics}


def discover_cases(root) -> List[Path]:
    """Walk a directory tree and return all manifest.yaml paths."""
    return sorted(Path(root).rglob("manifest.yaml"))
