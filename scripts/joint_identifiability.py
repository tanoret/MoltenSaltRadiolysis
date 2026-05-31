#!/usr/bin/env python3
"""Joint identifiability analysis across all reactions in the Cr+Zn chloride kernel
that have any rate-constant data.

This extends scripts/run_tier1_identifiability.py to the **joint** problem: all
reactions inferred simultaneously from the union of available rate observations.
The goal is to identify which parameters are practically identifiable from data
alone and which require informative priors (literature Arrhenius from
arrhenius_parameters.csv) under HBMAE.

The 8 reactions analysed:
  1. e_s- + Cr2+ -> Cr+               | Iwamatsu 2026 Fig. 2B  (20 points × T, [M])
  2. e_s- + Cr3+ -> Cr2+              | Iwamatsu 2026 Fig. 3B  (17 points × T, [M])
  3. e_s- + Zn2+ -> Zn+               | Iwamatsu 2022 Fig. 4B  (20 points × T, [M])
  4. Cl2•- + Cl2•- -> Cl3- + Cl-      | Iwamatsu 2022 Eq.4:  k(400)=2.2e9, Ea=26 kJ/mol
  5. Cl2•- + Cr2+ -> Cr3+ + 2 Cl-     | Iwamatsu 2026 Eq.9:  k(400)=7.2e9, NO Ea
  6. Cl2•- + Cr3+ -> products         | Iwamatsu 2026 Eq.10: k(400)=1.4e9, NO Ea
  7. Cl2•- + Zn+ -> 2 Cl- + Zn2+      | Iwamatsu 2022 Eq.15a: k(400)~2e10, Ea=25 kJ/mol
  8. e_s- impurity decay              | k_bg = 1.0e7 s^-1 at 400°C (background)

For each reaction we report:
  - Number of data points and their constraint type (transient vs scalar vs Arrhenius)
  - MLE (log A, Ea)
  - 95% profile-likelihood confidence intervals
  - Identifiability classification:
       FULL    = both log A and Ea finite CI
       RIDGE   = log A and Ea coupled but joint k(T_avg) determined
       PARTIAL = only one of A or Ea has finite CI
       PRIOR   = neither identifiable without literature prior
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from msr_radiolysis.validation.identifiability import (
    ParameterSpec,
    pseudo_first_order_nll,
    profile_likelihood,
)

R_GAS = 8.314462618


@dataclass
class ReactionData:
    """All available data for one reaction."""
    name: str
    pseudo_data: list   # list of (T_K, [M] mol/L, k_obs s^-1, sigma_log_k)
    point_data: list    # list of (T_K, log_k, sigma_log_k) — scalar reports
    arrhenius_prior_A: float | None  # literature A (M^-1 s^-1)
    arrhenius_prior_Ea: float | None  # literature Ea (J/mol)
    arrhenius_sigma_A: float | None
    arrhenius_sigma_Ea: float | None
    note: str = ""


def joint_nll(theta, rdata: ReactionData):
    """Combined NLL: pseudo-1st-order data + scalar rate data."""
    log_A, Ea = theta
    total = 0.0
    for T, M, k_obs, sigma in rdata.pseudo_data:
        log_k_pred = log_A - Ea / (R_GAS * T) + np.log(M)
        total += 0.5 * ((np.log(k_obs) - log_k_pred) / sigma) ** 2
    for T, log_k_obs, sigma in rdata.point_data:
        log_k_pred = log_A - Ea / (R_GAS * T)
        total += 0.5 * ((log_k_obs - log_k_pred) / sigma) ** 2
    return total


def load_pseudo_data(path: Path, conc_col: str):
    tuples = []
    with path.open() as f:
        cleaned = (line for line in f if not line.lstrip().startswith("#"))
        rows = list(csv.DictReader(cleaned))
    for r in rows:
        T_K = float(r["T_C"]) + 273.15
        M = float(r[conc_col]) * 1e-3
        k = float(r["k_obs_1e8_per_s"]) * 1e8
        tuples.append((T_K, M, k, 0.10))   # sigma will be recalibrated to empirical
    return tuples


def classify_identifiability(profile_log_A, profile_Ea) -> str:
    """Classify identifiability based on profile-likelihood CIs."""
    log_A_finite = profile_log_A.lower_ci is not None and profile_log_A.upper_ci is not None
    Ea_finite = profile_Ea.lower_ci is not None and profile_Ea.upper_ci is not None
    if log_A_finite and Ea_finite:
        return "FULL"
    if not log_A_finite and not Ea_finite:
        return "RIDGE"     # both at grid boundary => ridge-degenerate
    return "PARTIAL"


def fit_one_reaction(rdata: ReactionData, log_A_init: float, Ea_init: float, verbose=True):
    """MLE + profile likelihood for one reaction's joint dataset."""
    n_total = len(rdata.pseudo_data) + len(rdata.point_data)
    if n_total < 2:
        return None  # can't profile with <2 obs

    nll_raw = lambda theta: joint_nll(theta, rdata)
    res = minimize(nll_raw, x0=np.array([log_A_init, Ea_init]),
                   method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 5000})
    mle = res.x

    # Empirical sigma calibration (Bates-Watts §2.2.4)
    dof = n_total - 2
    sigma_assumed = 0.10
    if dof > 0:
        sigma_emp = float(np.sqrt(2.0 * res.fun * sigma_assumed**2 / dof))
    else:
        sigma_emp = sigma_assumed
    scale = (sigma_assumed / max(sigma_emp, 1e-6)) ** 2
    nll = lambda theta: scale * nll_raw(theta)

    params = [
        ParameterSpec(name="log_A", initial=log_A_init),
        ParameterSpec(name="Ea", initial=Ea_init, lower=-2e4, upper=2e5),
    ]
    prof_A = profile_likelihood(nll, params, mle, target_idx=0, grid_span=3.0, n_grid=41)
    prof_E = profile_likelihood(nll, params, mle, target_idx=1, grid_span=2e4, n_grid=41)
    classification = classify_identifiability(prof_A, prof_E)

    if verbose:
        print(f"  {rdata.name}")
        print(f"    n_data    : pseudo={len(rdata.pseudo_data)}, point={len(rdata.point_data)}")
        print(f"    MLE       : log A = {mle[0]:.3f} (A = {np.exp(mle[0]):.3e}), Ea = {mle[1]/1e3:.2f} kJ/mol")
        print(f"    sigma_emp : {sigma_emp:.4f}")
        print(f"    log A CI  : [{prof_A.lower_ci}, {prof_A.upper_ci}]")
        print(f"    Ea CI     : [{prof_E.lower_ci}, {prof_E.upper_ci}]")
        print(f"    class     : {classification}")
        if rdata.arrhenius_prior_A is not None:
            ratio_A = np.exp(mle[0]) / rdata.arrhenius_prior_A
            offset_Ea = (mle[1] - rdata.arrhenius_prior_Ea) / 1e3
            print(f"    vs prior  : A_MLE/A_lit = {ratio_A:.3f}, Ea_MLE − Ea_lit = {offset_Ea:+.2f} kJ/mol")

    return mle, prof_A, prof_E, classification, sigma_emp


def main():
    print("=" * 76)
    print("Joint identifiability for the Cr+Zn chloride kernel")
    print("Profile likelihood + empirical sigma + identifiability classification")
    print("=" * 76)

    val_root = REPO / "validation"

    # --- Reaction inventory with all available rate data ---
    reactions = []

    # 1. e_s- + Cr2+ -> Cr+
    reactions.append(ReactionData(
        name="e_s- + Cr2+ -> Cr+",
        pseudo_data=load_pseudo_data(val_root / "cr_licl_kcl/iwamatsu_2026_pccp/data/vision_fig2B_Cr2_pseudo1st_order.csv", "Cr2_mM"),
        point_data=[],
        arrhenius_prior_A=1.7e13, arrhenius_prior_Ea=33.5e3,
        arrhenius_sigma_A=0.2e13, arrhenius_sigma_Ea=0.6e3,
    ))

    # 2. e_s- + Cr3+ -> Cr2+
    reactions.append(ReactionData(
        name="e_s- + Cr3+ -> Cr2+",
        pseudo_data=load_pseudo_data(val_root / "cr_licl_kcl/iwamatsu_2026_pccp/data/vision_fig3B_Cr3_pseudo1st_order.csv", "Cr3_mM"),
        point_data=[],
        arrhenius_prior_A=2.0e13, arrhenius_prior_Ea=31.8e3,
        arrhenius_sigma_A=0.5e13, arrhenius_sigma_Ea=0.5e3,
    ))

    # 3. e_s- + Zn2+ -> Zn+
    reactions.append(ReactionData(
        name="e_s- + Zn2+ -> Zn+",
        pseudo_data=load_pseudo_data(val_root / "zn_licl_kcl/horne_2022_pccp/data/vision_fig4B_Zn2_pseudo1st_order.csv", "Zn2_mM"),
        point_data=[],
        arrhenius_prior_A=2.4e13, arrhenius_prior_Ea=35.6e3,
        arrhenius_sigma_A=0.5e13, arrhenius_sigma_Ea=1.2e3,
    ))

    # 4. Cl2•- + Cl2•- -> Cl3- + Cl-  (single-T k + reported Ea)
    # k(400°C) = 2.2e9 ± 2.0e9 M^-1 s^-1; Ea = 26 ± 2 kJ/mol
    # Build effective scalar observations at 400 and 600 °C using Arrhenius extrapolation
    # so the joint fit can resolve A and Ea. Use the published Ea-σ as the data σ.
    T1 = 673.15
    log_k_400 = np.log(2.2e9)
    # Derive log_k at 600°C from k(400) and Ea=26 kJ/mol:
    log_k_600 = log_k_400 - 26e3 / R_GAS * (1/873.15 - 1/T1)
    reactions.append(ReactionData(
        name="Cl2•- + Cl2•- -> Cl3- + Cl-",
        pseudo_data=[],
        point_data=[(T1, log_k_400, 0.5), (873.15, log_k_600, 0.5)],   # σ=0.5 reflects Ea±2
        arrhenius_prior_A=2.2e9 * np.exp(26e3 / (R_GAS * T1)),
        arrhenius_prior_Ea=26e3,
        arrhenius_sigma_A=None, arrhenius_sigma_Ea=2e3,
        note="reconstructed from k(400) and Ea=26 ± 2 kJ/mol",
    ))

    # 5. Cl2•- + Cr2+ -> Cr3+ + 2 Cl-  (single-T point only, no Ea)
    reactions.append(ReactionData(
        name="Cl2•- + Cr2+ -> Cr3+ + 2 Cl-",
        pseudo_data=[],
        point_data=[(673.15, np.log(7.2e9), 0.05)],
        arrhenius_prior_A=None, arrhenius_prior_Ea=None,
        arrhenius_sigma_A=None, arrhenius_sigma_Ea=None,
        note="only k(400°C) reported; A and Ea both unidentifiable from data alone",
    ))

    # 6. Cl2•- + Cr3+ -> Cr2+ + Cl2_diss
    reactions.append(ReactionData(
        name="Cl2•- + Cr3+ -> Cr2+ + Cl2_diss",
        pseudo_data=[],
        point_data=[(673.15, np.log(1.4e9), 0.07)],
        arrhenius_prior_A=None, arrhenius_prior_Ea=None,
        arrhenius_sigma_A=None, arrhenius_sigma_Ea=None,
        note="only k(400°C) reported; mechanism ambiguous (Cr4+ vs Cr2+)",
    ))

    # 7. Cl2•- + Zn+ -> 2 Cl- + Zn2+  (Ea reported, k at 400°C order-of-magnitude)
    # k(400°C) ~ 2e10, Ea = 25 ± 3 kJ/mol; reconstruct point at 600°C
    log_k_400_eq15 = np.log(2.0e10)
    log_k_600_eq15 = log_k_400_eq15 - 25e3 / R_GAS * (1/873.15 - 1/T1)
    reactions.append(ReactionData(
        name="Cl2•- + Zn+ -> 2 Cl- + Zn2+",
        pseudo_data=[],
        point_data=[(T1, log_k_400_eq15, 0.5), (873.15, log_k_600_eq15, 0.5)],
        arrhenius_prior_A=2.0e10 * np.exp(25e3 / (R_GAS * T1)),
        arrhenius_prior_Ea=25e3,
        arrhenius_sigma_A=None, arrhenius_sigma_Ea=3e3,
    ))

    # --- Run all reactions ---
    print()
    results = []
    for r in reactions:
        if r.arrhenius_prior_A is not None:
            init = (np.log(r.arrhenius_prior_A), r.arrhenius_prior_Ea)
        else:
            init = (np.log(1e10), 30e3)
        try:
            out = fit_one_reaction(r, init[0], init[1])
        except Exception as e:
            print(f"  {r.name}: FAILED ({e})")
            continue
        if out is None:
            print(f"  {r.name}: skipped (< 2 data points; needs prior anchoring)")
            continue
        mle, prof_A, prof_E, classification, sigma_emp = out
        results.append((r, mle, prof_A, prof_E, classification, sigma_emp))

    # --- Summary table ---
    print()
    print("=" * 76)
    print("Joint identifiability summary")
    print("=" * 76)
    print(f"{'Reaction':<35s} {'n':>4s} {'class':>10s} {'A vs lit':>12s} {'Ea vs lit':>12s}")
    print("-" * 76)
    for r, mle, _, _, classification, _ in results:
        n_total = len(r.pseudo_data) + len(r.point_data)
        if r.arrhenius_prior_A is not None:
            ratio = f"{np.exp(mle[0]) / r.arrhenius_prior_A:>11.3f}"
            offset = f"{(mle[1] - r.arrhenius_prior_Ea)/1e3:>+11.2f}"
        else:
            ratio = "n/a"
            offset = "n/a"
        print(f"{r.name:<35s} {n_total:>4d} {classification:>10s} {ratio:>12s} {offset:>12s}")
    print()

    # --- Implication for HBMAE ---
    full = [c for *_, c, _ in results if c == "FULL"]
    ridge = [c for *_, c, _ in results if c == "RIDGE"]
    partial = [c for *_, c, _ in results if c == "PARTIAL"]
    print("Implications for HBMAE Tier 2 calibration:")
    print(f"  FULL    identifiable (no prior needed): {len(full)}")
    print(f"  RIDGE   degenerate (need joint k(T_ref) anchor): {len(ridge)}")
    print(f"  PARTIAL (need 1 prior dimension): {len(partial)}")
    n_no_data = sum(1 for r in reactions if (len(r.pseudo_data) + len(r.point_data)) < 2)
    print(f"  PRIOR-only (no data, full prior): {n_no_data}")
    print()
    print("Conclusion: hierarchical informative priors are necessary for at least")
    print(f"{len(ridge) + len(partial) + n_no_data}/{len(reactions)} reactions in the Cr+Zn chloride kernel.")

    return results


if __name__ == "__main__":
    main()
