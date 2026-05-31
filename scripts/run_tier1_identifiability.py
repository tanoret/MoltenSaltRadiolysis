#!/usr/bin/env python3
"""Tier 1 identifiability + consistency driver.

Runs profile-likelihood on the Cr+Zn pseudo-1st-order rate model using the
Iwamatsu 2026 and Iwamatsu 2022 digitized data, and the cross-paper Arrhenius
consistency check on the Zn²⁺ + e_s⁻ reaction across Pikaev 1982 / Iwamatsu 2022.

Usage:
    python scripts/run_tier1_identifiability.py
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from msr_radiolysis.validation.identifiability import (
    ParameterSpec,
    ProfileResult,
    RateObservation,
    consistency_check,
    identifiability_report,
    profile_likelihood,
    pseudo_first_order_nll,
)


def load_pseudo_rate_csv(path: Path, T_col: str = "T_C", c_col: str = None,
                          k_col: str = None) -> list:
    """Load pseudo-1st-order rate data into (T_K, [M] mol/L, k_obs, sigma_log_k) tuples.

    Expects a CSV with header containing T_C, the metal-concentration column, and
    a k_obs_1e8_per_s column.
    """
    tuples = []
    with path.open() as f:
        # skip leading '#'-comment lines
        cleaned = (line for line in f if not line.lstrip().startswith("#"))
        rows = list(csv.DictReader(cleaned))
    # auto-detect concentration column if not given
    if c_col is None:
        for k in rows[0].keys():
            if k.endswith("_mM"):
                c_col = k
                break
    if k_col is None:
        for k in rows[0].keys():
            if k.startswith("k_obs"):
                k_col = k
                break
    for r in rows:
        T_K = float(r[T_col]) + 273.15
        Mconc = float(r[c_col]) * 1e-3   # mM -> M (mol/L)
        # k_obs given as values × 10^8 s^-1 in our vision CSV
        k_obs = float(r[k_col]) * 1e8
        # assume 10% noise in the absorbance trace -> sigma_log_k ~ 0.1
        sigma_log_k = 0.10
        tuples.append((T_K, Mconc, k_obs, sigma_log_k))
    return tuples


def run_profile_for_reaction(name: str, data_tuples, log_A_init: float, Ea_init: float):
    """Profile-likelihood for a 2-parameter Arrhenius from pseudo-1st-order data.

    The σ_log_k assumed in load_pseudo_rate_csv is conservative (0.10); we recalibrate to
    the empirical MLE-residual scale, which is the methodologically correct σ when the
    quoted measurement uncertainty does not match the actual fit residuals (Bates &
    Watts, Nonlinear Regression Analysis, 1988, §2.2). This avoids over- or under-stating
    parameter precision purely from a mis-specified noise model.
    """
    nll_raw = lambda theta: pseudo_first_order_nll(theta, data_tuples)
    res = minimize(nll_raw, x0=np.array([log_A_init, Ea_init]), method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 5000})
    mle = res.x

    # Empirical σ from MLE residuals (Bates-Watts §2.2.4)
    # NLL_raw = 0.5 * Σ ((log k_obs - log k_pred) / sigma_assumed)²
    # If sigma_actual differs from sigma_assumed, the *true* residual-squared sum is
    # 2 * NLL_raw * sigma_assumed². Empirical σ² = (2 NLL_raw σ_assumed²) / (n - p).
    n = len(data_tuples)
    dof = n - 2
    sigma_assumed = 0.10
    sigma_emp = float(np.sqrt(2.0 * res.fun * sigma_assumed**2 / dof)) if dof > 0 else sigma_assumed
    # Use empirical σ for the profile (rescale NLL accordingly)
    scale = (sigma_assumed / sigma_emp)**2
    nll = lambda theta: scale * nll_raw(theta)

    print(f"  {name}: MLE log_A = {mle[0]:.4f} (A = {np.exp(mle[0]):.3e}), Ea = {mle[1]/1e3:.2f} kJ/mol")
    print(f"        empirical sigma_log_k = {sigma_emp:.4f} (vs assumed {sigma_assumed:.2f}); NLL_calibrated = {scale * res.fun:.4f}")

    params = [
        ParameterSpec(name=f"log_A_{name}", initial=log_A_init, transform="identity"),
        ParameterSpec(name=f"Ea_{name}", initial=Ea_init, lower=0.0, upper=2.0e5, transform="identity"),
    ]
    profiles = []
    for i, p in enumerate(params):
        grid_span = 2.0 if i == 0 else 1.5e4
        pr = profile_likelihood(nll, params, mle, target_idx=i, grid_span=grid_span, n_grid=41)
        profiles.append(pr)
    return profiles


def main():
    print("=" * 72)
    print("Tier 1 identifiability + consistency analysis")
    print("=" * 72)

    # ---------------------------------------------------------------------
    # 1. Profile likelihoods for Cr2+, Cr3+, Zn2+ vs e_s-
    # ---------------------------------------------------------------------
    val_root = REPO / "validation"
    print("\n--- Profile likelihood analyses ---")
    all_profiles = []

    cr2_csv = val_root / "cr_licl_kcl/iwamatsu_2026_pccp/data/vision_fig2B_Cr2_pseudo1st_order.csv"
    cr2_data = load_pseudo_rate_csv(cr2_csv, c_col="Cr2_mM")
    print(f"\nCr2+: {len(cr2_data)} (T, [M], k_obs) tuples")
    profiles_cr2 = run_profile_for_reaction("Cr2", cr2_data, log_A_init=np.log(1.7e13), Ea_init=33.5e3)
    all_profiles.extend(profiles_cr2)

    cr3_csv = val_root / "cr_licl_kcl/iwamatsu_2026_pccp/data/vision_fig3B_Cr3_pseudo1st_order.csv"
    cr3_data = load_pseudo_rate_csv(cr3_csv, c_col="Cr3_mM")
    print(f"\nCr3+: {len(cr3_data)} (T, [M], k_obs) tuples")
    profiles_cr3 = run_profile_for_reaction("Cr3", cr3_data, log_A_init=np.log(2.0e13), Ea_init=31.8e3)
    all_profiles.extend(profiles_cr3)

    zn2_csv = val_root / "zn_licl_kcl/horne_2022_pccp/data/vision_fig4B_Zn2_pseudo1st_order.csv"
    zn2_data = load_pseudo_rate_csv(zn2_csv, c_col="Zn2_mM")
    print(f"\nZn2+: {len(zn2_data)} (T, [M], k_obs) tuples")
    profiles_zn2 = run_profile_for_reaction("Zn2", zn2_data, log_A_init=np.log(2.4e13), Ea_init=35.6e3)
    all_profiles.extend(profiles_zn2)

    # ---------------------------------------------------------------------
    # 2. Cross-paper consistency check: Zn2+ across Pikaev 1982, Iwamatsu 2022
    # ---------------------------------------------------------------------
    print("\n--- Cross-paper Arrhenius consistency ---\n")
    all_consistency = []

    # For e_s- + Zn2+: combine Pikaev 1982 (NaCl, KCl) + Iwamatsu 2022 (LiCl-KCl extrapolated at 5 T)
    zn_obs = []
    # Pikaev 1982 - NaCl at 850 C: k = 1.7e9 M^-1 s^-1
    # use a conservative log-σ of 0.5 (Iwamatsu 2022 describes these as "best-guess estimates")
    zn_obs.append(RateObservation(paper="Pikaev1982", reaction_id="e_s-+Zn2+",
                                    T_K=1123.15, log_k=np.log(1.7e9), sigma_log_k=0.5,
                                    salt="NaCl"))
    zn_obs.append(RateObservation(paper="Pikaev1982", reaction_id="e_s-+Zn2+",
                                    T_K=1073.15, log_k=np.log(2.8e9), sigma_log_k=0.5,
                                    salt="KCl"))
    # Iwamatsu 2022 - LiCl-KCl at 5 T from the published density-corrected Arrhenius
    # A = 2.4e13, Ea = 35.6 kJ/mol; sigma_log_k from sigma_A and sigma_Ea propagation
    A, sA, Ea, sEa = 2.4e13, 0.5e13, 35.6e3, 1.2e3
    for T_C in [400, 450, 500, 550, 600]:
        T_K = T_C + 273.15
        log_k = np.log(A) - Ea / (8.314 * T_K)
        # propagate (rel A)^2 + (Ea/RT)^2 (rel Ea / Ea)^2
        sigma_log_k = np.sqrt((sA/A)**2 + (sEa/(8.314*T_K))**2)
        zn_obs.append(RateObservation(paper="Iwamatsu2022", reaction_id="e_s-+Zn2+",
                                        T_K=T_K, log_k=log_k, sigma_log_k=sigma_log_k,
                                        salt="LiCl-KCl"))
    cr_zn = consistency_check(zn_obs, reaction_id="e_s- + Zn2+ -> Zn+")
    print(f"  {cr_zn.note}")
    print(f"  Best A: {cr_zn.A_best:.3e} M^-1 s^-1")
    print(f"  Best Ea: {cr_zn.Ea_best/1e3:.2f} kJ/mol")
    if not cr_zn.feasible:
        print(f"  Residuals (>2σ flagged):")
        for k, v in cr_zn.residuals.items():
            print(f"    {k}: {v:+.2f}σ" + ("  *** OUTLIER" if abs(v) > 2 else ""))
    all_consistency.append(cr_zn)

    # For e_s- + Cr2+: only Iwamatsu 2026, single paper -> can't test cross-paper, skip
    # For Cl2•- + Cl2•-: Hagiwara 1987 vs Iwamatsu 2022 (both report Ea, but only Iwamatsu reports k)
    # Skip due to lack of multi-paper k data

    # ---------------------------------------------------------------------
    # 3. Write report
    # ---------------------------------------------------------------------
    output_path = REPO / "validation" / "TIER1_REPORT.md"
    identifiability_report(all_profiles, all_consistency, output_path)
    print(f"\n--- Report written to {output_path.relative_to(REPO)} ---")

    return all_profiles, all_consistency


if __name__ == "__main__":
    main()
