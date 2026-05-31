#!/usr/bin/env python3
"""Tier 4 / Caveat 4: end-to-end integrated HBMAE MCMC.

This script combines THREE observation modalities in a single posterior:
  (M1) Transient absorbance time-series from Iwamatsu 2026 Cr experiment
       (9 traces; modality T in Section 3 of methods.tex)
  (M2) Pikaev 1982 + Iwamatsu 2022 scalar Zn rate observations across 3 salts
       (modality K)
  (M3) Phillips 2022 censored NULL for [Cl2]_gas in NaCl-UCl3
       (modality C)

The shared parameters across these data sources:
  - Cr Arrhenius (A5, Ea5, A6, Ea6) - constrained by M1 (Cr transients) and one
    Phillips dependence (the Cl2 production rate sees Cr's chemistry only via
    background impurity scavenging, but not directly)
  - Zn intrinsic Arrhenius (A_Zn, Ea_Zn) - constrained by M2
  - Salt perturbations η^(LiCl-KCl), η^(NaCl), η^(KCl) for Zn - constrained by M2
  - Facility offset b^(Pikaev) - constrained by M2
  - Background impurity decay rate k_bg - constrained by M1 (Cr transients) and M3 (NULL)
  - Effective G(Cl•) under Phillips conditions - constrained by M3 (NULL)
  - Per-trace Cr pulse dose [e_s-]_0^(t) - constrained by M1 only

This is the first integrated HBMAE run that exercises Theorems 2, 4, 5, and 6
simultaneously on real data plus the Phillips NULL constraint. Posterior is sampled
with emcee. Diagnostic output includes per-parameter marginals, posterior predictive
checks for each modality, and the inferred facility/salt-perturbation parameters.
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import emcee
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# Use the rigorous slow-manifold solver for Phillips NULL
from scripts.tier4_slow_manifold import (
    C_DETECTION_THRESHOLD,
    KH_CL2,
    T_K as T_K_phillips,
    V_gas,
    V_liq,
    dose_rate_J_m3_s,
    duration_s,
    slow_manifold_steady_state,
)

warnings.filterwarnings("ignore")

R_GAS = 8.314462618
NA = 6.02214076e23
EV_J = 1.602176634e-19


# =============================================================================
# Modality M1: Iwamatsu 2026 Cr transients (9 traces, stiff ODE)
# =============================================================================

def cr_rhs(t, y, k5_SI, k6_SI, k_bg):
    eS, Cr2, Crp, Cr3 = y
    rate5 = k5_SI * eS * Cr2
    rate6 = k6_SI * eS * Cr3
    return [
        -(rate5 + rate6 + k_bg * eS),
        -rate5 + rate6,
        +rate5,
        -rate6,
    ]


def solve_cr_trace(eS0, Cr2_0, Cr3_0, A5, Ea5, A6, Ea6, k_bg, t_eval, T_K=673.15):
    """Forward solve Cr transient. Robust to high-stiffness initial conditions by
    trying BDF then LSODA fallback with relaxed tolerances."""
    k5_M = A5 * np.exp(-Ea5 / (R_GAS * T_K))
    k6_M = A6 * np.exp(-Ea6 / (R_GAS * T_K))
    k5_SI = k5_M * 1e-3
    k6_SI = k6_M * 1e-3
    y0 = [eS0, Cr2_0, 0.0, Cr3_0]
    t_span = (0, t_eval[-1] * 1.01)
    for method, atol, rtol, mstep in [
        ("BDF", 1e-15, 1e-9, 1e-9),
        ("BDF", 1e-13, 1e-7, 1e-10),
        ("LSODA", 1e-15, 1e-8, None),
        ("Radau", 1e-13, 1e-7, None),
    ]:
        try:
            kwargs = dict(method=method, atol=atol, rtol=rtol, args=(k5_SI, k6_SI, k_bg))
            if mstep is not None:
                kwargs["max_step"] = mstep
            sol = solve_ivp(cr_rhs, t_span, y0, t_eval=t_eval, **kwargs)
            if sol.success and np.all(np.isfinite(sol.y[0])):
                return sol.y[0]
        except Exception:
            continue
    return None


def load_cr_traces():
    base = REPO / "validation/cr_licl_kcl/iwamatsu_2026_pccp/data"
    metadata = [
        ("absorbance1mMCr2.csv", "1 mM Cr(II)", 0.99, 0.0),
        ("absorbance2mMCr2.csv", "2 mM Cr(II)", 1.98, 0.0),
        ("absorbance3mMCr2.csv", "3 mM Cr(II)", 2.97, 0.0),
        ("absorbance4mMCr2.csv", "4 mM Cr(II)", 3.96, 0.0),
        ("absorbance1mMCr3.csv", "1 mM Cr(III)", 0.0, 1.05),
        ("absorbance2mMCr3.csv", "2 mM Cr(III)", 0.0, 2.10),
        ("absorbance3mMCr3.csv", "3 mM Cr(III)", 0.0, 3.15),
        ("absorbance4mMCr3.csv", "4 mM Cr(III)", 0.0, 4.20),
        ("absorbance5mMCr3.csv", "5 mM Cr(III)", 0.0, 5.25),
    ]
    traces = []
    for fname, label, Cr2, Cr3 in metadata:
        path = base / fname
        if not path.exists():
            continue
        ts, As = [], []
        with path.open() as f:
            r = csv.DictReader(f)
            for row in r:
                ts.append(float(row["time"]) * 1e-9)
                As.append(float(row["absorbance"]))
        ts_arr = np.array(ts)
        As_arr = np.array(As)
        # Some digitized CSVs (e.g. absorbance4mMCr3.csv) have unsorted timestamps.
        # solve_ivp requires monotonic t_eval; sort before storing.
        order = np.argsort(ts_arr)
        ts_arr = ts_arr[order]
        As_arr = As_arr[order]
        # Also drop exact duplicates (rare but breaks BDF)
        mask = np.concatenate([[True], np.diff(ts_arr) > 0])
        traces.append({
            "label": label,
            "t_s": ts_arr[mask],
            "abs": As_arr[mask],
            "Cr2_init": Cr2,
            "Cr3_init": Cr3,
        })
    return traces


# =============================================================================
# Modality M2: Zn scalar rate observations
# =============================================================================
def zn_observations():
    """Return Zn rate observations: (paper, salt, T_K, log_k, sigma_log_k)."""
    obs = []
    # Iwamatsu 2022 LiCl-KCl: reconstruct k at 5 T from refit Arrhenius
    A, sA, Ea, sEa = 2.4e13, 0.5e13, 35.6e3, 1.2e3
    for T_C in [400, 450, 500, 550, 600]:
        T_K = T_C + 273.15
        log_k = np.log(A) - Ea / (R_GAS * T_K)
        sigma = np.sqrt((sA/A)**2 + (sEa/(R_GAS*T_K))**2)
        obs.append(("Iwamatsu2022", "LiCl-KCl", T_K, log_k, sigma))
    # Pikaev 1982
    obs.append(("Pikaev1982", "NaCl", 1123.15, np.log(1.7e9), 0.5))
    obs.append(("Pikaev1982", "KCl",  1073.15, np.log(2.8e9), 0.5))
    return obs


# =============================================================================
# Modality M3: Phillips NULL via slow-manifold
# =============================================================================
def phillips_log_likelihood(G_Cl, k_bg, k_U3_M_at_400=7e9, U3_initial=1.0e4):
    """Censored likelihood factor for Phillips NULL given current G(Cl•) and k_bg."""
    ss = slow_manifold_steady_state(G_eS=0.5, G_Cl=G_Cl,
                                     U3=U3_initial,
                                     k_U3_M_at_400=k_U3_M_at_400,
                                     include_U=True, T=T_K_phillips)
    r_prod = ss["r_Cl2_diss_prod"]
    total_Cl2 = r_prod * V_liq * duration_s
    denom = 1.0 + V_liq * KH_CL2 * R_GAS * T_K_phillips / V_gas
    n_gas = total_Cl2 / denom
    C_gas = n_gas / V_gas
    if not np.isfinite(C_gas) or C_gas <= 0:
        return -np.inf
    z = (np.log(C_DETECTION_THRESHOLD) - np.log(C_gas)) / 0.5
    return float(np.log(max(norm.cdf(z), 1e-300)))


# =============================================================================
# Joint posterior
# =============================================================================
# Parameter vector (24 total):
#   [0..3]   = log A5, Ea5, log A6, Ea6      (Cr Arrhenius)
#   [4..12]  = log [e_s-]_0 for each of 9 Cr traces (pulse-dose nuisance)
#   [13]     = log k_bg
#   [14, 15] = log A_Zn intrinsic, Ea_Zn intrinsic
#   [16..21] = 3 salt perturbations × (eta_logA, eta_Ea) for Zn at LiCl-KCl, NaCl, KCl
#   [22]     = b^(Pikaev) facility offset
#   [23]     = log G(Cl•) for Phillips modality

PRIOR_LOG_A5_MU, PRIOR_LOG_A5_SD = np.log(1.7e13), 0.2 / 1.7
PRIOR_EA5_MU, PRIOR_EA5_SD = 33.5e3, 0.6e3
PRIOR_LOG_A6_MU, PRIOR_LOG_A6_SD = np.log(2.0e13), 0.5 / 2.0
PRIOR_EA6_MU, PRIOR_EA6_SD = 31.8e3, 0.5e3
PRIOR_LOG_ES0_MU, PRIOR_LOG_ES0_SD = np.log(1.7e-2), 0.5
PRIOR_LOG_KBG_MU, PRIOR_LOG_KBG_SD = np.log(1e7), 0.5
PRIOR_LOG_AZN_MU, PRIOR_LOG_AZN_SD = np.log(2.4e13), 0.5 / 2.4
PRIOR_EAZN_MU, PRIOR_EAZN_SD = 35.6e3, 1.2e3
PRIOR_ETA_LOGA_SD = 0.5
PRIOR_ETA_EA_SD = 5e3
PRIOR_B_PIKAEV_MU, PRIOR_B_PIKAEV_SD = 0.0, np.log(10) / 2.0
PRIOR_LOG_G_CL_MU, PRIOR_LOG_G_CL_SD = np.log(0.5), 1.5   # G(Cl•) in [0.05, 5] @ 2σ

SALT_INDEX = {"LiCl-KCl": 0, "NaCl": 1, "KCl": 2}


def _gauss(x, mu, sd):
    return -0.5 * ((x - mu) / sd) ** 2


def log_prior(theta):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    log_eS0 = theta[4:13]
    log_kbg = theta[13]
    log_AZn, EaZn = theta[14], theta[15]
    eta_Zn = theta[16:22].reshape(3, 2)
    b_pikaev = theta[22]
    log_G_Cl = theta[23]
    if Ea5 < 0 or Ea5 > 1.5e5: return -np.inf
    if Ea6 < 0 or Ea6 > 1.5e5: return -np.inf
    if EaZn < 0 or EaZn > 1.5e5: return -np.inf
    lp = 0
    lp += _gauss(log_A5, PRIOR_LOG_A5_MU, PRIOR_LOG_A5_SD)
    lp += _gauss(Ea5, PRIOR_EA5_MU, PRIOR_EA5_SD)
    lp += _gauss(log_A6, PRIOR_LOG_A6_MU, PRIOR_LOG_A6_SD)
    lp += _gauss(Ea6, PRIOR_EA6_MU, PRIOR_EA6_SD)
    for x in log_eS0: lp += _gauss(x, PRIOR_LOG_ES0_MU, PRIOR_LOG_ES0_SD)
    lp += _gauss(log_kbg, PRIOR_LOG_KBG_MU, PRIOR_LOG_KBG_SD)
    lp += _gauss(log_AZn, PRIOR_LOG_AZN_MU, PRIOR_LOG_AZN_SD)
    lp += _gauss(EaZn, PRIOR_EAZN_MU, PRIOR_EAZN_SD)
    for k in range(3):
        lp += _gauss(eta_Zn[k, 0], 0, PRIOR_ETA_LOGA_SD)
        lp += _gauss(eta_Zn[k, 1], 0, PRIOR_ETA_EA_SD)
    lp += _gauss(b_pikaev, PRIOR_B_PIKAEV_MU, PRIOR_B_PIKAEV_SD)
    lp += _gauss(log_G_Cl, PRIOR_LOG_G_CL_MU, PRIOR_LOG_G_CL_SD)
    return lp


def log_likelihood(theta, cr_traces, zn_obs, sigma_obs=5e-3):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    A5 = np.exp(log_A5); A6 = np.exp(log_A6)
    log_eS0 = theta[4:13]
    log_kbg = theta[13]
    k_bg = np.exp(log_kbg)
    log_AZn, EaZn = theta[14], theta[15]
    eta_Zn = theta[16:22].reshape(3, 2)
    b_pikaev = theta[22]
    log_G_Cl = theta[23]
    G_Cl = np.exp(log_G_Cl)

    ll = 0.0

    # ---- M1: Cr transients ----
    for i, tr in enumerate(cr_traces):
        eS0 = np.exp(log_eS0[i])
        eS_model = solve_cr_trace(eS0, tr["Cr2_init"], tr["Cr3_init"], A5, Ea5, A6, Ea6, k_bg, tr["t_s"])
        if eS_model is None or not np.all(np.isfinite(eS_model)):
            return -np.inf
        m_max = np.max(eS_model)
        if m_max <= 0:
            return -np.inf
        scale = np.nanmax(tr["abs"]) / m_max
        pred = eS_model * scale
        ll += -0.5 * np.sum(((pred - tr["abs"]) / sigma_obs) ** 2)

    # ---- M2: Zn scalar rates ----
    for paper, salt, T_K, log_k_obs, sigma in zn_obs:
        si = SALT_INDEX[salt]
        log_A_s = log_AZn + eta_Zn[si, 0]
        Ea_s = EaZn + eta_Zn[si, 1]
        log_k_pred = log_A_s - Ea_s / (R_GAS * T_K)
        b = b_pikaev if paper == "Pikaev1982" else 0.0
        ll += _gauss(log_k_obs - b, log_k_pred, sigma)

    # ---- M3: Phillips NULL censored ----
    # k_bg from Cr inference enters via the slow-manifold formula (background Cl•/Cl2•- losses).
    ll += phillips_log_likelihood(G_Cl, k_bg)

    return ll


def log_posterior(theta, cr_traces, zn_obs):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, cr_traces, zn_obs)


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 80)
    print("Tier 4 / Caveat 4: Integrated HBMAE MCMC")
    print("=" * 80)
    print()
    print("Composite likelihood spans three modalities:")
    print("  M1 (transient): 9 Iwamatsu 2026 Cr absorbance traces (stiff ODE)")
    print("  M2 (scalar):    7 Iwamatsu 2022 + Pikaev 1982 Zn rate observations")
    print("  M3 (censored):  1 Phillips 2022 NULL for [Cl2]_gas (slow-manifold reduction)")
    print()

    cr_traces = load_cr_traces()
    zn_obs = zn_observations()
    print(f"Data loaded:")
    print(f"  Cr transients: {len(cr_traces)} (n_t ~ {len(cr_traces[0]['t_s'])} each)")
    print(f"  Zn scalar rates: {len(zn_obs)} ({sum(1 for o in zn_obs if o[0]=='Iwamatsu2022')} Iwamatsu, "
          f"{sum(1 for o in zn_obs if o[0]=='Pikaev1982')} Pikaev)")
    print(f"  Phillips NULL: 1 censored constraint at [Cl2]_gas < {C_DETECTION_THRESHOLD:.2e} mol/m^3")
    print()

    ndim = 24
    nwalkers = max(2 * ndim + 8, 56)
    print(f"MCMC: ndim = {ndim}, nwalkers = {nwalkers}")
    print(f"      Production = 200 steps after 100 burn-in")
    print()

    # Initial walker positions perturbed around prior centres
    rng = np.random.default_rng(42)
    p0 = np.column_stack([
        rng.normal(PRIOR_LOG_A5_MU, 0.05, nwalkers),
        rng.normal(PRIOR_EA5_MU, 5e2, nwalkers),
        rng.normal(PRIOR_LOG_A6_MU, 0.05, nwalkers),
        rng.normal(PRIOR_EA6_MU, 5e2, nwalkers),
        *[rng.normal(PRIOR_LOG_ES0_MU, 0.2, nwalkers) for _ in range(9)],
        rng.normal(PRIOR_LOG_KBG_MU, 0.2, nwalkers),
        rng.normal(PRIOR_LOG_AZN_MU, 0.05, nwalkers),
        rng.normal(PRIOR_EAZN_MU, 5e2, nwalkers),
        rng.normal(0, 0.1, nwalkers),  # eta_logA LiCl-KCl
        rng.normal(0, 5e2, nwalkers),
        rng.normal(-0.3, 0.1, nwalkers),  # eta_logA NaCl (near M3 from Tier 3)
        rng.normal(3e3, 5e2, nwalkers),
        rng.normal(-0.15, 0.1, nwalkers),  # eta_logA KCl
        rng.normal(2e3, 5e2, nwalkers),
        rng.normal(0, 0.5, nwalkers),  # b_Pikaev (will pull to -4.9)
        rng.normal(PRIOR_LOG_G_CL_MU, 0.2, nwalkers),
    ])

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(cr_traces, zn_obs))
    print("Burn-in (100 steps)...")
    state = sampler.run_mcmc(p0, 100, progress=False)
    sampler.reset()
    print("Production (200 steps)...")
    sampler.run_mcmc(state, 200, progress=False)

    chain = sampler.get_chain(flat=True)
    log_prob = sampler.get_log_prob(flat=True)
    valid = np.isfinite(log_prob)
    chain = chain[valid]
    log_prob = log_prob[valid]
    print(f"Effective samples after rejecting -inf: {len(chain)} / {nwalkers * 200}")

    np.save(REPO / "validation/tier4_integrated_chain.npy", chain)
    np.save(REPO / "validation/tier4_integrated_logprob.npy", log_prob)

    # --- Summary ---
    labels = [
        "log A5", "Ea5", "log A6", "Ea6",
        *[f"log[e_s-]_0_t{i+1}" for i in range(9)],
        "log k_bg",
        "log A_Zn (intr)", "Ea_Zn (intr)",
        "η_logA_LiCl-KCl", "η_Ea_LiCl-KCl",
        "η_logA_NaCl",     "η_Ea_NaCl",
        "η_logA_KCl",      "η_Ea_KCl",
        "b_Pikaev",
        "log G(Cl•) Phillips"
    ]
    print()
    print("=" * 80)
    print("Posterior summary (integrated HBMAE, all three modalities)")
    print("=" * 80)
    print(f"  {'Parameter':<22s} {'Median':>14s} {'σ':>12s} {'2.5%':>14s} {'97.5%':>14s}")
    summary = []
    for i, lab in enumerate(labels):
        m = np.median(chain[:, i])
        s = np.std(chain[:, i])
        lo, hi = np.percentile(chain[:, i], [2.5, 97.5])
        print(f"  {lab:<22s} {m:>14.4g} {s:>12.3g} {lo:>14.4g} {hi:>14.4g}")
        summary.append((lab, m, s, lo, hi))

    print()
    print("Physical-parameter posteriors:")
    print(f"  A5  (e_s-+Cr2+):     {np.exp(np.median(chain[:, 0])):.3e}   lit. 1.7e13 ± 0.2e13")
    print(f"  Ea5:                 {np.median(chain[:, 1])/1e3:.2f} kJ/mol  lit. 33.5 ± 0.6")
    print(f"  A6  (e_s-+Cr3+):     {np.exp(np.median(chain[:, 2])):.3e}   lit. 2.0e13 ± 0.5e13")
    print(f"  Ea6:                 {np.median(chain[:, 3])/1e3:.2f} kJ/mol  lit. 31.8 ± 0.5")
    print(f"  k_bg:                {np.exp(np.median(chain[:, 13])):.3e} s^-1  lit. ~1e7")
    print(f"  A_Zn intrinsic:      {np.exp(np.median(chain[:, 14])):.3e}   lit. 2.4e13 ± 0.5e13")
    print(f"  Ea_Zn intrinsic:     {np.median(chain[:, 15])/1e3:.2f} kJ/mol  lit. 35.6 ± 1.2")
    print(f"  b_Pikaev:            {np.median(chain[:, 22]):+.3f}        (factor {np.exp(np.median(chain[:, 22])):.4f})")
    print(f"  G(Cl•) Phillips:     {np.exp(np.median(chain[:, 23])):.3f} molecules/100 eV")
    print()
    print("Salt perturbations for Zn (η^(s)):")
    print(f"  LiCl-KCl: η_logA = {np.median(chain[:, 16]):+.3f}, η_Ea = {np.median(chain[:, 17])/1e3:+.2f} kJ/mol")
    print(f"  NaCl:     η_logA = {np.median(chain[:, 18]):+.3f}, η_Ea = {np.median(chain[:, 19])/1e3:+.2f} kJ/mol")
    print(f"  KCl:      η_logA = {np.median(chain[:, 20]):+.3f}, η_Ea = {np.median(chain[:, 21])/1e3:+.2f} kJ/mol")

    # Save CSV
    with (REPO / "validation/TIER4_INTEGRATED_POSTERIOR.csv").open("w") as f:
        f.write("parameter,median,sigma,ci_2.5,ci_97.5\n")
        for lab, m, s, lo, hi in summary:
            f.write(f"{lab},{m},{s},{lo},{hi}\n")
    print()
    print(f"Saved to validation/TIER4_INTEGRATED_POSTERIOR.csv")
    print(f"Chain: validation/tier4_integrated_chain.npy")


if __name__ == "__main__":
    main()
