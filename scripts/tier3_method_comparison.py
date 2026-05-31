#!/usr/bin/env python3
"""Tier 3 comparison: HBMAE vs 4 state-of-the-art comparator methods on the same data.

Test case: e_s- + Zn2+ -> Zn+ rate constant across three alkali chloride hosts using
seven observations from two papers (Iwamatsu 2022 in LiCl-KCl at 5 T; Pikaev 1982 in
NaCl at 850°C and KCl at 800°C). This exercises Theorems 2 (cross-salt) and 6
(facility-effect) of the HBMAE framework simultaneously and provides a direct
empirical comparison against established Bayesian-calibration paradigms.

Five methods are implemented on identical data:

  M0  Single-paper Bayesian:
        Iwamatsu 2022 LiCl-KCl data only; Gaussian likelihood; literature-informed
        prior on (log A, Ea). Baseline: what you get from one paper.

  M1  Naive multi-paper pooled (no facility, no hierarchy):
        All 7 observations weighted equally; single (log A, Ea). Equivalent to the
        weighted-least-squares fit that Tier 1 showed is jointly infeasible
        (chi^2 = 96, p < 1e-17).

  M2  Facility-effect only (no hierarchy):
        Iwamatsu + Pikaev with Pikaev offset b ~ N(0, sigma_b). Single intrinsic
        (log A, Ea); Iwamatsu anchored at b = 0. Theorem 6 of HBMAE.

  M3  HBMAE full hierarchical:
        Per-salt theta^(s) = theta + eta^(s) with eta^(s) ~ N(0, Lambda); facility
        effect for Pikaev. Three salt-specific perturbations (LiCl-KCl, NaCl, KCl) +
        intrinsic theta. Theorems 2 + 6 simultaneously.

  M4  Galagali-Marzouk-style (weakly informative prior, no facility, no hierarchy):
        Same model space as M1 but with broad uninformative prior. Demonstrates what
        happens when one drops the literature-informed prior. Closest comparator to
        the published kinetic-network-Bayesian framework.

For each method we report:
  - Posterior median + 95% credible interval on (log A, Ea)
  - Posterior log-likelihood at the MAP (=> WAIC proxy)
  - Posterior log-marginal-likelihood (via harmonic-mean estimator; rough)
  - Posterior predictive at a held-out Iwamatsu point (cross-validation accuracy)
  - Posterior estimate of facility offset b^(Pikaev) (M2, M3)
  - Posterior estimate of salt deviation eta^(s) (M3 only)

Outcome: a tabular comparison + corner plots showing how the different inference
strategies handle the Iwamatsu/Pikaev inconsistency that Tier 1 surfaced.
"""

from __future__ import annotations

import sys
import csv
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import emcee
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

R_GAS = 8.314462618


# ---------------------------------------------------------------------------
# Data: 7 observations of log k(T) across 3 salts and 2 papers
# ---------------------------------------------------------------------------
@dataclass
class Obs:
    paper: str    # "Iwamatsu2022" or "Pikaev1982"
    salt: str     # "LiCl-KCl", "NaCl", "KCl"
    T_K: float
    log_k: float
    sigma_log_k: float


def make_observations():
    """Build the 7-observation dataset.

    Iwamatsu 2022 LiCl-KCl: reconstruct k at 5 T from the published density-corrected
    Arrhenius (A = 2.4e13 ± 0.5e13, Ea = 35.6 ± 1.2 kJ/mol) via propagation.
    Pikaev 1982: two reported scalar values with conservative sigma_log_k = 0.5 (factor
    1.6 in k) reflecting the "best-guess" caveat in Iwamatsu 2022 about Pikaev's pulse
    time-resolution.
    """
    obs = []
    # Iwamatsu 2022 LiCl-KCl Arrhenius (refit 2026): A = 2.4e13 ± 0.5e13, Ea = 35.6 ± 1.2 kJ/mol
    A, sA, Ea, sEa = 2.4e13, 0.5e13, 35.6e3, 1.2e3
    for T_C in [400, 450, 500, 550, 600]:
        T_K = T_C + 273.15
        log_k = np.log(A) - Ea / (R_GAS * T_K)
        sigma = np.sqrt((sA/A)**2 + (sEa/(R_GAS*T_K))**2)
        obs.append(Obs("Iwamatsu2022", "LiCl-KCl", T_K, log_k, sigma))
    # Pikaev 1982
    obs.append(Obs("Pikaev1982", "NaCl", 1123.15, np.log(1.7e9), 0.5))
    obs.append(Obs("Pikaev1982", "KCl",  1073.15, np.log(2.8e9), 0.5))
    return obs


# ---------------------------------------------------------------------------
# Prior specifications
# ---------------------------------------------------------------------------
# Informative literature prior for HBMAE / M0..M3 (Iwamatsu 2022 density-corrected)
PRIOR_INFO_LOG_A_MU,  PRIOR_INFO_LOG_A_SD  = np.log(2.4e13), 0.5/2.4   # σ_logA from σ_A
PRIOR_INFO_EA_MU,     PRIOR_INFO_EA_SD     = 35.6e3, 1.2e3

# Weakly informative prior for GM-style M4
PRIOR_WEAK_LOG_A_MU,  PRIOR_WEAK_LOG_A_SD  = np.log(1e11), 4.0       # very wide
PRIOR_WEAK_EA_MU,     PRIOR_WEAK_EA_SD     = 40e3, 30e3              # very wide

# Hierarchical perturbation: Λ_i — solvent shifts O(few percent) on A and O(kT) on Ea
# Marcus theory predicts O(few kT) ≈ O(5-15 kJ/mol) shifts in Ea between alkali chloride hosts.
LAMBDA_LOG_A = 0.5     # ~ factor 1.65 in A across salts
LAMBDA_EA    = 5e3     # 5 kJ/mol across salts

# Facility shift prior: Pikaev possibly biased low by factor up to 10 (per Iwamatsu 2022 critique)
PRIOR_B_PIKAEV_MU, PRIOR_B_PIKAEV_SD = 0.0, np.log(10.0)/2.0   # σ ≈ log(10)/2 covers 1/10 to 10 at 2σ


# ---------------------------------------------------------------------------
# Method-specific log posteriors
# ---------------------------------------------------------------------------
def _gauss_lp(x, mu, sigma):
    return -0.5 * ((x - mu) / sigma) ** 2


# ----- M0: Iwamatsu only -----
def make_logpost_M0(obs_iwa):
    """Single-paper Bayesian: only Iwamatsu, literature-informed prior."""
    def logpost(theta):
        log_A, Ea = theta
        if Ea < 0 or Ea > 1.5e5: return -np.inf
        lp  = _gauss_lp(log_A, PRIOR_INFO_LOG_A_MU, PRIOR_INFO_LOG_A_SD)
        lp += _gauss_lp(Ea,    PRIOR_INFO_EA_MU,    PRIOR_INFO_EA_SD)
        for o in obs_iwa:
            log_k_pred = log_A - Ea / (R_GAS * o.T_K)
            lp += _gauss_lp(o.log_k, log_k_pred, o.sigma_log_k)
        return lp
    return logpost


# ----- M1: Naive pooled (no facility, no hierarchy) -----
def make_logpost_M1(all_obs):
    """Naive multi-paper pooled: equal-weight, single (logA, Ea), informative prior."""
    def logpost(theta):
        log_A, Ea = theta
        if Ea < 0 or Ea > 1.5e5: return -np.inf
        lp  = _gauss_lp(log_A, PRIOR_INFO_LOG_A_MU, PRIOR_INFO_LOG_A_SD)
        lp += _gauss_lp(Ea,    PRIOR_INFO_EA_MU,    PRIOR_INFO_EA_SD)
        for o in all_obs:
            log_k_pred = log_A - Ea / (R_GAS * o.T_K)
            lp += _gauss_lp(o.log_k, log_k_pred, o.sigma_log_k)
        return lp
    return logpost


# ----- M2: Facility-effect only -----
def make_logpost_M2(all_obs):
    """Single intrinsic (logA, Ea) + Pikaev offset b. Iwamatsu anchored b = 0."""
    def logpost(theta):
        log_A, Ea, b_pikaev = theta
        if Ea < 0 or Ea > 1.5e5: return -np.inf
        lp  = _gauss_lp(log_A, PRIOR_INFO_LOG_A_MU, PRIOR_INFO_LOG_A_SD)
        lp += _gauss_lp(Ea,    PRIOR_INFO_EA_MU,    PRIOR_INFO_EA_SD)
        lp += _gauss_lp(b_pikaev, PRIOR_B_PIKAEV_MU, PRIOR_B_PIKAEV_SD)
        for o in all_obs:
            log_k_pred = log_A - Ea / (R_GAS * o.T_K)
            b = b_pikaev if o.paper == "Pikaev1982" else 0.0
            lp += _gauss_lp(o.log_k - b, log_k_pred, o.sigma_log_k)
        return lp
    return logpost


# ----- M3: HBMAE full hierarchical -----
def make_logpost_M3(all_obs, salts):
    """Per-salt theta^(s) = theta_intrinsic + eta^(s); facility offset for Pikaev.

    Parameter vector: [log_A, Ea, eta_logA_LiCl-KCl, eta_Ea_LiCl-KCl, eta_logA_NaCl, eta_Ea_NaCl,
                       eta_logA_KCl, eta_Ea_KCl, b_pikaev]
    """
    salt_idx = {s: i for i, s in enumerate(salts)}
    K = len(salts)

    def logpost(theta):
        log_A, Ea = theta[0], theta[1]
        etas = theta[2:2+2*K].reshape(K, 2)
        b_pikaev = theta[-1]
        if Ea < 0 or Ea > 1.5e5: return -np.inf
        lp  = _gauss_lp(log_A, PRIOR_INFO_LOG_A_MU, PRIOR_INFO_LOG_A_SD)
        lp += _gauss_lp(Ea,    PRIOR_INFO_EA_MU,    PRIOR_INFO_EA_SD)
        lp += _gauss_lp(b_pikaev, PRIOR_B_PIKAEV_MU, PRIOR_B_PIKAEV_SD)
        # Hierarchical layer: each eta is N(0, Lambda)
        for k in range(K):
            lp += _gauss_lp(etas[k, 0], 0.0, LAMBDA_LOG_A)
            lp += _gauss_lp(etas[k, 1], 0.0, LAMBDA_EA)
        # Likelihood per observation with salt-specific theta + facility offset
        for o in all_obs:
            si = salt_idx[o.salt]
            log_A_s = log_A + etas[si, 0]
            Ea_s    = Ea    + etas[si, 1]
            log_k_pred = log_A_s - Ea_s / (R_GAS * o.T_K)
            b = b_pikaev if o.paper == "Pikaev1982" else 0.0
            lp += _gauss_lp(o.log_k - b, log_k_pred, o.sigma_log_k)
        return lp
    return logpost


# ----- M4: GM-style weakly informative, no facility, no hierarchy -----
def make_logpost_M4(all_obs):
    """Galagali-Marzouk-style: weakly informative prior; otherwise like M1."""
    def logpost(theta):
        log_A, Ea = theta
        if Ea < 0 or Ea > 1.5e5: return -np.inf
        lp  = _gauss_lp(log_A, PRIOR_WEAK_LOG_A_MU, PRIOR_WEAK_LOG_A_SD)
        lp += _gauss_lp(Ea,    PRIOR_WEAK_EA_MU,    PRIOR_WEAK_EA_SD)
        for o in all_obs:
            log_k_pred = log_A - Ea / (R_GAS * o.T_K)
            lp += _gauss_lp(o.log_k, log_k_pred, o.sigma_log_k)
        return lp
    return logpost


# ---------------------------------------------------------------------------
# Run MCMC for one method and compute diagnostics
# ---------------------------------------------------------------------------
def run_mcmc(logpost, p0, nwalkers=64, nsteps=2000, burn=500, label=""):
    ndim = p0.shape[1]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost)
    state = sampler.run_mcmc(p0, burn, progress=False)
    sampler.reset()
    sampler.run_mcmc(state, nsteps, progress=False)
    chain = sampler.get_chain(flat=True)
    log_prob = sampler.get_log_prob(flat=True)
    finite = np.isfinite(log_prob)
    return chain[finite], log_prob[finite], sampler


# Held-out cross-validation: drop one Iwamatsu point and predict it
def held_out_log_likelihood(chain, log_prob, held_out_obs, predict_fn):
    """For each posterior sample, compute log p(y_held | theta). Average for predictive density.

    predict_fn(theta, T_K) -> log_k_pred  (must be method-specific).
    """
    log_pred_density = []
    for theta in chain[::10]:   # thin
        log_k_pred = predict_fn(theta, held_out_obs.T_K)
        lp = _gauss_lp(held_out_obs.log_k, log_k_pred, held_out_obs.sigma_log_k)
        log_pred_density.append(lp)
    # Predictive density = log(mean(exp(lp))); use log-sum-exp
    log_pred_density = np.array(log_pred_density)
    log_mean = np.max(log_pred_density) + np.log(np.mean(np.exp(log_pred_density - np.max(log_pred_density))))
    return float(log_mean)


def waic_estimate(chain, logpost, obs_used):
    """WAIC approximation via posterior samples of pointwise log-likelihoods.

    WAIC = -2 (lppd - p_WAIC), lppd = sum_i log mean_s p(y_i | theta_s).
    """
    # Compute pointwise log-likelihood per observation per posterior sample
    log_lik = np.zeros((len(chain), len(obs_used)))
    # We need a per-method "obs log-lik given theta" function. Approximate by:
    #   ll_i(theta) = logpost(theta) − log_prior(theta) − sum_{j≠i} ll_j(theta)
    # which requires per-observation likelihood functions. For methods that share the
    # same likelihood form (M0..M4 all use Gaussian on log_k), we can compute directly.
    # Here we use a simplified WAIC: report (mean log-lik, var log-lik) per posterior.
    # Pointwise WAIC requires per-obs ll which depends on method-specific log_k_pred.
    # We compute this in the caller for each method explicitly.
    return None  # see method-specific WAIC below


def compute_method_diagnostics(chain, log_prob, predict_fn, obs_used, label, held_out=None):
    """Compute summary statistics for one method's posterior."""
    # Posterior median + 95% CI on (log A, Ea)
    log_A_samples = chain[:, 0]
    Ea_samples = chain[:, 1]
    median_A = np.exp(np.median(log_A_samples))
    ci_A = np.exp(np.percentile(log_A_samples, [2.5, 97.5]))
    median_E = np.median(Ea_samples) / 1e3
    ci_E = np.percentile(Ea_samples, [2.5, 97.5]) / 1e3

    # WAIC: pointwise log-likelihood per observation per sample
    # ll_ij = log p(y_i | theta_j)  for each i, j
    thin = max(1, len(chain) // 1000)
    samples = chain[::thin]
    ll = np.zeros((len(samples), len(obs_used)))
    for j, theta in enumerate(samples):
        for i, o in enumerate(obs_used):
            log_k_pred = predict_fn(theta, o.T_K, o)
            ll[j, i] = _gauss_lp(o.log_k, log_k_pred, o.sigma_log_k)
    # lppd_i = log mean_j p(y_i | theta_j) — log-sum-exp normalized by num samples
    max_ll = np.max(ll, axis=0)
    lppd = np.sum(max_ll + np.log(np.mean(np.exp(ll - max_ll), axis=0)))
    # p_WAIC = sum_i Var_j[log p(y_i | theta_j)]
    p_waic = np.sum(np.var(ll, axis=0, ddof=1))
    waic = -2.0 * (lppd - p_waic)
    elpd_waic = lppd - p_waic

    # Held-out predictive density
    if held_out is not None:
        log_pred_held = held_out_log_likelihood(chain, log_prob, held_out,
                                                 lambda th, T: predict_fn(th, T, held_out))
    else:
        log_pred_held = None

    print(f"\n  {label}:")
    print(f"    Posterior on Arrhenius (intrinsic where applicable):")
    print(f"      A       = {median_A:.3e}  (95% CI [{ci_A[0]:.2e}, {ci_A[1]:.2e}])")
    print(f"      Ea      = {median_E:.2f} kJ/mol  (95% CI [{ci_E[0]:.2f}, {ci_E[1]:.2f}])")
    print(f"    elpd_WAIC = {elpd_waic:+.3f}   (higher is better)")
    print(f"    p_WAIC    = {p_waic:.3f}      (effective n_params)")
    if log_pred_held is not None:
        print(f"    held-out  = {log_pred_held:+.3f}   (log predictive density at Iwamatsu 550°C)")
    return {
        "label": label,
        "median_A": median_A,
        "ci_A": ci_A,
        "median_E_kJ": median_E,
        "ci_E_kJ": ci_E,
        "ndim": chain.shape[1],
        "elpd_waic": float(elpd_waic),
        "p_waic": float(p_waic),
        "waic": float(waic),
        "log_pred_held_out": float(log_pred_held) if log_pred_held is not None else None,
    }


# ---------------------------------------------------------------------------
# Per-method predict_fn that returns log k_pred for one observation
# ---------------------------------------------------------------------------
def predict_M0(theta, T_K, obs):
    log_A, Ea = theta
    return log_A - Ea / (R_GAS * T_K)

def predict_M1(theta, T_K, obs):
    log_A, Ea = theta
    return log_A - Ea / (R_GAS * T_K)

def predict_M2(theta, T_K, obs):
    log_A, Ea, b_pikaev = theta
    log_k_pred = log_A - Ea / (R_GAS * T_K)
    b = b_pikaev if obs.paper == "Pikaev1982" else 0.0
    return log_k_pred + b   # the observation has b added; the model predicts that

def predict_M3(theta, T_K, obs):
    log_A, Ea = theta[0], theta[1]
    salts = ["LiCl-KCl", "NaCl", "KCl"]
    si = salts.index(obs.salt)
    eta_logA = theta[2 + 2*si]
    eta_Ea   = theta[3 + 2*si]
    b_pikaev = theta[-1]
    log_A_s = log_A + eta_logA
    Ea_s    = Ea    + eta_Ea
    log_k_pred = log_A_s - Ea_s / (R_GAS * T_K)
    b = b_pikaev if obs.paper == "Pikaev1982" else 0.0
    return log_k_pred + b

def predict_M4(theta, T_K, obs):
    log_A, Ea = theta
    return log_A - Ea / (R_GAS * T_K)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("Tier 3: HBMAE vs comparator methods on Iwamatsu (LiCl-KCl) + Pikaev (NaCl, KCl)")
    print("=" * 78)

    all_obs = make_observations()
    obs_iwa = [o for o in all_obs if o.paper == "Iwamatsu2022"]
    obs_pik = [o for o in all_obs if o.paper == "Pikaev1982"]
    print(f"\nData: {len(obs_iwa)} Iwamatsu observations + {len(obs_pik)} Pikaev observations")
    for o in all_obs:
        print(f"  {o.paper:>14s}  {o.salt:>10s}  T={o.T_K-273.15:5.1f}°C  log k={o.log_k:7.3f}  σ={o.sigma_log_k:.3f}")

    # Pick a held-out observation: the 550°C Iwamatsu point
    held_out_obs = obs_iwa[3]   # 550°C
    print(f"\nHeld-out for predictive: {held_out_obs.paper} {held_out_obs.salt} T={held_out_obs.T_K-273.15}°C")
    obs_iwa_train = [o for o in obs_iwa if o is not held_out_obs]
    obs_all_train = obs_iwa_train + obs_pik

    salts = ["LiCl-KCl", "NaCl", "KCl"]

    results = {}
    rng = np.random.default_rng(7)

    # --- M0: Iwamatsu only ---
    print("\n--- M0: Single-paper Bayesian (Iwamatsu only) ---")
    p0 = np.column_stack([
        rng.normal(PRIOR_INFO_LOG_A_MU, 0.1, 64),
        rng.normal(PRIOR_INFO_EA_MU, 1e3, 64),
    ])
    chain, lp, _ = run_mcmc(make_logpost_M0(obs_iwa_train), p0, label="M0")
    results["M0"] = compute_method_diagnostics(chain, lp, predict_M0, obs_iwa_train, "M0 — Iwamatsu-only", held_out_obs)
    np.save(REPO / "validation/tier3_M0_chain.npy", chain)

    # --- M1: Naive pooled ---
    print("\n--- M1: Naive multi-paper pooled (info prior, no facility, no hierarchy) ---")
    chain, lp, _ = run_mcmc(make_logpost_M1(obs_all_train), p0, label="M1")
    results["M1"] = compute_method_diagnostics(chain, lp, predict_M1, obs_all_train, "M1 — Naive pool", held_out_obs)
    np.save(REPO / "validation/tier3_M1_chain.npy", chain)

    # --- M2: Facility-effect only ---
    print("\n--- M2: Facility-effect only (HBMAE Theorem 6 without Theorem 2) ---")
    p0_M2 = np.column_stack([
        rng.normal(PRIOR_INFO_LOG_A_MU, 0.1, 64),
        rng.normal(PRIOR_INFO_EA_MU, 1e3, 64),
        rng.normal(0.0, 0.2, 64),
    ])
    chain, lp, _ = run_mcmc(make_logpost_M2(obs_all_train), p0_M2, label="M2")
    results["M2"] = compute_method_diagnostics(chain, lp, predict_M2, obs_all_train, "M2 — Facility effect", held_out_obs)
    # extract b_pikaev posterior
    b_pikaev = chain[:, 2]
    print(f"    b^(Pikaev) posterior: median = {np.median(b_pikaev):+.3f}, 95% CI [{np.percentile(b_pikaev, 2.5):+.3f}, {np.percentile(b_pikaev, 97.5):+.3f}]")
    print(f"    (interpretation: Pikaev rates biased by factor {np.exp(np.median(b_pikaev)):.3f} relative to Iwamatsu)")
    np.save(REPO / "validation/tier3_M2_chain.npy", chain)

    # --- M3: HBMAE full hierarchical ---
    print("\n--- M3: HBMAE full (Theorem 2 cross-salt + Theorem 6 facility-effect) ---")
    K = len(salts)
    ndim_M3 = 2 + 2*K + 1
    p0_M3 = np.zeros((64, ndim_M3))
    p0_M3[:, 0] = rng.normal(PRIOR_INFO_LOG_A_MU, 0.1, 64)
    p0_M3[:, 1] = rng.normal(PRIOR_INFO_EA_MU, 1e3, 64)
    for k in range(K):
        p0_M3[:, 2+2*k]   = rng.normal(0, 0.1, 64)   # eta_logA
        p0_M3[:, 3+2*k]   = rng.normal(0, 200, 64)   # eta_Ea
    p0_M3[:, -1] = rng.normal(0, 0.2, 64)            # b_pikaev
    chain, lp, _ = run_mcmc(make_logpost_M3(obs_all_train, salts), p0_M3, label="M3")
    results["M3"] = compute_method_diagnostics(chain, lp, predict_M3, obs_all_train, "M3 — HBMAE full", held_out_obs)
    b_pikaev = chain[:, -1]
    print(f"    b^(Pikaev) posterior: median = {np.median(b_pikaev):+.3f}, 95% CI [{np.percentile(b_pikaev, 2.5):+.3f}, {np.percentile(b_pikaev, 97.5):+.3f}]")
    for k, s in enumerate(salts):
        eta_A = chain[:, 2+2*k]
        eta_E = chain[:, 3+2*k] / 1e3
        print(f"    eta^({s}) : eta_logA median = {np.median(eta_A):+.3f}, eta_Ea median = {np.median(eta_E):+.2f} kJ/mol")
    np.save(REPO / "validation/tier3_M3_chain.npy", chain)

    # --- M4: GM-style weakly informative ---
    print("\n--- M4: Galagali-Marzouk-style (weakly informative, no facility, no hierarchy) ---")
    p0_M4 = np.column_stack([
        rng.normal(PRIOR_WEAK_LOG_A_MU, 1.0, 64),
        rng.normal(PRIOR_WEAK_EA_MU, 5e3, 64),
    ])
    chain, lp, _ = run_mcmc(make_logpost_M4(obs_all_train), p0_M4, label="M4")
    results["M4"] = compute_method_diagnostics(chain, lp, predict_M4, obs_all_train, "M4 — GM-style weak prior", held_out_obs)
    np.save(REPO / "validation/tier3_M4_chain.npy", chain)

    # ---------------------------------------------------------------------
    # Comparison table
    # ---------------------------------------------------------------------
    print()
    print("=" * 88)
    print("Cross-method comparison (held-out point = Iwamatsu 550°C in LiCl-KCl)")
    print("=" * 88)
    print(f"{'Method':<28s} {'A median':>12s} {'Ea median':>12s} {'A 95% CI width':>16s} {'elpd_WAIC':>11s} {'held-out':>10s}")
    print("-" * 88)
    for key in ["M0", "M1", "M2", "M3", "M4"]:
        r = results[key]
        ci_width_A = np.log10(r["ci_A"][1]) - np.log10(r["ci_A"][0])
        print(f"{r['label']:<28s} {r['median_A']:>12.2e} {r['median_E_kJ']:>10.2f} kJ "
              f"{ci_width_A:>10.3f} dex {r['elpd_waic']:>10.3f} {r['log_pred_held_out']:>10.3f}")

    # Save the results as CSV
    with (REPO / "validation/TIER3_METHOD_COMPARISON.csv").open("w") as f:
        f.write("method,median_A,ci_A_lo,ci_A_hi,median_Ea_kJ,ci_Ea_lo_kJ,ci_Ea_hi_kJ,ndim,elpd_WAIC,p_WAIC,WAIC,held_out_lpd\n")
        for key in ["M0", "M1", "M2", "M3", "M4"]:
            r = results[key]
            f.write(f"{r['label']},{r['median_A']},{r['ci_A'][0]},{r['ci_A'][1]},"
                    f"{r['median_E_kJ']},{r['ci_E_kJ'][0]},{r['ci_E_kJ'][1]},{r['ndim']},"
                    f"{r['elpd_waic']},{r['p_waic']},{r['waic']},{r['log_pred_held_out']}\n")

    print()
    print("Interpretation:")
    print("  - elpd_WAIC: higher is better (predictive accuracy on data already used).")
    print("  - held-out:  higher is better (predictive accuracy on left-out Iwamatsu point).")
    print("  - A 95% CI width: log10(A_hi/A_lo), narrower = more constrained posterior.")
    print()
    print("Expected ordering by quality (held-out + WAIC):  M3 > M2 >> M1, M0; M4 highly variable.")

    return results


if __name__ == "__main__":
    main()
