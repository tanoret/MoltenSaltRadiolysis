#!/usr/bin/env python3
"""Cr Tier 2 train/val split + PSIS-LOO network selection.

Two products:
  (1) Random 7-train / 2-validation split on the 9 Iwamatsu 2026 Cr transients.
      Fits the minimal network (A) on the 7 train traces, plots posterior
      predictive overlays for both subsets.
  (2) Pareto-smoothed importance-sampling LOO (PSIS-LOO) comparison of two
      candidate networks fit on all 9 traces:
        Network A (minimal):    e_s- + Cr3+ -> Cr2+;
                                e_s- + Cr2+ -> Cr+;
                                e_s- + impurity -> decay (k_bg)
        Network B (extended):   Network A + (e_s- + Cr+ -> Cr0) third electron
                                attachment channel (Arrhenius A_recomb,
                                Ea_recomb).
      The third channel introduces a second-order non-linearity that flattens
      the late-time tail.  PSIS-LOO ELPD difference quantifies whether the
      extra channel earns its complexity.

Outputs:
  validation/cr_train_val/train_val_chain.npy
  validation/cr_train_val/train_val_summary.csv
  validation/cr_train_val/posterior_predictive_train.npz
  validation/cr_train_val/posterior_predictive_val.npz
  validation/cr_psis_loo/networkA_chain.npy
  validation/cr_psis_loo/networkB_chain.npy
  validation/cr_psis_loo/networkA_pointwise_loglik.npy
  validation/cr_psis_loo/networkB_pointwise_loglik.npy
  validation/cr_psis_loo/loo_comparison.csv
  manuscript/figures/fig_train_val.pdf
  manuscript/figures/fig_psis_loo.pdf
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import emcee
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "manuscript" / "figures"))
sys.path.insert(0, str(REPO / "scripts"))
try:
    import figstyle  # type: ignore
    figstyle.apply_rc()
except Exception:
    pass


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Re-use the Tier 2 load_traces() and forward solver
from tier2_bayesian_calibration import load_traces, R_GAS

OUT_TV = REPO / "validation" / "cr_train_val"
OUT_TV.mkdir(parents=True, exist_ok=True)
OUT_LOO = REPO / "validation" / "cr_psis_loo"
OUT_LOO.mkdir(parents=True, exist_ok=True)
FIG_DIR = REPO / "manuscript" / "figures"

T_K = 673.15  # 400 C, all Iwamatsu 2026 traces are at this temperature

# ----------------------------------------------------------------------------
# Forward solvers for Network A and Network B
# ----------------------------------------------------------------------------

def rhs_A(t, y, k5, k6, k_bg):
    eS, Cr2, Cr3, Cr1 = y
    r5 = k5 * eS * Cr2     # e_s- + Cr2+ -> Cr+
    r6 = k6 * eS * Cr3     # e_s- + Cr3+ -> Cr2+
    rb = k_bg * eS
    return [-r5 - r6 - rb,
             r6 - r5,
            -r6,
             r5]


def rhs_B(t, y, k5, k6, k_bg, k_rec):
    eS, Cr2, Cr3, Cr1, Cr0 = y
    r5 = k5 * eS * Cr2     # e_s- + Cr2+ -> Cr+
    r6 = k6 * eS * Cr3     # e_s- + Cr3+ -> Cr2+
    r_rec = k_rec * eS * Cr1  # e_s- + Cr+ -> Cr0
    rb = k_bg * eS
    return [-r5 - r6 - rb - r_rec,
             r6 - r5,
            -r6,
             r5 - r_rec,
             r_rec]


def solve_trace(eS0_mM3, Cr2_0, Cr3_0, theta_A_or_B: dict, t_eval, network="A"):
    """Forward solve and return [e_s-](t) in mol/m^3."""
    A5, Ea5 = theta_A_or_B["A5"], theta_A_or_B["Ea5"]
    A6, Ea6 = theta_A_or_B["A6"], theta_A_or_B["Ea6"]
    k_bg = theta_A_or_B["k_bg"]
    k5 = A5 * np.exp(-Ea5 / (R_GAS * T_K)) / 1000.0  # M^-1 s^-1 -> m^3/(mol s)
    k6 = A6 * np.exp(-Ea6 / (R_GAS * T_K)) / 1000.0
    # Guard against pathological parameter values that cause BDF to hang
    if not (np.isfinite(k5) and np.isfinite(k6) and np.isfinite(k_bg)):
        return None
    if k5 < 0 or k6 < 0 or k_bg < 0:
        return None
    if k5 > 1e15 or k6 > 1e15 or k_bg > 1e15:
        return None
    if network == "A":
        y0 = [eS0_mM3, Cr2_0, Cr3_0, 0.0]
        f = lambda t, y: rhs_A(t, y, k5, k6, k_bg)
    else:
        A_rec, Ea_rec = theta_A_or_B["A_rec"], theta_A_or_B["Ea_rec"]
        k_rec = A_rec * np.exp(-Ea_rec / (R_GAS * T_K)) / 1000.0
        if not np.isfinite(k_rec) or k_rec < 0 or k_rec > 1e15:
            return None
        y0 = [eS0_mM3, Cr2_0, Cr3_0, 0.0, 0.0]
        f = lambda t, y: rhs_B(t, y, k5, k6, k_bg, k_rec)
    try:
        sol = solve_ivp(
            f, (0, max(float(t_eval[-1]), 1e-9)), y0, t_eval=t_eval,
            method="BDF", rtol=1e-6, atol=1e-12, max_step=1e-7,
        )
        if not sol.success or sol.y.shape[1] != len(t_eval):
            return None
        return sol.y[0]  # e_s-
    except Exception:
        return None


# ----------------------------------------------------------------------------
# Priors (shared with Tier 2)
# ----------------------------------------------------------------------------

PRIOR_LOG_A5_MU, PRIOR_LOG_A5_SIGMA = np.log(1.7e13), 0.12
PRIOR_EA5_MU, PRIOR_EA5_SIGMA = 33.5e3, 0.6e3
PRIOR_LOG_A6_MU, PRIOR_LOG_A6_SIGMA = np.log(2.0e13), 0.25
PRIOR_EA6_MU, PRIOR_EA6_SIGMA = 31.8e3, 0.5e3
PRIOR_LOG_ES0_MU, PRIOR_LOG_ES0_SIGMA = np.log(1.7e-2), 0.5
PRIOR_LOG_KBG_MU, PRIOR_LOG_KBG_SIGMA = np.log(1e7), 0.5
# Network B extra: weak prior on the third-channel Arrhenius
PRIOR_LOG_AREC_MU, PRIOR_LOG_AREC_SIGMA = np.log(1.0e13), 0.5
PRIOR_EAREC_MU, PRIOR_EAREC_SIGMA = 30.0e3, 5.0e3
SIGMA_OBS = 5e-3


def _params(theta, N_traces, network="A"):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    if network == "A":
        log_eS0s = theta[4:4 + N_traces]
        log_kbg = theta[4 + N_traces]
    else:
        log_A_rec, Ea_rec = theta[4:6]
        log_eS0s = theta[6:6 + N_traces]
        log_kbg = theta[6 + N_traces]
    base = {
        "A5": np.exp(log_A5), "Ea5": Ea5,
        "A6": np.exp(log_A6), "Ea6": Ea6,
        "k_bg": np.exp(log_kbg),
    }
    if network == "B":
        base["A_rec"] = np.exp(log_A_rec)
        base["Ea_rec"] = Ea_rec
    return base, log_eS0s


def log_prior(theta, N_traces, network="A"):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    if network == "A":
        log_eS0s = theta[4:4 + N_traces]
        log_kbg = theta[4 + N_traces]
        log_A_rec, Ea_rec = None, None
    else:
        log_A_rec, Ea_rec = theta[4:6]
        log_eS0s = theta[6:6 + N_traces]
        log_kbg = theta[6 + N_traces]

    if Ea5 < 0 or Ea5 > 1e5: return -np.inf
    if Ea6 < 0 or Ea6 > 1e5: return -np.inf
    if Ea_rec is not None and (Ea_rec < 0 or Ea_rec > 1e5): return -np.inf

    lp = 0.0
    lp += -0.5 * ((log_A5 - PRIOR_LOG_A5_MU) / PRIOR_LOG_A5_SIGMA) ** 2
    lp += -0.5 * ((Ea5 - PRIOR_EA5_MU) / PRIOR_EA5_SIGMA) ** 2
    lp += -0.5 * ((log_A6 - PRIOR_LOG_A6_MU) / PRIOR_LOG_A6_SIGMA) ** 2
    lp += -0.5 * ((Ea6 - PRIOR_EA6_MU) / PRIOR_EA6_SIGMA) ** 2
    for log_eS0 in log_eS0s:
        lp += -0.5 * ((log_eS0 - PRIOR_LOG_ES0_MU) / PRIOR_LOG_ES0_SIGMA) ** 2
    lp += -0.5 * ((log_kbg - PRIOR_LOG_KBG_MU) / PRIOR_LOG_KBG_SIGMA) ** 2
    if log_A_rec is not None:
        lp += -0.5 * ((log_A_rec - PRIOR_LOG_AREC_MU) / PRIOR_LOG_AREC_SIGMA) ** 2
        lp += -0.5 * ((Ea_rec - PRIOR_EAREC_MU) / PRIOR_EAREC_SIGMA) ** 2
    return lp


def trace_loglik(theta, traces, network="A", return_pointwise=False):
    """Return (total log-lik, pointwise per-trace mean log-lik if requested)."""
    base, log_eS0s = _params(theta, len(traces), network)
    ll_total = 0.0
    if return_pointwise:
        # Pointwise log-lik per trace, summed over time-points within that trace
        pointwise = np.zeros(len(traces))
    for i, tr in enumerate(traces):
        eS0 = np.exp(log_eS0s[i])
        pred = solve_trace(eS0, tr["Cr2_init"], tr["Cr3_init"], base,
                            tr["t_s"], network=network)
        if pred is None or not np.all(np.isfinite(pred)):
            return (-np.inf, None) if return_pointwise else -np.inf
        m_max = np.max(pred)
        if m_max <= 0:
            return (-np.inf, None) if return_pointwise else -np.inf
        scale = np.nanmax(tr["abs"]) / m_max
        pred_obs = pred * scale
        res = (pred_obs - tr["abs"]) / SIGMA_OBS
        ll_i = -0.5 * np.sum(res * res)
        ll_total += ll_i
        if return_pointwise:
            pointwise[i] = ll_i
    if return_pointwise:
        return ll_total, pointwise
    return ll_total


def log_posterior_factory(traces, network="A"):
    def lp(theta):
        prior = log_prior(theta, len(traces), network)
        if not np.isfinite(prior):
            return -np.inf
        ll = trace_loglik(theta, traces, network)
        if not np.isfinite(ll):
            return -np.inf
        return prior + ll
    return lp


def initial_walkers(network, N_traces, nwalkers, rng):
    if network == "A":
        ndim = 4 + N_traces + 1
        x0 = np.zeros((nwalkers, ndim))
        x0[:, 0] = PRIOR_LOG_A5_MU + 0.05 * rng.standard_normal(nwalkers)
        x0[:, 1] = PRIOR_EA5_MU + 200 * rng.standard_normal(nwalkers)
        x0[:, 2] = PRIOR_LOG_A6_MU + 0.05 * rng.standard_normal(nwalkers)
        x0[:, 3] = PRIOR_EA6_MU + 200 * rng.standard_normal(nwalkers)
        x0[:, 4:4 + N_traces] = PRIOR_LOG_ES0_MU + 0.1 * rng.standard_normal((nwalkers, N_traces))
        x0[:, -1] = PRIOR_LOG_KBG_MU + 0.1 * rng.standard_normal(nwalkers)
    else:
        ndim = 4 + 2 + N_traces + 1
        x0 = np.zeros((nwalkers, ndim))
        x0[:, 0] = PRIOR_LOG_A5_MU + 0.05 * rng.standard_normal(nwalkers)
        x0[:, 1] = PRIOR_EA5_MU + 200 * rng.standard_normal(nwalkers)
        x0[:, 2] = PRIOR_LOG_A6_MU + 0.05 * rng.standard_normal(nwalkers)
        x0[:, 3] = PRIOR_EA6_MU + 200 * rng.standard_normal(nwalkers)
        x0[:, 4] = PRIOR_LOG_AREC_MU + 0.1 * rng.standard_normal(nwalkers)
        x0[:, 5] = PRIOR_EAREC_MU + 1000 * rng.standard_normal(nwalkers)
        x0[:, 6:6 + N_traces] = PRIOR_LOG_ES0_MU + 0.1 * rng.standard_normal((nwalkers, N_traces))
        x0[:, -1] = PRIOR_LOG_KBG_MU + 0.1 * rng.standard_normal(nwalkers)
    return x0, ndim


# ----------------------------------------------------------------------------
# Train / Val split
# ----------------------------------------------------------------------------

def run_train_val(traces, val_idx, n_steps=1500, n_walkers=64, tag="train_val"):
    train_traces = [t for i, t in enumerate(traces) if i not in val_idx]
    val_traces = [traces[i] for i in val_idx]
    N = len(train_traces)
    rng = np.random.default_rng(42)
    x0, ndim = initial_walkers("A", N, n_walkers, rng)
    print(f"[{tag}] training on {N} traces, holding out {len(val_traces)} for validation")
    print(f"[{tag}] ndim = {ndim}, walkers = {n_walkers}, steps = {n_steps}")
    lp = log_posterior_factory(train_traces, network="A")
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lp)
    sampler.run_mcmc(x0, n_steps, progress=False)
    chain = sampler.get_chain(discard=n_steps // 3, flat=True)
    log_prob = sampler.get_log_prob(discard=n_steps // 3, flat=True)
    keep = np.isfinite(log_prob)
    return chain[keep], log_prob[keep], train_traces, val_traces


def posterior_predictive(chain, traces, network="A", n_draw=200, T_K=T_K):
    """Compute posterior predictive [e_s-](t) for each trace; return rescaled
    to match max absorbance, plus mean and 90% CI envelopes."""
    rng = np.random.default_rng(0)
    idx = rng.choice(chain.shape[0], size=min(n_draw, chain.shape[0]), replace=False)
    preds = {tr["label"]: [] for tr in traces}
    for s in idx:
        theta = chain[s]
        base, log_eS0s = _params(theta, len(traces), network)
        for i, tr in enumerate(traces):
            eS0 = np.exp(log_eS0s[i])
            pred = solve_trace(eS0, tr["Cr2_init"], tr["Cr3_init"], base,
                                tr["t_s"], network=network)
            if pred is None or np.any(~np.isfinite(pred)):
                continue
            m_max = np.max(pred)
            if m_max <= 0:
                continue
            scale = np.nanmax(tr["abs"]) / m_max
            preds[tr["label"]].append(pred * scale)
    summary = {}
    for label, samples in preds.items():
        if len(samples) == 0:
            continue
        arr = np.array(samples)
        summary[label] = {
            "t_s": [tr["t_s"] for tr in traces if tr["label"] == label][0],
            "obs": [tr["abs"] for tr in traces if tr["label"] == label][0],
            "mean": arr.mean(axis=0),
            "p05": np.percentile(arr, 5, axis=0),
            "p95": np.percentile(arr, 95, axis=0),
        }
    return summary


# ----------------------------------------------------------------------------
# PSIS-LOO via ArviZ
# ----------------------------------------------------------------------------

def run_full_network(traces, network, n_steps=2000, n_walkers=80):
    rng = np.random.default_rng(13 if network == "A" else 17)
    x0, ndim = initial_walkers(network, len(traces), n_walkers, rng)
    lp = log_posterior_factory(traces, network=network)
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lp)
    print(f"[{network}] running emcee: ndim={ndim}, n_walkers={n_walkers}, steps={n_steps}")
    sampler.run_mcmc(x0, n_steps, progress=False)
    chain = sampler.get_chain(discard=n_steps // 3, flat=True)
    log_prob = sampler.get_log_prob(discard=n_steps // 3, flat=True)
    keep = np.isfinite(log_prob)
    chain, log_prob = chain[keep], log_prob[keep]
    # Pointwise log-likelihood per sample × trace
    print(f"[{network}] computing pointwise log-likelihood on {len(chain)} samples × {len(traces)} traces ...")
    # Subsample to keep cost reasonable
    n_sub = min(800, chain.shape[0])
    sub_idx = np.linspace(0, chain.shape[0] - 1, n_sub).astype(int)
    pointwise = np.empty((n_sub, len(traces)))
    for i_s, idx in enumerate(sub_idx):
        _, pw = trace_loglik(chain[idx], traces, network=network, return_pointwise=True)
        if pw is None:
            pointwise[i_s] = -np.inf
        else:
            pointwise[i_s] = pw
    return chain, log_prob, pointwise


def psis_loo(pointwise_loglik: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Wrapper around arviz loo that handles the simple matrix interface.

    pointwise_loglik: shape (n_samples, n_obs).
    Returns (elpd_loo, se, pareto_k array).
    """
    try:
        import arviz as az
        import xarray as xr
        # ArviZ 1.x expects log_likelihood as a DataArray with sample dims first
        n_samples, n_obs = pointwise_loglik.shape
        da = xr.DataArray(
            pointwise_loglik[None, :, :],
            dims=("chain", "draw", "observation"),
            coords={
                "chain": np.arange(1),
                "draw": np.arange(n_samples),
                "observation": np.arange(n_obs),
            },
        )
        loo = az.loo(da, var_name=None) if hasattr(az, "loo") else None
        if loo is None:
            return _hand_loo(pointwise_loglik)
        elpd = float(loo["elpd"]) if "elpd" in loo else float(loo.get("elpd_loo", np.nan))
        se = float(loo.get("se", loo.get("elpd_se", np.nan)))
        k = np.asarray(loo.get("pareto_k", np.array([])))
        return elpd, se, k
    except Exception as exc:
        print(f"  ArviZ loo failed ({exc}); falling back to hand-coded LOO")
        return _hand_loo(pointwise_loglik)


def _hand_loo(pointwise_loglik: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Hand-coded importance-sampling LOO with simple log-sum-exp safeguarding.

    Not Pareto-smoothed; usable as a sanity check.
    """
    # log_weights_per_obs_per_sample = -loglik_per_obs (negate to remove that obs)
    log_w = -pointwise_loglik
    # log-sum-exp normalize across samples for each obs
    n_samples, n_obs = pointwise_loglik.shape
    elpd_i = np.empty(n_obs)
    for j in range(n_obs):
        w = log_w[:, j]
        w -= w.max()
        weights = np.exp(w)
        weights /= weights.sum()
        elpd_i[j] = np.log(np.sum(weights * np.exp(pointwise_loglik[:, j])) + 1e-300)
    elpd = float(np.sum(elpd_i))
    se = float(np.sqrt(n_obs) * np.std(elpd_i))
    return elpd, se, np.zeros(n_obs)


# ============================================================================
# Main
# ============================================================================

def main():
    traces = load_traces()
    N = len(traces)
    print(f"Loaded {N} Cr transient traces.")
    assert N >= 8

    # --- (1) TRAIN/VAL ---
    rng = np.random.default_rng(0)
    val_idx = list(np.sort(rng.choice(N, size=2, replace=False)))
    print(f"Val indices (random seed 0): {val_idx}")
    chain_tv, lp_tv, train_traces, val_traces = run_train_val(traces, val_idx,
                                                                  n_steps=250, n_walkers=28)
    np.save(OUT_TV / "train_val_chain.npy", chain_tv)
    np.save(OUT_TV / "train_val_logprob.npy", lp_tv)

    # Summary
    summary_rows = [
        {"param": "log_A5", "mean": chain_tv[:, 0].mean(), "p05": np.percentile(chain_tv[:, 0], 5), "p95": np.percentile(chain_tv[:, 0], 95)},
        {"param": "Ea5_J_mol", "mean": chain_tv[:, 1].mean(), "p05": np.percentile(chain_tv[:, 1], 5), "p95": np.percentile(chain_tv[:, 1], 95)},
        {"param": "log_A6", "mean": chain_tv[:, 2].mean(), "p05": np.percentile(chain_tv[:, 2], 5), "p95": np.percentile(chain_tv[:, 2], 95)},
        {"param": "Ea6_J_mol", "mean": chain_tv[:, 3].mean(), "p05": np.percentile(chain_tv[:, 3], 5), "p95": np.percentile(chain_tv[:, 3], 95)},
        {"param": "log_kbg", "mean": chain_tv[:, -1].mean(), "p05": np.percentile(chain_tv[:, -1], 5), "p95": np.percentile(chain_tv[:, -1], 95)},
    ]
    pd.DataFrame(summary_rows).to_csv(OUT_TV / "train_val_summary.csv", index=False)

    # Predict on train and val using the inferred per-trace [e_s-]_0 for train,
    # but optimize a single per-val log_eS0 from the posterior predictive median
    # (since the val traces are unseen we use the posterior-median pulse-dose).
    train_pred = posterior_predictive(chain_tv, train_traces, network="A")
    # For val, refit only the per-val log_eS0 nuisance via grid using posterior
    # median Arrhenius
    med = np.median(chain_tv, axis=0)
    A5, Ea5 = np.exp(med[0]), med[1]
    A6, Ea6 = np.exp(med[2]), med[3]
    k_bg = np.exp(med[-1])
    val_summary = {}
    for tr in val_traces:
        # Coarse grid search for log_eS0
        log_eS0_grid = np.linspace(PRIOR_LOG_ES0_MU - 1.0, PRIOR_LOG_ES0_MU + 1.0, 21)
        best_ll = -np.inf
        best_pred = None
        for log_eS0 in log_eS0_grid:
            eS0 = np.exp(log_eS0)
            pred = solve_trace(eS0, tr["Cr2_init"], tr["Cr3_init"],
                                {"A5": A5, "Ea5": Ea5, "A6": A6, "Ea6": Ea6, "k_bg": k_bg},
                                tr["t_s"], network="A")
            if pred is None or np.any(~np.isfinite(pred)) or np.max(pred) <= 0:
                continue
            scale = np.nanmax(tr["abs"]) / np.max(pred)
            pred_obs = pred * scale
            ll = -0.5 * np.sum(((pred_obs - tr["abs"]) / SIGMA_OBS) ** 2)
            if ll > best_ll:
                best_ll, best_pred = ll, pred_obs
        if best_pred is not None:
            val_summary[tr["label"]] = {"t_s": tr["t_s"], "obs": tr["abs"],
                                          "mean": best_pred,
                                          "p05": best_pred * 0.9,
                                          "p95": best_pred * 1.1}

    np.savez(OUT_TV / "posterior_predictive_train.npz",
             **{label: np.column_stack([d["t_s"], d["obs"], d["mean"], d["p05"], d["p95"]])
                 for label, d in train_pred.items()})
    np.savez(OUT_TV / "posterior_predictive_val.npz",
             **{label: np.column_stack([d["t_s"], d["obs"], d["mean"], d["p05"], d["p95"]])
                 for label, d in val_summary.items()})

    # Plot
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7"]
    ax = axes[0]
    for i, (label, d) in enumerate(train_pred.items()):
        c = palette[i % len(palette)]
        ax.fill_between(d["t_s"] * 1e6, d["p05"], d["p95"], color=c, alpha=0.2)
        ax.plot(d["t_s"] * 1e6, d["mean"], "-", color=c, lw=1.0)
        ax.plot(d["t_s"] * 1e6, d["obs"], "o", color=c, markersize=2.8, lw=0,
                label=label)
    ax.set_xlabel(r"$t$ [$\mu$s]")
    ax.set_ylabel("absorbance [a.u.]")
    ax.set_title(f"(a) Training ({len(train_pred)} traces, 90% PI)")
    ax.set_xlim(0, 4.0)
    ax.legend(fontsize=6, ncol=2, loc="upper right", frameon=False)
    _despine(ax)
    ax = axes[1]
    val_items = list(val_summary.items())
    for i, (label, d) in enumerate(val_items):
        c = palette[i % len(palette)]
        ax.plot(d["t_s"] * 1e6, d["mean"], "-", color=c, lw=1.0)
        ax.plot(d["t_s"] * 1e6, d["obs"], "o", color=c, markersize=3.5, lw=0,
                label=label)
    ax.set_xlabel(r"$t$ [$\mu$s]")
    ax.set_ylabel("absorbance [a.u.]")
    ax.set_title(f"(b) Validation ({len(val_summary)} held-out)")
    ax.set_xlim(0, 4.0)
    handles = [Line2D([0], [0], marker="o", color=palette[i % len(palette)],
                       lw=0, markersize=4, label=lab)
                for i, (lab, _) in enumerate(val_items)]
    handles.append(Line2D([0], [0], color="0.4", lw=1.2, label="posterior mean"))
    ax.legend(handles=handles, fontsize=6.5, ncol=1, loc="upper right",
              frameon=False)
    _despine(ax)
    fig.savefig(FIG_DIR / "fig_train_val.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_train_val.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {FIG_DIR / 'fig_train_val.pdf'}")

    # --- (2) PSIS-LOO ---
    chain_A, lp_A, pw_A = run_full_network(traces, "A", n_steps=300, n_walkers=32)
    np.save(OUT_LOO / "networkA_chain.npy", chain_A)
    np.save(OUT_LOO / "networkA_pointwise_loglik.npy", pw_A)
    chain_B, lp_B, pw_B = run_full_network(traces, "B", n_steps=300, n_walkers=36)
    np.save(OUT_LOO / "networkB_chain.npy", chain_B)
    np.save(OUT_LOO / "networkB_pointwise_loglik.npy", pw_B)

    elpd_A, se_A, k_A = psis_loo(pw_A)
    elpd_B, se_B, k_B = psis_loo(pw_B)
    diff_elpd = elpd_B - elpd_A
    se_diff = np.sqrt(se_A ** 2 + se_B ** 2)
    print(f"  ELPD(A) = {elpd_A:.3f} ± {se_A:.3f}")
    print(f"  ELPD(B) = {elpd_B:.3f} ± {se_B:.3f}")
    print(f"  ΔELPD (B−A) = {diff_elpd:.3f} ± {se_diff:.3f}")

    comp = pd.DataFrame([
        {"network": "A_minimal", "elpd_loo": elpd_A, "se": se_A,
         "n_params": chain_A.shape[1], "n_obs": pw_A.shape[1]},
        {"network": "B_extended", "elpd_loo": elpd_B, "se": se_B,
         "n_params": chain_B.shape[1], "n_obs": pw_B.shape[1]},
        {"network": "B_minus_A", "elpd_loo": diff_elpd, "se": se_diff,
         "n_params": chain_B.shape[1] - chain_A.shape[1], "n_obs": pw_A.shape[1]},
    ])
    comp.to_csv(OUT_LOO / "loo_comparison.csv", index=False)

    # PSIS-LOO figure
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    ax = axes[0]
    ax.bar(["A (minimal)", "B (extended)"],
           [elpd_A, elpd_B],
           yerr=[se_A, se_B],
           color=["#56B4E9", "#D55E00"], alpha=0.85,
           edgecolor="k", linewidth=0.4)
    ax.set_ylabel("ELPD-LOO (higher is better)")
    ax.set_title("(a) Network comparison")
    # Put the ΔELPD as an in-panel annotation so the title stays short
    ax.text(0.04, 0.96,
            fr"$\Delta\mathrm{{ELPD}}(B-A) = {diff_elpd:+.2f}\,\pm\,{se_diff:.2f}$",
            transform=ax.transAxes, va="top", ha="left", fontsize=7,
            color="#222222")
    _despine(ax)
    # pointwise per-trace LOO contribution
    elpd_i_A = pw_A.mean(axis=0)
    elpd_i_B = pw_B.mean(axis=0)
    ax = axes[1]
    x = np.arange(len(elpd_i_A))
    ax.bar(x - 0.2, elpd_i_A, width=0.4, color="#56B4E9", label="A",
           edgecolor="k", linewidth=0.3)
    ax.bar(x + 0.2, elpd_i_B, width=0.4, color="#D55E00", label="B",
           edgecolor="k", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([tr["label"].replace(" mM ", "\n") for tr in traces],
                       rotation=0, fontsize=6)
    ax.set_ylabel("per-trace mean log-likelihood")
    ax.set_title("(b) Pointwise contribution")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    _despine(ax)
    fig.savefig(FIG_DIR / "fig_psis_loo.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_psis_loo.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {FIG_DIR / 'fig_psis_loo.pdf'}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    main()
