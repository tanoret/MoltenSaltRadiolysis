#!/usr/bin/env python3
"""Posterior visualization for Tier 2 calibration.

Loads the saved chain from scripts/tier2_bayesian_calibration.py and produces:
  (a) corner plot of the 4 Arrhenius parameters
  (b) posterior predictive overlay on the 9 transient traces
  (c) trace plots for chain convergence diagnostics
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import corner

from scripts.tier2_bayesian_calibration import load_traces, solve_trace


def main():
    out = REPO / "validation"
    chain = np.load(out / "tier2_chain.npy")          # shape (n_steps, n_walkers, ndim)
    log_prob = np.load(out / "tier2_log_prob.npy")
    print(f"Loaded chain: shape {chain.shape}, log_prob shape {log_prob.shape}")
    n_steps, n_walkers, ndim = chain.shape

    flat = chain.reshape(-1, ndim)

    # ---------- (a) corner of Arrhenius -----------
    arrh_idx = [0, 1, 2, 3]
    arrh_labels = ["log A5", "Ea5 (J/mol)", "log A6", "Ea6 (J/mol)"]
    truths = [np.log(1.7e13), 33.5e3, np.log(2.0e13), 31.8e3]
    fig = corner.corner(flat[:, arrh_idx], labels=arrh_labels, truths=truths,
                        truth_color="red", show_titles=True, quantiles=[0.025, 0.5, 0.975])
    fig.suptitle("HBMAE Tier 2 posterior — Arrhenius parameters (red = literature)", y=1.02, fontsize=12)
    fig.savefig(out / "tier2_corner_arrhenius.png", bbox_inches="tight", dpi=120)
    plt.close(fig)
    print(f"Wrote {(out / 'tier2_corner_arrhenius.png').relative_to(REPO)}")

    # ---------- (b) posterior predictive overlay -----------
    traces = load_traces()
    N = len(traces)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    rng = np.random.default_rng(42)
    n_pp = 200
    pp_idx = rng.integers(0, flat.shape[0], size=n_pp)
    for tr_idx, tr in enumerate(traces):
        ax = axes[0] if "Cr(II)" in tr["label"] else axes[1]
        # Plot experimental
        ax.plot(tr["t_s"] * 1e9, tr["abs"], "o", ms=3, alpha=0.4, label=f"exp {tr['label']}")
        # Posterior predictive band: 200 forward solves on subsampled posterior
        t_eval = np.linspace(0, 25e-9, 200)
        pp_curves = []
        for pi in pp_idx[:50]:   # 50 PP draws per trace for speed
            sample = flat[pi]
            log_A5, Ea5, log_A6, Ea6 = sample[0:4]
            log_eS0 = sample[4 + tr_idx]
            log_kbg = sample[-1]
            eS = solve_trace(np.exp(log_eS0), tr["Cr2_init"], tr["Cr3_init"],
                              np.exp(log_A5), Ea5, np.exp(log_A6), Ea6,
                              np.exp(log_kbg), t_eval)
            if eS is not None:
                scale = np.nanmax(tr["abs"]) / max(np.max(eS), 1e-12)
                pp_curves.append(eS * scale)
        if pp_curves:
            pp_arr = np.array(pp_curves)
            med = np.median(pp_arr, axis=0)
            lo, hi = np.percentile(pp_arr, [2.5, 97.5], axis=0)
            ax.plot(t_eval * 1e9, med, "-", lw=1.2, label=f"model {tr['label']}")
            ax.fill_between(t_eval * 1e9, lo, hi, alpha=0.15)
    for ax, title in zip(axes, ["Cr(II) initial", "Cr(III) initial"]):
        ax.set_xlabel("Time [ns]")
        ax.set_title(title)
        ax.set_xlim(0, 25)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=7, loc="upper right")
    axes[0].set_ylabel("Absorbance (model rescaled to exp amplitude)")
    plt.suptitle("HBMAE Tier 2 posterior predictive overlay (95% credible band)")
    plt.tight_layout()
    fig.savefig(out / "tier2_posterior_predictive.png", dpi=120)
    plt.close(fig)
    print(f"Wrote {(out / 'tier2_posterior_predictive.png').relative_to(REPO)}")

    # ---------- (c) trace plots for the 4 Arrhenius parameters -----------
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    for i, (idx, lab, truth) in enumerate(zip(arrh_idx, arrh_labels, truths)):
        for w in range(min(n_walkers, 16)):
            axes[i].plot(chain[:, w, idx], "k-", alpha=0.25, lw=0.5)
        axes[i].axhline(truth, color="red", lw=1.0, label=f"literature {truth:.3g}")
        axes[i].set_ylabel(lab)
        axes[i].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("MCMC step")
    plt.suptitle("HBMAE Tier 2 chain trace (16 walkers shown)")
    plt.tight_layout()
    fig.savefig(out / "tier2_chain_trace.png", dpi=120)
    plt.close(fig)
    print(f"Wrote {(out / 'tier2_chain_trace.png').relative_to(REPO)}")


if __name__ == "__main__":
    main()
