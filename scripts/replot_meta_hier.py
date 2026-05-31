#!/usr/bin/env python3
"""Re-render fig_meta_hier.pdf using the saved meta-hierarchical chain + LOMO results.

This avoids re-running the MCMC; it just rebuilds the 4-panel figure with the
publication style.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "manuscript" / "figures"))
sys.path.insert(0, str(REPO / "scripts"))
try:
    import figstyle  # type: ignore
    figstyle.apply_rc()
except Exception:
    pass

# Bring back the indexing logic from the upstream module by importing it.
from meta_hierarchical_layer import (
    build_dataset, build_chem_features, make_logpost, CHEM_FEATURES,
)

OUT = REPO / "validation" / "meta_hier"
FIG_DIR = REPO / "manuscript" / "figures"


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    chain = np.load(OUT / "full_posterior_chain.npy")
    lomo_df = pd.read_csv(OUT / "lomo_results.csv")
    df = build_dataset()

    metals_all = sorted(df["metal"].unique())
    hosts_all = sorted(df["host"].unique())
    X_all, feat_names = build_chem_features(metals_all)
    _, _, _, meta = make_logpost(df, metals_all, hosts_all, X_all, feat_names)
    indices = meta["indices"]
    metals_sorted = metals_all

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), constrained_layout=True)
    palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999", "#882255"]

    # (a) raw data scatter
    ax = axes[0, 0]
    host_marker = {"LiCl-KCl": "o", "NaCl": "s", "KCl": "^", "KBr": "D",
                    "NaBr": "P", "KI": "X", "H2O": "*"}
    for i, m in enumerate(metals_sorted):
        sub = df[df["metal"] == m]
        for host_, grp in sub.groupby("host"):
            ax.errorbar(1000.0 / grp["T_K"], grp["log10_k"], yerr=grp["sigma_log10_k"],
                         fmt=host_marker.get(host_, "o"),
                         color=palette[i % len(palette)], alpha=0.85,
                         markersize=4.5, lw=0.7,
                         label=f"{m} / {host_}" if host_ == "LiCl-KCl" or m == "Cf3+" else None)
    ax.set_xlabel(r"$1000/T$ [K$^{-1}$]")
    ax.set_ylabel(r"$\log_{10} k$ [M$^{-1}$ s$^{-1}$]")
    ax.set_title("(a) Multi-metal data")
    ax.legend(fontsize=6, ncol=2, loc="lower left", frameon=False)
    _despine(ax)

    # (b) Chemistry-feature loadings
    ax = axes[0, 1]
    betaA_samp = chain[:, indices["betaA"]]
    betaE_samp = chain[:, indices["betaE"]]
    n_feat = betaA_samp.shape[1]
    y = np.arange(n_feat)
    means_A = betaA_samp.mean(axis=0)
    err_A = np.array([means_A - np.percentile(betaA_samp, 5, axis=0),
                       np.percentile(betaA_samp, 95, axis=0) - means_A])
    means_E = betaE_samp.mean(axis=0)
    err_E = np.array([means_E - np.percentile(betaE_samp, 5, axis=0),
                       np.percentile(betaE_samp, 95, axis=0) - means_E])
    ax.errorbar(means_A, y - 0.15, xerr=err_A, fmt="o", color="#0072B2",
                markersize=5, label=r"$\log_{10} A$")
    ax.errorbar(means_E, y + 0.15, xerr=err_E, fmt="s", color="#D55E00",
                markersize=5, label=r"$E_a$ [kJ mol$^{-1}$]")
    ax.axvline(0, color="#bbbbbb", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(meta["feature_names"], fontsize=7)
    ax.set_xlabel("posterior coefficient (mean, 90% CI)")
    ax.set_title("(b) Feature loadings")
    ax.legend(fontsize=7, frameon=False)
    _despine(ax)

    # (c) LOMO scatter
    ax = axes[1, 0]
    if len(lomo_df) > 0:
        cmap = {m: palette[i % len(palette)] for i, m in enumerate(metals_sorted)}
        for m, grp in lomo_df.groupby("metal"):
            yerr_pred = np.array([
                grp["log10_k_pred_mean"] - grp["log10_k_pred_p05"],
                grp["log10_k_pred_p95"] - grp["log10_k_pred_mean"],
            ])
            ax.errorbar(grp["log10_k_obs"], grp["log10_k_pred_mean"],
                         yerr=yerr_pred, xerr=grp["sigma_log10_k_obs"],
                         fmt="o", color=cmap.get(m, "k"), alpha=0.85,
                         markersize=5, lw=0.7, label=m)
        lo = min(lomo_df["log10_k_obs"].min(), lomo_df["log10_k_pred_p05"].min())
        hi = max(lomo_df["log10_k_obs"].max(), lomo_df["log10_k_pred_p95"].max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.7)
    ax.set_xlabel(r"observed $\log_{10} k$")
    ax.set_ylabel(r"LOMO predicted $\log_{10} k$")
    ax.set_title("(c) LOMO predictions")
    ax.legend(fontsize=6, ncol=2, frameon=False)
    _despine(ax)

    # (d) Host effects
    ax = axes[1, 1]
    if len(meta["non_ref_hosts"]) > 0:
        b_samp = chain[:, indices["B"]]
        means = b_samp.mean(axis=0)
        err = np.array([means - np.percentile(b_samp, 5, axis=0),
                         np.percentile(b_samp, 95, axis=0) - means])
        y = np.arange(len(meta["non_ref_hosts"]))
        ax.errorbar(means, y, xerr=err, fmt="s", color="#009E73", markersize=6)
        ax.axvline(0, color="#bbbbbb", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(meta["non_ref_hosts"], fontsize=7)
        ax.set_xlabel(r"$b_h$ on $\log_{10} k$ (ref " + str(meta["ref_host"]) + ")")
        ax.set_title("(d) Host effects")
        _despine(ax)

    fig.savefig(FIG_DIR / "fig_meta_hier.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_meta_hier.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {FIG_DIR / 'fig_meta_hier.pdf'}")


if __name__ == "__main__":
    main()
