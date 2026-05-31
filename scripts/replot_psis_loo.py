#!/usr/bin/env python3
"""Re-render fig_psis_loo.pdf with publication style using saved pointwise log-likelihoods."""
from __future__ import annotations

import sys
import warnings
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

OUT_LOO = REPO / "validation" / "cr_psis_loo"
FIG_DIR = REPO / "manuscript" / "figures"


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _hand_loo(pw):
    log_w = -pw
    n_samples, n_obs = pw.shape
    elpd_i = np.empty(n_obs)
    for j in range(n_obs):
        w = log_w[:, j].copy()
        w -= w.max()
        weights = np.exp(w)
        weights /= weights.sum()
        elpd_i[j] = np.log(np.sum(weights * np.exp(pw[:, j])) + 1e-300)
    elpd = float(np.sum(elpd_i))
    se = float(np.sqrt(n_obs) * np.std(elpd_i))
    return elpd, se


def main():
    pw_A = np.load(OUT_LOO / "networkA_pointwise_loglik.npy")
    pw_B = np.load(OUT_LOO / "networkB_pointwise_loglik.npy")
    elpd_A, se_A = _hand_loo(pw_A)
    elpd_B, se_B = _hand_loo(pw_B)
    diff_elpd = elpd_B - elpd_A
    se_diff = np.sqrt(se_A ** 2 + se_B ** 2)

    # Try to use the values already in the comparison CSV (preserves figure paper text)
    try:
        comp = pd.read_csv(OUT_LOO / "loo_comparison.csv")
        elpd_A_csv = float(comp[comp["network"] == "A_minimal"]["elpd_loo"].iloc[0])
        elpd_B_csv = float(comp[comp["network"] == "B_extended"]["elpd_loo"].iloc[0])
        se_A_csv = float(comp[comp["network"] == "A_minimal"]["se"].iloc[0])
        se_B_csv = float(comp[comp["network"] == "B_extended"]["se"].iloc[0])
        diff_csv = float(comp[comp["network"] == "B_minus_A"]["elpd_loo"].iloc[0])
        se_diff_csv = float(comp[comp["network"] == "B_minus_A"]["se"].iloc[0])
        elpd_A, elpd_B, se_A, se_B = elpd_A_csv, elpd_B_csv, se_A_csv, se_B_csv
        diff_elpd, se_diff = diff_csv, se_diff_csv
    except Exception:
        pass

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    ax = axes[0]
    ax.bar(["A (minimal)", "B (extended)"],
           [elpd_A, elpd_B], yerr=[se_A, se_B],
           color=["#56B4E9", "#D55E00"], alpha=0.85,
           edgecolor="k", linewidth=0.4)
    ax.set_ylabel("ELPD-LOO (higher is better)")
    ax.set_title("(a) Network comparison")
    ax.text(0.04, 0.96,
            fr"$\Delta\mathrm{{ELPD}}(B-A) = {diff_elpd:+.2f}\,\pm\,{se_diff:.2f}$",
            transform=ax.transAxes, va="top", ha="left", fontsize=7,
            color="#222222")
    _despine(ax)

    elpd_i_A = pw_A.mean(axis=0)
    elpd_i_B = pw_B.mean(axis=0)
    ax = axes[1]
    x = np.arange(len(elpd_i_A))
    ax.bar(x - 0.2, elpd_i_A, width=0.4, color="#56B4E9", label="A",
           edgecolor="k", linewidth=0.3)
    ax.bar(x + 0.2, elpd_i_B, width=0.4, color="#D55E00", label="B",
           edgecolor="k", linewidth=0.3)
    ax.set_xticks(x)
    try:
        from tier2_bayesian_calibration import load_traces
        traces = load_traces()
        labels = [tr["label"].replace(" mM ", "\n") for tr in traces[:len(elpd_i_A)]]
    except Exception:
        labels = [f"trace {i+1}" for i in range(len(elpd_i_A))]
    ax.set_xticklabels(labels, rotation=0, fontsize=6)
    ax.set_ylabel("per-trace mean log-likelihood")
    ax.set_title("(b) Pointwise contribution")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    _despine(ax)

    fig.savefig(FIG_DIR / "fig_psis_loo.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_psis_loo.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {FIG_DIR / 'fig_psis_loo.pdf'}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
