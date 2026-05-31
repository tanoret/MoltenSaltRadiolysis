#!/usr/bin/env python3
"""Visualize Tier 3 method comparison.

Produces a side-by-side comparison plot of the (log A, Ea) posteriors from the
five methods, plus a bar chart of WAIC and held-out predictive accuracy.
"""

from __future__ import annotations

import sys
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent

METHODS = [
    ("M0", "M0 — Iwamatsu only",      "C0"),
    ("M1", "M1 — Naive pool",          "C3"),
    ("M2", "M2 — Facility effect",    "C2"),
    ("M3", "M3 — HBMAE full",          "C4"),
    ("M4", "M4 — GM-style weak prior", "C1"),
]
LIT_A, LIT_EA_KJ = 2.4e13, 35.6


def main():
    # ---------- Joint posterior in (log10 A, Ea) ----------
    fig, ax = plt.subplots(figsize=(10, 6))
    for key, label, color in METHODS:
        path = REPO / f"validation/tier3_{key}_chain.npy"
        if not path.exists():
            continue
        chain = np.load(path)
        log_A_samples = chain[:, 0]
        Ea_samples = chain[:, 1] / 1e3
        # Thin
        idx = np.random.default_rng(0).integers(0, chain.shape[0], size=min(800, chain.shape[0]))
        ax.scatter(log_A_samples[idx] / np.log(10), Ea_samples[idx],
                   s=3, alpha=0.35, color=color, label=label)
    ax.axhline(LIT_EA_KJ, color="black", lw=0.7, ls=":", label=f"Lit. Ea = {LIT_EA_KJ} kJ/mol")
    ax.axvline(np.log10(LIT_A), color="black", lw=0.7, ls="--", label=f"Lit. log10 A = {np.log10(LIT_A):.2f}")
    ax.set_xlabel("log10 A   (M⁻¹ s⁻¹)")
    ax.set_ylabel("Eₐ (kJ/mol)")
    ax.set_title("Tier 3: posterior on intrinsic Arrhenius across 5 inference methods")
    ax.legend(loc="lower right", markerscale=4)
    ax.grid(True, alpha=0.3)
    out1 = REPO / "validation/tier3_posterior_comparison.png"
    plt.tight_layout()
    fig.savefig(out1, dpi=120)
    plt.close(fig)
    print(f"Wrote {out1.relative_to(REPO)}")

    # ---------- Metrics bar chart ----------
    summary = {}
    with (REPO / "validation/TIER3_METHOD_COMPARISON.csv").open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            summary[row["method"][:2]] = row

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    keys = ["M0", "M1", "M2", "M3", "M4"]
    labels = ["M0\nIwamatsu", "M1\nNaive\npool", "M2\nFacility", "M3\nHBMAE\nfull", "M4\nGM weak\nprior"]
    colors = ["C0", "C3", "C2", "C4", "C1"]

    # 1. elpd_WAIC
    elpd = [float(summary[k]["elpd_WAIC"]) for k in keys]
    axes[0].bar(labels, elpd, color=colors)
    axes[0].set_ylabel("elpd_WAIC (higher = better predictive)")
    axes[0].set_title("In-sample predictive accuracy")
    axes[0].axhline(0, color="black", lw=0.5)
    axes[0].grid(True, alpha=0.3)

    # 2. held-out lpd
    held = [float(summary[k]["held_out_lpd"]) for k in keys]
    axes[1].bar(labels, held, color=colors)
    axes[1].set_ylabel("log p(y_held | data, M)")
    axes[1].set_title("Out-of-sample predictive (Iwamatsu 550°C held out)")
    axes[1].axhline(0, color="black", lw=0.5)
    axes[1].grid(True, alpha=0.3)

    # 3. bias in log10(A) and Ea
    log_A_bias = [np.log10(float(summary[k]["median_A"]) / LIT_A) for k in keys]
    Ea_bias = [(float(summary[k]["median_Ea_kJ"]) - LIT_EA_KJ) for k in keys]
    x = np.arange(len(keys))
    width = 0.35
    axes[2].bar(x - width/2, log_A_bias, width, color=colors, label="log10(A_post/A_lit)", alpha=0.7)
    ax2b = axes[2].twinx()
    ax2b.bar(x + width/2, Ea_bias, width, color=colors, label="Ea_post − Ea_lit (kJ/mol)", alpha=0.4, hatch="//")
    axes[2].set_xticks(x); axes[2].set_xticklabels(labels)
    axes[2].set_ylabel("log10(A_post / A_lit)", color="black")
    ax2b.set_ylabel("Eₐ_post − Eₐ_lit (kJ/mol)", color="gray")
    axes[2].set_title("Bias vs literature values")
    axes[2].axhline(0, color="black", lw=0.5)
    axes[2].grid(True, alpha=0.3)

    plt.suptitle("HBMAE Tier 3 method comparison — Iwamatsu (LiCl-KCl) + Pikaev (NaCl, KCl) data", fontsize=11)
    plt.tight_layout()
    out2 = REPO / "validation/tier3_metrics_comparison.png"
    fig.savefig(out2, dpi=120)
    plt.close(fig)
    print(f"Wrote {out2.relative_to(REPO)}")


if __name__ == "__main__":
    main()
