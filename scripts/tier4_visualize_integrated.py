#!/usr/bin/env python3
"""Visualize the integrated HBMAE posterior from scripts/tier4_integrated_hbmae.py."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import corner

REPO = Path(__file__).resolve().parent.parent

LABELS = [
    "log A5", "Ea5", "log A6", "Ea6",
    *[f"log[eS]0_{i+1}" for i in range(9)],
    "log k_bg",
    "log A_Zn", "Ea_Zn",
    "η_logA_LiCl-KCl", "η_Ea_LiCl-KCl",
    "η_logA_NaCl",     "η_Ea_NaCl",
    "η_logA_KCl",      "η_Ea_KCl",
    "b_Pikaev",
    "log G(Cl•)"
]

LITERATURE = {
    "log A5": np.log(1.7e13),
    "Ea5":    33.5e3,
    "log A6": np.log(2.0e13),
    "Ea6":    31.8e3,
    "log k_bg": np.log(1e7),
    "log A_Zn": np.log(2.4e13),
    "Ea_Zn":    35.6e3,
    "b_Pikaev": -4.89,    # from Tier-3 M2/M3
}


def main():
    chain = np.load(REPO / "validation/tier4_integrated_chain.npy")
    print(f"Loaded integrated chain: shape {chain.shape}")
    n_samples, ndim = chain.shape

    # ---------- Corner plot for the "core 6" Arrhenius parameters ----------
    core_indices = [0, 1, 2, 3, 14, 15]  # Cr5, Cr6, Zn intrinsic
    core_labels = [LABELS[i] for i in core_indices]
    core_truths = [LITERATURE.get(LABELS[i]) for i in core_indices]
    fig = corner.corner(chain[:, core_indices], labels=core_labels,
                        truths=core_truths, truth_color="red",
                        quantiles=[0.025, 0.5, 0.975], show_titles=True)
    fig.suptitle("HBMAE integrated MCMC — Arrhenius parameters", y=1.02, fontsize=12)
    fig.savefig(REPO / "validation/tier4_integrated_corner_arrhenius.png",
                bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("Wrote tier4_integrated_corner_arrhenius.png")

    # ---------- Salt perturbation and facility/Phillips parameters ----------
    side_indices = [16, 17, 18, 19, 20, 21, 22, 23]
    side_labels = [LABELS[i] for i in side_indices]
    fig = corner.corner(chain[:, side_indices], labels=side_labels,
                        quantiles=[0.025, 0.5, 0.975], show_titles=True)
    fig.suptitle("HBMAE integrated — salt perturbations, facility, Phillips G(Cl•)",
                 y=1.02, fontsize=11)
    fig.savefig(REPO / "validation/tier4_integrated_corner_hierarchy.png",
                bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("Wrote tier4_integrated_corner_hierarchy.png")

    # ---------- Summary marginals as bar plot ----------
    fig, ax = plt.subplots(figsize=(12, 6))
    medians = np.median(chain, axis=0)
    lo = np.percentile(chain, 2.5, axis=0)
    hi = np.percentile(chain, 97.5, axis=0)
    err_lo = medians - lo
    err_hi = hi - medians

    # Normalize each parameter to its prior centre for visualization
    # Just plot raw values for now
    x = np.arange(ndim)
    ax.errorbar(x, medians, yerr=[err_lo, err_hi], fmt="o", color="C4",
                ms=5, capsize=3, label="posterior 95% CI")
    # Mark literature
    for i, lab in enumerate(LABELS):
        if lab in LITERATURE:
            ax.scatter([i], [LITERATURE[lab]], marker="x", s=80, color="red", zorder=10)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Parameter value")
    ax.set_title("HBMAE integrated posterior medians + 95% CI (red x = literature)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    fig.savefig(REPO / "validation/tier4_integrated_marginals.png", dpi=110)
    plt.close(fig)
    print("Wrote tier4_integrated_marginals.png")


if __name__ == "__main__":
    main()
