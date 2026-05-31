#!/usr/bin/env python3
"""Figure: HBMAE vs four state-of-the-art comparator methods on the multi-paper Zn case.

Produces a 3-panel double-column figure for the Comparison section of the article.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from manuscript.figures.figstyle import (
    apply_rc, double_col, savefig, despine, WONG,
)

LIT_A_LOG = np.log10(2.4e13)
LIT_EA_KJ = 35.6


def load_summary():
    rows = {}
    with (REPO / "validation/TIER3_METHOD_COMPARISON.csv").open() as f:
        for row in csv.DictReader(f):
            rows[row["method"][:2]] = row
    return rows


def main():
    apply_rc()
    summary = load_summary()
    keys = ["M0", "M1", "M2", "M3", "M4"]
    labels = [r"M$_0$", r"M$_1$", r"M$_2$", r"M$_3$", r"M$_4$"]

    elpd = [float(summary[k]["elpd_WAIC"]) for k in keys]
    held = [float(summary[k]["held_out_lpd"]) for k in keys]
    log_A_bias = [np.log10(float(summary[k]["median_A"])) - LIT_A_LOG for k in keys]
    Ea_bias = [float(summary[k]["median_Ea_kJ"]) - LIT_EA_KJ for k in keys]

    # Use a consistent Wong subset
    colors = [WONG[2], WONG[6], WONG[3], WONG[5], WONG[1]]

    fig, axes = plt.subplots(1, 3, figsize=(double_col(), 2.5),
                              constrained_layout=True)
    x = np.arange(len(keys))

    # Panel (a) in-sample WAIC
    ax = axes[0]
    ax.bar(x, elpd, color=colors, width=0.7, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(r"$\mathrm{elpd}_{\mathrm{WAIC}}$ [nats]")
    ax.set_title("(a) In-sample elpd")
    ax.axhline(0, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="x", length=0)
    despine(ax)

    # Panel (b) held-out predictive
    ax = axes[1]
    ax.bar(x, held, color=colors, width=0.7, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(r"$\log p(y_{\mathrm{held}}\mid \mathcal{D}, M)$")
    ax.set_title("(b) Held-out elpd")
    ax.axhline(0, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="x", length=0)
    despine(ax)

    # Panel (c) bias against literature (twin axes for two scales)
    ax = axes[2]
    w = 0.35
    ax.bar(x - w/2, log_A_bias, w, color=colors, edgecolor="black",
           linewidth=0.4, alpha=0.85)
    ax.set_ylabel(r"$\log_{10}(A_{\mathrm{post}} / A_{\mathrm{lit}})$")
    ax2 = ax.twinx()
    ax2.bar(x + w/2, Ea_bias, w, color=colors, edgecolor="black",
            linewidth=0.4, alpha=0.5, hatch="///")
    ax2.set_ylabel(r"$E_{a,\mathrm{post}}-E_{a,\mathrm{lit}}$ [kJ mol$^{-1}$]")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_title("(c) Bias vs literature")
    ax.axhline(0, color="#bbbbbb", lw=0.5)
    ax.tick_params(axis="x", length=0)
    # On a twin-axis panel keep only the top spine off; right spine is needed
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    from matplotlib.patches import Patch
    legend_patches = [
        Patch(facecolor="0.85", edgecolor="black", linewidth=0.4,
              label=r"$\log_{10}(A_{\mathrm{post}}/A_{\mathrm{lit}})$"),
        Patch(facecolor="0.85", edgecolor="black", linewidth=0.4,
              hatch="///", alpha=0.5,
              label=r"$\Delta E_{a}$"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=6.5,
              frameon=False)

    savefig(fig, REPO / "manuscript/figures/fig_method_comparison.pdf")
    print("Wrote manuscript/figures/fig_method_comparison.{pdf,png}")


if __name__ == "__main__":
    main()
