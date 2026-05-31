#!/usr/bin/env python3
"""Re-render fig_train_val.pdf with correct ns x-axis range.

The Cr transient time-series is on a 0-50 ns scale, not 0-4 us. This is a
post-hoc fix for the figure produced by cr_train_val_and_psis_loo.py.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "manuscript" / "figures"))
try:
    import figstyle  # type: ignore
    figstyle.apply_rc()
except Exception:
    pass

OUT_TV = REPO / "validation" / "cr_train_val"
FIG_DIR = REPO / "manuscript" / "figures"

palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7"]


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_npz_curves(path: Path) -> dict:
    """Each key in the npz contains a (n_t, 5) array: [t_s, obs, mean, p05, p95]."""
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def main():
    train = load_npz_curves(OUT_TV / "posterior_predictive_train.npz")
    val = load_npz_curves(OUT_TV / "posterior_predictive_val.npz")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)

    # (a) Training set
    ax = axes[0]
    max_t_ns = 0.0
    for i, (label, arr) in enumerate(train.items()):
        c = palette[i % len(palette)]
        t_ns = arr[:, 0] * 1e9
        obs, mean, p05, p95 = arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
        ax.fill_between(t_ns, p05, p95, color=c, alpha=0.2)
        ax.plot(t_ns, mean, "-", color=c, lw=1.0)
        ax.plot(t_ns, obs, "o", color=c, markersize=2.8, lw=0, label=label)
        if len(t_ns) > 0:
            max_t_ns = max(max_t_ns, float(t_ns.max()))
    ax.set_xlabel(r"$t$ [ns]")
    ax.set_ylabel("absorbance [a.u.]")
    ax.set_title(f"(a) Training ({len(train)} traces, 90% PI)")
    ax.set_xlim(0, max_t_ns * 1.05 if max_t_ns > 0 else 1.0)
    ax.legend(fontsize=6, ncol=2, loc="upper right", frameon=False)
    _despine(ax)

    # (b) Validation set
    ax = axes[1]
    max_t_ns = 0.0
    val_items = list(val.items())
    for i, (label, arr) in enumerate(val_items):
        c = palette[i % len(palette)]
        t_ns = arr[:, 0] * 1e9
        obs, mean = arr[:, 1], arr[:, 2]
        p05, p95 = arr[:, 3], arr[:, 4]
        ax.fill_between(t_ns, p05, p95, color=c, alpha=0.2)
        ax.plot(t_ns, mean, "-", color=c, lw=1.0)
        ax.plot(t_ns, obs, "o", color=c, markersize=3.5, lw=0, label=label)
        if len(t_ns) > 0:
            max_t_ns = max(max_t_ns, float(t_ns.max()))
    ax.set_xlabel(r"$t$ [ns]")
    ax.set_ylabel("absorbance [a.u.]")
    ax.set_title(f"(b) Validation ({len(val)} held-out)")
    ax.set_xlim(0, max_t_ns * 1.05 if max_t_ns > 0 else 1.0)
    # Build a single legend that names each trace once and adds one "model" entry
    handles = [Line2D([0], [0], marker="o", color=palette[i % len(palette)],
                       lw=0, markersize=4, label=lab)
                for i, (lab, _) in enumerate(val_items)]
    handles.append(Line2D([0], [0], color="0.4", lw=1.2, label="posterior mean"))
    ax.legend(handles=handles, fontsize=6.5, ncol=1, loc="upper right",
              frameon=False)
    _despine(ax)

    fig.savefig(FIG_DIR / "fig_train_val.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_train_val.png", dpi=300, bbox_inches="tight")
    print(f"Re-rendered {FIG_DIR / 'fig_train_val.pdf'}")


if __name__ == "__main__":
    main()
