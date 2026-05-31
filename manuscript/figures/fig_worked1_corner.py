#!/usr/bin/env python3
"""Figure: joint posterior on (log A_5, E_{a,5}) for Worked Example I (Cr2+ + e_s-)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
# Avoid corner's optional arviz/jax dependency on machines where jaxlib has AVX issues
sys.modules.setdefault("jax", None)
import corner

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from manuscript.figures.figstyle import apply_rc, savefig, single_col

import matplotlib.pyplot as plt


def main():
    apply_rc()
    chain = np.load(REPO / "validation/tier2_chain.npy")
    flat = chain.reshape(-1, chain.shape[-1])
    arrh = flat[:, [0, 1]]                                # (logA5, Ea5)
    arrh[:, 1] /= 1e3                                     # kJ/mol

    fig = corner.corner(
        arrh,
        labels=[r"$\log A_5$", r"$E_{a,5}$ [kJ mol$^{-1}$]"],
        truths=[np.log(1.7e13), 33.5],
        truth_color="#0072B2",
        quantiles=[0.025, 0.5, 0.975],
        show_titles=True,
        title_fmt=".3g",
        title_kwargs={"fontsize": 8},
        label_kwargs={"fontsize": 9},
        hist_kwargs={"linewidth": 1.0, "color": "k"},
        contour_kwargs={"linewidths": 0.7},
        plot_density=True,
        plot_contours=True,
    )
    fig.set_size_inches(single_col(), single_col())
    savefig(fig, REPO / "manuscript/figures/fig_worked1_corner.pdf")
    print("Wrote fig_worked1_corner.{pdf,png}")


if __name__ == "__main__":
    main()
