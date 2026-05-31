"""Common helper for publication-quality figures.

Usage:
    from manuscript.figures.figstyle import setup, single_col, double_col

    fig, ax = setup(double_col(), height_in=3.5)
    ...
    fig.savefig("manuscript/figures/fig_X.pdf")
"""
from __future__ import annotations

from pathlib import Path
import matplotlib

matplotlib.use("Agg")  # interactive backend disabled (Mac display issue)
import matplotlib.pyplot as plt


STYLE_PATH = Path(__file__).resolve().parent / "mpl_style_nse.mplstyle"

# Journal column widths (inches).
SINGLE_COL_WIDTH_IN = 3.5
DOUBLE_COL_WIDTH_IN = 7.0

# Wong colorblind palette (Bang Wong, Nature Methods 2011).
WONG = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]


def apply_rc():
    """Apply the publication rcParams (also activated by plt.style.use)."""
    plt.style.use(str(STYLE_PATH))
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times", "Computer Modern Roman"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.4,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "figure.dpi": 110,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def setup(figsize=(SINGLE_COL_WIDTH_IN, 2.6), n_rows=1, n_cols=1, constrained=True):
    """Create a figure + axes preconfigured for publication."""
    apply_rc()
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize,
                              constrained_layout=constrained)
    return fig, axes


def single_col():
    return SINGLE_COL_WIDTH_IN


def double_col():
    return DOUBLE_COL_WIDTH_IN


def despine(ax, top=True, right=True):
    """Drop the top and right spines from an axes."""
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)


def panel_label(ax, label, x=0.02, y=0.97):
    """Place a panel label (e.g. 'a') in the upper-left inside the axes."""
    ax.text(x, y, label, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, fontweight="bold")


def savefig(fig, path, dpi=600):
    """Save with publication defaults. Also writes a PNG sibling at 300 dpi."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    if p.suffix.lower() == ".pdf":
        png = p.with_suffix(".png")
        fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.05)
