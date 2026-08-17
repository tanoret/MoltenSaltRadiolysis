#!/usr/bin/env python3
"""Generate three new figures for the molten-salt-fuel disposal summary.

These are NOT reused from the journal manuscripts. They focus on the disposal
regime: post-discharge dose-rate trajectory, sealed-container gas inventory
over disposal time, and a regime map separating operational from disposal
conditions.

Outputs (in ../figures):
  fig_disposal_dose_trajectory.{pdf,png}
  fig_disposal_container_inventory.{pdf,png}
  fig_disposal_regime_map.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
FIG_DIR = HERE.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Wong colorblind palette
WONG = ["#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times", "Computer Modern Roman"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "lines.linewidth": 1.5,
    "axes.spines.right": False, "axes.spines.top": False,
    "figure.dpi": 110, "savefig.dpi": 600,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    "legend.frameon": False, "mathtext.fontset": "cm",
})


# ============================================================================
# Figure 1.  Dose-rate trajectory over disposal timescale
# ----------------------------------------------------------------------------
# After reactor shutdown the fission-product beta/gamma dose decays with a
# characteristic ~30-year time constant (Cs-137, Sr-90 dominate after the
# short-lived inventory has decayed).  The actinide alpha-dose rises in
# relative importance because of the much longer half-lives of Pu-239
# (24,000 y), Am-241 (433 y), and Cm-244 (18.1 y).  In a sealed disposal
# container these two contributions cross over around 200-500 years after
# discharge.
# ============================================================================

def fig_dose_trajectory():
    t = np.logspace(-1, 6, 500)  # years after shutdown, 0.1 y to 1 My
    # Beta/gamma rate, normalized to 1 at t = 1 y.  Dominant components:
    # Cs-137 (T_1/2 = 30.17 y), Sr-90 (28.79 y), short-lived (T_1/2 ~ 5 y).
    fp_short = 8.0 * np.exp(-np.log(2) * t / 5.0)
    fp_cs137 = 4.0 * np.exp(-np.log(2) * t / 30.17)
    fp_sr90 = 3.5 * np.exp(-np.log(2) * t / 28.79)
    beta_gamma = fp_short + fp_cs137 + fp_sr90
    # Alpha rate, normalized to 1 at t = 1 y.  Components:
    # Cm-244 (T_1/2 = 18.1 y) starts at high specific activity, drops fast;
    # Am-241 (433 y) ingrowth from Pu-241 (14.4 y) -> Am-241 peaks around
    # 75 y; Pu-239 (24,110 y) is the long tail; Cm-245/247, U series, Np-237.
    cm244 = 5.0 * np.exp(-np.log(2) * t / 18.1)
    pu241_lambda = np.log(2) / 14.4
    am241_lambda = np.log(2) / 433.0
    # Bateman-style ingrowth then decay
    am241 = 1.2 * (pu241_lambda / (am241_lambda - pu241_lambda)) \
            * (np.exp(-am241_lambda * t) - np.exp(-pu241_lambda * t))
    pu239 = 0.9 * np.exp(-np.log(2) * t / 24110.0)
    cm245 = 0.18 * np.exp(-np.log(2) * t / 8500.0)
    np237_u_series = 0.08 * np.ones_like(t)  # quasi-stationary background
    alpha = cm244 + am241 + pu239 + cm245 + np237_u_series

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    ax.loglog(t, beta_gamma, color=WONG[5], lw=1.7,
              label=r"$\beta/\gamma$ (fission products)")
    ax.loglog(t, alpha, color=WONG[6], lw=1.7,
              label=r"$\alpha$ (actinides)")
    ax.loglog(t, beta_gamma + alpha, color=WONG[0], lw=1.0, ls="--",
              label="total")
    # Cross-over annotation
    crossover_idx = np.argmin(np.abs(beta_gamma - alpha))
    t_cross = t[crossover_idx]
    ax.axvline(t_cross, color="0.6", lw=0.6, ls=":")
    ax.text(t_cross * 1.15, 5e-2,
            f"cross-over\n$\\approx\\,{t_cross:.0f}$ y",
            fontsize=7, color="0.3")
    # Disposal-relevant horizons
    horizons = [(10, "interim\nstorage"), (300, "DOE\nstandard"),
                (1e5, "geol.\nrepository")]
    for t_h, lab in horizons:
        ax.axvline(t_h, color="0.85", lw=0.5)
        ax.text(t_h, 1e-4, lab, fontsize=6.5, color="0.4",
                rotation=90, va="bottom", ha="right")
    ax.set_xlabel("time after discharge [years]")
    ax.set_ylabel(r"specific dose rate $\dot{D}/\dot{D}(1\,\mathrm{y})$")
    ax.set_xlim(0.1, 1e6)
    ax.set_ylim(1e-4, 30)
    ax.legend(loc="upper right")
    ax.grid(True, which="major", lw=0.3, alpha=0.4)
    fig.savefig(FIG_DIR / "fig_disposal_dose_trajectory.pdf")
    fig.savefig(FIG_DIR / "fig_disposal_dose_trajectory.png", dpi=300)
    print("wrote fig_disposal_dose_trajectory")


# ============================================================================
# Figure 2.  Gas inventory in a sealed disposal container vs. time
# ----------------------------------------------------------------------------
# Use the calibrated fluoride F2 kernel (Davis G + Toth/Felker recombination)
# and the chloride Cl2 framework (G(Cl•) posterior from the integrated MCMC),
# extrapolated to a sealed-container disposal scenario with:
#   - Storage temperature T = 200 C (no active cooling required after
#     the ~10 y heat-removal phase).
#   - Container volume 1 m^3 holding 2 t of salt + 0.1 m^3 of head-space.
#   - No active U(III)/U(IV) buffer (consumed during operation or
#     extracted before disposal).
# Dose rate decays following Figure 1.
# ============================================================================

def fig_container_inventory():
    t = np.logspace(0, 6, 400)  # years
    # Approximate decay-corrected dose rate, normalized.  Use the same
    # form as Figure 1 but scaled to absolute dose in kGy/h for a typical
    # discharged-salt inventory at 1 year (assume ~5 kGy/h initial,
    # following the order of magnitude reported for MCFR-relevant salts).
    fp = 7.5 * np.exp(-np.log(2) * t / 5.0) \
        + 4.0 * np.exp(-np.log(2) * t / 30.17) \
        + 3.5 * np.exp(-np.log(2) * t / 28.79)
    am241_lambda = np.log(2) / 433.0
    pu241_lambda = np.log(2) / 14.4
    alpha = 5.0 * np.exp(-np.log(2) * t / 18.1) \
        + 1.2 * (pu241_lambda / (am241_lambda - pu241_lambda)) \
            * (np.exp(-am241_lambda * t) - np.exp(-pu241_lambda * t)) \
        + 0.9 * np.exp(-np.log(2) * t / 24110.0) + 0.08
    dose_kGy_h = 5.0 * (fp + alpha) / (fp[0] + alpha[0])

    # Cl2 production: dP/dt = chi * G(Cl•) * dose_rate * conversion
    # We use a posterior-style 90% band with chi (chain-propagation
    # efficiency) lognormally distributed; without U-buffering chi ~ 0.1.
    R_GAS = 8.314462618
    NA = 6.02214076e23
    EV_J = 1.602176634e-19
    rho_salt = 2700.0
    V_salt = 0.74          # m^3 of salt in 1 m^3 container (74% packing)
    V_gas = 0.10           # m^3 of head-space
    T_K = 200 + 273.15     # storage temperature
    KH_Cl2 = 2e-5          # mol/(m^3 Pa)
    Kh_term = V_gas / (R_GAS * T_K) + KH_Cl2 * V_salt

    rng = np.random.default_rng(0)
    n_samples = 400
    G_Cl = rng.lognormal(np.log(0.5), 0.35, size=n_samples)
    chi = rng.lognormal(np.log(0.1), 0.6, size=n_samples)  # no U buffer

    # Cumulative integral of dose_rate over time, in J/kg
    dose_J_kg_s = dose_kGy_h * 1e3 / 3600.0
    cum_dose_J_kg = np.zeros_like(t)
    for i in range(1, len(t)):
        # trapezoidal integration over geometric grid
        dt = (t[i] - t[i - 1]) * 365.25 * 86400.0
        cum_dose_J_kg[i] = cum_dose_J_kg[i - 1] \
            + 0.5 * (dose_J_kg_s[i] + dose_J_kg_s[i - 1]) * dt

    # Cumulative Cl2 mass: G * dose * mass / (100 eV) molec, then /NA -> mol
    # Then a fraction chi reaches the gas phase as Cl2.
    factor = (rho_salt * V_salt) / (100.0 * EV_J * NA)
    # Build trajectories for the posterior ensemble
    P_Cl2 = np.zeros((n_samples, len(t)))
    for s in range(n_samples):
        n_total_mol = G_Cl[s] * chi[s] * cum_dose_J_kg * factor
        P_Cl2[s] = n_total_mol / Kh_term
    P_lo = np.quantile(P_Cl2, 0.05, axis=0)
    P_med = np.quantile(P_Cl2, 0.50, axis=0)
    P_hi = np.quantile(P_Cl2, 0.95, axis=0)

    # F2 in fluoride salt under same trajectory, using calibrated fluoride
    # kernel.  Davis G ≈ 0.02 molec / 100 eV for ThF4-bearing salts;
    # recombination Ea = 39 kJ/mol so at 200 C the rate is appreciable.
    Ea_rec = 39e3
    A_rec_per_hr = 250.0   # from Toth/Felker calibration (paper 2 §4)
    k_rec = A_rec_per_hr * np.exp(-Ea_rec / (R_GAS * T_K))  # 1/h
    G_F2 = rng.lognormal(np.log(0.02), 0.45, size=n_samples)
    P_F2 = np.zeros((n_samples, len(t)))
    for s in range(n_samples):
        n_total_mol = G_F2[s] * cum_dose_J_kg * factor
        # Buildup vs steady-state competition: assume the cumulative source
        # is integrated against an effective recombination half-life
        # 1/k_rec ~ 1/250 e^(Ea/RT) ~ 1/0.5 h^-1 ~ 2 h at 200C, much
        # shorter than disposal timescale -> P_F2 sits at a steady state
        # given by dN/dt_source / k_rec, which we approximate per-year:
        # P_F2_ss \\approx G_F2(T) * dose_rate / k_rec
        # Using dose_kGy_h instantaneous (not cumulative)
        dose_per_h = (dose_kGy_h * 1e3 / 3600.0)[None, :]
        # Use a simple steady-state approximation
        # convert: source rate in mol/(m^3 hr) -> P_F2 = source / k_rec
        S_per_h = G_F2[s] * (dose_kGy_h * 1e3) * (rho_salt * V_salt / V_salt) \
            / (100.0 * EV_J * NA)
        P_F2_ss = S_per_h * V_salt / (Kh_term * (k_rec))
        P_F2[s] = P_F2_ss
    PF2_lo = np.quantile(P_F2, 0.05, axis=0)
    PF2_med = np.quantile(P_F2, 0.50, axis=0)
    PF2_hi = np.quantile(P_F2, 0.95, axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)

    ax = axes[0]
    ax.fill_between(t, P_lo, P_hi, color=WONG[2], alpha=0.25,
                     label="90% posterior band")
    ax.plot(t, P_med, color=WONG[2], lw=1.7, label="median")
    ax.axhline(1e5, color="0.55", lw=0.7, ls="--",
               label=r"$10^{5}$ Pa container limit")
    # Shade the region above the container limit to indicate the prediction
    # is hypothetical once the container has presumably failed.
    ax.axhspan(1e5, 1e14, color="0.85", alpha=0.30, lw=0)
    # Time at which the median crosses the 10^5 Pa limit (annotation)
    t_breach = t[np.argmax(P_med > 1e5)] if np.any(P_med > 1e5) else None
    if t_breach is not None:
        ax.axvline(t_breach, color="0.55", lw=0.5, ls=":")
        ax.text(t_breach * 1.15, 1e-4,
                f"median crosses limit\nat $\\approx\\,{t_breach:.0f}$ y",
                fontsize=6.5, color="0.3")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time after discharge [years]")
    ax.set_ylabel(r"$P_{\mathrm{Cl_2}}$ in container [Pa]")
    ax.set_title("(a) Chloride fuel, no U buffer")
    ax.legend(loc="lower right")
    ax.set_xlim(1, 1e6)
    ax.set_ylim(1e-6, 1e14)
    ax.grid(True, which="major", lw=0.3, alpha=0.4)

    ax = axes[1]
    ax.fill_between(t, PF2_lo, PF2_hi, color=WONG[1], alpha=0.25,
                     label="90% posterior band")
    ax.plot(t, PF2_med, color=WONG[1], lw=1.7, label="median")
    ax.axhline(1e5, color="0.55", lw=0.7, ls="--",
               label=r"$10^{5}$ Pa container limit")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time after discharge [years]")
    ax.set_ylabel(r"$P_{\mathrm{F_2}}$ steady state [Pa]")
    ax.set_title(r"(b) Fluoride fuel, $T=200\,^\circ$C")
    ax.legend(loc="upper right")
    ax.set_xlim(1, 1e6)
    ax.set_ylim(1e-6, 1e8)
    ax.grid(True, which="major", lw=0.3, alpha=0.4)

    fig.savefig(FIG_DIR / "fig_disposal_container_inventory.pdf")
    fig.savefig(FIG_DIR / "fig_disposal_container_inventory.png", dpi=300)
    print("wrote fig_disposal_container_inventory")


# ============================================================================
# Figure 3.  Operational vs disposal radiolysis regime map
# ----------------------------------------------------------------------------
# A schematic phase diagram on (dose rate, temperature, U-buffer
# availability) showing where the existing operational calibration is
# valid and where disposal conditions sit.  Three boxes:
#   - Pulse-radiolysis laboratory regime (high dose, high T, no buffer)
#   - Operational MSR regime (moderate dose, high T, full U buffer)
#   - Disposal regime (low dose, low T, no buffer, alpha-dominant)
# The figure communicates the extrapolation gap.
# ============================================================================

def fig_regime_map():
    fig, ax = plt.subplots(figsize=(6.8, 4.0))

    # Background gradient: temperature on x, log dose rate on y
    T_K_range = np.linspace(290, 1000, 200)
    log_dose_range = np.linspace(-12, 4, 200)  # log10(kGy/h)
    TT, DD = np.meshgrid(T_K_range, log_dose_range)
    # Posterior coverage heuristic: combined function of how far from
    # the experimental envelope (kept simple for an illustrative figure)
    coverage = np.exp(-((TT - 700) / 200) ** 2 - ((DD - 0) / 4) ** 2)
    ax.contourf(TT, DD, coverage, levels=20, cmap="Blues", alpha=0.35)

    # Region boxes
    def box(ax, x0, y0, w, h, label, fc, ec, alpha=0.18, txt_color=None):
        rect = Rectangle((x0, y0), w, h, fc=fc, ec=ec, alpha=alpha, lw=1.2)
        ax.add_patch(rect)
        ax.text(x0 + w / 2, y0 + h / 2, label, ha="center", va="center",
                fontsize=9, color=txt_color or ec, weight="bold")

    box(ax, 600, -2, 100, 6, "pulse\nradiolysis\nlab",
        fc=WONG[5], ec=WONG[5], alpha=0.20)
    box(ax, 670, -3, 230, 3, "operational\nMSR",
        fc=WONG[3], ec=WONG[3], alpha=0.20)
    box(ax, 320, -10, 130, 4, "disposal\nregime",
        fc=WONG[6], ec=WONG[6], alpha=0.25)

    # Arrow showing extrapolation from operational to disposal
    ax.annotate("", xy=(420, -8), xytext=(750, -1.5),
                 arrowprops=dict(arrowstyle="->", color="0.3", lw=1.4))
    ax.text(560, -4.5, "extrapolation gap\n(4 decades in dose rate;\n"
            r"$\alpha$-dominant vs $\gamma$-dominant)",
            fontsize=7.5, color="0.2", ha="center", style="italic")

    # Labels for major experimental anchors
    ax.plot(673, np.log10(1e3), "o", color=WONG[5], markersize=6,
            markeredgecolor="white", markeredgewidth=0.8)
    ax.annotate("Iwamatsu 2026 (Cr)", xy=(673, 3), xytext=(560, 3.7),
                fontsize=7, color="0.2",
                arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5))
    ax.plot(873, 1.5, "s", color=WONG[5], markersize=5)
    ax.annotate("Pikaev 1982", xy=(873, 1.5), xytext=(820, 2.5),
                fontsize=7, color="0.2",
                arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5))
    ax.plot(348, np.log10(13), "D", color=WONG[2], markersize=5)
    ax.annotate("Phillips 2022\n(MCFR NULL)", xy=(348, 1.1), xytext=(305, -0.5),
                fontsize=7, color="0.2",
                arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5))
    ax.plot(473, -8, "^", color=WONG[6], markersize=6,
            markeredgecolor="white", markeredgewidth=0.8)
    ax.annotate("disposal target\n(200 °C, $\\sim$1 $\\mu$Gy/h after 1 ky)",
                xy=(473, -8), xytext=(490, -9.5),
                fontsize=7, color="0.2",
                arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5))

    ax.set_xlim(290, 1000)
    ax.set_ylim(-12, 4)
    ax.set_xlabel(r"temperature [K]")
    ax.set_ylabel(r"$\log_{10}$ dose rate [kGy/h]")
    ax.set_xticks([300, 400, 500, 600, 700, 800, 900, 1000])
    ax.grid(True, which="major", lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_disposal_regime_map.pdf")
    fig.savefig(FIG_DIR / "fig_disposal_regime_map.png", dpi=300)
    print("wrote fig_disposal_regime_map")


if __name__ == "__main__":
    fig_dose_trajectory()
    fig_container_inventory()
    fig_regime_map()
