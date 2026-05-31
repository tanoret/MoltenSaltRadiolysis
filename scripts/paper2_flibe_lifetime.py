#!/usr/bin/env python3
"""Paper 2 §8(b): FLiBe-UF4 MSR F2 inventory over reactor lifetime.

Predicts the cover-gas F2 partial pressure for a fluoride-fueled MSR over a 60-year
operational lifetime under realistic dose rate.

Methodology:
  - Use the static F2-production kernel from the calibrated fluoride layer:
        G(F2) from Davis 2022 (FLiBe-UF4: 0.005-0.007 molec/100eV)
        recombination Arrhenius A_rec, Ea_rec from Toth-Felker 1990
  - Production rate: r_F2 = G(F2)/100/eV * dose_rate / N_A
  - Loss rate: k_rec(T) * F2_dissolved (1st order in dissolved F2)
  - Steady-state P(F2) governed by the balance equation; transient approach is
    fast (hours) compared to the 60-y reactor lifetime, so steady-state is
    achieved early.
  - Compare to a notional design-basis safety limit (100 Pa F2 in cover gas).

Outputs:
  manuscript/figures/fig_flibe_f2_lifetime.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
VAL = REPO / "validation"
FIGS = REPO / "manuscript" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

WONG = [
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
]
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times", "Computer Modern Roman"],
    "font.size": 9,
    "axes.labelsize": 9, "axes.titlesize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "lines.linewidth": 1.4,
    "axes.spines.right": False, "axes.spines.top": False,
    "figure.dpi": 110, "savefig.dpi": 600,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    "legend.frameon": False, "mathtext.fontset": "cm",
})


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

R_GAS = 8.314462618
NA = 6.02214076e23
EV_J = 1.602176634e-19
KH_F2 = 1e-5  # mol/(m^3·Pa) Henry constant for F2 in fluoride melt (estimate)

# =============================================================================
# FLiBe-UF4 MSR plant model
# Reference Davis 2022 used FLiBe-UF4 at HFR dose ~89 kGy/h to measure G(F2);
# Toth-Felker MSRE-composition tests used in-pile gamma at ~similar order.
# For an operational MSR we adopt 10 kGy/h as a representative spectrum-averaged
# dose rate to the salt (lower than research-reactor values because most fission
# energy escapes as fast neutrons absorbed elsewhere).
# =============================================================================
T_K = 873.15            # 600 °C operating temperature
T_C = T_K - 273.15
V_salt = 50.0           # salt loop volume (m^3)
V_cover = 5.0           # cover-gas volume (m^3)
mass_salt = 2200.0 * V_salt   # kg (FLiBe ρ ≈ 2200 kg/m^3 at 600°C)
dose_rate_kGy_h = 10.0  # representative dose rate
DAVIS_REF_DOSE_kGy_h = 89.0  # FLiBe-UF4 row in Davis Table III
scale_dose = dose_rate_kGy_h / DAVIS_REF_DOSE_kGy_h
dose_rate_J_kg_s = dose_rate_kGy_h * 1e3 / 3600.0
dose_rate_J_m3_s = dose_rate_J_kg_s * 2200.0
safety_limit_Pa = 100.0

# =============================================================================
# G(F2) posterior — from Davis 2022 FLiBe-UF4 row + max-slope row
# =============================================================================
davis = pd.read_csv(VAL / "eflibe_f2_yield/davis_2022_nse/data/G_values_table_III.csv",
                    comment="#")
flibe_rows = davis[davis["salt"].str.startswith("FLiBe_UF4")].dropna(subset=["G_F2"])
G_F2_mu = flibe_rows["G_F2"].astype(float).values
G_F2_sigma = flibe_rows["sigma"].astype(float).values
# Use combined Gaussian mixture: G ~ N(0.005, 0.001) ∪ N(0.007, 0.001)
n_samples = 500
rng = np.random.default_rng(2026)
mix = rng.integers(0, len(G_F2_mu), size=n_samples)
G_F2 = rng.normal(G_F2_mu[mix], G_F2_sigma[mix])
G_F2 = np.clip(G_F2, 0.002, 0.020)

# =============================================================================
# Toth-Felker recombination kernel
# =============================================================================
fkernel = pd.read_csv(VAL / "fluoride_kernel/calibrated_parameters.csv")
A_rec = float(fkernel[fkernel["parameter"] == "A_rec"]["value"].iloc[0])  # 1/h
Ea_rec_J_mol = float(fkernel[fkernel["parameter"] == "Ea_rec_J_mol"]["value"].iloc[0])
# Convert to SI: rate constant for first-order F2 loss = k_rec(T) (1/s) on dissolved F2
def k_rec(T_K_):
    """Return k_rec in 1/s; A_rec given in 1/h."""
    return A_rec / 3600.0 * np.exp(-Ea_rec_J_mol / (R_GAS * T_K_))

# =============================================================================
# Steady-state F2 partial pressure (production = recombination)
# d[F2_diss]/dt = r_F2 - k_rec * [F2_diss] ≈ 0 in steady state
# We use the pre-calibrated fluoride kernel which absorbs the Toth-Felker
# reference into A_rec; the kernel's predicted P_ss(FLiBe-UF4) at the
# Davis dose rate is read directly and scaled by the operational/reference
# dose-rate ratio.
#
# r_F2 = G_F2 × dose × Avogadro factor; P_ss ∝ r_F2 / k_rec(T)
# =============================================================================
factor = dose_rate_J_m3_s / (100.0 * EV_J) / NA
# Use kernel-based P_ss at this T, then scale by (G_used / G_kernel_anchor)
# G(FLiBe-UF4) median from Davis Table III = 0.005, max-slope = 0.007 (mean 0.006)
G_kernel_anchor = 0.006
fl_curves = pd.read_csv(VAL / "fluoride_kernel/predicted_buildup_curves.csv")
# Interpolate kernel P_ss(T) at operational T
P_ss_kernel_at_T = np.interp(T_K, fl_curves["T_K"], fl_curves["P_ss_Pa_FLiBe_UF4"])
# Scale by both dose-rate ratio and G(F2) sample / kernel-anchor G
P_F2_ss = P_ss_kernel_at_T * scale_dose * (G_F2 / G_kernel_anchor)
# (No Henry split needed: the kernel reports gas-phase partial pressure directly.)

# =============================================================================
# Transient approach to steady state (relaxation time = 1/k_rec)
# Hours scale — fast compared to 60-y lifetime. Plot only the lifetime-scale
# steady-state inventory.
# =============================================================================
T_years = np.linspace(0.0, 60.0, 121)
P_F2_t = np.zeros((n_samples, len(T_years)))
for i in range(n_samples):
    # Steady state reached within hours; constant over reactor life
    P_F2_t[i, :] = P_F2_ss[i]

# Build a posterior band
P_lo = np.quantile(P_F2_t, 0.05, axis=0)
P_med = np.quantile(P_F2_t, 0.50, axis=0)
P_hi = np.quantile(P_F2_t, 0.95, axis=0)

# Temperature-dependent steady state for the right panel
T_axis_C = np.linspace(40, 700, 60)
T_axis_K = T_axis_C + 273.15
P_kernel_T = np.interp(T_axis_K, fl_curves["T_K"], fl_curves["P_ss_Pa_FLiBe_UF4"])
P_ss_T = np.zeros((n_samples, len(T_axis_C)))
for i, g in enumerate(G_F2):
    P_ss_T[i, :] = P_kernel_T * scale_dose * (g / G_kernel_anchor)

PT_lo = np.quantile(P_ss_T, 0.05, axis=0)
PT_med = np.quantile(P_ss_T, 0.50, axis=0)
PT_hi = np.quantile(P_ss_T, 0.95, axis=0)

# =============================================================================
# Figure
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
ax, ax2 = axes

ax.fill_between(T_years, P_lo, P_hi, color=WONG[3], alpha=0.30,
                label="HBMAE 90% band")
ax.plot(T_years, P_med, color=WONG[3], lw=1.6, label="HBMAE median")
ax.axhline(safety_limit_Pa, color="#bbbbbb", lw=1, ls=":",
           label="100 Pa limit")
ax.set_yscale("log")
ax.set_xlabel("operating years")
ax.set_ylabel(r"$P_{\mathrm{F}_2}$ cover gas [Pa]")
ax.set_title(r"(a) Cover-gas inventory @ 600 $^\circ$C")
ax.set_xlim(0, 60)
ax.set_ylim(max(1e-12, P_med[0]/100), max(safety_limit_Pa*10, P_hi[-1]*10))
ax.legend(loc="upper right", fontsize=7, frameon=False)
_despine(ax)

# Right panel: P(F2) vs T (recombination kernel)
ax2.fill_between(T_axis_C, PT_lo, PT_hi, color=WONG[6], alpha=0.30,
                 label="HBMAE 90% band")
ax2.plot(T_axis_C, PT_med, color=WONG[6], lw=1.6, label="HBMAE median")
ax2.axhline(safety_limit_Pa, color="#bbbbbb", lw=1, ls=":",
            label="100 Pa limit")
ax2.axvline(150.0, color=WONG[1], lw=1, ls="--",
            label=r"$T_{\mathrm{bal}}=150\,^\circ$C")
ax2.set_yscale("log")
ax2.set_xlabel(r"$T$ [$^\circ$C]")
ax2.set_ylabel(r"$P_{\mathrm{F}_2,\mathrm{ss}}$ [Pa]")
ax2.set_title(r"(b) $P(\mathrm{F}_2)_{\mathrm{ss}}$ vs $T$")
ax2.legend(loc="upper right", fontsize=7, frameon=False)
_despine(ax2)

fig.savefig(FIGS / "fig_flibe_f2_lifetime.pdf", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "fig_flibe_f2_lifetime.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("FLiBe-UF4 F2 cover-gas pressure (steady state, T=873 K):")
print(f"  5%   : {P_lo[-1]:.3e} Pa")
print(f"  50%  : {P_med[-1]:.3e} Pa")
print(f"  95%  : {P_hi[-1]:.3e} Pa")
print(f"P(exceed limit): {(P_F2_t[:, -1] > safety_limit_Pa).mean():.3f}")
print(f"figures: {FIGS / 'fig_flibe_f2_lifetime.pdf'}")
