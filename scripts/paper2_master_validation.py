#!/usr/bin/env python3
"""Paper 2 master validation figure.

Builds a 4x4 multi-panel figure overlaying the HBMAE-calibrated model
predictions on every digitized experimental data set in the chloride +
fluoride database.  Outputs:

    manuscript/figures/fig_master_validation.{pdf,png}
    manuscript/figures/fig_residuals_all.{pdf,png}
    manuscript/figures/fig_data_landscape.{pdf,png}
"""
from __future__ import annotations

import csv
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

R_GAS = 8.314462618

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

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times", "Computer Modern Roman"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "lines.linewidth": 1.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "legend.frameon": False,
    "mathtext.fontset": "cm",
})


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def read_csv(path):
    return pd.read_csv(path, comment="#")


# ----- Posterior chains -----
T4 = np.load(VAL / "tier4_integrated_chain.npy")
# Columns: 0 logA5, 1 Ea5, 2 logA6, 3 Ea6, 4..12 log[e_s]0_t1..t9,
# 13 log k_bg, 14 logA_Zn, 15 Ea_Zn, 16 eta_logA_LK, 17 eta_Ea_LK,
# 18 eta_logA_Na, 19 eta_Ea_Na, 20 eta_logA_K, 21 eta_Ea_K,
# 22 b_Pikaev, 23 logG_Cl_Phillips
logA5 = T4[:, 0]; Ea5 = T4[:, 1]
logA6 = T4[:, 2]; Ea6 = T4[:, 3]
logA_Zn_intr = T4[:, 14]; Ea_Zn_intr = T4[:, 15]
eta_logA_LK = T4[:, 16]; eta_Ea_LK = T4[:, 17]
eta_logA_Na = T4[:, 18]; eta_Ea_Na = T4[:, 19]
eta_logA_K = T4[:, 20]; eta_Ea_K = T4[:, 21]
b_Pikaev = T4[:, 22]


def k_arr(logA, Ea, T_K):
    return np.exp(logA - Ea / (R_GAS * T_K))


# =============================================================================
fig, axes = plt.subplots(4, 4, figsize=(13.0, 13.0), constrained_layout=True)
axes = axes.ravel()

# Common axis labels
LBL_T = r"$T$ [K]"
LBL_K = r"$k$ [M$^{-1}$ s$^{-1}$]"
LBL_LOGK = r"$\log_{10} k$ [M$^{-1}$ s$^{-1}$]"
LBL_INVT = r"$1000/T$ [K$^{-1}$]"

# --- Panel a: Cr2+ + e_s- k(T)
ax = axes[0]
cr_kT = read_csv(VAL / "cr_licl_kcl/iwamatsu_2026_pccp/data/k_vs_T_from_arrhenius.csv")
cr2 = cr_kT[cr_kT["reaction_id"] == "eq5"]
cr3 = cr_kT[cr_kT["reaction_id"] == "eq6"]
T_grid = np.linspace(660, 880, 60)
ks = np.array([k_arr(logA5, Ea5, T) for T in T_grid])
q = np.quantile(ks, [0.05, 0.5, 0.95], axis=1)
ax.fill_between(T_grid, q[0], q[2], color=WONG[2], alpha=0.25, label="HBMAE 90%")
ax.plot(T_grid, q[1], color=WONG[2], lw=1.4, label="HBMAE median")
ax.errorbar(cr2["T_K"], cr2["k_M_inv_s"], yerr=cr2["sigma_k_M_inv_s"],
            fmt="o", color=WONG[6], mec="k", mew=0.4, ms=5, label="Iwamatsu 2026")
ax.set_yscale("log"); ax.set_xlabel(LBL_T); ax.set_ylabel(LBL_K)
ax.set_title(r"(a) $e_s^- + \mathrm{Cr}^{2+} \to \mathrm{Cr}^+$")
ax.legend(loc="lower right", frameon=False); _despine(ax)

# --- Panel b: Cr3+ + e_s- k(T)
ax = axes[1]
ks = np.array([k_arr(logA6, Ea6, T) for T in T_grid])
q = np.quantile(ks, [0.05, 0.5, 0.95], axis=1)
ax.fill_between(T_grid, q[0], q[2], color=WONG[3], alpha=0.25, label="HBMAE 90%")
ax.plot(T_grid, q[1], color=WONG[3], lw=1.4, label="HBMAE median")
ax.errorbar(cr3["T_K"], cr3["k_M_inv_s"], yerr=cr3["sigma_k_M_inv_s"],
            fmt="s", color=WONG[6], mec="k", mew=0.4, ms=5, label="Iwamatsu 2026")
ax.set_yscale("log"); ax.set_xlabel(LBL_T); ax.set_ylabel(LBL_K)
ax.set_title(r"(b) $e_s^- + \mathrm{Cr}^{3+} \to \mathrm{Cr}^{2+}$")
ax.legend(loc="lower right", frameon=False); _despine(ax)

# --- Panel c: Zn2+ + e_s- k(T) in LiCl-KCl
ax = axes[2]
zn = read_csv(VAL / "zn_licl_kcl/horne_2022_pccp/data/k_vs_T_from_arrhenius.csv")
zn_refit = zn[zn["reaction_id"] == "eq5_refit_2026"]
ks = np.array([k_arr(logA_Zn_intr + eta_logA_LK,
                     Ea_Zn_intr + eta_Ea_LK, T) for T in T_grid])
q = np.quantile(ks, [0.05, 0.5, 0.95], axis=1)
ax.fill_between(T_grid, q[0], q[2], color=WONG[1], alpha=0.25, label="HBMAE 90%")
ax.plot(T_grid, q[1], color=WONG[1], lw=1.4, label="HBMAE median")
ax.errorbar(zn_refit["T_K"], zn_refit["k_M_inv_s"],
            yerr=zn_refit["sigma_k_M_inv_s"], fmt="o",
            color=WONG[6], mec="k", mew=0.4, ms=5, label="Iwamatsu 2026 refit")
ax.set_yscale("log"); ax.set_xlabel(LBL_T); ax.set_ylabel(LBL_K)
ax.set_title(r"(c) $e_s^- + \mathrm{Zn}^{2+}$ LiCl-KCl")
ax.legend(loc="lower right", frameon=False); _despine(ax)

# --- Panel d: Pikaev/Iwamatsu Zn cross-paper resolution (NaCl + KCl)
ax = axes[3]
T_pik = 1123.15; T_kcl = 1073.15
ks_intr_Na = k_arr(logA_Zn_intr + eta_logA_Na, Ea_Zn_intr + eta_Ea_Na, T_pik)
ks_intr_K  = k_arr(logA_Zn_intr + eta_logA_K, Ea_Zn_intr + eta_Ea_K, T_kcl)
ks_obs_Na = 10**(np.log10(ks_intr_Na) + b_Pikaev)
ks_obs_K  = 10**(np.log10(ks_intr_K)  + b_Pikaev)
q_intr = np.quantile(np.log10(ks_intr_Na), [0.05, 0.5, 0.95])
q_obs  = np.quantile(np.log10(ks_obs_Na),  [0.05, 0.5, 0.95])
q_intrK = np.quantile(np.log10(ks_intr_K), [0.05, 0.5, 0.95])
q_obsK  = np.quantile(np.log10(ks_obs_K),  [0.05, 0.5, 0.95])

ax.errorbar([0.7], [q_intr[1]], yerr=[[q_intr[1]-q_intr[0]], [q_intr[2]-q_intr[1]]],
            fmt="D", color=WONG[2], mec="k", mew=0.4, ms=6, label="intr (NaCl)")
ax.errorbar([1.3], [q_obs[1]], yerr=[[q_obs[1]-q_obs[0]], [q_obs[2]-q_obs[1]]],
            fmt="s", color=WONG[1], mec="k", mew=0.4, ms=6,
            label=r"intr+$b^{\mathrm{Pik}}$ (NaCl)")
ax.errorbar([2.0], [np.log10(1.7e9)], yerr=[[0.13], [0.13]], fmt="o",
            color=WONG[6], mec="k", mew=0.4, ms=6, label="Pikaev (NaCl)")
ax.errorbar([2.7], [q_intrK[1]],
            yerr=[[q_intrK[1]-q_intrK[0]], [q_intrK[2]-q_intrK[1]]],
            fmt="D", color=WONG[3], mec="k", mew=0.4, ms=6, label="intr (KCl)")
ax.errorbar([3.3], [q_obsK[1]],
            yerr=[[q_obsK[1]-q_obsK[0]], [q_obsK[2]-q_obsK[1]]],
            fmt="s", color=WONG[7], mec="k", mew=0.4, ms=6,
            label=r"intr+$b^{\mathrm{Pik}}$ (KCl)")
ax.errorbar([4.0], [np.log10(2.8e9)], yerr=[[0.13], [0.13]], fmt="o",
            color=WONG[6], mec="k", mew=0.4, ms=6, label="Pikaev (KCl)")
ax.set_xticks([0.7, 1.3, 2.0, 2.7, 3.3, 4.0])
ax.set_xticklabels(["intr", "bias", "obs", "intr", "bias", "obs"], fontsize=7)
ax.set_ylabel(LBL_LOGK)
ax.set_title(r"(d) Zn$^{2+}$ Pikaev resolution")
ax.legend(loc="lower left", fontsize=6, ncol=1, frameon=False); _despine(ax)

# --- Panel e: Nd3+ Arrhenius
ax = axes[4]
nd = read_csv(VAL / "nd_licl_kcl/castro_baldivieso_2026_ic/data/k_vs_T.csv")
lomo = read_csv(VAL / "meta_hier/lomo_results.csv")
nd_lomo = lomo[lomo["metal"] == "Nd3+"].sort_values("T_K")
ax.fill_between(nd_lomo["T_K"], 10**nd_lomo["log10_k_pred_p05"],
                10**nd_lomo["log10_k_pred_p95"], color=WONG[7], alpha=0.18,
                label="meta-hier 90% (LOMO)")
ax.plot(nd_lomo["T_K"], 10**nd_lomo["log10_k_pred_mean"], color=WONG[7],
        lw=1.2, ls="--", label="meta-hier mean")
ax.errorbar(nd["T_K"], nd["k_M_inv_s"], yerr=nd["sigma_k_M_inv_s"],
            fmt="o", color=WONG[6], mec="k", mew=0.4, ms=5,
            label="Castro-Baldivieso 2026")
ax.set_yscale("log"); ax.set_xlabel(LBL_T); ax.set_ylabel(LBL_K)
ax.set_title(r"(e) $e_s^- + \mathrm{Nd}^{3+}$ LiCl-KCl")
ax.legend(loc="lower right", frameon=False); _despine(ax)

# --- Panel f: Cf3+ at 295 K (aqueous baseline)
ax = axes[5]
cf_lomo = lomo[lomo["metal"] == "Cf3+"].iloc[0]
ax.errorbar([1], [cf_lomo["log10_k_pred_mean"]],
            yerr=[[cf_lomo["log10_k_pred_mean"]-cf_lomo["log10_k_pred_p05"]],
                  [cf_lomo["log10_k_pred_p95"]-cf_lomo["log10_k_pred_mean"]]],
            fmt="D", color=WONG[2], mec="k", mew=0.4, ms=6,
            label="meta-hier LOMO")
ax.errorbar([2], [cf_lomo["log10_k_obs"]], yerr=[[0.02], [0.02]],
            fmt="o", color=WONG[6], mec="k", mew=0.4, ms=6,
            label="Rotermund 2024")
ax.set_xticks([1, 2]); ax.set_xticklabels(["pred", "obs"])
ax.set_xlim(0.5, 2.5)
ax.set_ylabel(LBL_LOGK)
ax.set_title(r"(f) Cf$^{3+}+e_{aq}^-$ @ 295 K")
ax.legend(loc="lower right", frameon=False); _despine(ax)

# --- Panel g: Pikaev multi-metal e_s- rates at hot temperatures
ax = axes[6]
pik = read_csv(VAL / "alkali_halide_baseline/pikaev_1982_rpc/data/es_rate_constants.csv")
hosts_color = {"NaCl": WONG[1], "KCl": WONG[2], "KBr": WONG[3]}
z_map = {"Na+": 1, "K+": 1, "Ag+": 1, "Tl+": 1, "Ba2+": 2, "Sr2+": 2,
         "Ca2+": 2, "Cd2+": 2, "Zn2+": 2}
rng = np.random.RandomState(1)
for _, r in pik.iterrows():
    if r["solute"] in z_map:
        x = z_map[r["solute"]] + rng.uniform(-0.1, 0.1)
        y = np.log10(float(r["k_M1_s1"]))
        c = hosts_color.get(r["melt"], WONG[0])
        ax.scatter(x, y, c=c, s=28, edgecolors="k", linewidths=0.3)
ax.set_xticks([1, 2]); ax.set_xlabel(r"metal charge $z$")
ax.set_ylabel(LBL_LOGK)
ax.set_title(r"(g) Pikaev 1982 $e_s^-$ scavengers")
for h, c in hosts_color.items():
    ax.scatter([], [], c=c, s=28, edgecolors="k", linewidths=0.3, label=h)
ax.legend(loc="lower right", title="host", frameon=False); _despine(ax)

# --- Panel h: Hagiwara Cl2- Arrhenius
ax = axes[7]
hag = read_csv(VAL / "cl2m_licl_kcl/hagiwara_1987_rpc/data/rate_constants.csv")
hag_pure = hag[hag["reaction"] == "2k_over_eps_l_second_order_decay"].copy()
ax.plot(1000.0/hag_pure["T_K"], np.log10(hag_pure["k_value"]),
        "o", color=WONG[6], mec="k", mew=0.4, ms=5, label="Hagiwara 1987")
A = np.vstack([1000.0/hag_pure["T_K"], np.ones(len(hag_pure))]).T
slope, intercept = np.linalg.lstsq(A, np.log10(hag_pure["k_value"]), rcond=None)[0]
Ea_h = -slope * np.log(10) * R_GAS
T_line = np.linspace(hag_pure["T_K"].min(), hag_pure["T_K"].max(), 50)
ax.plot(1000.0/T_line, slope*(1000.0/T_line) + intercept,
        "-", color=WONG[2], lw=1.2,
        label=fr"OLS $E_a$={Ea_h:.0f} J/mol")
ax.set_xlabel(LBL_INVT)
ax.set_ylabel(r"$\log_{10}(2k/\varepsilon\ell)$ [s$^{-1}$]")
ax.set_title(r"(h) Hagiwara Cl$_2^{\bullet-}$")
ax.legend(loc="upper right", fontsize=6, frameon=False); _despine(ax)

# --- Panel i: Conrad 2023 iodide perturbation
ax = axes[8]
conrad = [
    ("Cl$_2^{\\bullet-}$ Hagiwara 1987", 24.0, 0.0),
    ("Cl$_2^{\\bullet-}$ Iwamatsu 2022", 25.0, 1.0),
    ("ICl$^{\\bullet-}$ Conrad 2023", 23.4, 5.7),
]
yp = np.arange(len(conrad))[::-1]
for i, (label, ea, sig) in enumerate(conrad):
    ax.errorbar([ea], [yp[i]], xerr=[[sig], [sig]], fmt="o",
                color=WONG[1 + i], mec="k", mew=0.4, ms=6)
ax.set_yticks(yp); ax.set_yticklabels([c[0] for c in conrad], fontsize=6)
ax.set_xlabel(r"$E_a$ [kJ mol$^{-1}$]")
ax.set_xlim(15, 40)
ax.set_title(r"(i) Conrad iodide perturb.")
_despine(ax)

# --- Panel j: Davis 2022 G(F2) by salt
ax = axes[9]
davis = read_csv(VAL / "eflibe_f2_yield/davis_2022_nse/data/G_values_table_III.csv")
davis = davis.dropna(subset=["G_F2"]).head(10).reset_index(drop=True)
y = np.arange(len(davis))
sigs = davis["sigma"].fillna(0.001)
ax.barh(y, davis["G_F2"], xerr=sigs, color=WONG[3], edgecolor="k", linewidth=0.4)
ax.set_yticks(y); ax.set_yticklabels(davis["salt"], fontsize=6)
ax.set_xlabel(r"$G(\mathrm{F}_2)$ [molec / 100 eV]")
ax.set_title(r"(j) Davis fluoride $G$-values")
ax.invert_yaxis()
_despine(ax)

# --- Panel k: Toth-Felker FLiBe-UF4 buildup vs T
ax = axes[10]
fl = pd.read_csv(VAL / "fluoride_kernel/predicted_buildup_curves.csv")
ax.plot(fl["T_C"], fl["P_ss_Pa_FLiBe_UF4"], color=WONG[2], lw=1.4, label=r"FLiBe-UF$_4$")
ax.plot(fl["T_C"], fl["P_ss_Pa_BeF2"], color=WONG[6], lw=1.0, ls="--", label=r"BeF$_2$")
ax.plot(fl["T_C"], fl["P_ss_Pa_LiF"], color=WONG[3], lw=1.0, ls=":", label="LiF")
ax.plot(fl["T_C"], fl["P_ss_Pa_65LiF_29BeF2_5ZrF4_0.66UF4_Co60"],
        color=WONG[1], lw=1.0, ls="-.", label="MSRE Co-60")
ax.axvline(150.0, color="#bbbbbb", ls="--", lw=0.8,
           label=r"$T_{\mathrm{bal}}=150\,^\circ\mathrm{C}$")
ax.set_yscale("log"); ax.set_xlabel(r"$T$ [$^\circ$C]")
ax.set_ylabel(r"$P_{\mathrm{F}_2,\mathrm{ss}}$ [Pa]")
ax.set_title(r"(k) Steady-state $P(\mathrm{F}_2)$")
ax.legend(loc="upper right", fontsize=6, frameon=False); _despine(ax)

# --- Panel l: Makarov X2- molar absorptivity vs salt
ax = axes[11]
mak = read_csv(VAL / "oxidants_halide_melts/makarov_1982_bull/data/x2m_molar_absorptivity.csv")
mak = mak[(mak["band"] == "I") & mak["eps_M_inv_cm_inv"].notna()]
salts_x = list(mak["system"].unique())
xpos = {s: i for i, s in enumerate(salts_x)}
spec_color = {"Cl2-": WONG[6], "Br2-": WONG[2], "I2-": WONG[3]}
for _, r in mak.iterrows():
    ax.scatter(xpos[r["system"]], r["eps_M_inv_cm_inv"],
               color=spec_color.get(r["species"], WONG[0]),
               edgecolors="k", linewidths=0.3, s=36)
ax.set_xticks(list(xpos.values()))
ax.set_xticklabels(list(xpos.keys()), rotation=55, fontsize=5)
ax.set_ylabel(r"$\varepsilon_{\max}$ [M$^{-1}$ cm$^{-1}$]")
ax.set_title(r"(l) Makarov X$_2^{\bullet-}$ abs.")
for sp, c in spec_color.items():
    ax.scatter([], [], color=c, edgecolors="k", linewidths=0.3, s=36, label=sp)
ax.legend(loc="upper right", fontsize=6, frameon=False); _despine(ax)

# --- Panel m: Phillips 2022 NULL benchmark
ax = axes[12]
phil = read_csv(VAL / "TIER3_PHILLIPS_NULL.csv")
phil_num = phil[phil["G_Cl"].apply(lambda x: str(x).replace(".", "", 1).lstrip("-").replace("e", "").replace("E", "").isdigit())]
phil_num = phil_num.astype({"G_Cl": float, "Cl2_gas_mol_m3": float, "ratio_to_threshold": float})
ax.loglog(phil_num["G_Cl"], phil_num["Cl2_gas_mol_m3"], "o-",
          color=WONG[2], mec="k", mew=0.4, ms=5, label="with U(III/IV)")
without_U = phil[phil["G_Cl"] == "GAMMA_B_WITHOUT_U"]
if len(without_U) > 0:
    p_noU = float(without_U.iloc[0]["Cl2_gas_mol_m3"])
    ax.axhline(p_noU, color=WONG[6], ls="--", lw=1,
               label=fr"w/o U: {p_noU:.1e} mol m$^{{-3}}$")
ax.axhline(1.38e-2, color="#bbbbbb", ls=":", lw=1, label="Phillips limit")
ax.set_xlabel(r"$G(\mathrm{Cl}^\bullet)$ [molec / 100 eV]")
ax.set_ylabel(r"$[\mathrm{Cl}_2]_{\mathrm{gas}}$ [mol m$^{-3}$]")
ax.set_title(r"(m) Phillips NULL benchmark")
ax.legend(loc="lower right", fontsize=6, frameon=False); _despine(ax)

# --- Panel n: Cr3+ transient
ax = axes[13]
tr = pd.read_csv(VAL / "cr_licl_kcl/iwamatsu_2026_pccp/data/absorbance3mMCr3.csv")
t_us = tr["time"].values
abs_obs = tr["absorbance"].values
ax.plot(t_us, abs_obs, "o", color=WONG[6], mec="k", mew=0.3, ms=4,
        label=r"3 mM Cr$^{3+}$ @ 400 $^\circ$C")
k6_samp = k_arr(logA6, Ea6, 673.15)
k_obs_samp = k6_samp * 3e-3
t_grid = np.logspace(np.log10(t_us.min()), np.log10(t_us.max()), 100)
y_samp = np.array([abs_obs.max() * np.exp(-k * t_grid * 1e-6)
                   for k in k_obs_samp[:500]])
q = np.quantile(y_samp, [0.05, 0.5, 0.95], axis=0)
ax.fill_between(t_grid, q[0], q[2], color=WONG[2], alpha=0.25, label="HBMAE 90%")
ax.plot(t_grid, q[1], color=WONG[2], lw=1.2, label="HBMAE median")
ax.set_xscale("log"); ax.set_xlabel(r"$t$ [$\mu$s]")
ax.set_ylabel("absorbance [a.u.]")
ax.set_title(r"(n) Cr$^{3+}$ transient")
ax.legend(loc="upper right", fontsize=6, frameon=False); _despine(ax)

# --- Panel o: Kristoffersen-Metiu IE
ax = axes[14]
km = read_csv(VAL / "es_theory/kristoffersen_metiu_2018_jpcc/data/theory_vs_experiment.csv")
km_ie = km[km["quantity"].str.startswith("IE")]
labels = []
vals = []
sigs = []
for _, r in km_ie.iterrows():
    s = str(r["theory_value"])
    try:
        if "+/-" in s:
            v, sg = s.split("+/-")
            vals.append(float(v.strip()))
            sigs.append(float(sg.strip()))
        else:
            vals.append(float(s.strip()))
            sigs.append(0.5)
        labels.append(r["salt_composition"][:14])
    except ValueError:
        continue
y = np.arange(len(labels))
ax.errorbar(vals, y, xerr=sigs, fmt="o", color=WONG[2],
            mec="k", mew=0.4, ms=5, label="KM18 AIMD")
ax.axvline(2.9, color=WONG[6], ls="--", lw=1, label=r"Li WF 2.9 eV")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6)
ax.set_xlabel("IE [eV]")
ax.set_title(r"(o) KM18 bipolaron IE")
ax.legend(loc="upper right", fontsize=6, frameon=False); _despine(ax)

# --- Panel p: residual scatter (mini)
ax = axes[15]
metals_color = {
    "Cr3+": WONG[2], "Zn2+": WONG[1], "Nd3+": WONG[3], "Cf3+": WONG[6],
    "Cd2+": WONG[5], "Tl+": WONG[7], "Ag+": WONG[4], "Ba2+": WONG[0],
    "Sr2+": WONG[0], "Ca2+": WONG[0],
}
host_marker = {"LiCl-KCl": "o", "NaCl": "s", "KCl": "D",
               "KBr": "^", "H2O": "v"}
for _, r in lomo.iterrows():
    obs = float(r["log10_k_obs"]); pred = float(r["log10_k_pred_mean"])
    ax.scatter(pred, obs - pred,
               color=metals_color.get(r["metal"], WONG[0]),
               marker=host_marker.get(r["host"], "x"),
               s=26, edgecolors="k", linewidths=0.3)
ax.axhline(0, color="k", lw=0.6)
ax.axhline(0.3, color=WONG[6], lw=0.4, ls=":")
ax.axhline(-0.3, color=WONG[6], lw=0.4, ls=":")
ax.set_xlabel(r"predicted $\log_{10} k$")
ax.set_ylabel(r"residual $\log_{10} k$")
ax.set_title("(p) LOMO residuals")
_despine(ax)

fig.savefig(FIGS / "fig_master_validation.pdf", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "fig_master_validation.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Residuals all
# =============================================================================
fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
for _, r in lomo.iterrows():
    obs = float(r["log10_k_obs"]); pred = float(r["log10_k_pred_mean"])
    p5 = float(r["log10_k_pred_p05"]); p95 = float(r["log10_k_pred_p95"])
    color = metals_color.get(r["metal"], WONG[0])
    ax.errorbar(pred, obs - pred,
                yerr=[[max(pred-p5, 0)], [max(p95-pred, 0)]],
                fmt=host_marker.get(r["host"], "x"),
                color=color, ecolor=color,
                mec="k", mew=0.3, ms=6, alpha=0.85,
                elinewidth=0.5, capsize=2)
ax.axhline(0, color="k", lw=0.6)
ax.axhline(0.3, color=WONG[6], lw=0.6, ls="--", label=r"$\pm 0.3$ log units")
ax.axhline(-0.3, color=WONG[6], lw=0.6, ls="--")
ax.set_xlabel(r"predicted $\log_{10} k$ [M$^{-1}$ s$^{-1}$]")
ax.set_ylabel(r"$\log_{10} k_{\mathrm{obs}} - \log_{10} k_{\mathrm{pred}}$")
for m, c in metals_color.items():
    ax.scatter([], [], color=c, edgecolors="k", linewidths=0.3, s=30, label=m)
for h, mk in host_marker.items():
    ax.scatter([], [], color="grey", edgecolors="k", linewidths=0.3, s=30,
               marker=mk, label=h)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7,
          frameon=False, ncol=1)
_despine(ax)
fig.savefig(FIGS / "fig_residuals_all.pdf", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "fig_residuals_all.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# =============================================================================
# Data landscape
# =============================================================================
fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
for _, r in lomo.iterrows():
    obs = float(r["log10_k_obs"]); T = float(r["T_K"])
    ax.scatter(1000.0/T, obs, color=metals_color.get(r["metal"], WONG[0]),
               marker=host_marker.get(r["host"], "x"),
               s=36, edgecolors="k", linewidths=0.3, alpha=0.85)
ax.set_xlabel(r"$1000/T$ [K$^{-1}$]")
ax.set_ylabel(r"$\log_{10} k$ [M$^{-1}$ s$^{-1}$]")
for m, c in metals_color.items():
    ax.scatter([], [], color=c, edgecolors="k", linewidths=0.3, s=30, label=m)
for h, mk in host_marker.items():
    ax.scatter([], [], color="grey", edgecolors="k", linewidths=0.3, s=30,
               marker=mk, label=h)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7,
          frameon=False, ncol=1)
_despine(ax)
fig.savefig(FIGS / "fig_data_landscape.pdf", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "fig_data_landscape.png", dpi=300, bbox_inches="tight")
plt.close(fig)


print("Done; figures written to", FIGS)
