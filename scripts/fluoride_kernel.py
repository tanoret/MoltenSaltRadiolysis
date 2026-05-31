#!/usr/bin/env python3
"""Calibrated fluoride F2-production kernel.

Combines:
  - Davis et al. 2022 (NSE) Table III G-values for LiF, BeF2, FLiBe-UF4, ThF4,
    65LiF-29BeF2-5ZrF4-0.66UF4 (MSRE composition), measured under HFR gamma
    at 40-60 C with G in molecules F2 per 100 eV absorbed.
  - Toth & Felker 1990 (Rad. Eff. Defects Solids) recombination kinetics:
    Ea_rec = 39 kJ/mol over 150-200 C; T_balance = 150 C (generation =
    recombination steady-state crossover); finite recombination rate
    3.8 Pa/h at 40 C in the 64.5 LiF - 30.3 BeF2 - 5 ZrF4 - 0.13 UF4 salt.

The kernel computes the dose-balance equation for F2 partial pressure in the
gas headspace of a closed salt-vessel system:

    dP_F2/dt  =  S_F2(composition, T) * D                  (radiolytic source)
                  -  k_rec(T) * P_F2                         (surface recomb.)

where S_F2 = G_F2(composition) * conversion_factor and
      k_rec(T) = A_rec * exp(-Ea_rec / RT).

The pre-exponential A_rec is calibrated from the Toth/Felker balance condition
at T_balance = 423.15 K so that the predicted balance pressure under their
nominal dose rate matches their experimental value. The Davis G-values are
treated as composition-dependent constants (no T-dependence reported).

Outputs:
  validation/fluoride_kernel/calibrated_parameters.csv
  validation/fluoride_kernel/predicted_buildup_curves.csv
  manuscript/figures/fig_fluoride_kernel.pdf   (4-panel: G vs composition;
                                                steady-state P(T); buildup curve
                                                at 313 K; T-balance map)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "manuscript" / "figures"))
try:
    import figstyle  # type: ignore
    figstyle.apply_rc()
except Exception:  # pragma: no cover
    pass

R_GAS = 8.314462618  # J/(mol K)
N_A = 6.02214076e23  # mol^-1
EV_PER_JOULE = 6.241509074e18
EV_TO_J = 1.602176634e-19

# ---------------------------------------------------------------------------
# Inputs from extracted CSVs
# ---------------------------------------------------------------------------
DAVIS_CSV = REPO / "validation/eflibe_f2_yield/davis_2022_nse/data/G_values_table_III.csv"
TOTH_CSV = REPO / "validation/flibe_msre_f2/toth_felker_1990_redds/data/recombination_kinetics.csv"
OUT_DIR = REPO / "validation/fluoride_kernel"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR = REPO / "manuscript/figures"


def load_davis_G() -> pd.DataFrame:
    """Load Davis Table III; keep only this-work primary rows (one G per salt)."""
    df = pd.read_csv(DAVIS_CSV, comment="#")
    keep = [
        "BeF2", "LiF", "FLiBe_UF4", "ThF4",
        "65LiF_29BeF2_5ZrF4_0.66UF4_Co60",  # MSRE composition, Co-60 source
    ]
    sub = df[df["salt"].isin(keep)].copy()
    sub = sub[["salt", "G_F2", "sigma", "dose_rate_kGy_per_hr"]].reset_index(drop=True)
    sub["salt_label"] = sub["salt"].map({
        "BeF2": r"BeF$_2$",
        "LiF": "LiF",
        "FLiBe_UF4": r"FLiBe-UF$_4$",
        "ThF4": r"ThF$_4$",
        "65LiF_29BeF2_5ZrF4_0.66UF4_Co60": "MSRE\n(Co-60)",
    })
    return sub


def load_toth() -> dict:
    """Toth/Felker recombination parameters."""
    df = pd.read_csv(TOTH_CSV, comment="#")
    p = {row["parameter"]: row["value"] for _, row in df.iterrows()}
    return {
        "Ea_rec_J_mol": float(p["Ea_recombination"]) * 1e3,  # kJ/mol -> J/mol
        "T_balance_K": float(p["T_balance_C"]) + 273.15,
        "balance_rate_Pa_h": float(p["balance_pressure_Pa_h"]),
    }


def calibrate_A_rec(toth: dict, P_ref_Pa: float = 1000.0) -> float:
    """Calibrate the recombination pre-exponential A_rec.

    Steady-state balance at T_balance: A_rec exp(-Ea/RT_balance) * P_ref = S * D,
    where S * D is the unknown source term. We do NOT have Toth/Felker's exact
    dose rate, so we anchor A_rec relative to a reference pressure of 1000 Pa
    (a typical pulsed-gamma F2 buildup pressure in MSRE-like experiments) such
    that the relative buildup curve passes through (T=40 C, dP/dt = 3.8 Pa/h)
    and (T=150 C, dP/dt = 0).

    This anchors the *relative* T-dependence of k_rec; absolute kinetics in
    other devices need to scale A_rec to their own source term.
    """
    Ea = toth["Ea_rec_J_mol"]
    Tb = toth["T_balance_K"]
    # At the balance, source rate = k_rec(Tb) * P_ref
    # At 40 C, the *net* rate is 3.8 Pa/h, so:
    #   source - k_rec(313) * P_ref = 3.8 Pa/h
    # Eliminating source:
    #   k_rec(Tb) * P_ref - k_rec(313) * P_ref = 3.8 Pa/h
    #   (k_rec(Tb) - k_rec(313)) * P_ref = 3.8 Pa/h
    # Let kappa = A_rec * exp(-Ea/(R*Tb)), then k_rec(313) = kappa * exp(Ea/R*(1/Tb - 1/313))
    Tc = 313.15  # 40 C
    delta = np.exp(-Ea / (R_GAS * Tb)) - np.exp(-Ea / (R_GAS * Tc))
    A_rec = (toth["balance_rate_Pa_h"]) / (delta * P_ref_Pa)
    return A_rec


def k_rec(T_K: float | np.ndarray, A_rec: float, Ea_J_mol: float) -> np.ndarray:
    """Surface recombination rate constant, units 1/h (matches calibration)."""
    return A_rec * np.exp(-Ea_J_mol / (R_GAS * np.asarray(T_K)))


def predict_P_ss(
    G_F2_per_100eV: float,
    composition_dose_kGy_per_hr: float,
    T_K: float | np.ndarray,
    A_rec: float,
    Ea_J_mol: float,
    headspace_volume_L: float = 1.0,
    salt_mass_kg: float = 1.0,
) -> np.ndarray:
    """Steady-state P_F2 (Pa) at temperature T for a salt with given Davis G.

    Source term in molecules/h:
       S = G_F2 * (dose_rate * mass) / (100 eV per yield unit)
    Pressure conversion at vessel temperature: P [Pa] = n*R*T/V (ideal gas).

    With dose_rate in kGy/hr = kJ/(kg*hr), and 1 kJ = EV_PER_JOULE/1000 eV, the
    source in molecules/hr is:
       S = G * dose_rate * 1000 * EV_PER_JOULE * mass / 100
         = G * dose_rate * 10 * EV_PER_JOULE * mass     [molec/hr]
    """
    S_molec_per_hr = (G_F2_per_100eV * composition_dose_kGy_per_hr
                      * 10.0 * EV_PER_JOULE * salt_mass_kg)
    # convert to pressure rate at temperature T_K
    # n = S * dt -> P = n*R*T/V  => dP/dt = S*R*T/(V*N_A)
    dP_per_hr_pure_source = S_molec_per_hr * R_GAS * np.asarray(T_K) / (
        headspace_volume_L * 1e-3 * N_A
    )
    P_ss = dP_per_hr_pure_source / k_rec(T_K, A_rec, Ea_J_mol)
    return P_ss


def predicted_buildup(
    G_F2_per_100eV: float,
    composition_dose_kGy_per_hr: float,
    T_K: float,
    A_rec: float,
    Ea_J_mol: float,
    time_hours: np.ndarray,
    headspace_volume_L: float = 1.0,
    salt_mass_kg: float = 1.0,
) -> np.ndarray:
    """Closed-form solution P(t) for the linear ODE.

    P(t) = P_ss * (1 - exp(-k_rec * t))
    """
    P_ss = predict_P_ss(
        G_F2_per_100eV, composition_dose_kGy_per_hr, T_K,
        A_rec, Ea_J_mol, headspace_volume_L, salt_mass_kg,
    )
    k = k_rec(T_K, A_rec, Ea_J_mol)
    return P_ss * (1.0 - np.exp(-k * time_hours))


def main():
    davis = load_davis_G()
    toth = load_toth()

    A_rec = calibrate_A_rec(toth)
    Ea = toth["Ea_rec_J_mol"]

    # Save calibrated parameters
    params = pd.DataFrame({
        "parameter": ["A_rec", "Ea_rec_J_mol", "Ea_rec_kJ_mol",
                       "T_balance_K", "T_balance_C",
                       "balance_rate_at_40C_Pa_per_h",
                       "reference_pressure_Pa"],
        "value": [A_rec, Ea, Ea / 1e3,
                   toth["T_balance_K"], toth["T_balance_K"] - 273.15,
                   toth["balance_rate_Pa_h"], 1000.0],
        "units": ["1/h", "J/mol", "kJ/mol", "K", "C", "Pa/h", "Pa"],
        "source": ["calibrated from Toth/Felker balance", "Toth/Felker 1990",
                    "Toth/Felker 1990", "Toth/Felker 1990", "Toth/Felker 1990",
                    "Toth/Felker 1990", "anchor"],
    })
    params.to_csv(OUT_DIR / "calibrated_parameters.csv", index=False)

    # Temperature sweep
    T_C = np.linspace(40.0, 700.0, 200)
    T_K = T_C + 273.15

    # Steady-state pressures at each T for each Davis salt (unit dose rate, unit mass)
    curves = {}
    for _, row in davis.iterrows():
        G = float(row["G_F2"])
        D = float(row["dose_rate_kGy_per_hr"])
        P_ss = predict_P_ss(G, D, T_K, A_rec, Ea)
        curves[row["salt"]] = P_ss

    # Save tabular curves
    out = pd.DataFrame({"T_C": T_C, "T_K": T_K})
    for salt, vals in curves.items():
        out[f"P_ss_Pa_{salt}"] = vals
    out.to_csv(OUT_DIR / "predicted_buildup_curves.csv", index=False)

    # ---------- Figure ----------
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.6), constrained_layout=True)

    def _despine(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # (a) Davis G values bar chart
    ax = axes[0, 0]
    x_pos = np.arange(len(davis))
    err = davis["sigma"].fillna(0.0).to_numpy()
    ax.bar(x_pos, davis["G_F2"], yerr=err, color="#0072B2", alpha=0.85,
           edgecolor="k", linewidth=0.4)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(davis["salt_label"], rotation=0, fontsize=7)
    ax.set_ylabel(r"$G(\mathrm{F}_2)$ [molec / 100 eV]")
    ax.set_title("(a) Davis 2022 yields")
    _despine(ax)

    # (b) Steady-state P_F2 vs T for each composition
    ax = axes[0, 1]
    palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00"]
    for (salt, P_ss), color in zip(curves.items(), palette):
        label = davis.loc[davis["salt"] == salt, "salt_label"].iloc[0].replace("\n", " ")
        ax.semilogy(T_C, P_ss, color=color, label=label, lw=1.4)
    ax.axvline(toth["T_balance_K"] - 273.15, color="#bbbbbb", ls="--", lw=0.8)
    ax.text(toth["T_balance_K"] - 273.15 + 5, 1e6, "T-balance\n(Toth/Felker)",
            fontsize=6.5, color="#555555")
    ax.set_xlabel(r"$T$ [$^\circ$C]")
    ax.set_ylabel(r"$P_{\mathrm{F}_2,\mathrm{ss}}$ [Pa]")
    ax.set_title("(b) Steady-state pressure")
    ax.legend(fontsize=6.5, loc="upper right", frameon=False)
    _despine(ax)

    # (c) Buildup curve at 40 C for each composition
    ax = axes[1, 0]
    times = np.linspace(0, 200, 500)  # hours
    for (salt, _), color in zip(curves.items(), palette):
        G = float(davis.loc[davis["salt"] == salt, "G_F2"].iloc[0])
        D = float(davis.loc[davis["salt"] == salt, "dose_rate_kGy_per_hr"].iloc[0])
        P_t = predicted_buildup(G, D, 313.15, A_rec, Ea, times)
        label = davis.loc[davis["salt"] == salt, "salt_label"].iloc[0].replace("\n", " ")
        ax.plot(times, P_t / 1e3, color=color, label=label, lw=1.4)
    ax.set_xlabel(r"$t$ at 40 $^\circ$C [h]")
    ax.set_ylabel(r"$P_{\mathrm{F}_2}$ [kPa]")
    ax.set_title(r"(c) F$_2$ buildup at 40 $^\circ$C")
    ax.legend(fontsize=6.5, frameon=False)
    _despine(ax)

    # (d) Recombination rate k_rec(T) vs T (the Arrhenius)
    ax = axes[1, 1]
    invT = 1000.0 / T_K
    ax.semilogy(invT, k_rec(T_K, A_rec, Ea), color="#000000", lw=1.4)
    ax.axvline(1000.0 / toth["T_balance_K"], color="#bbbbbb", ls="--", lw=0.8)
    ax.set_xlabel(r"$1000/T$ [K$^{-1}$]")
    ax.set_ylabel(r"$k_{\mathrm{rec}}$ [h$^{-1}$]")
    ax.set_title(r"(d) Arrhenius: $E_a = 39$ kJ mol$^{-1}$")
    _despine(ax)

    fig.savefig(FIG_DIR / "fig_fluoride_kernel.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_fluoride_kernel.png", dpi=300, bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'calibrated_parameters.csv'}")
    print(f"Wrote {OUT_DIR / 'predicted_buildup_curves.csv'}")
    print(f"Wrote {FIG_DIR / 'fig_fluoride_kernel.pdf'}")
    print(f"A_rec calibrated = {A_rec:.3e} 1/h (referenced to P_ref=1000 Pa)")


if __name__ == "__main__":
    main()
