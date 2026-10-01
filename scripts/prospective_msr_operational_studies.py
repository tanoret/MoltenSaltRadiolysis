#!/usr/bin/env python3
"""Prospective operational screening studies for molten-salt radiolysis.

This script extends the static Article II calculations into transparent,
low-order operating scenarios.  The calculations are deliberately conditional:
none of the source, recovery, partition, removal, or redox-regeneration inputs
is promoted to a validated reactor-design parameter.

The script produces:
  * chloride cleanup-requirement and removal-outage studies;
  * chloride U(III) one-pass capacity and regeneration studies;
  * fluoride temperature/power operating envelopes and transients;
  * Monte Carlo sensitivity rankings that preserve the distinction between
    parametric uncertainty and structural scenario ranges.

All pressures are screening observables.  The 100 Pa line is a visualization
anchor, not a regulatory or design limit.
"""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr

R_GAS = 8.31446261815324  # J mol^-1 K^-1
N_A = 6.02214076e23       # mol^-1
E_100 = 100.0 * 1.602176634e-19  # J per 100 eV
HOURS_PER_YEAR = 365.25 * 24.0
SECONDS_PER_HOUR = 3600.0

# Article II / Phillips-screening anchors.
G_CL_REF = 0.2
ETA_TIMES_RECOVERY_REF = 1.670782136734233e-4
P_SCREEN = 100.0  # Pa; visualization anchor only

# Article II fluoride effective-pressure anchor.
T_FLIBE_REF_K = 873.15
T_FLIBE_REF_C = 600.0
P_FLIBE_REF = 31.0  # Pa
P_FLIBE_LOW = 19.0
P_FLIBE_HIGH = 44.0
K_REC_REF_H = 1.174  # h^-1 at 873.15 K
EA_REC_REF = 39.0e3  # J mol^-1


@dataclass(frozen=True)
class ChloridePlant:
    """Lumped chloride screening geometry and source assumptions."""

    rho_kg_m3: float = 2700.0
    liquid_volume_m3: float = 50.0
    headspace_volume_m3: float = 5.0
    temperature_K: float = 873.15
    henry_mol_m3_Pa: float = 2.0e-5
    dose_rate_kGy_h: float = 10.0
    G_cl_atoms_per_100eV: float = G_CL_REF

    @property
    def inventory_capacitance_mol_Pa(self) -> float:
        """Total equilibrium inventory per unit pressure.

        chi = V_g/(R T) + K_H V_l, so n_total = chi P.
        """
        return (
            self.headspace_volume_m3 / (R_GAS * self.temperature_K)
            + self.henry_mol_m3_Pa * self.liquid_volume_m3
        )

    @property
    def dose_rate_Gy_s(self) -> float:
        return self.dose_rate_kGy_h * 1.0e3 / SECONDS_PER_HOUR

    @property
    def primary_cl_atom_source_mol_h(self) -> float:
        source_mol_s = (
            self.G_cl_atoms_per_100eV
            * self.rho_kg_m3
            * self.liquid_volume_m3
            * self.dose_rate_Gy_s
            / (E_100 * N_A)
        )
        return source_mol_s * SECONDS_PER_HOUR


@dataclass(frozen=True)
class U3Inventory:
    """One-pass U(III) reducing-equivalent inventory used in Article II."""

    moles_U3: float = 312485.53307717235


def eta_upper_bound(
    recovery_fraction: float,
    G_cl: float = G_CL_REF,
    product_bound_at_ref_G: float = ETA_TIMES_RECOVERY_REF,
) -> float:
    """Conditional upper bound on stable-release fraction eta.

    The Phillips observation constrains eta * f_rec at a selected G value.
    Because that bound scales inversely with G, the released stable-product
    source is nearly independent of G while eta remains below one.
    """
    if not (0.0 < recovery_fraction <= 1.0):
        raise ValueError("recovery_fraction must be in (0, 1]")
    if G_cl <= 0.0:
        raise ValueError("G_cl must be positive")
    product_bound = product_bound_at_ref_G * (G_CL_REF / G_cl)
    return min(1.0, product_bound / recovery_fraction)


def chloride_stable_cl2_source_mol_h(
    plant: ChloridePlant,
    recovery_fraction: float,
) -> float:
    """Conditional stable Cl2 source at the Phillips-derived eta bound."""
    eta = eta_upper_bound(recovery_fraction, plant.G_cl_atoms_per_100eV)
    return 0.5 * eta * plant.primary_cl_atom_source_mol_h


def pressure_source_rate_Pa_h(
    plant: ChloridePlant,
    recovery_fraction: float,
) -> float:
    return (
        chloride_stable_cl2_source_mol_h(plant, recovery_fraction)
        / plant.inventory_capacitance_mol_Pa
    )


def required_removal_rate_h(
    pressure_source_Pa_h: float,
    target_pressure_Pa: float = P_SCREEN,
) -> float:
    if pressure_source_Pa_h < 0.0:
        raise ValueError("pressure_source_Pa_h must be nonnegative")
    if target_pressure_Pa <= 0.0:
        raise ValueError("target_pressure_Pa must be positive")
    return pressure_source_Pa_h / target_pressure_Pa


def maximum_half_life_h_for_target(
    pressure_source_Pa_h: float,
    target_pressure_Pa: float = P_SCREEN,
) -> float:
    k_req = required_removal_rate_h(pressure_source_Pa_h, target_pressure_Pa)
    return math.inf if k_req == 0.0 else math.log(2.0) / k_req


def simulate_linear_inventory(
    t_h: np.ndarray,
    source_mol_h: Callable[[float], float],
    removal_h: Callable[[float], float],
    n0_mol: float,
) -> np.ndarray:
    """Solve dn/dt = S(t) - k(t)n for a nonnegative inventory."""
    if t_h.ndim != 1 or len(t_h) < 2 or np.any(np.diff(t_h) <= 0.0):
        raise ValueError("t_h must be a strictly increasing one-dimensional grid")
    if n0_mol < 0.0:
        raise ValueError("n0_mol must be nonnegative")

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        n = max(float(y[0]), 0.0)
        return np.array([max(source_mol_h(t), 0.0) - max(removal_h(t), 0.0) * n])

    sol = solve_ivp(
        rhs,
        (float(t_h[0]), float(t_h[-1])),
        np.array([n0_mol]),
        t_eval=t_h,
        method="BDF",
        rtol=1.0e-9,
        atol=1.0e-13,
        max_step=max(0.05, (t_h[-1] - t_h[0]) / 1000.0),
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.maximum(sol.y[0], 0.0)


def periodic_power_fraction(t_h: float, on_h: float = 18.0, off_h: float = 6.0) -> float:
    period = on_h + off_h
    return 1.0 if (t_h % period) < on_h else 0.0


def k_rec_flibe_h(T_K: float, Ea_J_mol: float = EA_REC_REF) -> float:
    if T_K <= 0.0:
        raise ValueError("temperature must be positive")
    return K_REC_REF_H * math.exp(
        -Ea_J_mol / R_GAS * (1.0 / T_K - 1.0 / T_FLIBE_REF_K)
    )


def flibe_pressure_source_Pa_h(
    dose_fraction: float = 1.0,
    P_ref: float = P_FLIBE_REF,
) -> float:
    if dose_fraction < 0.0:
        raise ValueError("dose_fraction must be nonnegative")
    return K_REC_REF_H * P_ref * dose_fraction


def simulate_flibe_pressure(
    t_h: np.ndarray,
    dose_fraction: Callable[[float], float],
    temperature_K: Callable[[float], float],
    P0_Pa: float,
    P_ref: float = P_FLIBE_REF,
    Ea_J_mol: float = EA_REC_REF,
    extra_removal_h: Callable[[float], float] | None = None,
) -> np.ndarray:
    """Effective pressure balance dP/dt = source - removal*P.

    This is the dynamic counterpart of the Article II static closure.  It is
    not a gas-liquid mechanistic model; P is the effective pressure observable
    associated with the fixed reference geometry and partition assumptions.
    """
    if extra_removal_h is None:
        extra_removal_h = lambda _t: 0.0

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        P = max(float(y[0]), 0.0)
        source = flibe_pressure_source_Pa_h(max(dose_fraction(t), 0.0), P_ref)
        k_total = k_rec_flibe_h(temperature_K(t), Ea_J_mol) + max(extra_removal_h(t), 0.0)
        return np.array([source - k_total * P])

    sol = solve_ivp(
        rhs,
        (float(t_h[0]), float(t_h[-1])),
        np.array([max(P0_Pa, 0.0)]),
        t_eval=t_h,
        method="BDF",
        rtol=1.0e-9,
        atol=1.0e-11,
        max_step=max(0.01, (t_h[-1] - t_h[0]) / 2000.0),
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.maximum(sol.y[0], 0.0)


def time_to_U3_depletion_years(
    plant: ChloridePlant,
    inventory: U3Inventory,
    duty_fraction: float,
    regeneration_fraction_of_full_power_source: float,
) -> float:
    """Time to consume one-pass U(III) capacity in an equivalent balance.

    The source is expressed in chlorine-atom oxidizing equivalents.  The
    regeneration term is a prescribed fraction of the *full-power* source,
    which makes the roles of duty factor and continuous regeneration explicit.
    """
    if not (0.0 <= duty_fraction <= 1.0):
        raise ValueError("duty_fraction must be in [0, 1]")
    if regeneration_fraction_of_full_power_source < 0.0:
        raise ValueError("regeneration fraction must be nonnegative")
    source = plant.primary_cl_atom_source_mol_h
    net = duty_fraction * source - regeneration_fraction_of_full_power_source * source
    if net <= 0.0:
        return math.inf
    return inventory.moles_U3 / net / HOURS_PER_YEAR


def required_regeneration_fraction_for_horizon(
    plant: ChloridePlant,
    inventory: U3Inventory,
    horizon_years: float,
    duty_fraction: float = 1.0,
) -> float:
    if horizon_years <= 0.0:
        raise ValueError("horizon_years must be positive")
    allowable_net = inventory.moles_U3 / (horizon_years * HOURS_PER_YEAR)
    required = duty_fraction - allowable_net / plant.primary_cl_atom_source_mol_h
    return max(0.0, required)


def _loguniform(rng: np.random.Generator, low: float, high: float, size: int) -> np.ndarray:
    return np.exp(rng.uniform(np.log(low), np.log(high), size=size))


def _spearman_ranking(inputs: dict[str, np.ndarray], output: np.ndarray) -> list[tuple[str, float]]:
    ranking: list[tuple[str, float]] = []
    for name, values in inputs.items():
        rho, _ = spearmanr(values, output)
        ranking.append((name, float(rho)))
    return sorted(ranking, key=lambda item: abs(item[1]), reverse=True)


def chloride_uncertainty_study(n: int = 30000, seed: int = 20260818) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    recovery = _loguniform(rng, 0.01, 1.0, n)
    dose = _loguniform(rng, 3.0, 30.0, n)
    removal_half_life = _loguniform(rng, 1.0, 100.0, n)
    G = _loguniform(rng, 0.05, 1.0, n)
    KH = _loguniform(rng, 5.0e-6, 5.0e-5, n)
    headspace_ratio = rng.uniform(0.05, 0.20, n)
    temperature = rng.uniform(823.15, 973.15, n)

    rho_salt = 2700.0
    V_l = 50.0
    V_g = headspace_ratio * V_l
    dose_Gy_s = dose * 1.0e3 / SECONDS_PER_HOUR
    atom_source_h = G * rho_salt * V_l * dose_Gy_s / (E_100 * N_A) * SECONDS_PER_HOUR
    eta = np.minimum(1.0, ETA_TIMES_RECOVERY_REF * (G_CL_REF / G) / recovery)
    stable_cl2_source_h = 0.5 * eta * atom_source_h
    chi = V_g / (R_GAS * temperature) + KH * V_l
    k = np.log(2.0) / removal_half_life
    P_ss = stable_cl2_source_h / (chi * k)

    inputs = {
        "analytical recovery": recovery,
        "dose rate": dose,
        "cleanup half-life": removal_half_life,
        "G(Cl)": G,
        "Henry coefficient": KH,
        "headspace/liquid ratio": headspace_ratio,
        "temperature": temperature,
    }
    ranking = _spearman_ranking(inputs, np.log10(P_ss))
    return {
        "inputs": inputs,
        "P_ss": P_ss,
        "ranking": ranking,
        "quantiles": np.quantile(P_ss, [0.05, 0.5, 0.95]),
        "exceedance": float(np.mean(P_ss > P_SCREEN)),
    }


def flibe_uncertainty_study(n: int = 30000, seed: int = 20260819) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    P_ref = rng.uniform(P_FLIBE_LOW, P_FLIBE_HIGH, n)
    Ea = rng.uniform(25.0e3, 55.0e3, n)  # structural sweep, not posterior
    temperature_C = rng.uniform(500.0, 600.0, n)  # cool-side envelope; monotonic Arrhenius effect
    temperature_K = temperature_C + 273.15
    dose_fraction = _loguniform(rng, 0.1, 1.0, n)
    extra_half_life = _loguniform(rng, 0.5, 720.0, n)
    k_extra = np.log(2.0) / extra_half_life
    source = K_REC_REF_H * P_ref * dose_fraction
    k_rec = K_REC_REF_H * np.exp(-Ea / R_GAS * (1.0 / temperature_K - 1.0 / T_FLIBE_REF_K))
    P_ss = source / (k_rec + k_extra)

    inputs = {
        "600 C pressure anchor": P_ref,
        "removal activation energy": Ea,
        "temperature": temperature_C,
        "dose fraction": dose_fraction,
        "extra-cleanup half-life": extra_half_life,
    }
    ranking = _spearman_ranking(inputs, np.log10(P_ss))
    return {
        "inputs": inputs,
        "P_ss": P_ss,
        "ranking": ranking,
        "quantiles": np.quantile(P_ss, [0.05, 0.5, 0.95]),
        "exceedance": float(np.mean(P_ss > P_SCREEN)),
    }


def _write_rows(path: Path, fieldnames: Iterable[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def run_studies(output_root: Path, samples: int = 30000) -> dict[str, object]:
    fig_dir = output_root / "figures"
    result_dir = output_root / "results"
    fig_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    plant = ChloridePlant()
    inventory = U3Inventory()
    chi = plant.inventory_capacitance_mol_Pa

    # ------------------------------------------------------------------
    # Study 1: cleanup requirement and removal-loss transients.
    # ------------------------------------------------------------------
    recovery_grid = np.logspace(-2, 0, 300)
    dose_cases = [3.0, 10.0, 30.0]
    half_life_curves: dict[float, np.ndarray] = {}
    for dose in dose_cases:
        p = ChloridePlant(dose_rate_kGy_h=dose)
        half_life_curves[dose] = np.array([
            maximum_half_life_h_for_target(pressure_source_rate_Pa_h(p, f))
            for f in recovery_grid
        ])

    # Seven-day transient, initialized at the nominal steady state.  The cleanup
    # system is unavailable from day 2 through day 4 and then recovers.
    t_h = np.linspace(0.0, 7.0 * 24.0, 1601)
    nominal_half_life_h = 24.0
    k_nom = math.log(2.0) / nominal_half_life_h
    transient_pressures: dict[float, np.ndarray] = {}
    for f_rec in [1.0, 0.1, 0.01]:
        source = chloride_stable_cl2_source_mol_h(plant, f_rec)
        n0 = source / k_nom

        def source_fn(_t: float, source_value: float = source) -> float:
            return source_value

        def removal_fn(t: float) -> float:
            return 0.0 if 48.0 <= t < 96.0 else k_nom

        n = simulate_linear_inventory(t_h, source_fn, removal_fn, n0)
        transient_pressures[f_rec] = n / chi

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.45), constrained_layout=True)
    ax = axes[0]
    for dose in dose_cases:
        ax.plot(recovery_grid, half_life_curves[dose], label=f"{dose:g} kGy h$^{{-1}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("effective analytical recovery fraction")
    ax.set_ylabel("maximum cleanup half-life for 100 Pa [h]")
    ax.set_title("(a) Conditional cleanup requirement")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    for f_rec, P in transient_pressures.items():
        ax.plot(t_h / 24.0, P, label=fr"$f_{{\rm rec}}={f_rec:g}$")
    ax.axvspan(2.0, 4.0, alpha=0.12, label="cleanup unavailable")
    ax.axhline(P_SCREEN, ls=":", lw=1.2, label="100 Pa screening anchor")
    ax.set_yscale("log")
    ax.set_xlabel("time [days]")
    ax.set_ylabel(r"conditional $P_{\mathrm{Cl}_2}$ [Pa]")
    ax.set_title("(b) Forty-eight-hour cleanup interruption")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(frameon=False, fontsize=6.8)
    fig.savefig(fig_dir / "fig_prospective_chloride_cleanup.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_prospective_chloride_cleanup.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Study 2: duty cycle and U(III) regeneration requirements.
    # ------------------------------------------------------------------
    regen_grid = np.linspace(0.0, 1.05, 500)
    depletion_curves: dict[float, np.ndarray] = {}
    for dose in dose_cases:
        p = ChloridePlant(dose_rate_kGy_h=dose)
        vals = []
        for r in regen_grid:
            t_dep = time_to_U3_depletion_years(p, inventory, 1.0, r)
            vals.append(min(t_dep, 200.0) if math.isfinite(t_dep) else 200.0)
        depletion_curves[dose] = np.array(vals)

    duty_grid = np.linspace(0.1, 1.0, 181)
    regen2_grid = np.linspace(0.0, 1.1, 221)
    margin = np.empty((len(regen2_grid), len(duty_grid)))
    horizon = 60.0
    for i, regen in enumerate(regen2_grid):
        for j, duty in enumerate(duty_grid):
            t_dep = time_to_U3_depletion_years(plant, inventory, duty, regen)
            margin[i, j] = 2.0 if not math.isfinite(t_dep) else min(t_dep / horizon, 2.0)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.45), constrained_layout=True)
    ax = axes[0]
    for dose in dose_cases:
        ax.plot(regen_grid, depletion_curves[dose], label=f"{dose:g} kGy h$^{{-1}}$")
    ax.axhline(60.0, ls=":", lw=1.2, label="60-y horizon")
    ax.set_yscale("log")
    ax.set_ylim(0.3, 220.0)
    ax.set_xlabel("regeneration rate / full-power oxidant source")
    ax.set_ylabel("one-pass-equivalent depletion time [y]")
    ax.set_title("(a) Regeneration needed for long operation")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    mesh = ax.pcolormesh(duty_grid, regen2_grid, margin, shading="auto")
    contour = ax.contour(duty_grid, regen2_grid, margin, levels=[1.0], linewidths=1.4)
    ax.clabel(contour, fmt={1.0: "60-y boundary"}, fontsize=7)
    ax.plot(duty_grid, duty_grid, ls="--", lw=1.0, label="regeneration = average source")
    ax.set_xlabel("irradiation duty fraction")
    ax.set_ylabel("regeneration / full-power source")
    ax.set_title("(b) Sixty-year U(III)-capacity margin")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("depletion time / 60 y (capped at 2)")
    fig.savefig(fig_dir / "fig_prospective_chloride_buffer.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_prospective_chloride_buffer.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Study 3: FLiBe temperature/power envelope and prescribed transient.
    # ------------------------------------------------------------------
    T_C = np.linspace(500.0, 700.0, 241)
    dose_frac = np.linspace(0.05, 1.0, 191)
    TT, DD = np.meshgrid(T_C, dose_frac)
    k_grid = K_REC_REF_H * np.exp(
        -EA_REC_REF / R_GAS * (1.0 / (TT + 273.15) - 1.0 / T_FLIBE_REF_K)
    )
    P_grid = K_REC_REF_H * P_FLIBE_REF * DD / k_grid

    t_f = np.linspace(0.0, 48.0, 2401)

    def schedule_T(t: float) -> float:
        if t < 12.0:
            return 600.0 + 273.15
        if t < 18.0:
            return (600.0 - (100.0 / 6.0) * (t - 12.0)) + 273.15
        if t < 30.0:
            return 500.0 + 273.15
        if t < 36.0:
            return (500.0 + (100.0 / 6.0) * (t - 30.0)) + 273.15
        return 600.0 + 273.15

    def schedule_dose(t: float) -> float:
        if t < 12.0:
            return 1.0
        if t < 18.0:
            return 1.0 - (0.70 / 6.0) * (t - 12.0)
        if t < 30.0:
            return 0.30
        if t < 36.0:
            return 0.30 + (0.70 / 6.0) * (t - 30.0)
        return 1.0

    P_trans = simulate_flibe_pressure(t_f, schedule_dose, schedule_T, P_FLIBE_REF)
    P_const = np.full_like(t_f, P_FLIBE_REF)
    T_sched_C = np.array([schedule_T(t) - 273.15 for t in t_f])
    D_sched = np.array([schedule_dose(t) for t in t_f])

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.45), constrained_layout=True)
    ax = axes[0]
    mesh = ax.pcolormesh(T_C, dose_frac, P_grid, shading="auto")
    levels = [25.0, 50.0, 75.0, 100.0]
    valid_levels = [level for level in levels if P_grid.min() <= level <= P_grid.max()]
    if valid_levels:
        contours = ax.contour(T_C, dose_frac, P_grid, levels=valid_levels, linewidths=1.0)
        ax.clabel(contours, fmt="%g Pa", fontsize=7)
    ax.set_xlabel(r"salt temperature [$^\circ$C]")
    ax.set_ylabel("dose-rate fraction of reference")
    ax.set_title("(a) Conditional FLiBe pressure envelope")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"$P_{\mathrm{F}_2}^{\rm ss}$ [Pa]")

    ax = axes[1]
    ax.plot(t_f, P_trans, label="temperature-power maneuver")
    ax.plot(t_f, P_const, ls="--", label="600 C, full-power reference")
    ax.axhline(P_SCREEN, ls=":", lw=1.2, label="100 Pa screening anchor")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"conditional $P_{\mathrm{F}_2}$ [Pa]")
    ax.set_title("(b) Coupled cooldown and power reduction")
    ax.grid(True, alpha=0.2)
    ax2 = ax.twinx()
    ax2.plot(t_f, T_sched_C / 700.0, alpha=0.35, label="temperature / 700 C")
    ax2.plot(t_f, D_sched, alpha=0.35, ls=":", label="dose fraction")
    ax2.set_ylabel("normalized schedule")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, frameon=False, fontsize=6.3, loc="upper right")
    fig.savefig(fig_dir / "fig_prospective_flibe_operation.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_prospective_flibe_operation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Study 4: uncertainty/sensitivity hierarchy.
    # ------------------------------------------------------------------
    chloride_mc = chloride_uncertainty_study(samples)
    flibe_mc = flibe_uncertainty_study(samples)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), constrained_layout=True)
    for ax, title, study in [
        (axes[0], "(a) Chloride steady-pressure ranking", chloride_mc),
        (axes[1], "(b) Cool-side FLiBe pressure ranking", flibe_mc),
    ]:
        ranking = study["ranking"]
        labels = [name for name, _rho in ranking]
        values = [rho for _name, rho in ranking]
        y = np.arange(len(labels))
        ax.barh(y, values)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.axvline(0.0, lw=0.8)
        ax.set_xlim(-1.0, 1.0)
        ax.set_xlabel("Spearman rank correlation with log pressure")
        ax.set_title(title)
        ax.grid(True, axis="x", alpha=0.2)
    fig.savefig(fig_dir / "fig_prospective_uncertainty_ranking.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_prospective_uncertainty_ranking.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Auditable numerical summaries.
    # ------------------------------------------------------------------
    rows: list[dict[str, object]] = []
    for dose in dose_cases:
        p = ChloridePlant(dose_rate_kGy_h=dose)
        for f_rec in [1.0, 0.1, 0.01]:
            a = pressure_source_rate_Pa_h(p, f_rec)
            hmax = maximum_half_life_h_for_target(a)
            rows.append({
                "study": "chloride_cleanup_requirement",
                "case": f"dose={dose:g}_frec={f_rec:g}",
                "metric": "maximum_cleanup_half_life_for_100Pa",
                "value": hmax,
                "units_or_note": "h; conditional Phillips-bound scenario",
            })
    for f_rec, P in transient_pressures.items():
        rows.extend([
            {
                "study": "chloride_cleanup_outage",
                "case": f"frec={f_rec:g}",
                "metric": "peak_pressure",
                "value": float(np.max(P)),
                "units_or_note": "Pa; 48-h removal outage after nominal steady state",
            },
            {
                "study": "chloride_cleanup_outage",
                "case": f"frec={f_rec:g}",
                "metric": "pressure_at_day7",
                "value": float(P[-1]),
                "units_or_note": "Pa",
            },
        ])
    for regen in [0.0, 0.90, 0.95, 0.9788, 0.99, 1.0]:
        t_dep = time_to_U3_depletion_years(plant, inventory, 1.0, regen)
        rows.append({
            "study": "chloride_U3_capacity",
            "case": f"regeneration_fraction={regen:g}",
            "metric": "depletion_time",
            "value": t_dep if math.isfinite(t_dep) else "inf",
            "units_or_note": "y; oxidizing-equivalent balance",
        })
    rows.append({
        "study": "chloride_U3_capacity",
        "case": "60y_full_power",
        "metric": "minimum_regeneration_fraction",
        "value": required_regeneration_fraction_for_horizon(plant, inventory, 60.0, 1.0),
        "units_or_note": "fraction of full-power primary oxidant source",
    })
    rows.extend([
        {
            "study": "flibe_operating_envelope",
            "case": "temperature_power_maneuver",
            "metric": "peak_pressure",
            "value": float(np.max(P_trans)),
            "units_or_note": "Pa",
        },
        {
            "study": "flibe_operating_envelope",
            "case": "temperature_power_maneuver",
            "metric": "minimum_pressure",
            "value": float(np.min(P_trans)),
            "units_or_note": "Pa",
        },
        {
            "study": "flibe_operating_envelope",
            "case": "500C_full_reference_dose",
            "metric": "steady_pressure",
            "value": float(P_FLIBE_REF * K_REC_REF_H / k_rec_flibe_h(773.15)),
            "units_or_note": "Pa; conditional effective model",
        },
        {
            "study": "flibe_operating_envelope",
            "case": "700C_full_reference_dose",
            "metric": "steady_pressure",
            "value": float(P_FLIBE_REF * K_REC_REF_H / k_rec_flibe_h(973.15)),
            "units_or_note": "Pa; conditional effective model",
        },
    ])
    for label, study in [("chloride", chloride_mc), ("flibe", flibe_mc)]:
        q05, q50, q95 = study["quantiles"]
        for qname, qval in [("P05", q05), ("P50", q50), ("P95", q95)]:
            rows.append({
                "study": f"{label}_uncertainty_sweep",
                "case": f"n={samples}",
                "metric": qname,
                "value": float(qval),
                "units_or_note": "Pa; scenario sweep, not posterior interval",
            })
        rows.append({
            "study": f"{label}_uncertainty_sweep",
            "case": f"n={samples}",
            "metric": "fraction_above_100Pa",
            "value": float(study["exceedance"]),
            "units_or_note": "scenario fraction",
        })
    _write_rows(
        result_dir / "prospective_operational_case_summary.csv",
        ["study", "case", "metric", "value", "units_or_note"],
        rows,
    )

    ranking_rows: list[dict[str, object]] = []
    for label, study in [("chloride", chloride_mc), ("flibe", flibe_mc)]:
        for rank, (name, rho) in enumerate(study["ranking"], start=1):
            ranking_rows.append({
                "system": label,
                "rank": rank,
                "input": name,
                "spearman_rho": rho,
                "note": "correlation with log10 steady pressure in scenario ensemble",
            })
    _write_rows(
        result_dir / "prospective_uncertainty_ranking.csv",
        ["system", "rank", "input", "spearman_rho", "note"],
        ranking_rows,
    )

    envelope_rows: list[dict[str, object]] = []
    for T in [500.0, 550.0, 600.0, 650.0, 700.0]:
        for d in [0.1, 0.3, 0.5, 1.0]:
            pss = flibe_pressure_source_Pa_h(d) / k_rec_flibe_h(T + 273.15)
            envelope_rows.append({
                "temperature_C": T,
                "dose_fraction": d,
                "steady_pressure_Pa": pss,
                "note": "conditional effective production-removal model",
            })
    _write_rows(
        result_dir / "prospective_flibe_operating_envelope.csv",
        ["temperature_C", "dose_fraction", "steady_pressure_Pa", "note"],
        envelope_rows,
    )

    return {
        "plant": plant,
        "inventory": inventory,
        "transient_pressures": transient_pressures,
        "flibe_transient": P_trans,
        "chloride_mc": chloride_mc,
        "flibe_mc": flibe_mc,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="directory containing figures/ and results/ subdirectories",
    )
    parser.add_argument("--samples", type=int, default=30000)
    args = parser.parse_args()
    if args.samples < 1000:
        raise SystemExit("--samples must be at least 1000")

    result = run_studies(args.output_root, samples=args.samples)
    plant: ChloridePlant = result["plant"]
    print(f"Chloride inventory capacitance: {plant.inventory_capacitance_mol_Pa:.6e} mol/Pa")
    for f_rec in [1.0, 0.1, 0.01]:
        a = pressure_source_rate_Pa_h(plant, f_rec)
        hmax = maximum_half_life_h_for_target(a)
        peak = float(np.max(result["transient_pressures"][f_rec]))
        print(
            f"f_rec={f_rec:g}: source slope={a:.6g} Pa/h; "
            f"max half-life for 100 Pa={hmax:.4g} h; outage peak={peak:.6g} Pa"
        )
    req = required_regeneration_fraction_for_horizon(plant, result["inventory"], 60.0, 1.0)
    print(f"Minimum full-power-equivalent regeneration fraction for 60 y: {req:.6f}")
    print(
        "FLiBe transient: "
        f"min={np.min(result['flibe_transient']):.6g} Pa, "
        f"max={np.max(result['flibe_transient']):.6g} Pa"
    )
    for label in ["chloride_mc", "flibe_mc"]:
        study = result[label]
        q = study["quantiles"]
        print(
            f"{label}: P05/P50/P95={q[0]:.6g}/{q[1]:.6g}/{q[2]:.6g} Pa; "
            f"scenario fraction >100 Pa={study['exceedance']:.4f}"
        )
        print("  ranking: " + ", ".join(f"{n}={r:+.3f}" for n, r in study["ranking"]))


if __name__ == "__main__":
    main()
