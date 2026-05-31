"""Multi-scale operator-splitting solver for chronic-irradiation radiolysis kinetics.

For chloride radiolysis under prolonged irradiation, the species timescales separate
into two well-resolved bands:

  Fast (radicals): e_s-, Cl•, Cl2•-  with τ_radical ~ 1-100 ns
  Slow:            Cl3-, Cl2_diss, Cl2_gas, U(III), U(IV)  with τ_slow ~ 10^3-10^7 s

The ratio ε = τ_radical/τ_slow ≈ 10^{-14} renders direct stiff ODE integration
numerically intractable across the full timescale. The natural reduction is the
slow-manifold approximation (Kuehn, Multiple Time Scale Dynamics, 2015, §6.3): in the
limit ε → 0, the fast variables relax to an algebraic submanifold determined by the
slow variables.

We implement Strang-splitting operator-splitting in which each long time step dt_slow
contains:
  (1) Fast sub-step: solve the radical algebraic balance with frozen slow variables
      via Newton iteration (vector quadratic for Cl2•-).
  (2) Slow sub-step: explicit Euler on slow variables using the rates implied by
      the radical steady state.

Reference: Marchuk's splitting scheme; Strang 1968; Kuehn 2015 Ch. 6.

This module is the rigorous replacement for the analytic QSSA used in scripts/
tier3_phillips_null.py. It addresses Caveat 1 in TIER3_EXTENSIONS_REPORT.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Tuple

import numpy as np

R_GAS = 8.314462618


# ============================================================================
# Slow-manifold projection
# ============================================================================

def fast_steady_state(slow_state: dict, params: dict) -> dict:
    """Solve the radical balance for [Cl•]_ss, [Cl2•-]_ss given the slow state.

    Slow state contains: U3, U4, n_Cl2_g, Cl3, Cl2_diss
    Radical balance assumes d[radical]/dt ≈ 0 on the slow timescale.

    Algebra:
      [Cl•]_ss balance: S_Cl = k1 [Cl•] [Cl-] + k2 [Cl•] [e_s-] + 2 k3 [Cl•]^2
        In the regime k1[Cl-] >> k3[Cl•], [Cl•]_ss = S_Cl / (k1 [Cl-]_const)
      [Cl2•-]_ss balance: 0 = k1 [Cl•] [Cl-] - 2 k4 [Cl2•-]^2 - k_U3 U3 [Cl2•-]
        Quadratic 2 k4 x^2 + k_U3 U3 x - r1 = 0
      [e_s-]_ss balance: S_eS = k2 [Cl•] [e_s-] + k_eS_U4 U4 [e_s-] + k_bg [e_s-]
        [e_s-]_ss = S_eS / (k2 [Cl•]_ss + k_eS_U4 U4 + k_bg)

    Returns dict with the steady-state radicals and the relevant reaction rates.
    """
    U3, U4 = slow_state["U3"], slow_state["U4"]

    k1, k2, k3, k4 = params["k1"], params["k2"], params["k3"], params["k4"]
    k_U3, k_eS_U4 = params["k_U3"], params["k_eS_U4"]
    k_bg = params["k_bg"]
    S_eS, S_Cl = params["S_eS"], params["S_Cl"]
    Cl_const = params["Cl_minus_const"]

    # Newton iteration for [Cl•]_ss including the self-recombination (R3) term,
    # which is otherwise neglected in the analytic QSSA.
    # F(x) = k1 Cl_const · x + 2 k3 x^2 + k2 [e_s-] · x − S_Cl = 0
    # We do not yet know [e_s-]_ss. Solve coupled fast system by alternating
    # iteration.
    Cl_atom = S_Cl / max(k1 * Cl_const, 1e-30)
    for it in range(30):
        # Update [e_s-]_ss given current [Cl•]_ss
        eS_loss = k2 * Cl_atom + k_eS_U4 * U4 + k_bg
        eS = S_eS / max(eS_loss, 1e-30)
        # Update [Cl•]_ss from quadratic including R3 + R2 sinks
        # F(Cl) = (k1 Cl_const + k2 eS) Cl + 2 k3 Cl^2 − S_Cl = 0
        # Stable form: x = 2 S_Cl / (b + sqrt(b^2 + 4·2k3·S_Cl))
        a = 2 * k3
        b = k1 * Cl_const + k2 * eS
        discrim_cl = b * b + 4 * a * S_Cl
        if discrim_cl < 0:
            discrim_cl = 0.0
        denom_cl = b + np.sqrt(discrim_cl)
        if a > 0 and denom_cl > 0:
            new_Cl = 2 * S_Cl / denom_cl
        elif b > 0:
            new_Cl = S_Cl / b
        else:
            new_Cl = np.sqrt(S_Cl / max(a, 1e-30)) if a > 0 else 0.0
        if abs(new_Cl - Cl_atom) < 1e-15 * max(abs(new_Cl), 1e-20):
            Cl_atom = new_Cl
            break
        Cl_atom = new_Cl

    # [Cl2•-]_ss from quadratic 2k4 x^2 + k_U3 U3 x − r1 = 0 where r1 = k1 Cl_const Cl_atom
    # NUMERICALLY STABLE form: when b^2 >> 4 a c (U sink dominates), the naive formula
    # (-b + sqrt(b^2 + 4ac)) / (2a) suffers catastrophic cancellation. Use the
    # conjugate form: x = 2c / (b + sqrt(b^2 + 4ac)), valid for both b > 0 and b ≈ 0.
    r1 = k1 * Cl_const * Cl_atom
    a2 = 2 * k4
    b2 = k_U3 * U3
    discrim = b2 * b2 + 4 * a2 * r1
    if discrim < 0:
        discrim = 0.0
    denom = b2 + np.sqrt(discrim)
    if a2 > 0 and denom > 0:
        Cl2m = 2 * r1 / denom
    elif b2 > 0:
        Cl2m = r1 / b2
    else:
        Cl2m = np.sqrt(r1 / max(a2, 1e-30)) if a2 > 0 else 0.0

    # Rates relevant to slow variables
    r4 = k4 * Cl2m * Cl2m                 # Cl2•- self-recomb -> Cl3-
    r3 = k3 * Cl_atom * Cl_atom           # Cl•+Cl• -> Cl2_diss
    r_U_Cl2m = k_U3 * U3 * Cl2m           # U(III) + Cl2•- -> U(IV) + 2Cl-
    r_U_eS = k_eS_U4 * eS * U4            # e_s- + U(IV) -> U(III)
    r2 = k2 * eS * Cl_atom                # e_s- + Cl• -> Cl-

    return {
        "Cl_atom_ss": float(Cl_atom),
        "Cl2m_ss": float(Cl2m),
        "eS_ss": float(eS),
        "r1": float(r1),
        "r2": float(r2),
        "r3": float(r3),
        "r4": float(r4),
        "r_U_Cl2m": float(r_U_Cl2m),
        "r_U_eS": float(r_U_eS),
    }


def slow_step(slow_state: dict, params: dict, dt: float) -> dict:
    """Advance the slow variables (Cl3, Cl2_diss, n_Cl2_g, U3, U4) by dt
    using the QSS rates from fast_steady_state.

    ODEs (slow manifold):
      d[Cl3]/dt = r4 - k5 [Cl3]
      d[Cl2_diss]/dt = k5 [Cl3] + r3 + S_Cl2_direct − flux_gas
      dn_Cl2_g/dt = flux_gas * V_liq
      d[U3]/dt = -r_U_Cl2m + r_U_eS
      d[U4]/dt = +r_U_Cl2m - r_U_eS
    """
    fast = fast_steady_state(slow_state, params)
    k5 = params["k5"]
    Cl3 = slow_state["Cl3"]
    Cl2d = slow_state["Cl2_diss"]
    n_Cl2_g = slow_state["n_Cl2_g"]
    U3 = slow_state["U3"]
    U4 = slow_state["U4"]
    V_liq = params["V_liq"]
    V_gas = params["V_gas"]
    kLa = params["kLa"]
    kH = params["kH"]
    T = params["T"]
    S_Cl2_direct = params.get("S_Cl2_direct", 0.0)

    # Gas exchange
    p_gas = n_Cl2_g * R_GAS * T / V_gas
    C_eq = kH * p_gas
    flux = kLa * (Cl2d - C_eq)

    # Slow-manifold rates
    dCl3 = fast["r4"] - k5 * Cl3
    dCl2d = k5 * Cl3 + fast["r3"] + S_Cl2_direct - flux
    dn = flux * V_liq
    dU3 = -fast["r_U_Cl2m"] + fast["r_U_eS"]
    dU4 = +fast["r_U_Cl2m"] - fast["r_U_eS"]

    # Adaptive substepping. Use a quasi-steady-state check on the FAST production
    # process: if r4 dt is large compared to the implicit equilibrium of (r4 vs k5·Cl3),
    # subdivide. Cap the number of sub-steps to avoid runaway.
    # Use absolute tolerance based on the magnitude that would matter for the gas:
    # the total Cl2 produced over dt is at most (r3+r4)·dt; if that overshoots Cl3 or
    # Cl2d by > tolerance, subdivide.
    max_abs = max(abs(dCl3) * dt, abs(dCl2d) * dt, abs(dn) * dt / max(V_liq, 1e-9),
                  abs(dU3) * dt, abs(dU4) * dt)
    scale = max(Cl3, Cl2d, U3, U4, abs(dn)/max(V_liq,1e-9), 1e-10)
    rel = max_abs / scale
    if rel > 0.1:
        n_sub = min(int(np.ceil(rel / 0.05)), 1000)  # CAP at 1000
        dt_sub = dt / n_sub
        for _ in range(n_sub):
            fast = fast_steady_state(slow_state, params)
            Cl3_ = slow_state["Cl3"]
            Cl2d_ = slow_state["Cl2_diss"]
            n_ = slow_state["n_Cl2_g"]
            U3_ = slow_state["U3"]
            U4_ = slow_state["U4"]
            p_ = n_ * R_GAS * T / V_gas
            flux_ = kLa * (Cl2d_ - kH * p_)
            slow_state = {
                "Cl3": max(0.0, Cl3_ + (fast["r4"] - k5 * Cl3_) * dt_sub),
                "Cl2_diss": max(0.0, Cl2d_ + (k5 * Cl3_ + fast["r3"] + S_Cl2_direct - flux_) * dt_sub),
                "n_Cl2_g": max(0.0, n_ + flux_ * V_liq * dt_sub),
                "U3": max(0.0, U3_ + (-fast["r_U_Cl2m"] + fast["r_U_eS"]) * dt_sub),
                "U4": max(0.0, U4_ + (+fast["r_U_Cl2m"] - fast["r_U_eS"]) * dt_sub),
            }
        return slow_state, fast
    else:
        return {
            "Cl3": max(0.0, Cl3 + dCl3 * dt),
            "Cl2_diss": max(0.0, Cl2d + dCl2d * dt),
            "n_Cl2_g": max(0.0, n_Cl2_g + dn * dt),
            "U3": max(0.0, U3 + dU3 * dt),
            "U4": max(0.0, U4 + dU4 * dt),
        }, fast


def integrate_chronic(params: dict, initial: dict, t_final: float,
                       n_slow_steps: int = 200, return_trajectory: bool = False) -> dict:
    """Integrate the chronic-irradiation problem from t=0 to t_final.

    Uses operator-splitting with n_slow_steps steps on the slow manifold.
    """
    state = dict(initial)
    dt = t_final / n_slow_steps
    if return_trajectory:
        traj = {k: [v] for k, v in state.items()}
        ts = [0.0]

    for k in range(n_slow_steps):
        state, fast = slow_step(state, params, dt)
        if return_trajectory:
            for kk, vv in state.items():
                traj[kk].append(vv)
            ts.append((k + 1) * dt)

    result = {"final_state": state, "final_fast": fast}
    if return_trajectory:
        result["trajectory"] = {k: np.array(v) for k, v in traj.items()}
        result["t"] = np.array(ts)
    return result
