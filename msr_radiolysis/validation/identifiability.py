"""Tier 1 identifiability analysis: profile likelihood + cross-paper consistency check.

This module implements the Raue–Timmer (2009) profile-likelihood machinery for ODE-based
kinetic models, plus a B2BDC-style cross-paper consistency check. It runs BEFORE any
Bayesian inference to answer:

  (Q1) Practical identifiability: for each parameter θ_i, can the data alone determine
       a finite-width credible interval? (profile-likelihood threshold test)

  (Q2) Structural identifiability: is the parameter map θ ↦ observable injective?
       (profile-likelihood: flat profile signals non-identifiability)

  (Q3) Cross-paper consistency: is there a single Arrhenius (A, Ea) consistent with
       all reported per-T rate constants from different papers within their σ?

References:
- Raue, A.; Kreutz, C.; Maiwald, T. et al. "Structural and practical identifiability
  analysis of partially observed dynamical models by exploiting the profile likelihood."
  Bioinformatics 25(15), 1923-1929 (2009).
- Frenklach, M. et al. "Bound-to-Bound Data Collaboration" methodology.

Usage:
    from msr_radiolysis.validation.identifiability import (
        profile_likelihood, consistency_check, identifiability_report,
    )

    # See cli.py for end-to-end driver.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2, norm

R_GAS = 8.314462618  # J / (mol K)

# ---------------------------------------------------------------------------
# Profile-likelihood machinery
# ---------------------------------------------------------------------------

@dataclass
class ParameterSpec:
    """One free parameter for identifiability analysis.

    Attributes:
        name: human-readable identifier (e.g., 'log_A_eq5').
        initial: starting value for fits (units: native parameter units).
        lower: lower bound (None for unconstrained).
        upper: upper bound (None for unconstrained).
        transform: 'identity' or 'log' (log-space parameter).
    """
    name: str
    initial: float
    lower: Optional[float] = None
    upper: Optional[float] = None
    transform: str = "identity"

    def to_internal(self, x: float) -> float:
        return np.log(x) if self.transform == "log" else x

    def from_internal(self, z: float) -> float:
        return float(np.exp(z)) if self.transform == "log" else float(z)


@dataclass
class ProfileResult:
    """Profile-likelihood output for one parameter."""
    parameter: ParameterSpec
    grid_values: np.ndarray         # values of θ_i (native units) at which profile evaluated
    profile_nll: np.ndarray         # NLL at each grid value with other params re-optimized
    mle_value: float                # full MLE of this parameter
    mle_nll: float                  # NLL at full MLE
    lower_ci: Optional[float]       # 95% CI lower bound (None if hits grid edge)
    upper_ci: Optional[float]       # 95% CI upper bound (None if hits grid edge)
    flat_profile: bool              # True if max(profile) - min(profile) < threshold
    note: str = ""


def profile_likelihood(
    nll: Callable[[np.ndarray], float],
    parameters: Sequence[ParameterSpec],
    mle: np.ndarray,
    target_idx: int,
    grid_span: float = 3.0,
    n_grid: int = 21,
    confidence_level: float = 0.95,
) -> ProfileResult:
    """Compute the profile NLL for parameter ``target_idx``.

    The profile is defined as
        PL(θ_i = v) = min over θ_{-i} of  NLL(θ_i = v, θ_{-i})
    We sweep v over a grid and re-optimize the remaining parameters at each grid point.

    The 1D profile threshold for a (1-α) confidence interval is
        Δ = 0.5 * χ²_1(1-α)   (with the factor of 1/2 since NLL = -log L)
    e.g., Δ ≈ 1.92 for α = 0.05.

    Args:
        nll: function taking a parameter vector (internal-space) and returning -log L.
        parameters: spec for each parameter.
        mle: full MLE (internal-space) — the starting point and reference NLL value.
        target_idx: index of the parameter to profile.
        grid_span: half-width of the grid in internal-space units around MLE.
        n_grid: number of grid points.
        confidence_level: confidence level for the CI threshold.

    Returns:
        ProfileResult with the swept profile, MLE, and 95% CIs.
    """
    spec = parameters[target_idx]
    threshold = 0.5 * chi2.ppf(confidence_level, df=1)

    grid_internal = np.linspace(mle[target_idx] - grid_span, mle[target_idx] + grid_span, n_grid)

    profile_nll = np.full(n_grid, np.nan)
    for i, target_val in enumerate(grid_internal):
        free_idx = [j for j in range(len(parameters)) if j != target_idx]

        def neg_logL_other(theta_other):
            theta = mle.copy()
            theta[target_idx] = target_val
            theta[free_idx] = theta_other
            return nll(theta)

        x0_other = mle[free_idx].copy()
        bounds_other = []
        for j in free_idx:
            p = parameters[j]
            lo = p.to_internal(p.lower) if (p.lower is not None) else None
            hi = p.to_internal(p.upper) if (p.upper is not None) else None
            bounds_other.append((lo, hi))

        try:
            res = minimize(neg_logL_other, x0_other, method="L-BFGS-B", bounds=bounds_other,
                           options={"maxiter": 200, "ftol": 1e-10})
            profile_nll[i] = res.fun if res.success else np.nan
        except Exception as e:
            profile_nll[i] = np.nan

    grid_native = np.array([spec.from_internal(z) for z in grid_internal])
    mle_native = spec.from_internal(mle[target_idx])

    # Find CIs by interpolating threshold crossings
    delta = profile_nll - nll(mle)
    lower_ci, upper_ci = None, None
    valid = ~np.isnan(delta)
    if valid.any():
        below = delta < threshold
        if below.any() and (~below).any():
            # find leftmost and rightmost crossings
            idx_lt = np.where(below)[0]
            i_min = idx_lt.min()
            i_max = idx_lt.max()
            if i_min > 0:
                lower_ci = float(np.interp(threshold, [delta[i_min-1], delta[i_min]], [grid_native[i_min-1], grid_native[i_min]]))
            if i_max < n_grid - 1:
                upper_ci = float(np.interp(threshold, [delta[i_max+1], delta[i_max]], [grid_native[i_max+1], grid_native[i_max]]))

    flat = (np.nanmax(delta) - np.nanmin(delta)) < 0.5 * threshold

    return ProfileResult(
        parameter=spec,
        grid_values=grid_native,
        profile_nll=profile_nll,
        mle_value=mle_native,
        mle_nll=nll(mle),
        lower_ci=lower_ci,
        upper_ci=upper_ci,
        flat_profile=bool(flat),
        note=("flat profile -> structurally non-identifiable" if flat else
              ("no CI bound -> practically non-identifiable on one side" if (lower_ci is None or upper_ci is None) else "")),
    )


# ---------------------------------------------------------------------------
# Cross-paper consistency check (B2BDC-style)
# ---------------------------------------------------------------------------

@dataclass
class RateObservation:
    """One rate-constant observation from one paper at one T."""
    paper: str
    reaction_id: str
    T_K: float
    log_k: float          # log of measured rate constant
    sigma_log_k: float    # 1σ uncertainty in log space
    salt: str             # salt system identifier


@dataclass
class ConsistencyResult:
    """Output of the cross-paper Arrhenius consistency check."""
    reaction_id: str
    n_observations: int
    n_papers: int
    feasible: bool                  # True if a single (A, Ea) is consistent with all obs
    A_best: float                   # best-fit A
    Ea_best: float                  # best-fit Ea (J/mol)
    chi2_min: float                 # minimum χ²
    chi2_critical: float            # χ²_{n - 2, 0.95}
    p_value: float                  # p-value of the χ² test
    residuals: Dict[str, float]     # per-observation standardized residuals
    note: str = ""


def consistency_check(observations: Sequence[RateObservation],
                      reaction_id: str,
                      alpha: float = 0.05) -> ConsistencyResult:
    """Check whether all observations are consistent with a single Arrhenius (A, Ea).

    Fits log k_i = log A − Ea / (R T_i) by weighted least squares; tests the χ² goodness
    of fit at level α. The test is the standard F-statistic for nested linear models,
    appropriate when measurement σ is well-calibrated.

    Args:
        observations: list of (paper, T, log_k, sigma_log_k) tuples for ONE reaction.
        reaction_id: identifier for the reaction (for output).
        alpha: significance level (default 0.05).

    Returns:
        ConsistencyResult with the fit and feasibility assessment.
    """
    obs = list(observations)
    n = len(obs)
    if n < 2:
        raise ValueError(f"Need at least 2 observations for consistency check; got {n}")

    # Design matrix for log k_i = log A − Ea/(R T_i)
    T = np.array([o.T_K for o in obs])
    y = np.array([o.log_k for o in obs])
    sigma = np.array([o.sigma_log_k for o in obs])
    X = np.column_stack([np.ones(n), -1.0 / (R_GAS * T)])  # columns: [log A coef, Ea coef]

    # Weighted least squares
    W = np.diag(1.0 / sigma**2)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    try:
        beta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return ConsistencyResult(
            reaction_id=reaction_id, n_observations=n, n_papers=len(set(o.paper for o in obs)),
            feasible=False, A_best=np.nan, Ea_best=np.nan, chi2_min=np.inf,
            chi2_critical=np.nan, p_value=0.0, residuals={},
            note="design matrix singular - one-T observations cannot identify (A, Ea)",
        )

    log_A_hat, Ea_hat = beta
    y_pred = X @ beta
    residuals = (y - y_pred) / sigma
    chi2_min = float(np.sum(residuals**2))
    dof = n - 2
    chi2_critical = float(chi2.ppf(1 - alpha, df=dof)) if dof > 0 else np.nan
    p_value = float(1.0 - chi2.cdf(chi2_min, df=dof)) if dof > 0 else np.nan
    feasible = (chi2_min <= chi2_critical) if dof > 0 else True

    return ConsistencyResult(
        reaction_id=reaction_id,
        n_observations=n,
        n_papers=len(set(o.paper for o in obs)),
        feasible=feasible,
        A_best=float(np.exp(log_A_hat)),
        Ea_best=float(Ea_hat),
        chi2_min=chi2_min,
        chi2_critical=chi2_critical,
        p_value=p_value,
        residuals={f"{o.paper}_T{int(o.T_K-273.15)}C": float(r) for o, r in zip(obs, residuals)},
        note=(f"PASS: {n} obs from {len(set(o.paper for o in obs))} papers consistent" if feasible
              else f"FAIL: χ²={chi2_min:.2f} > χ²_crit={chi2_critical:.2f} at α={alpha}; check residuals for outliers"),
    )


# ---------------------------------------------------------------------------
# Driver: pseudo-1st-order rate as identifiability target
# ---------------------------------------------------------------------------

def pseudo_first_order_nll(theta: np.ndarray,
                            data_tuples: Sequence[Tuple[float, float, float, float]]) -> float:
    """Negative log-likelihood under the pseudo-1st-order kinetic model.

    Model: k_obs = A · exp(-Ea/(R T)) · [M]  (units: 1/s)
    Observations: (T, [M], k_obs_measured, sigma_log_k)

    θ = [log A, Ea].

    Returns sum of squared standardized residuals in log space (= -2 log L for Gaussian
    log-noise model, up to constants).
    """
    log_A, Ea = theta
    total = 0.0
    for T, Mconc, k_obs, sigma in data_tuples:
        log_k_pred = log_A - Ea / (R_GAS * T) + np.log(Mconc)
        residual = (np.log(k_obs) - log_k_pred) / sigma
        total += 0.5 * residual**2
    return total


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def identifiability_report(profile_results: Sequence[ProfileResult],
                            consistency_results: Sequence[ConsistencyResult],
                            output_path: Path) -> None:
    """Write a markdown identifiability report."""
    lines = [
        "# Tier 1 identifiability + consistency report",
        "",
        "Auto-generated by `msr_radiolysis.validation.identifiability`.",
        "",
        "## Profile-likelihood results",
        "",
        "Threshold for 95% CI: Δ_NLL ≈ 1.92 (χ²_1 = 3.84).",
        "",
        "| Parameter | MLE | 95% CI | Identifiable? | Note |",
        "|---|---|---|---|---|",
    ]
    for pr in profile_results:
        ci_str = (
            f"[{pr.lower_ci:.3g}, {pr.upper_ci:.3g}]" if (pr.lower_ci is not None and pr.upper_ci is not None)
            else f"[{pr.lower_ci or '−∞'}, {pr.upper_ci or '+∞'}]"
        )
        status = "✓" if (pr.lower_ci is not None and pr.upper_ci is not None and not pr.flat_profile) else "✗"
        lines.append(f"| {pr.parameter.name} | {pr.mle_value:.4g} | {ci_str} | {status} | {pr.note} |")

    lines += [
        "",
        "## Cross-paper Arrhenius consistency",
        "",
        "Reactions with multi-paper rate data tested against a single (A, Ea). Significance α = 0.05.",
        "",
        "| Reaction | n_obs | n_papers | Best A | Best Ea (kJ/mol) | χ² / χ²_crit | p-value | Feasible? | Note |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for cr in consistency_results:
        chi_str = f"{cr.chi2_min:.2f} / {cr.chi2_critical:.2f}" if np.isfinite(cr.chi2_critical) else f"{cr.chi2_min:.2f} / —"
        status = "✓" if cr.feasible else "✗"
        Ea_kJ = cr.Ea_best / 1e3 if np.isfinite(cr.Ea_best) else np.nan
        A_str = f"{cr.A_best:.2e}" if np.isfinite(cr.A_best) else "—"
        lines.append(f"| {cr.reaction_id} | {cr.n_observations} | {cr.n_papers} | {A_str} | {Ea_kJ:.2f} | {chi_str} | {cr.p_value:.4f} | {status} | {cr.note} |")
        if not cr.feasible:
            lines.append("")
            lines.append(f"  Standardized residuals for {cr.reaction_id}:")
            for k, v in cr.residuals.items():
                tag = " ←  OUTLIER" if abs(v) > 2.0 else ""
                lines.append(f"  - {k}: {v:+.2f} σ{tag}")
            lines.append("")

    output_path.write_text("\n".join(lines))
