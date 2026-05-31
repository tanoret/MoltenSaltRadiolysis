#!/usr/bin/env python3
"""Tier 2 Bayesian calibration: minimal HBMAE implementation on Iwamatsu 2026 Cr transients.

This is the first full Bayesian calibration: 9 transient absorbance traces (Iwamatsu 2026
Figs. 2A and 3A digitized as CSVs) are fit through the chloride+Cr kernel via stiff
ODE forward solves. The posterior is sampled with the affine-invariant ensemble MCMC of
Goodman & Weare (2010) as implemented in emcee \citep{ForemanMackey2013}.

We adopt the HBMAE prior structure with one simplification: only a single salt
(LiCl-KCl) is used, so the hierarchical layer degenerates to the within-salt prior.
The Pikaev/Iwamatsu cross-salt analysis would require Tier 2 on at least two salts
simultaneously (the next step).

What is calibrated:
  - Arrhenius (A, Ea) for e_s- + Cr2+ -> Cr+        (Iwamatsu Eq. 5)
  - Arrhenius (A, Ea) for e_s- + Cr3+ -> Cr2+        (Iwamatsu Eq. 6)
  - Pulse-dose nuisance: initial [e_s-]_0 per trace  (informative log-normal prior)
  - Background e_s- impurity decay rate k_bg          (weakly informative)

Priors:
  - Arrhenius parameters: literature-informed Gaussian from arrhenius_parameters.csv
  - [e_s-]_0: LogNormal(μ=log(1.7e-2 mol/m³), σ=0.5) — covers 5e-3 to 5e-2 mol/m³ at 1σ
  - k_bg: LogNormal(μ=log(1e7 s^-1), σ=0.5)

Likelihood:
  - Absorbance traces compared to model e_s- prediction in scale-free observable mode
    (kind = "absorbance_scale_free" in HBMAE terms): the model trace is rescaled
    so its max matches the experimental amplitude, and the residuals are computed
    on the rescaled curve. This avoids needing the unknown ε_λ at 700 nm.
  - Per-point Gaussian noise σ_obs = 0.005 (calibrated to experimental absorbance scale).
"""

from __future__ import annotations

import csv
import os
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# Silence harmless warnings from stiff-ODE BDF on edge cases
warnings.filterwarnings("ignore", category=UserWarning)
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import emcee

R_GAS = 8.314462618

# ============================================================================
# Reduced ODE model: e_s- + (Cr2+, Cr3+) kinetics
# ============================================================================
# State vector: y = [e_s-, Cr2+, Cr+, Cr3+]   (units: mol/m^3)
# Other species (Cl-, Cl•, Cl2•-) are pseudo-stationary in the timescale studied (≤ 25 ns).
#
# Reactions:
#   r5: e_s- + Cr2+ -> Cr+       (Eq. 5, k5 = A5 exp(-Ea5/RT))
#   r6: e_s- + Cr3+ -> Cr2+      (Eq. 6, k6 = A6 exp(-Ea6/RT))
#   r_bg: e_s- -> (impurity)     (k_bg = first-order background decay, 1e7 s^-1)
#
# ODEs (per 1 m^3 solution, concentrations in mol/m^3, rates k in m^3/(mol s)):
#   d[e_s-]/dt = -(k5 [Cr2+] + k6 [Cr3+] + k_bg) [e_s-]
#   d[Cr2+]/dt = -k5 [e_s-][Cr2+] + k6 [e_s-][Cr3+]
#   d[Cr+]/dt  = +k5 [e_s-][Cr2+]
#   d[Cr3+]/dt = -k6 [e_s-][Cr3+]

def rhs(t, y, k5_SI, k6_SI, k_bg):
    eS, Cr2, Crp, Cr3 = y
    rate5 = k5_SI * eS * Cr2
    rate6 = k6_SI * eS * Cr3
    deS_dt = -(rate5 + rate6 + k_bg * eS)
    dCr2_dt = -rate5 + rate6
    dCrp_dt = +rate5
    dCr3_dt = -rate6
    return [deS_dt, dCr2_dt, dCrp_dt, dCr3_dt]


# ============================================================================
# Forward solve for one trace
# ============================================================================
def solve_trace(eS0, Cr2_0, Cr3_0, A5, Ea5, A6, Ea6, k_bg, t_eval, T_K=673.15):
    """Forward ODE solve. Returns e_s-(t_eval) array. A in M^-1 s^-1, Ea in J/mol."""
    k5_M = A5 * np.exp(-Ea5 / (R_GAS * T_K))  # M^-1 s^-1
    k6_M = A6 * np.exp(-Ea6 / (R_GAS * T_K))
    k5_SI = k5_M * 1e-3   # m^3 / (mol s)
    k6_SI = k6_M * 1e-3
    y0 = [eS0, Cr2_0, 0.0, Cr3_0]
    try:
        sol = solve_ivp(rhs, (0, t_eval[-1] * 1.01), y0, t_eval=t_eval,
                        args=(k5_SI, k6_SI, k_bg),
                        method="BDF", atol=1e-15, rtol=1e-9, max_step=1e-9)
        if not sol.success:
            return None
        return sol.y[0]   # e_s-(t)
    except Exception:
        return None


# ============================================================================
# Load digitized transients
# ============================================================================
def load_traces():
    """Load the 9 Iwamatsu 2026 absorbance traces.

    Returns list of dicts with keys: label, t_s, abs, Cr2_init, Cr3_init.
    """
    traces = []
    base = REPO / "validation/cr_licl_kcl/iwamatsu_2026_pccp/data"
    metadata = [
        ("absorbance1mMCr2.csv", "1 mM Cr(II)", 0.99, 0.0),
        ("absorbance2mMCr2.csv", "2 mM Cr(II)", 1.98, 0.0),
        ("absorbance3mMCr2.csv", "3 mM Cr(II)", 2.97, 0.0),
        ("absorbance4mMCr2.csv", "4 mM Cr(II)", 3.96, 0.0),
        ("absorbance1mMCr3.csv", "1 mM Cr(III)", 0.0, 1.05),
        ("absorbance2mMCr3.csv", "2 mM Cr(III)", 0.0, 2.10),
        ("absorbance3mMCr3.csv", "3 mM Cr(III)", 0.0, 3.15),
        ("absorbance4mMCr3.csv", "4 mM Cr(III)", 0.0, 4.20),
        ("absorbance5mMCr3.csv", "5 mM Cr(III)", 0.0, 5.25),
    ]
    for fname, label, Cr2, Cr3 in metadata:
        path = base / fname
        if not path.exists():
            continue
        ts, As = [], []
        with path.open() as f:
            r = csv.DictReader(f)
            for row in r:
                t_ns = float(row["time"])
                a = float(row["absorbance"])
                ts.append(t_ns * 1e-9)
                As.append(a)
        ts_arr = np.array(ts)
        As_arr = np.array(As)
        # Some digitized CSVs (e.g. absorbance4mMCr3.csv) have unsorted timestamps from
        # the WebPlotDigitizer extraction. solve_ivp requires monotonic t_eval; sort
        # and de-duplicate before storing. (Bug discovered in Tier 4; documented.)
        order = np.argsort(ts_arr)
        ts_arr, As_arr = ts_arr[order], As_arr[order]
        mask = np.concatenate([[True], np.diff(ts_arr) > 0])
        traces.append({
            "label": label,
            "t_s": ts_arr[mask],
            "abs": As_arr[mask],
            "Cr2_init": Cr2,   # mol/m³
            "Cr3_init": Cr3,
        })
    return traces


# ============================================================================
# Priors
# ============================================================================
# Prior layout (free parameters, in this order):
#   theta[0] = log A5
#   theta[1] = Ea5  (J/mol)
#   theta[2] = log A6
#   theta[3] = Ea6  (J/mol)
#   theta[4..4+N-1] = log [e_s-]_0 for each of N traces (per-trace pulse-dose nuisance)
#   theta[-1] = log k_bg (s^-1)

# Literature priors (Iwamatsu 2026, density-corrected)
PRIOR_LOG_A5_MU, PRIOR_LOG_A5_SIGMA = np.log(1.7e13), 0.2 / 1.7   # log(A ± 0.2e13/A)
PRIOR_EA5_MU,    PRIOR_EA5_SIGMA = 33.5e3, 0.6e3
PRIOR_LOG_A6_MU, PRIOR_LOG_A6_SIGMA = np.log(2.0e13), 0.5 / 2.0
PRIOR_EA6_MU,    PRIOR_EA6_SIGMA = 31.8e3, 0.5e3

PRIOR_LOG_ES0_MU, PRIOR_LOG_ES0_SIGMA = np.log(1.7e-2), 0.5  # mol/m^3, covers 5e-3 to 5e-2
PRIOR_LOG_KBG_MU, PRIOR_LOG_KBG_SIGMA = np.log(1e7), 0.5


def log_prior(theta, N_traces):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    log_eS0s = theta[4:4 + N_traces]
    log_kbg = theta[-1]

    if Ea5 < 0 or Ea5 > 1e5: return -np.inf
    if Ea6 < 0 or Ea6 > 1e5: return -np.inf

    lp = 0.0
    lp += -0.5 * ((log_A5 - PRIOR_LOG_A5_MU) / PRIOR_LOG_A5_SIGMA) ** 2
    lp += -0.5 * ((Ea5 - PRIOR_EA5_MU) / PRIOR_EA5_SIGMA) ** 2
    lp += -0.5 * ((log_A6 - PRIOR_LOG_A6_MU) / PRIOR_LOG_A6_SIGMA) ** 2
    lp += -0.5 * ((Ea6 - PRIOR_EA6_MU) / PRIOR_EA6_SIGMA) ** 2
    for log_eS0 in log_eS0s:
        lp += -0.5 * ((log_eS0 - PRIOR_LOG_ES0_MU) / PRIOR_LOG_ES0_SIGMA) ** 2
    lp += -0.5 * ((log_kbg - PRIOR_LOG_KBG_MU) / PRIOR_LOG_KBG_SIGMA) ** 2
    return lp


def log_likelihood(theta, traces, sigma_obs=5e-3):
    log_A5, Ea5, log_A6, Ea6 = theta[0:4]
    A5 = np.exp(log_A5)
    A6 = np.exp(log_A6)
    log_eS0s = theta[4:4 + len(traces)]
    log_kbg = theta[-1]
    k_bg = np.exp(log_kbg)

    ll = 0.0
    for i, tr in enumerate(traces):
        eS0 = np.exp(log_eS0s[i])
        # Forward solve at the trace's experimental time grid
        eS_model = solve_trace(eS0, tr["Cr2_init"], tr["Cr3_init"], A5, Ea5, A6, Ea6, k_bg, tr["t_s"])
        if eS_model is None or not np.all(np.isfinite(eS_model)):
            return -np.inf
        # Scale-free observable: rescale model to match max(abs)
        m_max = np.max(eS_model)
        if m_max <= 0:
            return -np.inf
        scale = np.nanmax(tr["abs"]) / m_max
        pred = eS_model * scale
        # Per-point Gaussian residual
        res = (pred - tr["abs"]) / sigma_obs
        ll += -0.5 * np.sum(res * res)
    return ll


def log_posterior(theta, traces):
    lp = log_prior(theta, len(traces))
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, traces)


# ============================================================================
# Run MCMC
# ============================================================================
def main():
    print("=" * 76)
    print("Tier 2: Bayesian calibration of Cr+e_s- kinetics in LiCl-KCl")
    print("Forward model: stiff ODE; sampler: emcee affine-invariant ensemble")
    print("=" * 76)

    traces = load_traces()
    N = len(traces)
    print(f"Loaded {N} traces from Iwamatsu 2026 Figs. 2A/3A")
    for tr in traces:
        print(f"  {tr['label']:>18s} : n_t={len(tr['t_s'])}, t∈[{tr['t_s'][0]*1e9:.2f}, {tr['t_s'][-1]*1e9:.2f}] ns, "
              f"Cr2+={tr['Cr2_init']:.2f}, Cr3+={tr['Cr3_init']:.2f}")

    # Parameter dim: 4 Arrhenius + N pulse-dose + 1 k_bg
    ndim = 4 + N + 1
    nwalkers = 2 * ndim + 8

    # Initial positions: small perturbations of prior means
    rng = np.random.default_rng(0)
    initial = np.array([
        PRIOR_LOG_A5_MU, PRIOR_EA5_MU,
        PRIOR_LOG_A6_MU, PRIOR_EA6_MU,
        *([PRIOR_LOG_ES0_MU] * N),
        PRIOR_LOG_KBG_MU,
    ])
    scales = np.array([
        0.05, 0.5e3, 0.05, 0.5e3,
        *([0.2] * N),
        0.2,
    ])
    p0 = initial + scales * rng.normal(size=(nwalkers, ndim))

    print(f"\nMCMC config:")
    print(f"  ndim       = {ndim}  ({4} Arrhenius + {N} pulse-dose + 1 k_bg)")
    print(f"  nwalkers   = {nwalkers}")
    print(f"  steps      = 400 (burn-in 100 + production 300)")

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(traces,))

    print("\nBurn-in (100 steps)...")
    state = sampler.run_mcmc(p0, 100, progress=False)
    sampler.reset()

    print("Production (300 steps)...")
    sampler.run_mcmc(state, 300, progress=False)

    # Diagnostics
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        eff_n = sampler.iteration * nwalkers / np.mean(tau)
    except Exception:
        tau = None
        eff_n = sampler.iteration * nwalkers / 50  # rough

    print()
    print("=" * 76)
    print("Posterior summary")
    print("=" * 76)
    chain = sampler.get_chain(flat=True)
    log_prob = sampler.get_log_prob(flat=True)
    print(f"  Total samples: {chain.shape[0]}; effective: ~{eff_n:.0f}")
    print(f"  log_prob max:  {np.max(log_prob):.2f}")

    labels = ["log A5", "Ea5", "log A6", "Ea6"] + [f"log[e_s]0_t{i+1}" for i in range(N)] + ["log k_bg"]
    print()
    print(f"{'Parameter':<18s} {'Mean':>12s} {'σ':>10s} {'2.5%':>12s} {'97.5%':>12s}")
    print("-" * 70)
    summary = []
    for i, lab in enumerate(labels):
        m = np.mean(chain[:, i])
        s = np.std(chain[:, i])
        lo, hi = np.percentile(chain[:, i], [2.5, 97.5])
        print(f"{lab:<18s} {m:>12.4g} {s:>10.3g} {lo:>12.4g} {hi:>12.4g}")
        summary.append((lab, m, s, lo, hi))

    # Convert log_A and Ea to A and check against literature
    print()
    print("Posterior on physical parameters:")
    A5_samples = np.exp(chain[:, 0])
    Ea5_samples = chain[:, 1]
    A6_samples = np.exp(chain[:, 2])
    Ea6_samples = chain[:, 3]
    print(f"  A5 (e_s- + Cr2+):   posterior median = {np.median(A5_samples):.3e}  (literature 1.7e13 ± 0.2e13)")
    print(f"  Ea5:                posterior median = {np.median(Ea5_samples)/1e3:.2f} kJ/mol  (literature 33.5 ± 0.6)")
    print(f"  A6 (e_s- + Cr3+):   posterior median = {np.median(A6_samples):.3e}  (literature 2.0e13 ± 0.5e13)")
    print(f"  Ea6:                posterior median = {np.median(Ea6_samples)/1e3:.2f} kJ/mol  (literature 31.8 ± 0.5)")

    # Save chain and summary
    out_dir = REPO / "validation"
    np.save(out_dir / "tier2_chain.npy", sampler.get_chain())
    np.save(out_dir / "tier2_log_prob.npy", sampler.get_log_prob())
    with (out_dir / "TIER2_POSTERIOR_SUMMARY.csv").open("w") as f:
        f.write("parameter,mean,sigma,ci_2.5,ci_97.5\n")
        for lab, m, s, lo, hi in summary:
            f.write(f"{lab},{m},{s},{lo},{hi}\n")
    print(f"\nSaved chain to {out_dir.relative_to(REPO)}/tier2_chain.npy and posterior summary CSV.")
    return sampler, summary


if __name__ == "__main__":
    main()
