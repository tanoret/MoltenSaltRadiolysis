"""
Quantitative comparison between experimental (paper) and simulated (CRN)
absorbance vs. time data for the COMPOSITE Cl2•- + eS- signal at 400 nm
or 700 nm in molten LiCl-KCl eutectic containing CrCl2, as studied by
Iwamatsu et al. 2026.

FIT_MODE controls how ε values are fitted cl2 ε are exact:

  "cl2_only"   → Hold ε(eS-,400nm) fixed at its spectral estimate and fit
                  only ε(Cl2•-,400nm). Recommended: fewest assumptions about
                  the unknown ε, and eS- estimate is better constrained than
                  Cl2•- at this wavelength since we know the band is broad.

  "both"       → Fit both ε(Cl2•-) and ε(eS-) simultaneously using bounded
                  2D minimisation. More flexible but requires caution: the two
                  parameters are partially correlated at early times.

  "fixed"      → Use literature/spectral estimates for both; no fitting.
                  Useful for forward prediction or sensitivity checks.


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import pearsonr

# ─── CONFIG ──────────────────────────────────────────────────────────────────
EXP_CSV = "absorbance1mMCr3.csv"
SIM_CSV = "species_concentration.csv"   # two-species CSV from integration script

# Column indices in CSVs
EXP_TIME_COL = 0
EXP_ABS_COL  = 1
# SIM columns are read by header name (see load_sim_csv)

# ── Measurement conditions ────────────────────────────────────────────────────
PATH_LENGTH   = 0.5    # cm
WAVELENGTH_NM = 700    # nm

# ── Extinction coefficient priors and bounds ──────────────────────────────────
#
# Cl2•- at 700 nm
EPS_CL2_PRIOR = 0          # M⁻¹cm⁻¹ (For 700 nm, set to zero)
EPS_CL2_BOUNDS = (0, 0)  # M⁻¹cm⁻¹  fit bounds (exact)

# eS- at 700 nm
# Estimated from Fig. 1(A) Iwamatsu 2022:
#   Gε(400nm) ≈ 28×10³  (reading broad eS- tail at 5ns)
#   G(eS-)    ≈ 2.8     → ε ≈ 10,000 M⁻¹cm⁻¹
# Range: 8,900–12,500 M⁻¹cm⁻¹ (G = 2.5–3.1)
EPS_ES_PRIOR  = 15957         # M⁻¹cm⁻¹  central spectral estimate
EPS_ES_BOUNDS = (13062, 18852)  # M⁻¹cm⁻¹  fit bounds

# ── Fitting mode ──────────────────────────────────────────────────────────────
# "cl2_only" : fit ε(Cl2•-), hold ε(eS-) at EPS_ES_PRIOR  [recommended]
# "both"     : fit both ε(Cl2•-) and ε(eS-)
# "fixed"    : use both priors, no fitting
FIT_MODE = "both"

# ── Time axis ─────────────────────────────────────────────────────────────────
TIME_SCALE_FACTOR = 1e9    # multiply sim time [s] → ns
Sim_Time_Offset = 7
HAS_EXP_HEADER    = True
HAS_SIM_HEADER    = True
SPARSE_GAP_FRACTION = 0.05
TIME_LABEL = "Time (ns)"
# ─────────────────────────────────────────────────────────────────────────────


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_exp_csv(path: str):
    """Load experimental absorbance CSV. Returns (time_ns, absorbance)."""
    header = 0 if HAS_EXP_HEADER else None
    df = pd.read_csv(path, header=header)
    x = df.iloc[:, EXP_TIME_COL].to_numpy(dtype=float)
    y = df.iloc[:, EXP_ABS_COL].to_numpy(dtype=float)
    order = np.argsort(x)
    return x[order], y[order]


def load_sim_csv(path: str):
    """
    Load two-species simulation CSV.
    Expected columns: t_s, Cl2rad- (or Cl2•-), e_s-
    Returns (time_s, conc_Cl2, conc_eS) as numpy arrays.
    Raises KeyError with helpful message if columns are missing.
    """
    header = 0 if HAS_SIM_HEADER else None
    df = pd.read_csv(path, header=header)
    df.columns = [c.strip() for c in df.columns]

    # Find time column (first column or named t_s)
    t_col = df.columns[0]
    t = df[t_col].to_numpy(dtype=float)

    # Find Cl2•- column — try several possible sanitised names
    cl2_candidates = ["Cl2rad-", "Cl2•-", "Cl2_rad-", "Cl2.-"]
    cl2_col = None
    for cand in cl2_candidates:
        if cand in df.columns:
            cl2_col = cand
            break
    if cl2_col is None:
        raise KeyError(
            f"Cannot find Cl2•- column in '{path}'.\n"
            f"Available columns: {list(df.columns)}\n"
            f"Tried: {cl2_candidates}\n"
            f"Check that export_species_csv was called with "
            f"species_names=['Cl2•-', 'e_s-']."
        )

    # Find eS- column
    es_candidates = ["e_s-", "e_srad-", "eS-", "e_s"]
    es_col = None
    for cand in es_candidates:
        if cand in df.columns:
            es_col = cand
            break
    if es_col is None:
        raise KeyError(
            f"Cannot find eS- column in '{path}'.\n"
            f"Available columns: {list(df.columns)}\n"
            f"Tried: {es_candidates}\n"
            f"Check that export_species_csv was called with "
            f"species_names=['Cl2•-', 'e_s-']."
        )

    order = np.argsort(t)
    return (t[order],
            df[cl2_col].to_numpy(dtype=float)[order],
            df[es_col].to_numpy(dtype=float)[order])


# ── Spline + coverage  ─────────────────────────────────────────────────

def make_splines(sim_t, sim_cl2, sim_es):
    """Fit cubic splines to both species. Returns (cs_cl2, cs_es)."""
    cs_cl2 = CubicSpline(sim_t, sim_cl2)
    cs_es  = CubicSpline(sim_t, sim_es)
    for name, arr, cs in [("Cl2•-", sim_cl2, cs_cl2), ("eS-", sim_es, cs_es)]:
        resid = np.abs(cs(sim_t) - arr)
        print(f"  Spline [{name}]: max |residual| = {resid.max():.3e}, "
              f"mean = {resid.mean():.3e}")
        if resid.max() > 1e-9:
            print(f"  ⚠  Non-negligible spline residuals for {name}.")
    return cs_cl2, cs_es


def check_coverage(exp_t, sim_t):
    """Return boolean mask of exp points inside sim time domain."""
    in_domain = (exp_t >= sim_t.min()) & (exp_t <= sim_t.max())
    n_out = (~in_domain).sum()
    if n_out:
        print(f"  ⚠  {n_out} experimental point(s) outside simulation range — excluded.")
        print(f"     Exp: [{exp_t.min():.4g}, {exp_t.max():.4g}] ns  "
              f"Sim: [{sim_t.min():.4g}, {sim_t.max():.4g}] ns")
        print(f"     → Increase t_final in the simulation (≥ 200 ns recommended).")
    gap_thresh = SPARSE_GAP_FRACTION * (sim_t.max() - sim_t.min())
    for i, et in enumerate(exp_t):
        if not in_domain[i]:
            continue
        idx = np.clip(np.searchsorted(sim_t, et), 1, len(sim_t) - 1)
        if sim_t[idx] - sim_t[idx - 1] > gap_thresh:
            print(f"  ⚠  Exp point t={et:.4g} ns in sparse simulation gap.")
    return in_domain


# ── Composite absorbance + fitting ────────────────────────────────────────────

def composite_absorbance(eps_cl2: float, eps_es: float,
                          c_cl2: np.ndarray, c_es: np.ndarray) -> np.ndarray:
    """
    A_pred = ε(Cl2•-) * l * C(Cl2•-) + ε(eS-) * l * C(eS-)
    Units: M⁻¹cm⁻¹ × cm × mol/L = dimensionless (AU)
    """
    return PATH_LENGTH * (eps_cl2 * c_cl2 + eps_es * c_es)


def fit_eps_cl2_only(exp_abs, c_cl2_interp, c_es_interp, eps_es_fixed):
    """
    Hold ε(eS-) fixed; fit ε(Cl2•-) by minimising RMSE.
    Analytical solution exists:
        d/dε Σ[ε*l*C_Cl2 + l*ε_es*C_eS - A_exp]² = 0
        → ε* = Σ[C_Cl2*(A_exp - l*ε_es*C_eS)] / (l * Σ[C_Cl2²])
    """
    l = PATH_LENGTH
    residual_from_es = exp_abs - l * eps_es_fixed * c_es_interp
    eps_cl2_opt = (np.sum(c_cl2_interp * residual_from_es)
                   / (l * np.sum(c_cl2_interp ** 2)))
    eps_cl2_opt = float(np.clip(eps_cl2_opt, *EPS_CL2_BOUNDS))
    return eps_cl2_opt


def fit_both_eps(exp_abs, c_cl2_interp, c_es_interp):
    """
    Fit both ε(Cl2•-) and ε(eS-) simultaneously by minimising RMSE.
    Uses scipy.optimize.minimize with L-BFGS-B and bounds.
    Returns (eps_cl2_opt, eps_es_opt).
    """
    def rmse(params):
        e_cl2, e_es = params
        pred = composite_absorbance(e_cl2, e_es, c_cl2_interp, c_es_interp)
        return np.sqrt(np.mean((pred - exp_abs) ** 2))

    result = minimize(
        rmse,
        x0=[EPS_CL2_PRIOR, EPS_ES_PRIOR],
        bounds=[EPS_CL2_BOUNDS, EPS_ES_BOUNDS],
        method="L-BFGS-B",
        options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": 2000},
    )
    eps_cl2_opt, eps_es_opt = result.x
    return eps_cl2_opt, eps_es_opt


# ── Metrics ───────────────────────────────────────────────────────────────────

def compute_metrics(exp_abs, pred_abs):
    residuals = pred_abs - exp_abs
    mae  = np.mean(np.abs(residuals))
    rmse = np.sqrt(np.mean(residuals ** 2))
    nonzero = exp_abs != 0
    mape = np.mean(np.abs(residuals[nonzero] / exp_abs[nonzero])) * 100
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((exp_abs - exp_abs.mean()) ** 2)
    r2   = 1 - ss_res / ss_tot
    r, _ = pearsonr(exp_abs, pred_abs)
    mu_e, mu_p = exp_abs.mean(), pred_abs.mean()
    var_e, var_p = exp_abs.var(), pred_abs.var()
    cov  = np.cov(exp_abs, pred_abs, ddof=0)[0, 1]
    ccc  = (2 * cov) / (var_e + var_p + (mu_e - mu_p) ** 2)
    return dict(MAE=mae, RMSE=rmse, MAPE_pct=mape, R2=r2,
                Pearson_r=r, Lins_CCC=ccc, residuals=residuals)


def print_results(eps_cl2, eps_es, metrics, fit_mode):
    fmts   = {"MAE": ".4e", "RMSE": ".4e", "MAPE_pct": ".3f",
              "R2": ".6f", "Pearson_r": ".6f", "Lins_CCC": ".6f"}
    labels = {"MAE": "MAE", "RMSE": "RMSE", "MAPE_pct": "MAPE (%)",
              "R2": "R²", "Pearson_r": "Pearson r", "Lins_CCC": "Lin's CCC"}
    mode_str = {"cl2_only": "Fit ε(Cl₂˙⁻), fixed ε(eS⁻)",
                "both":     "Fit both ε(Cl₂˙⁻) and ε(eS⁻)",
                "fixed":    "All ε fixed (no fitting)"}[fit_mode]
    print(f"\n{'='*62}")
    print(f"  Composite Cl₂˙⁻ + eS⁻ Absorbance at {WAVELENGTH_NM} nm")
    print(f"{'='*62}")
    print(f"  Fitting mode      : {mode_str}")
    print(f"  ε(Cl₂˙⁻, {WAVELENGTH_NM} nm): {eps_cl2:.0f} M⁻¹cm⁻¹"
          f"  (bounds {EPS_CL2_BOUNDS[0]:.0f}–{EPS_CL2_BOUNDS[1]:.0f})")
    print(f"  ε(eS⁻,   {WAVELENGTH_NM} nm): {eps_es:.0f} M⁻¹cm⁻¹"
          f"  (bounds {EPS_ES_BOUNDS[0]:.0f}–{EPS_ES_BOUNDS[1]:.0f})")
    print(f"  Path length       : {PATH_LENGTH} cm")
    print(f"\n  {'Metric':<20}  {'Value':>14}")
    print(f"  {'-'*36}")
    for k, v in metrics.items():
        if k == "residuals":
            continue
        print(f"  {labels[k]:<20}  {v:{fmts[k]}}")
    print(f"{'='*62}\n")


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_results(exp_t, exp_abs,
                 sim_t, sim_cl2, sim_es,
                 c_cl2_interp, c_es_interp,
                 pred_abs, eps_cl2, eps_es,
                 metrics, fit_mode):

    diffs     = metrics["residuals"]
    means     = (exp_abs + pred_abs) / 2
    mean_diff = diffs.mean()
    sd_diff   = diffs.std(ddof=1)
    loa_hi    = mean_diff + 1.96 * sd_diff
    loa_lo    = mean_diff - 1.96 * sd_diff

    # Determine contribution fractions at each experimental time point
    A_cl2_part = PATH_LENGTH * eps_cl2 * c_cl2_interp
    A_es_part  = PATH_LENGTH * eps_es  * c_es_interp
    frac_es    = np.where(pred_abs > 0, A_es_part / pred_abs * 100, 0.0)

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    mode_str = {"cl2_only": f"Fit ε(Cl₂˙⁻), ε(eS⁻)={eps_es:.0f} fixed",
                "both":     "Both ε fitted",
                "fixed":    "All ε fixed"}[fit_mode]
    fig.suptitle(
        f"Composite Cl₂˙⁻ + eS⁻ at {WAVELENGTH_NM} nm  |  {mode_str}\n"
        f"ε(Cl₂˙⁻)={eps_cl2:.0f} M⁻¹cm⁻¹,  ε(eS⁻)={eps_es:.0f} M⁻¹cm⁻¹,  "
        f"l={PATH_LENGTH} cm",
        fontsize=10, fontweight="bold"
    )
    panel = 0

    # ── Panel 1: Full overlay — concentration (dual y) ────────────────────────
    ax = axes[panel]; panel += 1
    ax2 = ax.twinx()
    ax.plot(sim_t, sim_cl2, lw=1.2, color="steelblue",  alpha=0.6, label="Sim [Cl₂˙⁻]")
    ax.plot(sim_t, sim_es,  lw=1.2, color="darkorange",  alpha=0.6, label="Sim [eS⁻]")
    ax2.scatter(exp_t, exp_abs, s=30, zorder=4, color="tomato", label="Exp absorbance")
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel("Concentration (mol/L)", color="steelblue")
    ax2.set_ylabel(f"Absorbance ({WAVELENGTH_NM} nm)", color="tomato")
    ax.set_title("Simulated Concentrations + Exp Signal")
    ax.tick_params(axis="y", labelcolor="steelblue")
    ax2.tick_params(axis="y", labelcolor="tomato")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 2: Absorbance comparison with component breakdown ───────────────
    ax = axes[panel]; panel += 1
    ax.scatter(exp_t, exp_abs,  s=40, zorder=5, color="tomato",
               label="Experiment")
    ax.scatter(exp_t, pred_abs, s=40, zorder=5, color="black", marker="^",
               label=f"Sim total (ε_Cl2={eps_cl2:.0f}, ε_eS={eps_es:.0f})")
    ax.fill_between(exp_t, 0, A_cl2_part, alpha=0.25, color="steelblue",
                    label=f"Cl₂˙⁻ contribution")
    ax.fill_between(exp_t, A_cl2_part, pred_abs, alpha=0.25, color="darkorange",
                    label=f"eS⁻ contribution")
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel(f"Absorbance ({WAVELENGTH_NM} nm)")
    ax.set_title("Absorbance Comparison\n(stacked contributions)")
    ax.legend(fontsize=7)
    ax.grid(True, ls="--", alpha=0.4)
    txt = (f"R² = {metrics['R2']:.4f}\n"
           f"Lin's CCC = {metrics['Lins_CCC']:.4f}\n"
           f"RMSE = {metrics['RMSE']:.3e}\n"
           f"MAPE = {metrics['MAPE_pct']:.2f} %")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=8,
            va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", alpha=0.85))

    # ── Panel 3: eS- fractional contribution over time ────────────────────────
    ax = axes[panel]; panel += 1
    ax.plot(exp_t, frac_es, "o-", color="darkorange", lw=1.5, ms=5)
    ax.axhline(50, color="gray", ls=":", lw=0.8, label="50% line")
    ax.axhline(10, color="gray", ls="--", lw=0.8, label="10% line")
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel("eS⁻ fraction of A_pred (%)")
    ax.set_title("eS⁻ Fractional Contribution\nto Predicted Absorbance")
    ax.set_ylim(-2, 102)
    ax.legend(fontsize=8)
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 4: Residuals vs time ────────────────────────────────────────────
    ax = axes[panel]; panel += 1
    ax.scatter(exp_t, diffs, s=35, color="steelblue", zorder=4)
    ax.axhline(0,         color="gray",      lw=0.8, ls=":")
    ax.axhline(mean_diff, color="steelblue", lw=1.5, ls="--",
               label=f"Mean = {mean_diff:.3e}")
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel("Sim − Exp (absorbance)")
    ax.set_title("Residuals vs Time")
    ax.legend(fontsize=8)
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 5: Bland-Altman ─────────────────────────────────────────────────
    ax = axes[panel]; panel += 1
    ax.scatter(means, diffs, s=35, color="steelblue", zorder=4)
    ax.axhline(mean_diff, color="steelblue", lw=1.8, ls="-",
               label=f"Bias = {mean_diff:.3e}")
    ax.axhline(loa_hi,    color="tomato",    lw=1.4, ls="--",
               label=f"+1.96 SD = {loa_hi:.3e}")
    ax.axhline(loa_lo,    color="tomato",    lw=1.4, ls="--",
               label=f"−1.96 SD = {loa_lo:.3e}")
    ax.axhline(0,         color="gray",      lw=0.8, ls=":")
    ax.set_xlabel("Mean Absorbance")
    ax.set_ylabel("Sim − Exp Absorbance")
    ax.set_title("Bland-Altman Plot")
    ax.legend(fontsize=7)
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 6: Derived [Cl2•-] concentration comparison ────────────────────
    ax = axes[panel]; panel += 1
    exp_abs_corrected = exp_abs - PATH_LENGTH * eps_es * c_es_interp
    exp_cl2_derived   = exp_abs_corrected / (eps_cl2 * PATH_LENGTH)
    ax.scatter(exp_t, exp_cl2_derived, s=40, zorder=4, color="tomato",
               label="Exp (eS⁻-corrected)")
    ax.scatter(exp_t, c_cl2_interp,   s=40, zorder=4, color="steelblue",
               marker="^", alpha=0.85, label="Sim [Cl₂˙⁻]")
    ax.plot(sim_t, sim_cl2, lw=1.0, color="steelblue", alpha=0.4)
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel("[Cl₂˙⁻] (mol/L)")
    ax.set_title("Derived [Cl₂˙⁻] Comparison\n(exp corrected for eS⁻)")
    ax.legend(fontsize=7)
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 7: Derived [eS-] concentration comparison ──────────────────────
    ax = axes[panel]; panel += 1
    exp_abs_cl2_removed = exp_abs - PATH_LENGTH * eps_cl2 * c_cl2_interp
    exp_es_derived      = exp_abs_cl2_removed / (eps_es * PATH_LENGTH)
    ax.scatter(exp_t, exp_es_derived, s=40, zorder=4, color="darkorange",
               label="Exp (Cl₂˙⁻-corrected)")
    ax.scatter(exp_t, c_es_interp,   s=40, zorder=4, color="purple",
               marker="^", alpha=0.85, label="Sim [eS⁻]")
    ax.plot(sim_t, sim_es, lw=1.0, color="purple", alpha=0.4)
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel("[eS⁻] (mol/L)")
    ax.set_title("Derived [eS⁻] Comparison\n(exp corrected for Cl₂˙⁻)")
    ax.legend(fontsize=7)
    ax.grid(True, ls="--", alpha=0.4)

    # ── Panel 8: Summary annotation ───────────────────────────────────────────
    ax = axes[panel]; panel += 1
    ax.axis("off")
    summary = (
        f"Composite Fit Summary\n"
        f"{'─'*32}\n"
        f"Measurement λ : {WAVELENGTH_NM} nm\n"
        f"Path length   : {PATH_LENGTH} cm\n"
        f"Fit mode      : {fit_mode}\n\n"
        f"ε(Cl₂˙⁻, {WAVELENGTH_NM} nm) = {eps_cl2:.0f} M⁻¹cm⁻¹\n"
        f"  Prior       : {EPS_CL2_PRIOR:.0f} M⁻¹cm⁻¹\n"
        f"  Bounds      : {EPS_CL2_BOUNDS[0]:.0f}–{EPS_CL2_BOUNDS[1]:.0f} M⁻¹cm⁻¹\n\n"
        f"ε(eS⁻, {WAVELENGTH_NM} nm)   = {eps_es:.0f} M⁻¹cm⁻¹\n"
        f"  Prior       : {EPS_ES_PRIOR:.0f} M⁻¹cm⁻¹\n"
        f"  Bounds      : {EPS_ES_BOUNDS[0]:.0f}–{EPS_ES_BOUNDS[1]:.0f} M⁻¹cm⁻¹\n\n"
        f"R²            = {metrics['R2']:.4f}\n"
        f"Lin's CCC     = {metrics['Lins_CCC']:.4f}\n"
        f"RMSE          = {metrics['RMSE']:.3e}\n"
        f"MAPE          = {metrics['MAPE_pct']:.2f} %\n"
        f"Pearson r     = {metrics['Pearson_r']:.4f}"
    )
    ax.text(0.05, 0.97, summary, transform=ax.transAxes, fontsize=9,
            va="top", ha="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", alpha=0.9))

    plt.tight_layout()
    plt.savefig("comparison_plots.png", dpi=150, bbox_inches="tight")
    print("Plot saved to comparison_plots.png")
    plt.show()


# ── Export ────────────────────────────────────────────────────────────────────

def export_results(exp_t, exp_abs, pred_abs, c_cl2_interp, c_es_interp,
                   eps_cl2, eps_es, metrics, fit_mode):
    """Write per-time-point results and summary metrics to CSV."""
    # Per-point CSV
    A_cl2_part = PATH_LENGTH * eps_cl2 * c_cl2_interp
    A_es_part  = PATH_LENGTH * eps_es  * c_es_interp
    df_pts = pd.DataFrame({
        "time_ns":          exp_t,
        "A_exp":            exp_abs,
        "A_pred_total":     pred_abs,
        "A_pred_Cl2rad":    A_cl2_part,
        "A_pred_eS":        A_es_part,
        "C_Cl2rad_mol_L":   c_cl2_interp,
        "C_eS_mol_L":       c_es_interp,
        "residual":         metrics["residuals"],
    })
    df_pts.to_csv("comparison_timeseries.csv", index=False)
    print("Time-series results saved to comparison_timeseries.csv")

    # Summary metrics CSV
    row = {
        "fit_mode":             fit_mode,
        "wavelength_nm":        WAVELENGTH_NM,
        "eps_Cl2_M-1cm-1":     eps_cl2,
        "eps_eS_M-1cm-1":      eps_es,
        "path_length_cm":       PATH_LENGTH,
        **{k: v for k, v in metrics.items() if k != "residuals"},
    }
    pd.DataFrame([row]).to_csv("metrics_summary.csv", index=False)
    print("Metrics saved to metrics_summary.csv")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 62)
    print(f"  Composite Cl₂˙⁻ + eS⁻ Comparison at {WAVELENGTH_NM} nm")
    print(f"  FIT_MODE = {FIT_MODE}")
    print("=" * 62 + "\n")

    # ── Load experimental data ────────────────────────────────────────────────
    print(f"Loading experimental data : {EXP_CSV}")
    exp_t, exp_abs = load_exp_csv(EXP_CSV)
    print(f"  {len(exp_t)} points | time [{exp_t.min():.4g}, {exp_t.max():.4g}] ns")

    # ── Load simulation data ──────────────────────────────────────────────────
    print(f"Loading simulation data   : {SIM_CSV}")
    sim_t_raw, sim_cl2, sim_es = load_sim_csv(SIM_CSV)
    sim_t = sim_t_raw * TIME_SCALE_FACTOR - Sim_Time_Offset
    print(f"  {len(sim_t)} points | time [{sim_t.min():.4g}, {sim_t.max():.4g}] ns")
    print(f"  [Cl₂˙⁻] range: [{sim_cl2.min():.3e}, {sim_cl2.max():.3e}] mol/L")
    print(f"  [eS⁻]   range: [{sim_es.min():.3e},  {sim_es.max():.3e}] mol/L")

    # ── Domain check ─────────────────────────────────────────────────────────
    overlap_lo = max(exp_t.min(), sim_t.min())
    overlap_hi = min(exp_t.max(), sim_t.max())
    if overlap_lo >= overlap_hi:
        print(f"\n  ✗ FATAL: No time overlap between experiment and simulation.")
        print(f"    → Increase t_final (recommend ≥ 200 ns).")
        raise SystemExit(1)
    n_in = int(np.sum((exp_t >= sim_t.min()) & (exp_t <= sim_t.max())))
    print(f"  Overlap: [{overlap_lo:.4g}, {overlap_hi:.4g}] ns "
          f"({n_in}/{len(exp_t)} exp points covered)")

    # ── Splines ───────────────────────────────────────────────────────────────
    print("\nFitting splines to simulation concentrations...")
    cs_cl2, cs_es = make_splines(sim_t, sim_cl2, sim_es)

    # ── Coverage mask + interpolation ─────────────────────────────────────────
    print("\nChecking experimental point coverage...")
    mask = check_coverage(exp_t, sim_t)
    exp_t_m      = exp_t[mask]
    exp_abs_m    = exp_abs[mask]
    c_cl2_interp = np.maximum(0.0, cs_cl2(exp_t_m))
    c_es_interp  = np.maximum(0.0, cs_es(exp_t_m))

    # ── Fit / set ε values ────────────────────────────────────────────────────
    if FIT_MODE == "cl2_only":
        eps_es = EPS_ES_PRIOR
        print(f"\nFitting ε(Cl₂˙⁻) with ε(eS⁻) fixed at {eps_es:.0f} M⁻¹cm⁻¹...")
        eps_cl2 = fit_eps_cl2_only(exp_abs_m, c_cl2_interp, c_es_interp, eps_es)
        print(f"  Best-fit ε(Cl₂˙⁻, {WAVELENGTH_NM} nm) = {eps_cl2:.0f} M⁻¹cm⁻¹")

    elif FIT_MODE == "both":
        print(f"\nFitting both ε(Cl₂˙⁻) and ε(eS⁻) simultaneously...")
        eps_cl2, eps_es = fit_both_eps(exp_abs_m, c_cl2_interp, c_es_interp)
        print(f"  Best-fit ε(Cl₂˙⁻, {WAVELENGTH_NM} nm) = {eps_cl2:.0f} M⁻¹cm⁻¹")
        print(f"  Best-fit ε(eS⁻,   {WAVELENGTH_NM} nm) = {eps_es:.0f} M⁻¹cm⁻¹")

    else:  # "fixed"
        eps_cl2 = EPS_CL2_PRIOR
        eps_es  = EPS_ES_PRIOR
        print(f"\nUsing fixed ε values: ε(Cl₂˙⁻)={eps_cl2:.0f}, ε(eS⁻)={eps_es:.0f} M⁻¹cm⁻¹")

    # ── Compute composite absorbance and metrics ───────────────────────────────
    pred_abs = composite_absorbance(eps_cl2, eps_es, c_cl2_interp, c_es_interp)
    metrics  = compute_metrics(exp_abs_m, pred_abs)
    print_results(eps_cl2, eps_es, metrics, FIT_MODE)

    # ── Export ────────────────────────────────────────────────────────────────
    export_results(exp_t_m, exp_abs_m, pred_abs, c_cl2_interp, c_es_interp,
                   eps_cl2, eps_es, metrics, FIT_MODE)

    # ── Plot ──────────────────────────────────────────────────────────────────
    plot_results(
        exp_t_m, exp_abs_m,
        sim_t, sim_cl2, sim_es,
        c_cl2_interp, c_es_interp,
        pred_abs, eps_cl2, eps_es,
        metrics, FIT_MODE,
    )


if __name__ == "__main__":
    main()