"""
Run all 19 Iwamatsu et al. (2026) traces with direct_integration.py and compare
predicted vs measured absorbance.

Observation model:  A(t) = l * [eps_e(lambda) * C_e(t) + eps_Cl2(lambda) * C_Cl2(t)]
  l = 0.5 cm; eps_Cl2(700) = 0; eps_Cl2(400) = 3400 M^-1 cm^-1 (fixed, aqueous, Hug 1981)
  eps_e(lambda) is the only optical parameter estimated from data.

Two estimation modes:
  per-trace : eps_e fitted to each trace separately (as in the original manuscript)
  global    : one eps_e per wavelength, fitted on the CALIBRATION set (Cr ~1, ~3, ~5 mM,
              both oxidation states), frozen, then scored on the HELD-OUT set (~2, ~4 mM).
Since A is linear in eps_e, the least-squares estimate is closed form:
  eps_e = sum(C_e * (A - l*eps_Cl2*C_Cl2)) / (l * sum(C_e^2))
"""
import json, copy, sys
import numpy as np, pandas as pd
from scipy.interpolate import CubicSpline
import direct_integration as di

L_CM = 0.5
EPS_CL2 = {700: 0.0, 400: 3400.0}
DATA = "data/"

BASE = {
    "kernel": "chloride",
    "temperature_K": 673.15,
    "radiation": "gamma",
    "pulse_dose_Gy": 22.5,          # midpoint of 15-30 Gy (Iwamatsu 2026)
    "density_g_cm3": 1.677,         # ESI Table S1, 400 C
    "pulse_width_s": 0.0,           # LEAF target B: < 120 ps (Wishart et al. 2004) -> instantaneous
    "G_values_override": {"e_s-": 2.8, "Cl2•-": 2.8},
    "initial_concentrations": {"Cl-": 30.0},   # LiCl-KCl eutectic 58:42 mol% at 1.677 g/cm3
    "metals": {"chromium": {"Cr2+": 0.0, "Cr3+": 0.0, "Cr+": 0.0}},
}
# Background pseudo-first-order losses: intercepts of k_obs vs [Cr], ESI Tables S2-S4
BACKGROUND = {
    "Cr2+": {"e_s-": 2.2e7, "Cl2•-": 3.3e6},
    "Cr3+": {"e_s-": 0.0,   "Cl2•-": 1.9e6},
}

# (label, lambda, metal, [Cr] mM, file, role)
CASES = [
    ("700 Cr2+ 0.99", 700, "Cr2+", 0.99, "absorbance1mMCr2.csv", "cal"),
    ("700 Cr2+ 2.09", 700, "Cr2+", 2.09, "absorbance2mMCr2.csv", "hold"),
    ("700 Cr2+ 3.03", 700, "Cr2+", 3.03, "absorbance3mMCr2.csv", "cal"),
    ("700 Cr2+ 4.00", 700, "Cr2+", 4.00, "absorbance4mMCr2.csv", "hold"),
    ("700 Cr3+ 1.05", 700, "Cr3+", 1.05, "absorbance1mMCr3.csv", "cal"),
    ("700 Cr3+ 2.02", 700, "Cr3+", 2.02, "absorbance2mMCr3.csv", "hold"),
    ("700 Cr3+ 2.86", 700, "Cr3+", 2.86, "absorbance3mMCr3.csv", "cal"),
    ("700 Cr3+ 3.81", 700, "Cr3+", 3.81, "absorbance4mMCr3.csv", "hold"),
    ("700 Cr3+ 5.00", 700, "Cr3+", 5.00, "absorbance5mMCr3.csv", "cal"),
    ("400 Cr2+ 0.99", 400, "Cr2+", 0.99, "Cl2Absorbance1mMCr2.csv", "cal"),
    ("400 Cr2+ 2.01", 400, "Cr2+", 2.01, "Cl2Absorbance2mMCr2.csv", "hold"),
    ("400 Cr2+ 3.01", 400, "Cr2+", 3.01, "Cl2Absorbance3mMCr2.csv", "cal"),
    ("400 Cr2+ 4.03", 400, "Cr2+", 4.03, "Cl2Absorbance4mMCr2.csv", "hold"),
    ("400 Cr2+ 4.99", 400, "Cr2+", 4.99, "Cl2Absorbance5mMCr2.csv", "cal"),
    ("400 Cr3+ 1.00", 400, "Cr3+", 1.00, "400nmAbsorbance1mMCr3.csv", "cal"),
    ("400 Cr3+ 2.02", 400, "Cr3+", 2.02, "400nmAbsorbance2mMCr3.csv", "hold"),
    ("400 Cr3+ 2.98", 400, "Cr3+", 2.98, "400nmAbsorbance3mMCr3.csv", "cal"),
    ("400 Cr3+ 3.91", 400, "Cr3+", 3.91, "400nmAbsorbance4mMCr3.csv", "hold"),
    ("400 Cr3+ 4.97", 400, "Cr3+", 4.97, "400nmAbsorbance5mMCr3.csv", "cal"),
]

def load_exp(f):
    d = pd.read_csv(DATA + f, skipinitialspace=True)
    x, y = d.iloc[:, 0].to_numpy(float), d.iloc[:, 1].to_numpy(float)
    o = np.argsort(x)
    return x[o], y[o]

def case_config(metal, conc_mM, base=None, background=None):
    cfg = copy.deepcopy(base or BASE)
    cfg["metals"]["chromium"][metal] = conc_mM * 1e-3
    bg = (background or BACKGROUND)[metal]
    cfg["background_first_order_s^-1"] = dict(bg)
    return cfg

def simulate_at(cfg, exp_t_ns, n_out=20000):
    """Concentrations of e_s- and Cl2•- at experimental times. Exp t=0 is the end of the pulse."""
    sysm = di.build_system(cfg)
    tau_ns = cfg.get("pulse_width_s", 0.0) * 1e9
    t, Y, info = di.integrate_system(sysm, (exp_t_ns.max() + tau_ns + 1) * 1e-9, n_out=n_out)
    ts = t * 1e9 - tau_ns
    i = sysm.species_index
    ce = np.maximum(0, CubicSpline(ts, Y[:, i["e_s-"]])(exp_t_ns))
    cc = np.maximum(0, CubicSpline(ts, Y[:, i["Cl2•-"]])(exp_t_ns))
    return ce, cc, (t, Y, sysm, info)

def metrics(a, p):
    r = p - a
    cov = np.cov(a, p, ddof=0)[0, 1]
    return dict(n=len(a), RMSE=float(np.sqrt(np.mean(r**2))),
                R2=float(1 - np.sum(r**2) / np.sum((a - a.mean())**2)),
                r=float(np.corrcoef(a, p)[0, 1]),
                CCC=float(2 * cov / (a.var() + p.var() + (a.mean() - p.mean())**2)),
                bias=float(r.mean()))

def ls_eps(A, ce, cc, lam):
    y = A - L_CM * EPS_CL2[lam] * cc
    return float(ce @ y / (L_CM * (ce @ ce)))

def run(base=None, background=None, verbose=True):
    sims = []
    for lab, lam, metal, c, f, role in CASES:
        et, ea = load_exp(f)
        ce, cc, _ = simulate_at(case_config(metal, c, base, background), et)
        sims.append(dict(label=lab, lam=lam, metal=metal, conc=c, role=role, t=et, A=ea, ce=ce, cc=cc))
    glob = {}
    for lam in (700, 400):
        cal = [s for s in sims if s["lam"] == lam and s["role"] == "cal"]
        A = np.concatenate([s["A"] for s in cal]); ce = np.concatenate([s["ce"] for s in cal])
        cc = np.concatenate([s["cc"] for s in cal])
        glob[lam] = ls_eps(A, ce, cc, lam)
    rows = []
    for s in sims:
        lam = s["lam"]
        e_own = ls_eps(s["A"], s["ce"], s["cc"], lam)
        pred_g = L_CM * (glob[lam] * s["ce"] + EPS_CL2[lam] * s["cc"])
        pred_o = L_CM * (e_own * s["ce"] + EPS_CL2[lam] * s["cc"])
        mg, mo = metrics(s["A"], pred_g), metrics(s["A"], pred_o)
        s["pred_global"] = pred_g
        rows.append(dict(case=s["label"], role=s["role"], eps_own=round(e_own), R2_own=mo["R2"],
                         eps_global=round(glob[lam]), **{f"{k}_global": v for k, v in mg.items()}))
    df = pd.DataFrame(rows)
    pooled = {}
    for lam in (700, 400):
        for role in ("cal", "hold"):
            ss = [s for s in sims if s["lam"] == lam and s["role"] == role]
            pooled[(lam, role)] = metrics(np.concatenate([s["A"] for s in ss]),
                                          np.concatenate([s["pred_global"] for s in ss]))
    if verbose:
        pd.set_option("display.width", 200)
        print(df.to_string(float_format=lambda x: f"{x:.4g}"))
        print("\nGlobal eps_e:", {k: round(v) for k, v in glob.items()})
        for k, v in pooled.items():
            print("pooled", k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    return df, glob, pooled, sims

if __name__ == "__main__":
    df, glob, pooled, sims = run()
    df.to_csv("results_revision2.csv", index=False)
