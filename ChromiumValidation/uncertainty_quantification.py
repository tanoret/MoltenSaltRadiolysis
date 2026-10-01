"""
Monte Carlo uncertainty propagation and group-wise sensitivity for the model.

Every MC sample:
  1. draws kinetic / source inputs, integrates all 19 traces (stores C_e, C_Cl2 at exp times),
  2. draws optical inputs (eps_Cl2(400), l),
  3. RE-CALIBRATES eps_e(lambda) on the calibration set (closed-form LS), exactly as in the
     nominal analysis (calibration inside the loop), and
  4. predicts absorbance for all traces with the frozen eps_e.
Prediction bands therefore include the eps_e re-calibration response to each input.

Input distributions (1 sigma unless stated) — see UQ_INPUTS below and revision_log.md.
"""
import sys, json, copy, time
import numpy as np
from multiprocessing import Pool
import run_cases as rc

rng_master = np.random.default_rng(20260924)

K_NOM = {"R1": 1e10, "R2": 1e10, "R3": 1e10, "R4": 1e10, "R5": 2.2e9, "R6": 1e4, "R7": 1e3,
         "R8": 7.2e9, "R9": 1.4e9, "R10": 1.7e10, "R11": 4.1e10, "R12": 6.1e10}

UQ_INPUTS = {
    # group: description
    "dose":        "D ~ U(15, 30) Gy, shared by all traces (Iwamatsu 2026: 15-30 Gy per pulse)",
    "G_level":     "G(e_s-) = 2.8 x Gε reading (U +/-11%: half a 5e3 minor tick on a 20e3-major axis, read value 22.5e3) x eps(Cl2•-,340) (logN 20%: 8000 vs 8800 aq. + melt transfer)",
    "G_ratio":     "G(Cl2•-)/G(e_s-) ~ logN(median 1, 10%)",
    "k_measured":  "R8, R9, R11, R12 ~ N(reported value, reported 1 sigma), Iwamatsu 2026",
    "k_R5":        "R5 ~ logN(median 2.2e9, sigma_ln 0.78)  (reported 2.2 +/- 2.0e9)",
    "k_R10":       "R10 ~ log-U(1.7e10/3, 1.7e10*3)  (estimate, no sigma reported)",
    "k_closure":   "R1-R4, R6, R7 ~ log-U(k/10, k*10)  (order-of-magnitude estimates)",
    "background":  "k_bg ~ N(Table 1 value, regression SE); negative draws set to 0 (median = nominal)",
    "eps_Cl2_400": "eps(Cl2•-, 400 nm) ~ logN(median 3400, 15%)  (aqueous, Hug 1981)",
    "path_length": "l ~ N(0.5, 0.01) cm",
}
# Nominal background rates = the Table 1 values used everywhere else (rc.BACKGROUND);
# 1-sigma = standard error of the ESI intercept regressions.
BG_SD = {"Cr2+": {"e_s-": 4.0e6, "Cl2•-": 9.9e5},
         "Cr3+": {"e_s-": 1.6e7, "Cl2•-": 3.6e5}}
KINETIC_GROUPS = ["dose", "G_level", "G_ratio", "k_measured", "k_R5", "k_R10", "k_closure", "background"]
OPTICAL_GROUPS = ["eps_Cl2_400", "path_length"]

def _truncnorm(rng, mu, sd, lo=0.0):
    while True:
        x = rng.normal(mu, sd)
        if x >= lo:
            return x

def draw_kinetic(rng, groups):
    """Return (base_cfg_changes, rate_overrides, background) for the requested groups; others nominal."""
    D, G, ratio = 22.5, 2.8, 1.0
    k = dict(K_NOM)
    bg = {m: dict(d) for m, d in rc.BACKGROUND.items()}
    if "dose" in groups:
        D = rng.uniform(15, 30)
    if "G_level" in groups:
        G = 2.8 * rng.uniform(0.89, 1.11) * rng.lognormal(0, 0.20)
    if "G_ratio" in groups:
        ratio = rng.lognormal(0, 0.10)
    if "k_measured" in groups:
        for t, sd in (("R8", 0.3e9), ("R9", 0.1e9), ("R11", 0.2e10), ("R12", 0.3e10)):
            k[t] = _truncnorm(rng, K_NOM[t], sd, 0.0)
    if "k_R5" in groups:
        k["R5"] = 2.2e9 * rng.lognormal(0, 0.78)   # sigma_ln matching a 91% relative SD
    if "k_R10" in groups:
        k["R10"] = 1.7e10 * 3 ** rng.uniform(-1, 1)
    if "k_closure" in groups:
        for t in ("R1", "R2", "R3", "R4", "R6", "R7"):
            k[t] = K_NOM[t] * 10 ** rng.uniform(-1, 1)
    if "background" in groups:
        # Normal centred on the nominal value; a negative draw is set to 0 (a rate cannot be negative).
        # This keeps the median of every distribution equal to its nominal value, including the
        # Cr3+ e- background whose nominal is 0.
        bg = {m: {s: max(0.0, rng.normal(rc.BACKGROUND[m][s], BG_SD[m][s])) for s in d}
              for m, d in BG_SD.items()}
    base = copy.deepcopy(rc.BASE)
    base["pulse_dose_Gy"] = D
    base["G_values_override"] = {"e_s-": G, "Cl2•-": G * ratio}
    base["rate_overrides"] = k
    return base, bg, dict(D=D, G=G, ratio=ratio, **k)

def draw_optical(rng, groups):
    eps_cl2 = 3400 * rng.lognormal(0, 0.15) if "eps_Cl2_400" in groups else 3400.0
    l = rng.normal(0.5, 0.01) if "path_length" in groups else 0.5
    return eps_cl2, l

EXP = [rc.load_exp(f) for (_, _, _, _, f, _) in rc.CASES]

def simulate_all(args):
    base, bg = args
    ce_all, cc_all = [], []
    for (lab, lam, metal, c, f, role), (et, ea) in zip(rc.CASES, EXP):
        ce, cc, _ = rc.simulate_at(rc.case_config(metal, c, base, bg), et)   # 20,000 output points, as in run_cases
        ce_all.append(ce); cc_all.append(cc)
    return ce_all, cc_all

def calibrate_predict(ce_all, cc_all, eps_cl2_400, l):
    """Closed-form eps_e per wavelength on the calibration set, then predict every trace."""
    epsc = {700: 0.0, 400: eps_cl2_400}
    eps_e = {}
    for lam in (700, 400):
        num = den = 0.0
        for i, (lab, lm, metal, c, f, role) in enumerate(rc.CASES):
            if lm == lam and role == "cal":
                y = EXP[i][1] - l * epsc[lam] * cc_all[i]
                num += ce_all[i] @ y; den += l * (ce_all[i] @ ce_all[i])
        eps_e[lam] = max(0.0, num / den)   # physical constraint eps >= 0 (1-parameter NNLS)
    preds = [l * (eps_e[lm] * ce_all[i] + epsc[lm] * cc_all[i])
             for i, (lab, lm, *_ ) in enumerate(rc.CASES)]
    return eps_e, preds

def run_mc(groups, n, seed, procs=2):
    rng = np.random.default_rng(seed)
    kin_groups = [g for g in groups if g in KINETIC_GROUPS]
    draws = [draw_kinetic(rng, kin_groups) for _ in range(n)]
    opt = [draw_optical(rng, groups) for _ in range(n)]
    if kin_groups:
        with Pool(procs) as p:
            sims = p.map(simulate_all, [(b, bg) for b, bg, _ in draws], chunksize=4)
    else:
        b0, bg0, _ = draw_kinetic(rng, [])
        nom = simulate_all((b0, bg0))
        sims = [nom] * n
    eps = np.zeros((n, 2)); P = [[] for _ in rc.CASES]
    for j, ((ce, cc), (ec, l)) in enumerate(zip(sims, opt)):
        e, preds = calibrate_predict(ce, cc, ec, l)
        eps[j] = e[700], e[400]
        for i, p in enumerate(preds):
            P[i].append(p)
    P = [np.array(p) for p in P]
    params = [d[2] | {"eps_Cl2_400": o[0], "l": o[1]} for d, o in zip(draws, opt)]
    run_mc.last_sims = sims
    return eps, P, params

if __name__ == "__main__":
    mode = sys.argv[1]
    n = int(sys.argv[2])
    t0 = time.time()
    if mode == "joint":
        eps, P, params = run_mc(KINETIC_GROUPS + OPTICAL_GROUPS, n, 1)
        CE = [np.array([s[0][i] for s in run_mc.last_sims]) for i in range(len(rc.CASES))]
        CC = [np.array([s[1][i] for s in run_mc.last_sims]) for i in range(len(rc.CASES))]
        np.savez("uq_joint.npz", eps=eps, params=json.dumps(params),
                 **{f"CE{i}": x for i, x in enumerate(CE)}, **{f"CC{i}": x for i, x in enumerate(CC)},
                 **{f"P{i}": p for i, p in enumerate(P)})
    elif mode == "onegroup":
        g = sys.argv[3]; gi = (KINETIC_GROUPS + OPTICAL_GROUPS).index(g)
        eps, P, _ = run_mc([g], n, 100 + gi)
        old = dict(np.load("uq_groups.npz"))
        old[f"{g}__eps"] = eps
        for i, p in enumerate(P):
            old[f"{g}__P{i}"] = p
        np.savez("uq_groups.npz", **old)
    elif mode == "groups":
        out = {}
        for gi, g in enumerate(KINETIC_GROUPS + OPTICAL_GROUPS):
            eps, P, _ = run_mc([g], n, 100 + gi)
            out[f"{g}__eps"] = eps
            for i, p in enumerate(P):
                out[f"{g}__P{i}"] = p
            print(g, "done", round(time.time() - t0), "s", flush=True)
        np.savez("uq_groups.npz", **out)
    print("finished", mode, n, round(time.time() - t0), "s")
