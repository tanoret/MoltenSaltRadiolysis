"""Post-process MC results: bands, residual metrics, source contributions, figures."""
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
import run_cases as rc
import uq
 
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
MODEL = "#2a78d6"           # categorical slot 1
CL2 = "#eb6834"             # slot 2 (Cl2•- contribution, 400 nm only)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
 
def obs_sigma(t, a):
    """Robust point-scatter estimate: deviation of each interior point from the line through its neighbours."""
    t, idx = np.unique(t, return_index=True); a = a[idx]
    w = (t[1:-1] - t[:-2]) / (t[2:] - t[:-2])
    d = a[1:-1] - ((1 - w) * a[:-2] + w * a[2:])
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(1 + w**2 + (1 - w)**2).mean()
 
nom_eps, nom_P, _ = uq.run_mc([], 1, 0)
nom_P = [p[0] for p in nom_P]
J = np.load("uq_joint.npz")
eps_j = J["eps"]
Pj = [J[f"P{i}"] for i in range(len(rc.CASES))]
G = np.load("uq_groups.npz")
groups = uq.KINETIC_GROUPS + uq.OPTICAL_GROUPS
 
# ---------------- per-trace metrics and bands
rows, bands = [], []
for i, (lab, lam, metal, c, f, role) in enumerate(rc.CASES):
    t, A = uq.EXP[i]
    pn = nom_P[i]
    lo, hi = np.percentile(Pj[i], [2.5, 97.5], axis=0)
    sd_model = Pj[i].std(axis=0)
    s_obs = obs_sigma(t, A)
    sd_pred = np.sqrt(sd_model**2 + s_obs**2)
    res = A - pn
    m = rc.metrics(A, pn)
    cover = np.mean(np.abs(res) <= 1.96 * sd_pred)
    rows.append(dict(case=lab, role=role, n=m["n"], RMSE=m["RMSE"], NRMSE_pct=100 * m["RMSE"] / A.max(),
                     R2=m["R2"], CCC=m["CCC"], bias=m["bias"], sigma_obs=s_obs,
                     mean_model_sd_pct=100 * sd_model.mean() / A.max(),
                     chi2_red=float(np.mean((res / sd_pred) ** 2)), coverage95=cover))
    bands.append((t, A, pn, lo, hi, sd_pred, res, s_obs))
met = pd.DataFrame(rows)
met.to_csv("uq_trace_metrics.csv", index=False)
 
pooled = []
for lam in (700, 400):
    for role in ("cal", "hold"):
        ii = [i for i, cse in enumerate(rc.CASES) if cse[1] == lam and cse[5] == role]
        A = np.concatenate([bands[i][1] for i in ii]); P = np.concatenate([bands[i][2] for i in ii])
        sd = np.concatenate([bands[i][5] for i in ii])
        m = rc.metrics(A, P)
        pooled.append(dict(wavelength=lam, set=role, **m,
                           coverage95=float(np.mean(np.abs(A - P) <= 1.96 * sd))))
pooled = pd.DataFrame(pooled)
pooled.to_csv("uq_pooled_metrics.csv", index=False)
 
# ---------------- source contributions (group-wise OAT MC)
contrib = []
for lam in (700, 400):
    ii = [i for i, cse in enumerate(rc.CASES) if cse[1] == lam]
    Amax = {i: uq.EXP[i][1].max() for i in ii}
    for g in groups:
        V = np.concatenate([(G[f"{g}__P{i}"] / Amax[i]).var(axis=0) for i in ii])
        e = G[f"{g}__eps"][:, 0 if lam == 700 else 1]
        contrib.append(dict(wavelength=lam, group=g, rel_sd_pred_pct=100 * np.sqrt(V.mean()),
                            eps_e_sd_pct=100 * e.std() / e.mean()))
    Vj = np.concatenate([(Pj[i] / Amax[i]).var(axis=0) for i in ii])
    e = eps_j[:, 0 if lam == 700 else 1]
    contrib.append(dict(wavelength=lam, group="JOINT", rel_sd_pred_pct=100 * np.sqrt(Vj.mean()),
                        eps_e_sd_pct=100 * e.std() / e.mean()))
contrib = pd.DataFrame(contrib)
for lam in (700, 400):
    sub = contrib[(contrib.wavelength == lam) & (contrib.group != "JOINT")]
    tot = (sub.rel_sd_pred_pct ** 2).sum()
    contrib.loc[sub.index, "variance_share_pct"] = 100 * sub.rel_sd_pred_pct ** 2 / tot
contrib.to_csv("uq_source_contributions.csv", index=False)
 
eps_summary = {lam: dict(nominal=float(nom_eps[0][k]), mean=float(eps_j[:, k].mean()),
                         sd=float(eps_j[:, k].std()),
                         p2_5=float(np.percentile(eps_j[:, k], 2.5)),
                         p97_5=float(np.percentile(eps_j[:, k], 97.5)))
               for k, lam in enumerate((700, 400))}
params = json.loads(str(J["params"]))
GD = np.array([p["G"] * p["D"] for p in params])
eps_summary["G*eps_700"] = dict(mean=float((np.array([p["G"] for p in params]) * eps_j[:, 0]).mean()),
                                sd=float((np.array([p["G"] for p in params]) * eps_j[:, 0]).std()))
eps_summary["corr(eps700, G*D)"] = float(np.corrcoef(eps_j[:, 0], GD)[0, 1])
eps_summary["corr(eps400, G*D)"] = float(np.corrcoef(eps_j[:, 1], GD)[0, 1])
json.dump(eps_summary, open("uq_eps_summary.json", "w"), indent=2)
 
# ---------------- figures: predicted vs measured absorbance + residuals, one figure per wavelength
def trace_grid(lam, fname):
    series = ["Cr2+", "Cr3+"]
    ii = {s: [i for i, cse in enumerate(rc.CASES) if cse[1] == lam and cse[2] == s] for s in series}
    ncol = max(len(v) for v in ii.values())
    fig = plt.figure(figsize=(2.3 * ncol, 6.6))
    outer = fig.add_gridspec(2, ncol, hspace=0.45, wspace=0.28)
    ymax = max(max(bands[i][1].max(), bands[i][2].max()) for v in ii.values() for i in v) * 1.08
    rmax = max(np.abs(bands[i][6]).max() for v in ii.values() for i in v) * 1.15
    for r, s in enumerate(series):
        for cidx, i in enumerate(ii[s]):
            lab, _, _, conc, _, role = rc.CASES[i]
            t, A, pn, lo, hi, sdp, res, s_obs = bands[i]
            sub = outer[r, cidx].subgridspec(2, 1, height_ratios=[3, 1.2], hspace=0.18)
            ax = fig.add_subplot(sub[0]); axr = fig.add_subplot(sub[1], sharex=ax)
            ax.fill_between(t, lo, hi, color=MODEL, alpha=0.18, lw=0, label="model 95% band")
            ax.plot(t, pn, color=MODEL, lw=2, label="model (frozen ε)")
            ax.plot(t, A, "o", ms=3.2, mfc="white", mec=INK, mew=0.9, label="measured")
            if lam == 400:
                ce, cc, _ = rc.simulate_at(rc.case_config(s, conc), t)
                ax.plot(t, rc.L_CM * 3400 * cc, color=CL2, lw=1.2, ls="--", label="Cl₂•⁻ part")
            ax.set_ylim(0, ymax); ax.grid(axis="y", color=GRID, lw=0.6)
            m = met.iloc[i]
            ax.set_title(f"{conc:.2f} mM {s.replace('+','⁺').replace('2','²').replace('3','³')}"
                         f"  ·  {'held-out' if role=='hold' else 'calibration'}", fontsize=8.5, color=INK)
            ax.text(0.97, 0.95, f"R² {m.R2:.3f}\nNRMSE {m.NRMSE_pct:.1f}%", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5, color=INK2)
            plt.setp(ax.get_xticklabels(), visible=False)
            axr.axhspan(-1.96 * s_obs, 1.96 * s_obs, color=MUTED, alpha=0.18, lw=0)
            axr.axhline(0, color=MUTED, lw=0.8)
            axr.plot(t, res, "o", ms=2.6, color=INK)
            axr.set_ylim(-rmax, rmax); axr.set_xlabel("Time after pulse (ns)")
            if cidx == 0:
                ax.set_ylabel("Absorbance"); axr.set_ylabel("Resid.")
            else:
                plt.setp(ax.get_yticklabels(), visible=False); plt.setp(axr.get_yticklabels(), visible=False)
            if r == 0 and cidx == 0:
                h, l = ax.get_legend_handles_labels()
    h.append(matplotlib.patches.Patch(color=MUTED, alpha=0.3)); l.append("residual = measured − model; grey = ±1.96σ data scatter")
    fig.legend(h, l, loc="upper center", ncol=len(l), frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    ep = eps_summary[lam]
    fig.suptitle(f"{lam} nm: measured vs predicted absorbance, ε(eₛ⁻) = {ep['nominal']:.0f} M⁻¹cm⁻¹ "
                 f"(one value, calibrated on ~1/3/5 mM, frozen)", y=1.035, fontsize=10, color=INK)
    fig.savefig(fname, dpi=200, bbox_inches="tight"); plt.close(fig)
 
trace_grid(700, "fig_700nm_pred_vs_meas.png")
trace_grid(400, "fig_400nm_pred_vs_meas.png")
 
# ---------------- source-contribution figure
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
order = groups
names = {"dose": "Pulse dose (15–30 Gy)", "G_level": "G-value (Gε reading + ε transfer)",
         "G_ratio": "G(Cl₂•⁻)/G(eₛ⁻) ratio", "k_measured": "Measured k (R8, R9, R11, R12)",
         "k_R5": "k R5 (Cl₂•⁻ disproportionation)", "k_R10": "k R10 (comproportionation)",
         "k_closure": "Estimated k (R1–R4, R6, R7)", "background": "Background decay rates",
         "eps_Cl2_400": "ε(Cl₂•⁻, 400 nm)", "path_length": "Path length"}
for ax, lam in zip(axes, (700, 400)):
    sub = contrib[(contrib.wavelength == lam) & (contrib.group != "JOINT")].set_index("group").loc[order]
    y = np.arange(len(order))[::-1]
    ax.barh(y, sub.rel_sd_pred_pct, color=MODEL, height=0.6)
    for yy, v in zip(y, sub.rel_sd_pred_pct):
        ax.text(v + 0.05, yy, f"{v:.2f}", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([names[g] for g in order])
    jt = contrib[(contrib.wavelength == lam) & (contrib.group == "JOINT")].rel_sd_pred_pct.iloc[0]
    ax.set_title(f"{lam} nm  (all sources jointly: {jt:.2f}%)", fontsize=9, color=INK)
    ax.set_xlabel("Prediction SD, % of trace peak"); ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_xlim(0, contrib.rel_sd_pred_pct.max() * 1.25)
fig.tight_layout(); fig.savefig("fig_uq_sources.png", dpi=200, bbox_inches="tight"); plt.close(fig)
 
pd.set_option("display.width", 220)
print(met.to_string(float_format=lambda x: f"{x:.4g}"))
print(pooled.to_string(float_format=lambda x: f"{x:.4g}"))
print(contrib.to_string(float_format=lambda x: f"{x:.3g}"))
print(json.dumps(eps_summary, indent=1))
