#!/usr/bin/env python3
"""Meta-hierarchical layer + leave-one-metal-out (LOMO) predictive demo.

Pulls together every metal-vs-(e_s-/e_aq-) second-order rate constant from the
digitized chloride/aqueous data:

  Iwamatsu 2026 Cr      (Cr3+ + e_s-,   LiCl-KCl, 5 temperatures)
  Iwamatsu 2022 Zn      (Zn2+ + e_s-,   LiCl-KCl, Arrhenius)
  Castro-Baldivieso Nd  (Nd3+ + e_s-,   LiCl-KCl, 5 temperatures)
  Rotermund 2024 Cf     (Cf3+ + e_aq-,  aqueous,  single temperature)
  Pikaev 1982            (M^(z+) + e_s-, NaCl/KCl/KBr at ~1073 K
                           for Zn, Cd, Tl, Ag, Ca, Sr, Ba)

Each observation is log10 k_obs.  A partial-pooling regression layer ties the
per-metal (log10 A_m, Ea_m) to chemistry features X_m = (charge z, ionic radius,
Pauling electronegativity, log10(z/r)) standardized across the metal set.

  log10 k_(m, h, T) = log10 A_m + b_h - (Ea_m / R) * (1/T) / ln(10)
  log10 A_m         = mu_A + X_m . beta_A + eps_A,m  ;  eps_A ~ N(0, sigma_A)
  Ea_m              = mu_E + X_m . beta_E + eps_E,m  ;  eps_E ~ N(0, sigma_E)
  b_h               ~ N(0, sigma_b)   (LiCl-KCl is the reference, b = 0)

Outputs:
  validation/meta_hier/metal_data.csv             (the assembled observation table)
  validation/meta_hier/chemistry_features.csv     (X_m, standardized)
  validation/meta_hier/full_posterior_chain.npy   (emcee chain, full fit)
  validation/meta_hier/posterior_summary.csv      (means / 90% CIs per param)
  validation/meta_hier/lomo_results.csv           (held-out vs predicted, per metal)
  manuscript/figures/fig_meta_hier.pdf            (4-panel: data, posterior, LOMO,
                                                   chemistry-feature loadings)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import emcee
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
except Exception:
    pass


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

R_GAS = 8.314462618  # J/(mol*K)
LN10 = np.log(10.0)

OUT = REPO / "validation" / "meta_hier"
OUT.mkdir(parents=True, exist_ok=True)
FIG_DIR = REPO / "manuscript" / "figures"

# ----------------------------------------------------------------------------
# 1.  Assemble the (metal, host, T, log10 k, sigma) observation table
# ----------------------------------------------------------------------------

def build_dataset() -> pd.DataFrame:
    rows = []

    # --- Iwamatsu 2026 Cr + refit Zn (5 T points from LiCl-KCl) -------------
    cr_path = REPO / "validation/cr_licl_kcl/iwamatsu_2026_pccp/data/k_vs_T_from_arrhenius.csv"
    if cr_path.exists():
        cr = pd.read_csv(cr_path, comment="#")
        for _, r in cr.iterrows():
            rxn = str(r.get("reaction", ""))
            k = float(r["k_M_inv_s"])
            sig_k = float(r["sigma_k_M_inv_s"])
            T = float(r["T_K"])
            if "Cr3+" in rxn:
                metal = "Cr3+"; src = "Iwamatsu 2026 PCCP eq6"
            elif "Zn2+" in rxn:
                metal = "Zn2+"; src = "Iwamatsu 2026 PCCP eq3 (refit)"
            else:
                continue  # skip Cr2+ rows (Cr2+ + e_s- -> Cr+, harder to compare)
            rows.append({
                "metal": metal, "host": "LiCl-KCl", "T_K": T,
                "log10_k": np.log10(k),
                "sigma_log10_k": max(0.02, sig_k / (k * LN10)),
                "source": src,
            })

    # --- Castro-Baldivieso 2026 Nd (5 T in LiCl-KCl) ------------------------
    nd_path = REPO / "validation/nd_licl_kcl/castro_baldivieso_2026_ic/data/k_vs_T.csv"
    if nd_path.exists():
        nd = pd.read_csv(nd_path, comment="#")
        for _, r in nd.iterrows():
            k = float(r["k_M_inv_s"])
            sig_k = float(r.get("sigma_k_M_inv_s", k * 0.1))
            rows.append({
                "metal": "Nd3+", "host": "LiCl-KCl", "T_K": float(r["T_K"]),
                "log10_k": np.log10(k),
                "sigma_log10_k": max(0.02, sig_k / (k * LN10)),
                "source": "Castro-Baldivieso 2026 Inorg. Chem. Table 1",
            })

    # --- Rotermund 2024 Cf (aqueous, single T) -------------------------------
    # Reference value from Eq. (2): k(Cf3+ + e_aq-) = (3.63 +/- 0.14)e10 M^-1 s^-1 at 295 K
    rows.append({
        "metal": "Cf3+", "host": "H2O", "T_K": 295.15,
        "log10_k": np.log10(3.63e10),
        "sigma_log10_k": 0.14 / (3.63 * LN10),
        "source": "Rotermund 2024 JPC A Eq. 2 (cited from Horne 2022)",
    })

    # --- Pikaev 1982 (multi-host, ~1073 K, multi-metal) ----------------------
    pk_path = REPO / "validation/alkali_halide_baseline/pikaev_1982_rpc/data/es_rate_constants.csv"
    if pk_path.exists():
        pk = pd.read_csv(pk_path, comment="#")
        for _, r in pk.iterrows():
            sol = str(r["solute"])
            # Skip matrix-cation reactions (Na+, K+) - those are the e_s- decay in
            # the pure melt and we don't model them in this hierarchical layer.
            if sol in ("Na+", "K+"):
                continue
            metal = sol
            host = str(r["melt"])
            T_K = float(r["T_K"])
            k = float(r["k_M1_s1"])
            log10k = np.log10(k)
            sig = float(r["k_sigma_rel"]) / LN10  # 0.30 -> ~0.13 in log10
            rows.append({
                "metal": metal, "host": host, "T_K": T_K,
                "log10_k": log10k, "sigma_log10_k": sig,
                "source": "Pikaev 1982 RPC Table 2",
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "metal_data.csv", index=False)
    return df


# ----------------------------------------------------------------------------
# 2.  Chemistry-feature table for each metal
# ----------------------------------------------------------------------------

CHEM_FEATURES = {
    # metal, z, r_pm, chi_Pauling, n_d (d-electrons in target oxidation state)
    "Cr3+":  (3,  62, 1.66, 3),
    "Zn2+":  (2,  74, 1.65, 10),
    "Nd3+":  (3,  98, 1.14, 0),      # f-electrons, treat n_d=0
    "Cf3+":  (3,  95, 1.30, 0),
    "Cd2+":  (2,  95, 1.69, 10),
    "Tl+":   (1, 150, 2.04, 10),
    "Ag+":   (1, 115, 1.93, 10),
    "Ca2+":  (2, 100, 1.04, 0),
    "Sr2+":  (2, 118, 0.95, 0),
    "Ba2+":  (2, 135, 0.89, 0),
}


def build_chem_features(metals: list[str]) -> tuple[np.ndarray, list[str]]:
    raw = np.array([
        [CHEM_FEATURES[m][0],                  # z (charge)
         CHEM_FEATURES[m][1],                  # r_pm
         CHEM_FEATURES[m][2],                  # chi
         np.log10(CHEM_FEATURES[m][0] / CHEM_FEATURES[m][1]),  # log10(z/r)
         ] for m in metals
    ])
    # Standardize per column
    mu = raw.mean(axis=0)
    sd = raw.std(axis=0) + 1e-9
    X = (raw - mu) / sd
    feature_names = ["z_std", "r_std", "chi_std", "log_z_over_r_std"]
    return X, feature_names


# ----------------------------------------------------------------------------
# 3.  Hierarchical Bayesian model + emcee sampler
# ----------------------------------------------------------------------------

def make_logpost(df: pd.DataFrame, metals_in_fit: list[str],
                  hosts_in_fit: list[str],
                  X_m: np.ndarray, feature_names: list[str]):
    """Build a closure log_posterior(theta) for emcee.

    Parameter layout (theta):
        [0]                       mu_A          (intercept of log10 A regression)
        [1 .. n_feat]             beta_A        (regression coeffs for log10 A)
        [1+n_feat]                sigma_A       (log10 A scatter across metals)
        [2+n_feat]                mu_E          (kJ/mol intercept for Ea)
        [3+n_feat .. 2+2*n_feat]  beta_E        (regression coeffs for Ea)
        [3+2*n_feat]              sigma_E       (kJ/mol scatter across metals)
        [4+2*n_feat .. 3+2*n_feat + n_m]    eps_A_m (one per metal)
        [4+2*n_feat + n_m .. 3+2*n_feat + 2*n_m]  eps_E_m (one per metal)
        [next]                    b_h for each non-reference host
        [next]                    sigma_obs
    Reference host: LiCl-KCl (b = 0 forced).
    """
    n_feat = X_m.shape[1]
    n_m = len(metals_in_fit)
    # Reference host is LiCl-KCl unless that's not present, then take the first
    if "LiCl-KCl" in hosts_in_fit:
        ref_host = "LiCl-KCl"
    else:
        ref_host = hosts_in_fit[0]
    non_ref_hosts = [h for h in hosts_in_fit if h != ref_host]
    n_h = len(non_ref_hosts)
    # Indices
    iA0 = 0
    iAbeta = slice(1, 1 + n_feat)
    isigmaA = 1 + n_feat
    iE0 = 2 + n_feat
    iEbeta = slice(3 + n_feat, 3 + 2 * n_feat)
    isigmaE = 3 + 2 * n_feat
    ieps_A = slice(4 + 2 * n_feat, 4 + 2 * n_feat + n_m)
    ieps_E = slice(4 + 2 * n_feat + n_m, 4 + 2 * n_feat + 2 * n_m)
    iB = slice(4 + 2 * n_feat + 2 * n_m, 4 + 2 * n_feat + 2 * n_m + n_h)
    isigmaO = 4 + 2 * n_feat + 2 * n_m + n_h
    ndim = isigmaO + 1

    # Pre-compute index arrays for the data
    m_idx = np.array([metals_in_fit.index(m) for m in df["metal"]])
    h_map = {h: (-1 if h == ref_host else non_ref_hosts.index(h)) for h in hosts_in_fit}
    h_idx = np.array([h_map[h] for h in df["host"]])
    T_K = df["T_K"].to_numpy()
    log10_k_obs = df["log10_k"].to_numpy()
    sigma_obs_data = df["sigma_log10_k"].to_numpy()

    def predict(theta: np.ndarray) -> np.ndarray:
        muA = theta[iA0]
        betaA = theta[iAbeta]
        muE = theta[iE0]
        betaE = theta[iEbeta]
        epsA = theta[ieps_A]
        epsE = theta[ieps_E]
        b_h = theta[iB] if n_h > 0 else np.array([])
        # Per-metal log10 A and Ea (kJ/mol -> J/mol later)
        logA_m = muA + X_m @ betaA + epsA
        Ea_m = muE + X_m @ betaE + epsE   # kJ/mol
        # Per-observation prediction
        logA_obs = logA_m[m_idx]
        Ea_obs = Ea_m[m_idx] * 1e3   # -> J/mol
        host_contrib = np.where(h_idx >= 0, b_h[np.clip(h_idx, 0, max(n_h - 1, 0))], 0.0)
        host_contrib = np.where(h_idx >= 0, host_contrib, 0.0)
        pred = logA_obs + host_contrib - (Ea_obs / R_GAS) * (1.0 / T_K) / LN10
        return pred

    def log_posterior(theta: np.ndarray) -> float:
        muA = theta[iA0]
        sigmaA = theta[isigmaA]
        muE = theta[iE0]
        sigmaE = theta[isigmaE]
        epsA = theta[ieps_A]
        epsE = theta[ieps_E]
        sigmaO = theta[isigmaO]
        # Constraints
        if sigmaA <= 0 or sigmaE <= 0 or sigmaO <= 0:
            return -np.inf
        if sigmaA > 3.0 or sigmaE > 50.0 or sigmaO > 1.0:
            return -np.inf
        # Priors
        lp = 0.0
        lp += -0.5 * ((muA - 10.5) / 2.0) ** 2  # weak prior, log10 A ~ 10.5
        lp += -0.5 * ((muE - 30.0) / 30.0) ** 2 # weak prior, Ea ~ 30 kJ/mol
        lp += -0.5 * np.sum((theta[iAbeta] / 1.0) ** 2)
        lp += -0.5 * np.sum((theta[iEbeta] / 20.0) ** 2)
        lp += -0.5 * np.sum((epsA / sigmaA) ** 2) - len(epsA) * np.log(sigmaA)
        lp += -0.5 * np.sum((epsE / sigmaE) ** 2) - len(epsE) * np.log(sigmaE)
        # half-normal priors on the hyper-scales
        lp += -0.5 * (sigmaA / 1.5) ** 2
        lp += -0.5 * (sigmaE / 20.0) ** 2
        lp += -0.5 * (sigmaO / 0.3) ** 2
        if n_h > 0:
            b_h = theta[iB]
            lp += -0.5 * np.sum((b_h / 1.5) ** 2)
        # Likelihood
        pred = predict(theta)
        resid = (log10_k_obs - pred) / np.hypot(sigma_obs_data, sigmaO)
        lp += -0.5 * np.sum(resid ** 2)
        lp += -np.sum(np.log(np.hypot(sigma_obs_data, sigmaO)))
        return lp

    return log_posterior, predict, ndim, {
        "metals_in_fit": metals_in_fit,
        "ref_host": ref_host,
        "non_ref_hosts": non_ref_hosts,
        "feature_names": feature_names,
        "indices": {
            "muA": iA0, "betaA": iAbeta, "sigmaA": isigmaA,
            "muE": iE0, "betaE": iEbeta, "sigmaE": isigmaE,
            "epsA": ieps_A, "epsE": ieps_E,
            "B": iB, "sigmaO": isigmaO,
        },
    }


def initial_walkers(ndim: int, nwalkers: int, n_feat: int, n_m: int, n_h: int) -> np.ndarray:
    rng = np.random.default_rng(0)
    x0 = np.zeros((nwalkers, ndim))
    x0[:, 0] = 10.5 + 0.1 * rng.standard_normal(nwalkers)            # muA
    x0[:, 1:1+n_feat] = 0.1 * rng.standard_normal((nwalkers, n_feat)) # betaA
    x0[:, 1+n_feat] = 0.6 + 0.1 * rng.uniform(size=nwalkers)         # sigmaA
    x0[:, 2+n_feat] = 30.0 + 1.0 * rng.standard_normal(nwalkers)     # muE
    x0[:, 3+n_feat:3+2*n_feat] = 1.0 * rng.standard_normal((nwalkers, n_feat))
    x0[:, 3+2*n_feat] = 10.0 + 1.0 * rng.uniform(size=nwalkers)       # sigmaE
    eps_start = 4 + 2 * n_feat
    x0[:, eps_start:eps_start+n_m] = 0.05 * rng.standard_normal((nwalkers, n_m))
    x0[:, eps_start+n_m:eps_start+2*n_m] = 0.5 * rng.standard_normal((nwalkers, n_m))
    iB_start = eps_start + 2 * n_m
    if n_h > 0:
        x0[:, iB_start:iB_start+n_h] = 0.1 * rng.standard_normal((nwalkers, n_h))
    x0[:, iB_start+n_h] = 0.05 + 0.01 * rng.uniform(size=nwalkers)   # sigmaO
    return x0


def run_mcmc(df: pd.DataFrame, n_steps: int = 4000, n_walkers: int = 96,
              tag: str = "full"):
    metals = sorted(df["metal"].unique())
    hosts = sorted(df["host"].unique())
    X_m, feature_names = build_chem_features(metals)
    log_post, predict, ndim, meta = make_logpost(df, metals, hosts, X_m, feature_names)
    n_feat = X_m.shape[1]
    n_m = len(metals)
    n_h = len(meta["non_ref_hosts"])
    x0 = initial_walkers(ndim, n_walkers, n_feat, n_m, n_h)

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_post)
    print(f"[{tag}] Starting emcee: ndim={ndim}, n_walkers={n_walkers}, n_steps={n_steps}")
    sampler.run_mcmc(x0, n_steps, progress=False)
    chain = sampler.get_chain(discard=n_steps // 3, flat=True)
    log_prob = sampler.get_log_prob(discard=n_steps // 3, flat=True)
    keep = np.isfinite(log_prob)
    chain = chain[keep]
    log_prob = log_prob[keep]
    print(f"[{tag}] Posterior samples kept: {chain.shape[0]} / {n_walkers * n_steps}")
    return chain, log_prob, meta, predict, X_m


def predict_held_out(chain: np.ndarray, meta: dict, X_held: np.ndarray,
                      host: str, T_K: float) -> np.ndarray:
    """Posterior-predictive log10 k for a metal characterised by chemistry-features
    X_held (1, n_feat) at given host and T."""
    n_feat = X_held.shape[1]
    iA0 = meta["indices"]["muA"]
    iAbeta = meta["indices"]["betaA"]
    isigmaA = meta["indices"]["sigmaA"]
    iE0 = meta["indices"]["muE"]
    iEbeta = meta["indices"]["betaE"]
    isigmaE = meta["indices"]["sigmaE"]
    iB = meta["indices"]["B"]
    non_ref = meta["non_ref_hosts"]
    ref = meta["ref_host"]
    if host == ref:
        b_h = np.zeros(chain.shape[0])
    elif host in non_ref:
        j = non_ref.index(host)
        b_h = chain[:, iB][:, j]
    else:
        b_h = np.random.default_rng(0).standard_normal(chain.shape[0]) * 0.5
    muA = chain[:, iA0]
    muE = chain[:, iE0]
    betaA = chain[:, iAbeta]
    betaE = chain[:, iEbeta]
    sigA = chain[:, isigmaA]
    sigE = chain[:, isigmaE]
    rng = np.random.default_rng(0)
    eps_A_draw = sigA * rng.standard_normal(chain.shape[0])
    eps_E_draw = sigE * rng.standard_normal(chain.shape[0])
    logA = muA + X_held[0] @ betaA.T + eps_A_draw
    Ea = (muE + X_held[0] @ betaE.T + eps_E_draw) * 1e3
    pred = logA + b_h - (Ea / R_GAS) * (1.0 / T_K) / LN10
    return pred


def main():
    df = build_dataset()
    print(f"Assembled {len(df)} rows across {df['metal'].nunique()} metals "
          f"and {df['host'].nunique()} hosts.")
    print(df.groupby("metal").size())

    # Save chemistry features in their unstandardized form for reference
    metals_all = sorted(df["metal"].unique())
    X_all, feat_names = build_chem_features(metals_all)
    raw_feat = np.array([[*CHEM_FEATURES[m]] for m in metals_all])
    chem_df = pd.DataFrame(raw_feat, columns=["z", "r_pm", "chi", "n_d"], index=metals_all)
    chem_df.index.name = "metal"
    for j, name in enumerate(feat_names):
        chem_df[name] = X_all[:, j]
    chem_df.to_csv(OUT / "chemistry_features.csv")

    # ---- Full fit ----
    chain, logp, meta, predict, X_m = run_mcmc(df, n_steps=3000, tag="full")
    np.save(OUT / "full_posterior_chain.npy", chain)
    np.save(OUT / "full_posterior_logprob.npy", logp)

    # Save posterior summary
    indices = meta["indices"]
    summary_rows = []
    summary_rows.append({"param": "mu_log10_A", "mean": np.mean(chain[:, indices["muA"]]),
                         "p05": np.percentile(chain[:, indices["muA"]], 5),
                         "p95": np.percentile(chain[:, indices["muA"]], 95)})
    for j, name in enumerate(meta["feature_names"]):
        col = chain[:, indices["betaA"]][:, j]
        summary_rows.append({"param": f"beta_logA_{name}", "mean": col.mean(),
                             "p05": np.percentile(col, 5), "p95": np.percentile(col, 95)})
    summary_rows.append({"param": "sigma_logA", "mean": np.mean(chain[:, indices["sigmaA"]]),
                         "p05": np.percentile(chain[:, indices["sigmaA"]], 5),
                         "p95": np.percentile(chain[:, indices["sigmaA"]], 95)})
    summary_rows.append({"param": "mu_Ea_kJmol", "mean": np.mean(chain[:, indices["muE"]]),
                         "p05": np.percentile(chain[:, indices["muE"]], 5),
                         "p95": np.percentile(chain[:, indices["muE"]], 95)})
    for j, name in enumerate(meta["feature_names"]):
        col = chain[:, indices["betaE"]][:, j]
        summary_rows.append({"param": f"beta_Ea_{name}", "mean": col.mean(),
                             "p05": np.percentile(col, 5), "p95": np.percentile(col, 95)})
    summary_rows.append({"param": "sigma_Ea_kJmol", "mean": np.mean(chain[:, indices["sigmaE"]]),
                         "p05": np.percentile(chain[:, indices["sigmaE"]], 5),
                         "p95": np.percentile(chain[:, indices["sigmaE"]], 95)})
    summary_rows.append({"param": "sigma_obs_log10", "mean": np.mean(chain[:, indices["sigmaO"]]),
                         "p05": np.percentile(chain[:, indices["sigmaO"]], 5),
                         "p95": np.percentile(chain[:, indices["sigmaO"]], 95)})
    for j, h in enumerate(meta["non_ref_hosts"]):
        col = chain[:, indices["B"]][:, j]
        summary_rows.append({"param": f"b_host_{h}", "mean": col.mean(),
                             "p05": np.percentile(col, 5), "p95": np.percentile(col, 95)})
    pd.DataFrame(summary_rows).to_csv(OUT / "posterior_summary.csv", index=False)

    # ---- LOMO ----
    lomo_rows = []
    metals_for_lomo = [m for m in metals_all
                      if (df["metal"] == m).sum() >= 1 and m in CHEM_FEATURES]
    for held in metals_for_lomo:
        df_train = df[df["metal"] != held].reset_index(drop=True)
        if df_train["metal"].nunique() < 3:
            continue
        print(f"  LOMO holding out {held} ({(df['metal']==held).sum()} obs)")
        try:
            chain_h, _, meta_h, _, X_h = run_mcmc(df_train, n_steps=2000, n_walkers=72, tag=f"LOMO_{held}")
        except Exception as e:
            print(f"  ! LOMO {held} failed: {e}")
            continue
        # Compute X_held for the held-out metal under the training-set's
        # standardization
        train_metals = sorted(df_train["metal"].unique())
        raw_train = np.array([[*CHEM_FEATURES[m]] for m in train_metals])
        mu_train = raw_train.mean(axis=0)
        sd_train = raw_train.std(axis=0) + 1e-9
        raw_held = np.array([*CHEM_FEATURES[held]])
        feats_held = np.array([
            CHEM_FEATURES[held][0],
            CHEM_FEATURES[held][1],
            CHEM_FEATURES[held][2],
            np.log10(CHEM_FEATURES[held][0] / CHEM_FEATURES[held][1]),
        ])
        feats_train_raw = np.array([
            [CHEM_FEATURES[m][0], CHEM_FEATURES[m][1], CHEM_FEATURES[m][2],
             np.log10(CHEM_FEATURES[m][0] / CHEM_FEATURES[m][1])]
            for m in train_metals
        ])
        mu_t = feats_train_raw.mean(axis=0)
        sd_t = feats_train_raw.std(axis=0) + 1e-9
        X_held_std = ((feats_held - mu_t) / sd_t)[None, :]
        # Predict at every (host, T) in the held-out data
        df_held = df[df["metal"] == held].reset_index(drop=True)
        for _, r in df_held.iterrows():
            host, T = r["host"], r["T_K"]
            pred = predict_held_out(chain_h, meta_h, X_held_std, host, T)
            lomo_rows.append({
                "metal": held, "host": host, "T_K": T,
                "log10_k_obs": r["log10_k"], "sigma_log10_k_obs": r["sigma_log10_k"],
                "log10_k_pred_mean": np.mean(pred),
                "log10_k_pred_p05": np.percentile(pred, 5),
                "log10_k_pred_p95": np.percentile(pred, 95),
            })
    lomo_df = pd.DataFrame(lomo_rows)
    lomo_df.to_csv(OUT / "lomo_results.csv", index=False)
    print(f"Saved {len(lomo_df)} LOMO predictions.")

    # ---------- Figure: 4 panels ----------
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), constrained_layout=True)
    palette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999", "#882255"]

    # (a) raw data scatter: log10 k vs 1000/T, color = metal, marker = host
    ax = axes[0, 0]
    host_marker = {"LiCl-KCl": "o", "NaCl": "s", "KCl": "^", "KBr": "D",
                    "NaBr": "P", "KI": "X", "H2O": "*"}
    metals_sorted = sorted(df["metal"].unique())
    for i, m in enumerate(metals_sorted):
        sub = df[df["metal"] == m]
        for host_, grp in sub.groupby("host"):
            ax.errorbar(1000.0 / grp["T_K"], grp["log10_k"], yerr=grp["sigma_log10_k"],
                         fmt=host_marker.get(host_, "o"),
                         color=palette[i % len(palette)], alpha=0.85,
                         markersize=4.5, lw=0.7,
                         label=f"{m} / {host_}" if host_ == "LiCl-KCl" or m == "Cf3+" else None)
    ax.set_xlabel(r"$1000/T$ [K$^{-1}$]")
    ax.set_ylabel(r"$\log_{10} k$ [M$^{-1}$ s$^{-1}$]")
    ax.set_title("(a) Multi-metal data")
    ax.legend(fontsize=6, ncol=2, loc="lower left", frameon=False)
    _despine(ax)

    # (b) Chemistry-feature loadings on log10 A and Ea
    ax = axes[0, 1]
    betaA_samp = chain[:, indices["betaA"]]
    betaE_samp = chain[:, indices["betaE"]]
    n_feat = betaA_samp.shape[1]
    y = np.arange(n_feat)
    means_A = betaA_samp.mean(axis=0)
    err_A = np.array([
        means_A - np.percentile(betaA_samp, 5, axis=0),
        np.percentile(betaA_samp, 95, axis=0) - means_A,
    ])
    means_E = betaE_samp.mean(axis=0)
    err_E = np.array([
        means_E - np.percentile(betaE_samp, 5, axis=0),
        np.percentile(betaE_samp, 95, axis=0) - means_E,
    ])
    ax.errorbar(means_A, y - 0.15, xerr=err_A, fmt="o", color="#0072B2",
                markersize=5, label=r"$\log_{10} A$")
    ax.errorbar(means_E, y + 0.15, xerr=err_E, fmt="s", color="#D55E00",
                markersize=5, label=r"$E_a$ [kJ mol$^{-1}$]")
    ax.axvline(0, color="#bbbbbb", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(meta["feature_names"], fontsize=7)
    ax.set_xlabel("posterior coefficient (mean, 90% CI)")
    ax.set_title("(b) Feature loadings")
    ax.legend(fontsize=7, frameon=False)
    _despine(ax)

    # (c) LOMO scatter: predicted vs observed log10 k
    ax = axes[1, 0]
    if len(lomo_df) > 0:
        cmap = {m: palette[i % len(palette)] for i, m in enumerate(metals_sorted)}
        for m, grp in lomo_df.groupby("metal"):
            yerr_pred = np.array([
                grp["log10_k_pred_mean"] - grp["log10_k_pred_p05"],
                grp["log10_k_pred_p95"] - grp["log10_k_pred_mean"],
            ])
            ax.errorbar(grp["log10_k_obs"], grp["log10_k_pred_mean"],
                         yerr=yerr_pred, xerr=grp["sigma_log10_k_obs"],
                         fmt="o", color=cmap.get(m, "k"), alpha=0.85,
                         markersize=5, lw=0.7, label=m)
        lo = min(lomo_df["log10_k_obs"].min(), lomo_df["log10_k_pred_p05"].min())
        hi = max(lomo_df["log10_k_obs"].max(), lomo_df["log10_k_pred_p95"].max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.7)
    ax.set_xlabel(r"observed $\log_{10} k$")
    ax.set_ylabel(r"LOMO predicted $\log_{10} k$")
    ax.set_title("(c) LOMO predictions")
    ax.legend(fontsize=6, ncol=2, frameon=False)
    _despine(ax)

    # (d) Host effects b_h posterior
    ax = axes[1, 1]
    if len(meta["non_ref_hosts"]) > 0:
        b_samp = chain[:, indices["B"]]
        means = b_samp.mean(axis=0)
        err = np.array([
            means - np.percentile(b_samp, 5, axis=0),
            np.percentile(b_samp, 95, axis=0) - means,
        ])
        y = np.arange(len(meta["non_ref_hosts"]))
        ax.errorbar(means, y, xerr=err, fmt="s", color="#009E73", markersize=6)
        ax.axvline(0, color="#bbbbbb", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(meta["non_ref_hosts"], fontsize=7)
        ax.set_xlabel(r"$b_h$ on $\log_{10} k$ (ref " + str(meta["ref_host"]) + ")")
        ax.set_title("(d) Host effects")
        _despine(ax)

    fig.savefig(FIG_DIR / "fig_meta_hier.pdf", dpi=600, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig_meta_hier.png", dpi=300, bbox_inches="tight")
    print(f"Wrote figure {FIG_DIR / 'fig_meta_hier.pdf'}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()
