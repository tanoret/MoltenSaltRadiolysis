# Model adequacy and calibration for molten-salt radiolysis networks

**Status:** methodological review — selects candidate methods for (1) building the most likely
radiolysis network for chloride/fluoride salts and (2) calibrating closure coefficients to data.

**Scope.** This review surveys the state-of-the-art in Bayesian and deterministic uncertainty
quantification (UQ) for chemical reaction networks, with emphasis on what is appropriate for
*stiff* ODE systems where data is sparse but accurate (the radiation-chemistry regime). It does
NOT attempt to cover every UQ technique — only those that (a) have been demonstrated on kinetic
networks, (b) handle the model-selection-plus-calibration coupling explicitly, and (c) have mature,
maintained software.

---

## 1. Problem statement

Let `θ ∈ R^p` be the vector of closure coefficients (Arrhenius `A_i`, `Ea_i`, G-values `G_j`, Henry
constants, mass-transfer coefficients). Let `M_k` index a finite or open-ended set of candidate
*networks* (subsets of plausible reactions). Each network defines a stiff ODE system
`dC/dt = f(C; θ, M_k)` on species concentrations `C ∈ R^{n_k}`. Observations `y` are noisy
functionals of the trajectory (transient absorbances at multiple `T`, [M] pairs; per-T 2nd-order
rates with σ; gas detection thresholds).

We want:

1. **Network selection / inference**: posterior `p(M_k | y)` over candidate networks, or
   posterior inclusion probabilities for individual reactions when the candidate set is
   combinatorially large.
2. **Closure-coefficient calibration**: posterior `p(θ | y, M_k)` for each retained network,
   with **propagated uncertainty** that includes both data noise and model inadequacy.
3. **Model adequacy diagnostics**: where (in `T`, [M], time-scale) does the calibrated model
   systematically disagree with data? Is the disagreement informative — i.e., does it point
   at a missing physical mechanism rather than just measurement noise?

The methodological literature for (1)–(3) is mature, primarily through 20 years of combustion-kinetics
UQ. Radiation chemistry has lagged this development by ~15 years; **the central scientific
contribution of this work can be the careful application of combustion-grade UQ to the radiolysis
network problem.**

---

## 2. Methodological landscape

### 2.1 Foundational paradigms

There are **three coexisting paradigms** for kinetic-model UQ, each with different epistemology:

| Paradigm | Stance | Output | Pros | Cons |
|---------|--------|--------|------|------|
| **Bayesian** (Marzouk, Najm, Berry, Le Maître) | Probabilistic inference with priors. | Posterior distributions over θ and M. | Naturally propagates uncertainty; supports model averaging and BOED. | Likelihood must be specified; computational cost; identifiability issues with imperfect models. |
| **Bound-to-Bound Data Collaboration (B2BDC)** (Frenklach, Russi, Packard) | Optimization on a feasible set defined by data + uncertainty bounds. | Feasible region in parameter space; consistency check on data. | Distribution-free; detects dataset inconsistencies; deterministic. | No notion of likelihood; less direct connection to predictive uncertainty. |
| **Sparse identification** (Brunton, Kutz, Hoffmann) | Discover the network as sparsest fit to time-series. | Network topology + coefficients (point estimates). | Discovers structure with minimal prior; works when network is unknown. | Less natural for sparse, multi-modality data; UQ usually retro-fitted. |

For the radiolysis problem, where the network *candidates* are known (the chemistry is constrained
by mass and charge balance and known precedents in combustion / aqueous radiolysis), the
**Bayesian paradigm is the strongest fit**: prior literature gives us informative priors,
identified networks are few, and we want explicit posterior predictive intervals to compare
against the Phillips 2022 NULL benchmark and the published Iwamatsu/Conrad rate constants.
B2BDC is complementary — useful for cross-paper consistency checks before any Bayesian inference.

### 2.2 The Galagali–Marzouk framework (closest precedent)

The single most relevant prior work for goal (1) is:

> **Galagali, N.; Marzouk, Y. M.** "Bayesian inference of chemical kinetic models from
> proposed reactions." *Chemical Engineering Science* **2015**, 123, 170–190.
> [doi:10.1016/j.ces.2014.11.030](https://doi.org/10.1016/j.ces.2014.11.030) (MIT preprint:
> http://web.mit.edu/ymarz/www/papers/GalagaliM_CES2015.pdf).

This paper does *exactly* the joint problem Mauricio is posing. The framework:

1. Enumerate a set of candidate elementary reactions whose union forms an over-complete network.
2. Place a **spike-and-slab** (or point-mass mixture) prior on each reaction's rate constant: with
   probability `π` the reaction is included with `log A ~ N(μ_A, σ_A^2)`; with probability `1-π`
   it is excluded (`log A = -∞`).
3. Run **reversible-jump MCMC** (RJMCMC) to jointly sample (a) which reactions are in the model
   and (b) their rate constants conditional on inclusion.
4. Apply an **adaptive independence MCMC** scheme with online EM updating of the proposal — this
   addresses the notorious tuning difficulty of vanilla RJMCMC.
5. Read off: posterior inclusion probability per reaction; posterior over network topologies
   (the most-probable mechanism); posterior over parameters with model-structure uncertainty
   marginalized out.

A 2019 follow-up extends this to large networks by exploiting reaction-network topology to design
better between-model proposals:

> **Galagali, N.; Marzouk, Y. M.** "Exploiting network topology for large-scale inference of
> nonlinear reaction models." *J. R. Soc. Interface* **2019**, 16(152), 20180766.
> [doi:10.1098/rsif.2018.0766](https://doi.org/10.1098/rsif.2018.0766).

**This is the methodological backbone we should adopt** for the chloride/fluoride-network selection
problem. The candidate-reaction list comes naturally from our [database.yaml](../msr_radiolysis/data/database.yaml) plus
literature candidates (Cl3⁻ vs Cl2 branching of Cl2•⁻ disproportionation; Cr3+ + Cl2•⁻ to Cr4+
or to Cr2+; F2 recombination; iodine extension reactions).

### 2.3 Closely related Bayesian SciML developments (2023–2025)

> **Li, Q.; Chen, H.; Koenig, B. C.; Deng, S.** "Bayesian chemical reaction neural network for
> autonomous kinetic uncertainty quantification." *Phys. Chem. Chem. Phys.* **2023**, 25,
> 3707–3717. [doi:10.1039/D2CP05083H](https://doi.org/10.1039/D2CP05083H).

Embeds the law of mass action and Arrhenius equation as architectural constraints in a neural
network, then runs MCMC OR variational inference (VI) over the network weights. Treats the
network topology and rate parameters jointly. Extension in:

> **Li, Q.; Chen, H.; Deng, S.** "Uncertainty quantified discovery of chemical reaction systems
> via Bayesian scientific machine learning." *Frontiers in Systems Biology* **2024**, 4, 1338518.

The CRNN line is interesting but has a key limitation for our problem: it discovers reactions
*from data alone*, throwing away the very rich prior information we have on candidate reactions.
For radiolysis we know the species; we don't want to discover spurious cross-couplings between
e.g. Cr²⁺ and Zn⁺. **The Galagali–Marzouk framework with informative spike-and-slab priors is
a better fit than CRNN.**

> **Yu, X.; Zhang, S.; Tao, F.; Liu, G.** "Uncertainty quantification and reduction for combustion
> kinetic modeling: A case study of NH3/H2 models." *Fuel* **2025**, 396, 135180.
> [Search link](https://www.sciencedirect.com/science/article/pii/S0016236125015352).

Recent example of the full UQ-on-combustion-mechanism workflow including reduction.

> **Chen, H.** "Bayesian inference and experimental design of combustion kinetic models." MIT MS
> thesis, 2023. [https://dspace.mit.edu/handle/1721.1/153682](https://dspace.mit.edu/handle/1721.1/153682).

Demonstrates an end-to-end pipeline: MCMC calibration (400k steps) + BIC-based selection +
optimal experimental design suggesting the next shock-tube experiment to run.

### 2.4 The Frenklach school (B2BDC + PrIMe) — complementary

> **Frenklach, M.** "Transforming data into knowledge — Process Informatics for combustion
> chemistry." *Proc. Combust. Inst.* **2007**, 31, 125–140.
> **Hegde, A.; Li, W.; Oreluk, J.; Packard, A.; Frenklach, M.** "Consistency analysis for
> massively inconsistent datasets in B2BDC." *Combust. Flame* **2018**, 196, 509–517.
> **Frenklach, M.; Packard, A.; Garcia-Donato, G.; Paulo, R.; Sacks, J.** "Comparison of
> statistical and deterministic frameworks of uncertainty quantification." *SIAM/ASA JUQ* **2016**.
> **Oreluk, J.; Packard, A.; Frenklach, M.** "Representing model discrepancy in B2BDC."
> *SIAM/ASA JUQ* **2021**, 9(2), 462–489.
> [doi:10.1137/19M1270185](https://doi.org/10.1137/19M1270185).

B2BDC defines a *feasible set* in parameter space that is consistent with all literature data
within reported uncertainty bounds. Key uses for us:

- **Consistency check before inference**: are Pikaev 1982 (k = 1.7e9 in NaCl @ 850 °C) and
  Iwamatsu 2022 (k = 16.1e10 in LiCl-KCl @ 600 °C extrapolated) **simultaneously consistent**
  with the same Arrhenius form? Iwamatsu argues no (Pikaev had inadequate time resolution).
  B2BDC formalizes this check.
- **Outlier detection**: which experimental constraints are jointly infeasible, before we attempt
  to calibrate a model to them?

This is *not* a substitute for Bayesian inference (it doesn't yield a posterior) but should be
run as a preprocessing step.

### 2.5 Model adequacy / discrepancy — the Kennedy–O'Hagan framework

> **Kennedy, M. C.; O'Hagan, A.** "Bayesian calibration of computer models." *J. R. Stat. Soc. B*
> **2001**, 63(3), 425–464. [doi:10.1111/1467-9868.00294](https://doi.org/10.1111/1467-9868.00294).

The canonical formulation:
```
  y(x) = ζ(x, θ*) + δ(x) + ε
       ────────  ──────  ────
       model     model    obs
       output    discrep. noise
```
with `δ(x)` a Gaussian process (GP) capturing unmodeled physics, `θ*` the "true" closure
coefficients, and `ε ~ N(0, σ²)`. Estimating `θ` and `δ` *jointly* often disagrees with
estimating `θ` against a *perfect* model.

But the K–O framework has a serious caveat:

> **Brynjarsdóttir, J.; O'Hagan, A.** "Learning about physical parameters: the importance of
> model discrepancy." *Inverse Problems* **2014**, 30, 114007.
> [doi:10.1088/0266-5611/30/11/114007](https://doi.org/10.1088/0266-5611/30/11/114007).

They show that a flexible (unconstrained) GP `δ` introduces **identifiability problems**: the
posterior on `θ*` can be biased even with infinite data, and the bias does not vanish.
Three mitigations:

1. **Physics-informed priors on δ**: bound the magnitude, smoothness, or asymptotic behavior of
   the discrepancy using known physics.
2. **Staged calibration**: estimate `θ` from data regimes where the model is known to be adequate,
   then estimate `δ` against remaining data.
3. **Constrained discrepancy**: replace the GP with a parametric form whose flexibility is
   bounded (e.g., linear in `T` and log `[M]` with bounded coefficients).

For radiolysis, the natural physics-informed prior is: `|δ|` must be small relative to the
observed absorbance amplitudes; `δ → 0` as dose → 0; `δ` should be smooth in `T` over the
experimental range. This is *implementable* in the K–O framework via informative GP hyperpriors.

### 2.6 Identifiability — the prerequisite for inference

Before fitting, ask: *can* the data identify the parameters at all?

> **Raue, A.; Kreutz, C.; Maiwald, T.; Bachmann, J.; Schilling, M.; Klingmüller, U.; Timmer, J.**
> "Structural and practical identifiability analysis of partially observed dynamical models by
> exploiting the profile likelihood." *Bioinformatics* **2009**, 25(15), 1923–1929.
> [doi:10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).

Profile likelihood is the most pragmatic tool: for each parameter `θ_i`, find the maximum
likelihood as a function of `θ_i` while re-optimizing the others. The shape of the profile
identifies:
- **Structural non-identifiability**: profile is flat → parameter is fundamentally undetermined
  by any amount of data (e.g., when only the ratio `A·exp(-Ea/RT)` is measurable at one T).
- **Practical non-identifiability**: profile is curved but does not reach the threshold → data
  amount/quality insufficient; informative prior or new experiment is needed.

> **Gutenkunst, R. N.; Waterfall, J. J.; Casey, F. P.; Brown, K. S.; Myers, C. R.; Sethna, J. P.**
> "Universally sloppy parameter sensitivities in systems biology models." *PLoS Comput. Biol.*
> **2007**, 3(10), e189. [doi:10.1371/journal.pcbi.0030189](https://doi.org/10.1371/journal.pcbi.0030189).

> **Apgar, J. F.; Witmer, D. K.; White, F. M.; Tidor, B.** "Sloppy models, parameter uncertainty,
> and the role of experimental design." *Mol. BioSyst.* **2010**, 6, 1890–1900.

Sethna et al. show that biological ODE models are universally *sloppy*: a few "stiff" parameter
combinations matter; many "sloppy" combinations don't affect any observable. Apgar et al. and the
2014 "Sloppy Models Can Be Identifiable" paper (Janzén et al., arXiv:1403.1417) show that
sloppiness *is reducible* given the right experiment.

For the radiolysis problem, run profile likelihood on the chloride kernel *before* committing to a
calibration strategy. Likely findings:
- `A_i × exp(-Ea_i / R T_ref)` at a single `T_ref` is well-identified, but `A_i` and `Ea_i`
  separately are not → need multi-T data (which Iwamatsu provides; Pikaev 1982 doesn't).
- G-values are coupled with initial pulse-dose calibration → need to fix one (or include the
  dose as a nuisance parameter with strong prior).

### 2.7 Model comparison (selecting among candidate networks)

For a *finite* set of candidate networks `M_1, ..., M_K`:

> **Vehtari, A.; Gelman, A.; Gabry, J.** "Practical Bayesian model evaluation using leave-one-out
> cross-validation and WAIC." *Stat. Comput.* **2017**, 27(5), 1413–1432.
> [doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4).
>
> **Vehtari, A.; Simpson, D.; Gelman, A.; Yao, Y.; Gabry, J.** "Pareto smoothed importance
> sampling." *J. Mach. Learn. Res.* **2024**, 25(72), 1–58. [arXiv:1507.02646](https://arxiv.org/abs/1507.02646).
>
> Living FAQ: [Vehtari Cross-Validation FAQ](https://users.aalto.fi/~ave/CV-FAQ.html).

**PSIS-LOO** (Pareto-smoothed importance sampling leave-one-out cross-validation) estimates
out-of-sample predictive accuracy. It is the canonical Bayesian model-comparison tool today.
Recommended workflow:

1. Fit each candidate network via MCMC.
2. Compute PSIS-LOO `elpd_loo` for each.
3. Compare via `loo_compare`: difference in `elpd` ± standard error.
4. **Caveat**: if `k̂` (Pareto shape diagnostic) exceeds 0.7 for some observations, the LOO
   estimate is unreliable; refit with K-fold CV instead.

Reject the Bayes Factor / marginal likelihood approach unless using bridge sampling or thermodynamic
integration — naive harmonic-mean estimators are notoriously unstable for ODE models.

### 2.8 Optimal experimental design (BOED)

Once we have a calibrated posterior, choose the next experiment to maximize expected information
gain:

> **Lindley, D. V.** "On a measure of the information provided by an experiment." *Ann. Math.
> Statist.* **1956**, 27, 986–1005.
>
> **Foster, A.; Ivanova, D. R.; Malik, I.; Rainforth, T.** "DEEP-Adaptive Design: Amortizing
> sequential Bayesian experimental design." *PMLR* **2021**, 139, 3384–3395.
> [arXiv:2103.02438](https://arxiv.org/abs/2103.02438).

For radiolysis, EIG could rank candidate next experiments:
- Pulse Iwamatsu-style measurement of a new metal in LiCl-KCl;
- Akiyama-style F⁻ + e_s⁻ in FLiNaK at 500 °C;
- High-dose chronic Phillips-style for U³⁺/U⁴⁺ in NaCl-UCl₃;
- Variation of pulse dose to reduce G-value vs k_2 correlation.

This is the right tool to motivate which paper Mauricio should prioritize digitizing or which
experiment INL/MSEE should fund next.

### 2.9 Surrogate-assisted inference for expensive forward solves

If each ODE solve is expensive (it is not — minutes at most for a stiff Cr+Zn chloride system),
the polynomial-chaos (PC) and Gaussian-process (GP) surrogate route accelerates Bayesian
inference dramatically. Key references:

> **Marzouk, Y. M.; Najm, H. N.** "Dimensionality reduction and polynomial chaos acceleration
> of Bayesian inference in inverse problems." *J. Comput. Phys.* **2009**, 228(6), 1862–1902.
> [doi:10.1016/j.jcp.2008.11.024](https://doi.org/10.1016/j.jcp.2008.11.024).
>
> **Wang, J.; Zabaras, N.** "Hierarchical Bayesian models for inverse problems in heat conduction."
> *Inverse Problems* **2005**.

For our problem, PC may not be needed initially — direct NUTS through `solve_ivp` with `BDF`
should be tractable. Hold this in reserve for the model-averaging step if it dominates compute.

### 2.10 Practical software ecosystem

Three viable stacks, ranked for the radiolysis problem:

| Rank | Stack | Strengths | Weaknesses |
|------|-------|-----------|------------|
| 1 | **Julia + Turing.jl + DifferentialEquations.jl + DiffEqBayes.jl + ArviZ.jl** | Best stiff-ODE solvers (Rodas, KenCarp); native ergonomic integration with NUTS/HMC; SciMLSensitivity for adjoint gradients; growing combustion-UQ user base. | Re-implement the model from Python; team must learn Julia. |
| 2 | **Python + PyMC + sunode (CVODE wrapper) + ArviZ** | Stays in current Python ecosystem; `sunode` provides backward-pass for NUTS; ArviZ is the canonical diagnostics package. | sunode is fragile for very stiff systems; less mature than DifferentialEquations.jl. |
| 3 | **Python + emcee + SciPy.integrate.solve_ivp** | Trivial to wire; works with the existing code as-is; affine-invariant ensemble sidesteps gradient computation. | Sampler doesn't scale beyond ~50 parameters; no PSIS-LOO support natively; no native model-comparison machinery. |

**Recommendation**: write the calibration code in **Python with sunode + PyMC** first (Tier 1
below), then *if* gradient/scaling pain becomes severe, port the ODE+inference to Julia.
The PyMC ecosystem now has mature PSIS-LOO via ArviZ and supports HMC/NUTS over ODEs.

For SBI (simulation-based inference) — only relevant if the likelihood becomes intractable:

> **Boelts, J. et al.** "sbi reloaded: a toolkit for simulation-based inference workflows."
> 2024. [arXiv:2411.17337](https://arxiv.org/abs/2411.17337). [https://sbi-dev.github.io/sbi/](https://sbi-dev.github.io/sbi/).

Hold in reserve; not needed for the explicit-likelihood radiolysis case.

---

## 3. Recommended methodology for the radiolysis problem

Based on the survey above, here is the concrete recommendation. **Three tiers**, each with explicit
deliverables. Tiers are *sequential* — completing each unlocks the next.

### Tier 1 — Identifiability and consistency (preprocessing, **no expensive MCMC**)

**Goal**: establish what's learnable from the data we have *before* committing to any inference
framework.

1. **Run profile-likelihood analysis** on the current chloride kernel (or the Cr-extended version)
   using only the digitized Iwamatsu 2026 transients. For each parameter `θ_i`, profile the
   negative-log-likelihood and check for structural vs practical non-identifiability.
   Deliverable: a table of identifiable / sloppy parameter combinations, e.g., "log(A) and Ea are
   only jointly identifiable via the linear combination log(A) − Ea/(R·T_avg)".
   Library: `pyABC`, `petab`, or hand-rolled `scipy.optimize` on a fixed grid.

2. **Run a B2BDC-style consistency check** across all reported rate constants (Pikaev 1982,
   Hagiwara 1987, Iwamatsu 2022, Iwamatsu 2026). Question: is there a single Arrhenius
   `(A, Ea)` that fits all reported values within their σ?
   Deliverable: pass/fail per reaction, with a list of inconsistent (paper, T) data points to flag.

3. **Sensitivity analysis** via Morris elementary effects or Sobol indices on the full Cr+Zn
   chloride model to rank parameter importance.
   Library: `SALib` (Python).

**Outcome of Tier 1**: a defensible model formulation (which parameters to fit jointly, which to
fix from literature, which to declare structurally indeterminate) and a vetted set of
consistent literature data. **This is where the Bayesian "enormous rigor" begins** — without this
step, the posterior is hard to interpret.

### Tier 2 — Calibration with informative priors and explicit discrepancy

**Goal**: posterior `p(θ | data, M_k)` for *each* candidate network with proper discrepancy modeling.

1. **Likelihood**: Gaussian on `log(absorbance)` for transient traces (multiplicative noise → log
   scale), and Gaussian on `log(k_obs)` for reported per-T rates. Use the σ values already in
   `arrhenius_parameters.csv` and `reported_rate_constants.csv`.

2. **Priors**:
   - For each Arrhenius `(A, Ea)` with published values: Gaussian on `(log A, Ea)` centered on
     the literature mean with the literature σ.
   - For unmeasured G-values: weakly informative log-normal `LogNormal(log(0.3), 1.0)` (covers
     `0.04` to `2.2` at 1σ), within the range of plausible G-values for halide melts.
   - For pulse-dose nuisance: `LogNormal(log(20 Gy), 0.3)` (covers 14–28 Gy per pulse, matching
     Iwamatsu's reported 15–30 Gy range).
   - For Henry constants and `k_L a`: bounded uniform on physical ranges.

3. **Inference**: PyMC with `sunode` (CVODE backend for the stiff ODE), NUTS sampler. Burn-in
   ~1000, posterior ~2000–5000 samples. Diagnose with R-hat, ESS, divergences.

4. **Model discrepancy**: include a parametric `δ(T, log [M])` term modeled as a low-order
   polynomial (NOT a full GP) with bounded coefficients, per the Brynjarsdóttir–O'Hagan
   identifiability guidance. Explicitly, for each observable type (e_s⁻ decay rate, Cl2•⁻
   self-rate, metal-redox rate):
   ```
   log k_pred(T, [M]) = log k_model(T, [M]; θ) + a + b·(1/T − 1/T_ref) + c·log [M]
   ```
   with `(a, b, c) ~ N(0, σ_disc)` and `σ_disc` set so `|δ| ≲ 20%` at the boundary of the data.

5. **Diagnostics**: posterior predictive checks (Gelman, Meng, Stern 1996); calibration plots
   (Cook, Gelman, Rubin 2006); test-statistic-based PPC (the Phillips 2022 NULL is the natural
   test: does the posterior predict `[Cl₂] < 1000 ppm` after 31 MGy at 600 °C?).

**Library**: PyMC ≥ 5.0 (5.x has improved ODE backends and ArviZ integration). ArviZ for
diagnostics + PSIS-LOO.

### Tier 3 — Network selection / inference

**Goal**: posterior over network topologies.

Two options depending on the size of the candidate set.

**Option A (small set, ≤ ~10 networks)**: fit each network in Tier 2, compute PSIS-LOO `elpd`,
compare via `loo_compare`. Report best model + competitive alternatives with `Δelpd ± SE` ratios.
This is the easy path and is enough for the *initial* paper.

**Option B (open-ended set, e.g., evaluating every plausible Cl-radical reaction)**: implement the
Galagali–Marzouk RJMCMC framework. Each reaction has a Bernoulli inclusion variable `γ_i ∈ {0, 1}`
with prior `P(γ_i = 1) = π_prior_i` (informative from literature precedence). Sample with
RJMCMC. Output: posterior inclusion probability per reaction, posterior over network topologies.

Recommend Option A for the first paper, with Option B held as a follow-on.

### Tier 4 (optional but powerful) — Optimal experimental design

**Goal**: rank candidate next experiments by expected information gain (EIG).

Use the Tier 2/3 posterior as the prior, define a parameterized space of candidate experiments
(e.g., `(salt, T, [M_dopant], pulse_dose, observable_λ)` for a pulse-rad experiment;
`(salt, T, dose_rate, duration)` for a steady-state experiment). Maximize EIG via Foster–Ivanova
variational BOED.

This is the right tool to motivate Mauricio's INL Akiyama-1994-paywall request or to argue for a
new FLiNaK experiment, in terms of *quantified* posterior shrinkage on the fluoride G-values.

---

## 4. Concrete candidate methods, with selection rationale

The user asked for **candidate methods** — here are the five I propose, ranked by how
defensible they are for the specific radiolysis problem:

| # | Method | Use for | Why | Cost |
|---|--------|---------|-----|------|
| 1 | **PyMC + sunode + NUTS** with informative literature priors | Tier 2 closure-coefficient calibration | Mature; ergonomic for ODE; canonical diagnostics. The default workhorse. | Minutes-to-hours per network. |
| 2 | **PSIS-LOO via ArviZ** | Tier 3 model selection among ≤ 10 candidate networks | Canonical Bayesian model comparison; gives uncertainty on the comparison; no extra fits needed. | Trivial post-MCMC. |
| 3 | **Profile likelihood (hand-rolled or PEtab)** | Tier 1 identifiability | The only tool that *honestly* answers "what does the data identify?" before inference. Mandatory under "enormous rigor". | Hours of CPU on a coarse parameter grid. |
| 4 | **Galagali–Marzouk RJMCMC** | Tier 3 network selection over an open-ended candidate set | The single most rigorous published approach to the joint network+parameter inference problem. | Days of CPU; significant implementation effort. |
| 5 | **Foster–Ivanova variational BOED** | Tier 4 experiment ranking | Justifies which paper to digitize next or which experiment to fund, in expected-information-gain units. | Modest, given a working Tier 2 posterior. |

**My recommendation**: start with Tier 1 (#3) and Tier 2 (#1 + parametric discrepancy from
Brynjarsdóttir–O'Hagan), then do Tier 3 with PSIS-LOO (#2) on 4–6 candidate chloride networks
(varying the Cl3⁻ branching, the Cr3+ + Cl2•⁻ pathway, and the Zn+ + Cl2•⁻ inclusion). The
Galagali–Marzouk machinery (#4) is the right tool *if* the PSIS-LOO comparison reveals that
the candidate set is too restricted — but I would not implement it speculatively.

---

## 5. Open questions / risks that need to be resolved before commitment

These are the points where I do not yet have enough information to make a recommendation with
confidence — flagged for the user's review:

1. **Computational budget**. The recommendations above assume ~CPU-hours per posterior sample
   for a Cr+Zn chloride network with NUTS. If the user wants posteriors over the *full*
   chloride+fluoride+metals network simultaneously, NUTS may not scale; we should benchmark
   the forward solve cost first.

2. **Likelihood specification for transients**. The Iwamatsu absorbance traces are noisy;
   choosing a per-point Gaussian likelihood is the default, but a *process-noise* formulation
   (where the noise is correlated in time, e.g., AR(1)) may be more honest. This is a defensible-
   methodology question; both can be reported.

3. **Cross-paper systematic biases**. Pikaev 1982 vs Iwamatsu 2022 for Zn²⁺ + e_s⁻ disagree by an
   order of magnitude. The Bayesian framework will either (a) downweight the older data via the
   discrepancy term, (b) declare them inconsistent (B2BDC consistency check), or (c) absorb the
   difference into a "facility effect" hierarchical model. **Which is correct is a substantive
   scientific question**, not a methodological one. The framework should let us *report* the
   ambiguity rather than hide it.

4. **Model discrepancy parameterization**. The B–O'H mitigation of physics-informed priors on
   `δ` is right in spirit but requires a concrete bound. For radiolysis, what is the *a priori*
   maximum allowed model–data discrepancy? 10%? 30%? This sets `σ_disc` in Tier 2.4 and needs
   to be defended from physical reasoning.

5. **Cross-salt transfer**. The fluoride data is sparse. Do we treat the chloride and fluoride
   kernels as *separate* models, or do we use a hierarchical Bayesian model that pools
   information across salts via shared hyperparameters (e.g., a common Eyring entropy of
   activation distribution)? This is a research design question: pooled models are more powerful
   but require defensible mechanistic justification.

---

## 6. Next steps (proposed)

1. **You decide** which of the candidate methods to commit to. I recommend the path: Tier 1
   profile-likelihood + B2BDC-style consistency check → Tier 2 PyMC/sunode NUTS with
   parametric discrepancy → Tier 3 PSIS-LOO over 4–6 candidate chloride networks → Tier 4
   variational BOED for the next-experiment question.

2. Once the path is agreed, I draft a **methods document** at the level of an article methods
   section: equations, software stack, prior justifications, likelihood specification, MCMC
   diagnostic thresholds, PSIS-LOO acceptance criteria, discrepancy parametrization. This is
   the document the manuscript will eventually quote.

3. Then we implement Tier 1. The current model is already in `msr_radiolysis/`; profile
   likelihood needs only the model plus a 1-D optimizer, so this is hours of work and
   immediately tells us whether the data identifies the parameters at all.

4. Then we proceed to Tier 2 + Tier 3 in PyMC. The validation manifests in [validation/](.)
   already hold the data structures the inference will use.

I am ready to draft the methods document the moment you confirm the path. If you would prefer a
different path — say, Julia / Turing.jl from the start, or a pure SINDy/CRNN discovery approach
rather than Bayesian inference on a known candidate set — say so before we commit, because the
implementation effort downstream is substantial.

---

## Appendix A — Key references in citation form

(Selected; full list in the paper draft.)

- Apgar, Witmer, White, Tidor. "Sloppy models, parameter uncertainty, and the role of
  experimental design." *Mol. BioSyst.* **2010**, 6, 1890–1900.
- Brynjarsdóttir, O'Hagan. "Learning about physical parameters: the importance of model
  discrepancy." *Inverse Problems* **2014**, 30, 114007. [doi:10.1088/0266-5611/30/11/114007](https://doi.org/10.1088/0266-5611/30/11/114007).
- Chen, H. "Bayesian inference and experimental design of combustion kinetic models." MIT MS
  thesis, **2023**. [https://dspace.mit.edu/handle/1721.1/153682](https://dspace.mit.edu/handle/1721.1/153682).
- Foster, A.; Ivanova, D. R.; Malik, I.; Rainforth, T. "Deep adaptive design: Amortizing
  sequential Bayesian experimental design." *PMLR* **2021**, 139, 3384–3395.
- Frenklach, M. "Transforming data into knowledge — Process Informatics for combustion chemistry."
  *Proc. Combust. Inst.* **2007**, 31, 125–140.
- Galagali, N.; Marzouk, Y. M. "Bayesian inference of chemical kinetic models from proposed
  reactions." *Chem. Eng. Sci.* **2015**, 123, 170–190. [doi:10.1016/j.ces.2014.11.030](https://doi.org/10.1016/j.ces.2014.11.030).
- Galagali, N.; Marzouk, Y. M. "Exploiting network topology for large-scale inference of nonlinear
  reaction models." *J. R. Soc. Interface* **2019**, 16(152), 20180766.
  [doi:10.1098/rsif.2018.0766](https://doi.org/10.1098/rsif.2018.0766).
- Gutenkunst, R. N.; Waterfall, J. J.; Casey, F. P.; Brown, K. S.; Myers, C. R.; Sethna, J. P.
  "Universally sloppy parameter sensitivities in systems biology models." *PLoS Comput. Biol.*
  **2007**, 3(10), e189. [doi:10.1371/journal.pcbi.0030189](https://doi.org/10.1371/journal.pcbi.0030189).
- Hegde, A.; Li, W.; Oreluk, J.; Packard, A.; Frenklach, M. "Consistency analysis for massively
  inconsistent datasets in B2BDC." *Combust. Flame* **2018**, 196, 509–517.
- Ji, W.; Qiu, W.; Shi, Z.; Pan, S.; Deng, S. "Stiff-PINN: Physics-informed neural network for
  stiff chemical kinetics." *J. Phys. Chem. A* **2021**, 125(36), 8098–8106.
- Kennedy, M. C.; O'Hagan, A. "Bayesian calibration of computer models." *J. R. Stat. Soc. B*
  **2001**, 63(3), 425–464.
- Li, Q.; Chen, H.; Koenig, B. C.; Deng, S. "Bayesian chemical reaction neural network for
  autonomous kinetic uncertainty quantification." *Phys. Chem. Chem. Phys.* **2023**, 25,
  3707–3717. [doi:10.1039/D2CP05083H](https://doi.org/10.1039/D2CP05083H).
- Marzouk, Y. M.; Najm, H. N. "Dimensionality reduction and polynomial chaos acceleration of
  Bayesian inference in inverse problems." *J. Comput. Phys.* **2009**, 228(6), 1862–1902.
  [doi:10.1016/j.jcp.2008.11.024](https://doi.org/10.1016/j.jcp.2008.11.024).
- Oreluk, J.; Packard, A.; Frenklach, M. "Representing model discrepancy in B2BDC."
  *SIAM/ASA JUQ* **2021**, 9(2), 462–489. [doi:10.1137/19M1270185](https://doi.org/10.1137/19M1270185).
- Raue, A.; Kreutz, C.; Maiwald, T.; Bachmann, J.; Schilling, M.; Klingmüller, U.; Timmer, J.
  "Structural and practical identifiability analysis of partially observed dynamical models by
  exploiting the profile likelihood." *Bioinformatics* **2009**, 25(15), 1923–1929.
  [doi:10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).
- Vehtari, A.; Gelman, A.; Gabry, J. "Practical Bayesian model evaluation using leave-one-out
  cross-validation and WAIC." *Stat. Comput.* **2017**, 27(5), 1413–1432.
  [doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4).
- Vehtari, A. "Cross-validation FAQ" (living document).
  [https://users.aalto.fi/~ave/CV-FAQ.html](https://users.aalto.fi/~ave/CV-FAQ.html).
- Yu, X.; Zhang, S.; Tao, F.; Liu, G. "Uncertainty quantification and reduction for combustion
  kinetic modeling: A case study of NH3/H2 models." *Fuel* **2025**, 396, 135180.
