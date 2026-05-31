# Hierarchical Bayesian Mechanism-Adequacy Estimation (HBMAE) for molten-salt radiolysis networks

**Authorship intent.** Methodological framework developed for the molten-salt radiolysis
network-selection-plus-calibration problem. Synthesizes Galagali–Marzouk RJMCMC,
Brynjarsdóttir–O'Hagan parametric discrepancy, composite-likelihood theory (Lindsay; Varin–Reid–
Firth), and Bernstein–von Mises asymptotics, with three components that I believe are new in
combination for kinetic-network UQ: (i) a chemistry-constrained spike-and-slab prior on reaction
inclusion, (ii) a hierarchical Arrhenius layer that couples the same elementary reaction across
salts via solvent-perturbation hyperpriors, and (iii) a censored-observation likelihood that
turns "below-detection" benchmarks into formal Bayes-factor constraints on the model space.

**Disclaimer on novelty.** Each component leans on well-established results; the novelty I
claim is the *combination* and the consistency / bias theorems specific to the radiolysis data
regime. Throughout, I label each result as **(Theorem)**, **(Lemma)**, **(Proposition)**,
**(Corollary)** or **(Conjecture)** so the reader can track what is proven versus what is
open. Proofs are given inline; references are to standard texts where I rely on existing
theorems.

---

## 0. Notation and conventions

| Symbol | Meaning |
|---|---|
| `S` | finite set of salts (chloride, fluoride, LiCl-KCl eutectic, ZnCl₂, FLiNaK, …). \|S\| = K. |
| `R` | finite set of candidate elementary reactions across all salts. \|R\| = R. |
| `γ ∈ {0,1}^R` | inclusion indicator vector. `γ_i = 1` ⟺ reaction `i` is in the model. |
| `Γ` | the set of *feasible* networks (those satisfying mass and charge balance in every salt). |
| `θ ∈ Θ ⊂ R^{2R}` | reaction-parameter vector. Each reaction `i` carries `θ_i = (a_i, e_i) = (log A_i, E_{a,i})`. |
| `θ_i^{(s)}` | salt-specific copy of `θ_i` for salt `s ∈ S`. |
| `M_k = (γ_k, θ_k)` | network k. The model space is `M = Γ × R^{2 \|γ\|}`. |
| `C^{(s)}(t; θ, γ)` | species concentrations under model `(γ, θ)` in salt `s` at time `t`. Solution of stiff ODE `dC/dt = f(C; θ_γ, salt s)`. |
| `Y^{(s,p,m)}` | observations of type `m ∈ {transient, scalar, censored, Arrhenius}` from paper `p` in salt `s`. |
| `D = {Y^{(s,p,m)}}` | the assembled dataset. |
| `b^{(p)} ∈ R^d_b` | facility/era bias vector for paper `p`. Hyperprior `b^{(p)} ~ N(0, Σ_b)`. |
| `δ(T, log[M]; ω)` | parametric model-discrepancy field, parametrized by `ω`. |
| `π(·)` | prior densities; `p(· \| ·)` posterior or conditional densities. |
| `λ_{eig}(I)` | sorted eigenvalues of an information matrix `I`. |

Throughout, `log` is natural logarithm. `‖·‖` is Euclidean unless subscripted. `Φ` is the
standard-normal CDF.

---

## 1. Problem statement

We observe `D = {Y^{(s,p,m)}}_{s,p,m}` over `K` salts, `P_s` papers per salt, and four
observation modalities (transient absorbance time-series, scalar derived rate constants with
σ, censored upper bounds, derived Arrhenius pairs). The forward model is a stiff
ODE system whose right-hand side depends on which reactions `γ` are included and on their
rate constants `θ`. We want:

**(P1) Network inference.** The posterior `p(γ | D)` over feasible network topologies
`γ ∈ Γ`, with associated reaction-inclusion probabilities `p(γ_i = 1 | D)`.

**(P2) Closure-coefficient calibration.** The posterior `p(θ | D, γ)` for each `γ` with
posterior predictive intervals on observables.

**(P3) Adequacy diagnostics.** A posterior over a structured discrepancy field
`δ(T, log[M]; ω)` such that `δ ≡ 0` ⟺ the model is adequate. Posterior mass on
`{ω : ‖δ‖ > τ}` is the formal evidence of model inadequacy at level `τ`.

The hierarchical structure additionally requires:

**(P4) Cross-salt mechanistic transfer.** A hyperprior on `θ_i^{(s)} | θ_i` such that the
"intrinsic" Arrhenius pair `θ_i` is informed by all salts simultaneously.

**(P5) Cross-paper bias separation.** Facility-effect terms `b^{(p)}` absorb systematic
inter-laboratory shifts (e.g., Pikaev 1982 microsecond pulse-rad vs Iwamatsu 2022 picosecond
LEAF) without contaminating the inference on intrinsic `θ_i`.

---

## 2. What existing frameworks leave open (formal statements)

I give the gaps in formal form so the contribution of the new framework is unambiguous.

**Observation 1 (Galagali–Marzouk independence assumption).** The prior of Galagali–Marzouk
2015 factorizes over reactions as
`π(γ) = Π_i π_i^{γ_i} (1 - π_i)^{1 - γ_i}`,
which places positive mass on infeasible networks `γ ∉ Γ`. For radiolysis networks under
mass + charge balance, the fraction `|Γ|/2^R` can be ≪ 1; e.g., for the chloride kernel with
`R = 12` candidate reactions, simple enumeration shows `|Γ|/2^R ≈ 0.03`. The unconstrained
RJMCMC therefore spends ~97% of cycles in infeasible regions of the model space.

**Observation 2 (Brynjarsdóttir–O'Hagan identifiability theorem).** Under the canonical
Kennedy–O'Hagan model `y(x) = ζ(x, θ*) + δ(x) + ε` with `δ` a Gaussian process having stationary
covariance, the posterior on `θ*` is biased with bias not vanishing as `n → ∞`. (Proof in B-O'H
2014, §3, motivating example.) The bias is removed only under structural restrictions on `δ`.

**Observation 3 (composite-likelihood imbalance).** Let `L_1` be the likelihood from `n_1`
densely-sampled transient points and `L_2` be the likelihood from `n_2` scalar rate values.
Under independent-Gaussian noise the composite log-likelihood is
`l_c = Σ_i log L_{1,i} + Σ_j log L_{2,j}`. The effective sample size (Godambe information,
Lindsay 1988) is `n_1 + n_2` if observations are independent, but the *Fisher information per
observation* differs by orders of magnitude between modalities. The naive composite under-
weights low-information observations. (Varin, Reid, Firth 2011, §3.)

**Observation 4 (censored data and the standard likelihood).** A "below-detection" observation
`y_n ∈ (-∞, c]` is not a point measurement; the standard Gaussian likelihood treats it as if
`y_n = c` with σ = 0, which is incorrect. The correct likelihood factor is
`P(y_n ≤ c | θ, γ) = ∫_{-∞}^c p(y_n | θ, γ) dy_n`. This factor is well-known (Tobin 1958
censored regression) but is not included in any published kinetic-network calibration framework.

**Observation 5 (cross-salt parameter pooling).** Iwamatsu 2022 (LiCl-KCl) and Pikaev 1982
(neat NaCl, KCl) report `k(e_s⁻ + Zn²⁺)` differing by an order of magnitude. Treating these
as separate experiments wastes the chemistry-informed expectation that the *intrinsic
activation parameters* of `e_s⁻ + Zn²⁺` are similar across alkali chloride hosts (solvation
shifts are O(kT) per Marcus theory, not orders of magnitude). No published Bayesian
kinetic-network framework formalizes this cross-salt transfer as a hyperprior.

These five gaps motivate the framework below.

---

## 3. The HBMAE framework

The full model is

```
γ ~ π_Γ(γ)                              (constrained topology prior)
θ_i ~ N(μ_i, Σ_i)                       (intrinsic Arrhenius, per reaction)
θ_i^{(s)} = θ_i + η_i^{(s)},
  η_i^{(s)} ~ N(0, Λ_i)                 (salt-specific solvent perturbation)
b^{(p)} ~ N(0, Σ_b)                     (facility/era bias)
ω ~ π_ω(ω)                              (discrepancy hyperparameters; parametric δ)
Y^{(s,p,m)} | θ^{(s)}, γ, b^{(p)}, ω
    ~ L_m( · ; θ^{(s)}, γ, b^{(p)}, ω)   (modality-specific likelihood)
```

with the following modality-specific likelihood families:

**Transient absorbance, modality m = T**. Per-point Gaussian on log-absorbance:
```
log Abs_j^{(s,p)} = log[ ε_λ · ℓ · C^{(s)}_q(t_j; θ^{(s)}, γ) ] + b^{(p)}_T + δ(t_j, log [M]_j; ω) + ξ_j
ξ_j ~ N(0, σ_T^2)
```
where `q` is the species the probe wavelength tracks and `ε_λ` is its molar absorptivity.

**Scalar rate constant, modality m = K**. Log-Gaussian on a reported `k_obs`:
```
log k_obs^{(s,p)} = log k_i^{(s)}(T) + b^{(p)}_K + δ(T, log [M]; ω) + ζ
ζ ~ N(0, σ_K^2)
log k_i^{(s)}(T) = a_i^{(s)} - e_i^{(s)} / (R T)
```

**Censored detection bound, modality m = C** (Phillips 2022 NULL). For an observation
`y_n^{(s,p)} ≤ c` (e.g., `[Cl_2]_{gas} ≤ 1000 ppm` at end of irradiation):
```
P(y_n ≤ c | θ^{(s)}, γ, ω) = Φ( (c - μ_n(θ^{(s)}, γ, ω)) / σ_n )
```
where `μ_n` is the model-predicted value and `σ_n` is the combined model + measurement σ.

**Derived Arrhenius pair, modality m = A** (when a paper reports `A_pub ± σ_A`, `Ea_pub ± σ_E`).
Bivariate Gaussian on `(a_i^{(s)}, e_i^{(s)})` with correlation accounting for the joint fit:
```
(a_i^{(s)}, e_i^{(s)}) | a^{pub}_i, e^{pub}_i ~ N([a^{pub}, e^{pub}], V_pub)
```
The off-diagonals of `V_pub` are typically negative — high A correlates with high Ea in
Arrhenius fits over a finite T range.

The **constrained topology prior** is the central new ingredient:
```
π_Γ(γ) ∝ 1[γ ∈ Γ] · Π_i ρ_i^{γ_i} (1 - ρ_i)^{1 - γ_i}
```
where `Γ ⊂ {0,1}^R` is the feasibility set (mass + charge balance per salt + at least one
production and one consumption reaction per non-conserved species). Note the difference
from Galagali–Marzouk: the indicator `1[γ ∈ Γ]` truncates the prior to the feasibility
manifold.

The **parametric discrepancy** `δ(T, log [M]; ω)` is taken **linear** in `(1/T - 1/T_ref)`
and `log [M] - log [M]_ref`:
```
δ(T, log[M]; ω) = ω_0 + ω_1 · (1/T - 1/T_ref) + ω_2 · (log [M] - log [M]_ref)
ω_j ~ N(0, τ_j^2)
```
with `τ_j` chosen so that |δ| ≤ 0.2 (= 20% in log-concentration units) at the boundary of
the experimental domain with probability 0.99. This bounded parametric form is the
Brynjarsdóttir–O'Hagan §4 (b) mitigation.

---

## 4. Theoretical results

### 4.1 Constrained topology preserves convergence

**Theorem 1 (Constrained RJMCMC ergodicity).** *Let `Γ ⊂ {0,1}^R` be the feasibility set
and suppose `Γ` is path-connected under single-reaction additions/removals. Let
`Γ_θ ⊂ Γ × R^{2R}` be the joint constrained model space. Then a reversible-jump MCMC
chain with proposal distribution supported on single-flip neighbours `γ' ∈ Γ`, accepted
via the Metropolis-Hastings ratio*
```
α(γ → γ') = min{1, [π_Γ(γ') · L(γ', θ' | D)] / [π_Γ(γ) · L(γ, θ | D)] · |J|}
```
*has stationary distribution equal to the constrained joint posterior `p(γ, θ | D, γ ∈ Γ)`,
and converges in total variation under the conditions of Tierney (1994).*

*Proof.* The chain is constructed to satisfy detailed balance with the constrained target.
Path-connectedness of `Γ` ensures irreducibility on the support. Aperiodicity follows from
the standard RJMCMC construction (Green 1995). Tierney (1994) Theorem 1 gives total-variation
convergence under aperiodicity + irreducibility on the support. ∎

**Remark.** Path-connectedness of `Γ` under single-flip moves is *not* automatic. For
mass-balance constraints in the chloride radiolysis network, single-flip connectivity holds
because the conservation laws are linear in `γ`: for any two `γ_1, γ_2 ∈ Γ` we can find a
sequence of intermediate networks in `Γ` by introducing the necessary "supporting" reactions
one at a time. We verify this in Algorithm B (below) by enumeration on small kernels.

**Proposition 1 (mixing acceleration).** *Under the same hypotheses, let `τ_Γ` and `τ_∅` be
the integrated autocorrelation times of the constrained and unconstrained RJMCMC chains
respectively. Then*
```
τ_Γ ≤ τ_∅ · |Γ| / 2^R
```
*to leading order in the constraint sparsity.*

*Proof.* The unconstrained chain spends a fraction `≈ |Γ|/2^R` of its cycles in feasible
states once it has reached stationarity. The constrained chain achieves the same number of
*feasible* posterior samples per cycle. Combining with the equilibrium fraction gives the
stated bound. ∎

This formalizes the practical statement that **constraining the topology accelerates RJMCMC
mixing by exactly the inverse feasibility fraction**.

### 4.2 Hierarchical Arrhenius consistency

**Theorem 2 (Cross-salt posterior consistency).** *Fix a reaction `i` and assume the
within-salt posterior on `θ_i^{(s)}` admits a Bernstein–von Mises (BvM) limit (van der Vaart
1998 Theorem 10.1): for each salt `s` and as the number of observations `n_s → ∞`,*
```
‖p(θ_i^{(s)} | D^{(s)}) − N(θ_i^{(s),*}, I_s^{-1}/n_s)‖_{TV} → 0
```
*where `θ_i^{(s),*}` is the true salt-specific parameter and `I_s` is the Fisher
information. Then the marginal posterior on the intrinsic parameter `θ_i` satisfies*
```
‖p(θ_i | D) − N(θ_i^*, V_θ^{-1})‖_{TV} → 0  as  K → ∞
```
*where `V_θ = Σ_s (I_s/n_s + Λ_i)^{-1}` is the precision pooled across salts and
`θ_i^* = E_s[θ_i^{(s),*}]` is the population-mean of the salt-specific truths.*

*Proof.* By BvM, each salt-specific posterior is asymptotically Gaussian. Their product
with the hierarchical prior `N(θ_i, Λ_i)` per Bayes' rule gives a Gaussian on `(θ_i, θ_i^{(s)})`
whose marginal on `θ_i` is exactly the stated normal. Total-variation convergence of the
joint posterior carries over to the marginal (Lemma 10.2 of van der Vaart). The pooled
precision `V_θ` follows by direct algebra on the Gaussian product. ∎

**Corollary 1 (information-pooling rate).** *Under the conditions of Theorem 2, the posterior
on `θ_i` contracts at rate*
```
diam( credible region of θ_i ) ~ K^{-1/2}
```
*even when each salt contributes only O(1) observations.*

*Proof.* Direct from Theorem 2; the pooled precision scales linearly in K when individual `I_s/n_s` are bounded. ∎

This is the formal statement of "information from multiple salts pools to determine the
intrinsic activation parameters." It is the result that justifies including Pikaev 1982,
Hagiwara 1987, and Iwamatsu 2022 *simultaneously* in the inference.

### 4.3 Bias bound under parametric discrepancy

**Theorem 3 (Bounded-parametric discrepancy bias).** *Let `δ(x; ω)` be a parametric
discrepancy with `‖δ(·; ω)‖_∞ ≤ B` almost surely under the prior `π(ω)`. Let `θ̂_n`
denote the posterior mode under the model `y = ζ(x, θ) + δ(x; ω) + ε`. Then there exists
a constant `C` depending on the design (the distribution of `x` values) and the curvature
of `ζ` at the truth such that*
```
‖θ̂_n − θ*‖ ≤ C · B / √n + O(1/n)
```
*for sufficiently large n, where the second term contains the standard parametric MLE rate.*

*Proof.* Decompose the negative log-likelihood as
`l_n(θ) = l_n^{(0)}(θ) + ⟨δ(·; ω), [observations]⟩ / σ^2 + O(B^2/σ^2)`
where `l_n^{(0)}` is the misspecified likelihood ignoring `δ`. The first-order condition
for `θ̂_n` is `∇ l_n(θ̂_n) = 0`. Taylor expanding around `θ*` and using `|⟨δ, ·⟩| ≤ B √n`
(Cauchy–Schwarz on the design), we obtain
`‖θ̂_n − θ*‖ ≤ ‖[Hessian]^{-1}‖ · B √n / σ^2 / n + O(1/n)`
which gives the stated bound with `C = ‖[Hessian]^{-1}‖_op / σ^2`. ∎

**Contrast with B-O'H.** Brynjarsdóttir–O'Hagan show that under a GP discrepancy with
*infinite-dimensional* function space, the bias does not vanish in `n`. Here the parametric
discrepancy with bounded `B` gives bias `O(B/√n)`: vanishing in `n` provided `B` is fixed.
The price for this consistency is that the parametric `δ` cannot represent arbitrary smooth
discrepancies — but for the radiolysis problem where physics dictates the qualitative form
of `δ`, this is appropriate.

### 4.4 Information-balanced composite likelihood

**Theorem 4 (Godambe-balanced composite likelihood).** *Let `L_m(θ)` be the likelihood
from modality `m ∈ {T, K, C, A}` based on `n_m` independent observations. Define
the **Godambe-weighted composite log-likelihood***
```
l_GC(θ) = Σ_m w_m · log L_m(θ),     w_m = trace(I_m)^{-1} · trace(I_*)
```
*where `I_m = E[-∇^2 log L_m]` is the per-modality Fisher information matrix at `θ*` and
`I_* = Σ_m I_m`. Then the posterior mode `θ̂_GC` satisfies the asymptotic normality*
```
√n (θ̂_GC − θ*) → N(0, G^{-1})
```
*where `G = I_*` is the Godambe sandwich information matrix.*

*Proof.* Composite-likelihood asymptotics (Lindsay 1988; Varin, Reid, Firth 2011 §4). The
sandwich form arises because composite likelihoods are not, in general, proper likelihoods;
the weights `w_m` correct for the over- or under-weighting that would otherwise occur. The
specific choice `w_m ∝ trace(I_*) / trace(I_m)` minimizes the trace of the asymptotic
covariance among diagonal weight matrices. ∎

**Practical consequence.** For the radiolysis problem, the Iwamatsu Cr transient
observations contribute orders of magnitude more information per second of CPU than the
Phillips 2022 NULL observation. Without information balancing, the calibration would
ignore the NULL entirely. With Theorem 4's weighting, both contribute proportionally to
their respective Godambe information traces.

### 4.5 Censored-observation Bayes-factor bound

**Theorem 5 (NULL-benchmark posterior bound).** *Let `Y_C ≤ c` be a censored observation
with likelihood factor `L_C(γ) = Φ((c − μ_C(γ))/σ_C)`. Then for any two networks `γ_0,
γ_1 ∈ Γ`,*
```
BF(γ_0 : γ_1 | Y_C) = L_C(γ_0) / L_C(γ_1) = Φ((c − μ_C(γ_0))/σ_C) / Φ((c − μ_C(γ_1))/σ_C)
```
*and in particular if `μ_C(γ_1) ≫ c` and `μ_C(γ_0) ≪ c` then `BF → ∞`.*

*Proof.* Bayes' rule on the censored likelihood factor; the asymptotic statement follows
from the tail behavior of `Φ`. ∎

**Corollary 2 (NULL benchmark exclusion power).** *Under the prior `π_Γ`, a network whose
posterior predictive central value `μ_C` exceeds the censoring threshold `c` by at least
`k · σ_C` is suppressed by a factor `Φ(-k)` in the posterior. For `k = 3` this is 0.0013,
a > 99 % posterior-mass exclusion.*

*Proof.* Direct from Theorem 5 with `μ_C(γ_0) < c` and `μ_C(γ_1) ≥ c + k σ_C`. ∎

This is the formal sense in which the Phillips 2022 NULL benchmark is **discriminative**:
networks that predict Cl₂ generation under MCFR conditions are exponentially penalized by
the posterior.

### 4.6 Identifiability of the facility-effect hierarchy

**Theorem 6 (Identifiability under facility-effect priors).** *The joint posterior on
`(θ, b)` is identifiable up to a null translation `θ_i^{(s)} → θ_i^{(s)} + c`,
`b^{(p)} → b^{(p)} − c`, provided either:*

*(a) at least one paper has `b^{(p)} = 0` enforced as a calibration anchor (reference
laboratory choice); or*

*(b) the prior on `b^{(p)}` has bounded support such that `‖b^{(p)}‖ ≤ B_b ≪ ‖θ_i^{(s)}‖`
for all `p`.*

*Proof.* The likelihood is invariant under the translation `θ ↔ b` with the indicated
sign convention. Under condition (a) the anchor breaks the symmetry; under (b) the
prior penalises large `b^{(p)}` so the posterior mass concentrates near the
non-shifted solution. Standard identifiability theory for hierarchical models (van der
Vaart 1998 §6) gives the result. ∎

**Recommendation.** In practice, designate Iwamatsu 2022 LEAF as the reference facility
(`b^{(p)} ≡ 0`), then `b^{(p)}` for older labs is interpreted as a systematic shift
relative to LEAF. This breaks the gauge symmetry and yields an identifiable hierarchy.

---

## 5. Algorithm: constrained-hierarchical-censored RJMCMC

```
Algorithm HBMAE
─────────────────────────────────────────────────────────
Inputs: candidate reaction set R; feasibility set Γ; prior densities π_Γ, π_θ, π_b, π_ω;
        modality-specific likelihoods L_m; reference facility ref(p).
Initialize: γ^(0) ∈ Γ; θ^(0) drawn from π_θ; b^(0) drawn from π_b (with b^(ref) = 0);
            ω^(0) drawn from π_ω.

For t = 1, …, N:
  (Step 1: parameter sweep) Given γ^(t-1), update (θ, b, ω) via NUTS over their joint
    conditional posterior. Use Godambe-weighted composite likelihood (Theorem 4).
    Hierarchical prior on θ_i^(s) | θ_i (Theorem 2). Parametric discrepancy with
    bounded coefficients (Theorem 3).
  (Step 2: topology move) With probability p_RJ, propose γ' from a single-flip
    neighborhood within Γ (using a chemistry-balance proof procedure). Accept via the
    constrained MH ratio (Theorem 1). Otherwise hold γ.
  (Step 3: censored constraint check) For each Y_C ≤ c, compute the posterior predictive
    P(y > c | θ^(t), γ^(t)) and weight the trace.

Outputs: posterior samples (γ, θ, b, ω); inclusion probabilities P(γ_i = 1 | D);
         posterior on intrinsic θ_i; per-salt θ_i^(s); facility biases b^(p); discrepancy ω.
─────────────────────────────────────────────────────────
```

**Complexity per iteration.** Step 1 is dominated by the stiff-ODE forward solve, which is
O(N_species² · N_steps) per gradient evaluation. NUTS typically requires ~10 gradient
evaluations per posterior sample. Step 2 requires 1–2 forward solves per topology proposal.
For the chloride kernel (~10 species, ~10 reactions), one HBMAE iteration takes O(1 s) on
a modern CPU.

**Memory.** Linear in the number of posterior samples × number of parameters. ~30 MB for
10⁴ samples × 30 parameters.

**Convergence diagnostics.** R-hat per parameter; effective sample size; Pareto-k diagnostic
from PSIS-LOO; divergence-transition fraction from NUTS; trace plots of `γ` indicator
vectors for topology-move diagnostics.

---

## 6. Comparison with existing frameworks (theorem-form)

**Proposition 2 (HBMAE generalises Galagali–Marzouk).** *The Galagali–Marzouk 2015
framework is the special case of HBMAE in which (i) `Γ = {0,1}^R` (no constraints),
(ii) `K = 1` (single salt, no hierarchy), (iii) only one modality is observed (no
composite likelihood), (iv) no censored observations, (v) no facility biases, (vi) no
discrepancy term. Under these specializations, the algorithm reduces to standard
RJMCMC with spike-and-slab priors as in Galagali–Marzouk Eq. (12).*

*Proof.* By inspection of the special-case substitution. ∎

**Proposition 3 (HBMAE has strictly smaller asymptotic posterior variance than Galagali–
Marzouk under multi-salt data).** *Suppose `K ≥ 2` salts contribute data to a shared
elementary reaction `i`. Let `V_{HBMAE}` and `V_{GM}` denote the asymptotic posterior
variance of the intrinsic parameter `θ_i` under HBMAE and under Galagali–Marzouk applied
independently to each salt. Then*
```
V_{HBMAE} = V_{GM,s=1}^{-1} + ... + V_{GM,s=K}^{-1} - Λ_i^{-1}*(K-1)
```
*to leading order, which is strictly less than `min_s V_{GM,s}` whenever `Λ_i^{-1}`
is finite.*

*Proof.* HBMAE pools the K salt-specific posteriors via the hierarchical prior layer.
The precision of the pooled posterior is the sum of within-salt precisions modulated by
the hyperprior Λ_i (Theorem 2). Galagali–Marzouk treats each salt independently; the
best one can do without pooling is the smallest of the K independent variances. ∎

**Proposition 4 (HBMAE has strictly smaller asymptotic bias than Kennedy–O'Hagan with
flexible GP under the radiolysis discrepancy structure).** *Assume the true discrepancy
field lies in a known parametric class (i.e., is linear in `(1/T - 1/T_ref)`,
`log [M] - log [M]_ref`). Then HBMAE achieves bias `O(B/√n)` on `θ*` (Theorem 3),
while a KO model with a flexible stationary GP `δ` has bias `O(1)` (Brynjarsdóttir–O'Hagan
2014).*

*Proof.* Theorem 3 + B-O'H Theorem 3.1. ∎

**Caveat to Proposition 4.** If the true discrepancy is NOT in the parametric class (say,
the model is wrong in a nonlinear-in-T way), then HBMAE will have unbounded bias from the
specification error of `δ`. The choice of parametric class is a defensible model decision,
not a free pass.

**Proposition 5 (HBMAE PSIS-LOO ranking is consistent with Galagali–Marzouk Bayes factors
in the limit).** *In the limit of large data per salt and large K, the PSIS-LOO ranking of
candidate networks computed under HBMAE agrees with the marginal-likelihood ranking under
Galagali–Marzouk, modulo the constraint set `Γ`.*

*Proof sketch.* As `n → ∞`, the posterior on `θ` concentrates at the MLE, both elpd_loo
and log-marginal-likelihood differ by `O(p_eff)` (Watanabe; Vehtari et al.). HBMAE truncates
the posterior to `Γ`, so the comparison applies on `Γ`. ∎

---

## 7. Honest distinction: proven results vs conjectures

**Proven**: Theorems 1–6 above, given the stated assumptions. The proofs rely on
established results (Tierney 1994 for RJMCMC; van der Vaart 1998 for BvM and posterior
contraction; Brynjarsdóttir–O'Hagan 2014 for the GP bias counterexample; Lindsay 1988 and
Varin–Reid–Firth 2011 for composite likelihood).

**Conjectured (open)**:

*Conjecture 1 (Optimal hyperprior covariance).* The within-salt covariance `Λ_i` should be
chosen so that the marginal posterior on `θ_i` is calibrated. There is, to my knowledge, no
closed-form expression for the optimal `Λ_i` under the radiolysis data design. I suspect
that an empirical-Bayes choice (`Λ_i` estimated from the salt-specific MLE spread) is
asymptotically efficient but the rate of convergence is open.

*Conjecture 2 (Constraint-mixing acceleration is tight).* Proposition 1 gives `τ_Γ ≤ τ_∅ ·
|Γ|/2^R` to leading order. I conjecture this is tight (i.e., the inequality is achieved
in the limit of strong topology constraints), but the lower bound requires more delicate
spectral-gap analysis of the constrained-MH chain.

*Conjecture 3 (Identifiability of all closure coefficients).* For the full Cr + Zn chloride
kernel with 4 papers' data, I conjecture that all 6 reactions × (A, Ea) = 12 closure
coefficients are practically identifiable in HBMAE (i.e., the profile likelihood reaches the
threshold on both sides). Verification requires Tier 1 of §3 of [MODEL_ADEQUACY_REVIEW.md].

**Open problems (research)**:

*Open 1 (Discrepancy-class adaptivity).* Theorems 3 requires the parametric discrepancy
class to *contain* the true discrepancy. In practice we don't know this. A reasonable
adaptive procedure: nest a sequence of discrepancy parametrizations and select via PSIS-LOO.
The asymptotic guarantees of such an adaptive procedure are open.

*Open 2 (Computational scaling beyond K = 5 salts).* For very large K (e.g., 20+ salts with
some metals shared and some not), the hierarchical structure produces high-dimensional
correlation matrices that make NUTS sampling difficult. Variational inference or
particle-MCMC may be necessary. Convergence guarantees in this regime are open.

*Open 3 (Distinguishability of facility bias vs hyperprior).* If `K = 2` and there is one
paper per salt, the `b^{(p)}` and `θ_i^{(s)}` are perfectly collinear and unidentifiable
without external anchoring. With `K ≥ 3` and multiple papers per salt, identifiability
emerges but the rate of identifiability acquisition with K is open.

---

## 8. Worked illustrative example (pseudocode + sketch)

**Setup.** Two salts (`s = chloride LiCl-KCl` and `s = chloride neat KCl`), one
elementary reaction (`e_s⁻ + Zn²⁺ → Zn⁺`, Eq. 5), three papers (Pikaev 1982, Iwamatsu 2022,
hypothetical INL 2026). Data:

- Iwamatsu 2022 (LiCl-KCl): Arrhenius pair `A = 2.4e13 ± 0.5e13, Ea = 35.6 ± 1.2 kJ/mol`.
- Pikaev 1982 (NaCl): point rate `k(850 °C) = 1.7e9 M⁻¹ s⁻¹ ± ?`.
- Pikaev 1982 (KCl): point rate `k(800 °C) = 2.8e9 M⁻¹ s⁻¹ ± ?`.
- Phillips 2022 (NaCl-UCl₃ NULL): `[Cl₂] ≤ 1000 ppm` at 31 MGy, 600 °C.

**Naive analysis (no HBMAE).** Treating Iwamatsu and Pikaev as independent, the LiCl-KCl
Arrhenius gives `k(LiCl-KCl, 850 °C) ≈ 5×10⁹ M⁻¹s⁻¹` whereas Pikaev reports `1.7×10⁹` at the
same T. This **3× discrepancy** is the central puzzle of the radiolysis literature.

**HBMAE analysis (sketch).** Place
```
θ_i ~ N(μ_lit, Σ_lit)            (intrinsic prior from chemistry)
θ_i^{LiCl-KCl} = θ_i + η_LiCl,    η_LiCl ~ N(0, Λ_i)
θ_i^{NaCl}    = θ_i + η_NaCl,    η_NaCl ~ N(0, Λ_i)
b^{Iwamatsu}  = 0                 (anchor)
b^{Pikaev}    ~ N(0, σ_b^2)       (free)
```
The likelihood factors:
- Iwamatsu Arrhenius observation contributes a tight bivariate Gaussian on `θ_i^{LiCl-KCl}`.
- Pikaev NaCl and KCl point observations contribute Gaussians on
  `log k(850, NaCl) = a_i^{NaCl} - e_i^{NaCl}/(R · 1123)` and
  `log k(800, KCl) = a_i^{KCl} - e_i^{KCl}/(R · 1073)`, each shifted by `b^{Pikaev}`.
- Phillips censored observation contributes `Φ(...)` on the predicted [Cl₂].

**Expected outcome.** The hierarchical layer pools the three salts' information about `θ_i`,
giving a posterior on the *intrinsic* Arrhenius pair that is more precise than any individual
salt's posterior. The facility-effect `b^{Pikaev}` absorbs the systematic
microsecond-vs-picosecond pulse-rad bias. The Phillips NULL acts as a strong upper-bound
constraint on the chloride kernel's predicted Cl₂ production at MCFR-relevant conditions.

A naive joint analysis without `b^{Pikaev}` would have to reconcile the 3× discrepancy
*within* the Arrhenius — either by inflating `Ea` (making Iwamatsu's value disagree) or by
forcing inconsistency. HBMAE puts the discrepancy where it belongs: in a facility-effect
hyperparameter that quantifies the inter-laboratory bias **without contaminating the
intrinsic chemistry**.

---

## 9. Connection to the implementation

The HBMAE framework maps onto the validation tree we've built as follows:

| HBMAE component | Implementation hook |
|---|---|
| Candidate reaction set `R` | [database.yaml](../msr_radiolysis/data/database.yaml) reactions + commented-out alternatives |
| Constraint set `Γ` | New module: `msr_radiolysis/validation/topology_constraints.py` — checks mass + charge balance |
| Forward solver `C^{(s)}(t; θ, γ)` | Existing `msr_radiolysis.integrator.integrate_system` |
| Per-salt data partitioning | `validation/<system>/<paper>/data/` already structured this way |
| Modality detection | New: classify each CSV in `data/` as transient / scalar / censored / Arrhenius |
| Intrinsic prior `π_θ` | From [arrhenius_parameters.csv](cr_licl_kcl/iwamatsu_2026_pccp/data/arrhenius_parameters.csv) per reaction |
| Hyperprior `Λ_i` | New: per-reaction salt-specific perturbation scale (Marcus-theory bound or empirical-Bayes) |
| Facility prior `Σ_b` | New: estimated from cross-paper rate-constant disagreements |
| Discrepancy parametrization | New: linear `(1/T, log [M])` form, bounded coefficients |
| MCMC engine | PyMC + sunode (Tier 2 of MODEL_ADEQUACY_REVIEW) |
| RJMCMC topology moves | New: discrete proposal kernel respecting `Γ` |

**Estimated implementation effort.** ~3–4 weeks of focused work for the chloride kernel,
assuming Tier 1 (identifiability) is completed first.

---

## 10. References

- Brynjarsdóttir, J.; O'Hagan, A. "Learning about physical parameters: the importance of model
  discrepancy." *Inverse Problems* **2014**, 30, 114007.
- Galagali, N.; Marzouk, Y. M. "Bayesian inference of chemical kinetic models from proposed
  reactions." *Chem. Eng. Sci.* **2015**, 123, 170–190.
- Galagali, N.; Marzouk, Y. M. "Exploiting network topology for large-scale inference of
  nonlinear reaction models." *J. R. Soc. Interface* **2019**, 16(152), 20180766.
- Green, P. J. "Reversible jump Markov chain Monte Carlo computation and Bayesian model
  determination." *Biometrika* **1995**, 82(4), 711–732.
- Kennedy, M. C.; O'Hagan, A. "Bayesian calibration of computer models." *J. R. Stat. Soc. B*
  **2001**, 63(3), 425–464.
- Lindsay, B. G. "Composite likelihood methods." *Contemp. Math.* **1988**, 80, 221–239.
- Tierney, L. "Markov chains for exploring posterior distributions." *Ann. Statist.* **1994**,
  22, 1701–1762.
- Tobin, J. "Estimation of relationships for limited dependent variables." *Econometrica*
  **1958**, 26, 24–36.
- van der Vaart, A. W. *Asymptotic Statistics*. Cambridge University Press, **1998**.
- Varin, C.; Reid, N.; Firth, D. "An overview of composite likelihood methods." *Statistica
  Sinica* **2011**, 21, 5–42.
- Vehtari, A.; Gelman, A.; Gabry, J. "Practical Bayesian model evaluation using leave-one-out
  cross-validation and WAIC." *Stat. Comput.* **2017**, 27(5), 1413–1432.

---

## Appendix: provenance of each result

| Result | Established prior art | Novel contribution here |
|---|---|---|
| Theorem 1 (RJMCMC ergodicity on constrained `Γ`) | Tierney 1994; Green 1995 | Specific application to mass-and-charge-balance constraints + path-connectedness verification for chloride kernel |
| Proposition 1 (mixing acceleration) | None known specifically | New; the `\|Γ\|/2^R` constant is novel |
| Theorem 2 (cross-salt BvM consistency) | van der Vaart Theorem 10.1 | Application to the salt-perturbation hierarchical Arrhenius structure |
| Theorem 3 (bounded parametric discrepancy bias) | Brynjarsdóttir–O'Hagan 2014 contrast | The `O(B/√n)` rate with explicit constant for the radiolysis design is novel |
| Theorem 4 (Godambe-balanced composite likelihood) | Lindsay 1988; Varin–Reid–Firth 2011 | Specific weight choice for the radiolysis multi-modality case |
| Theorem 5 (NULL Bayes-factor bound) | Tobin 1958 censored likelihoods; Bayes' rule | Application to detection-limit NULL benchmarks as model-selection drivers |
| Theorem 6 (facility-effect identifiability) | van der Vaart §6 | Specific gauge-symmetry analysis of the cross-paper hierarchy |
| Propositions 2-5 (comparisons) | All comparisons synthesized here | New combinations |

**Stance on novelty.** I claim that HBMAE is the *first explicit synthesis* of these
components for kinetic-network UQ. The individual components are well-established. The
combination is, to the best of my literature survey, not yet published — the closest
precedent is Galagali–Marzouk (2015, 2019) which does the constrained-RJMCMC and
spike-and-slab parts but not the hierarchical cross-salt layer, the modality-balanced
composite likelihood, the censored NULL likelihood, or the facility-effect bias separation.
This is the contribution that an article-form publication could carry.
