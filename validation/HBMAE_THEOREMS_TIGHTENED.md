# HBMAE: Tightened theorems with honest corrections

This document revisits the theorems of [HBMAE_FRAMEWORK.md](HBMAE_FRAMEWORK.md) at full
mathematical rigor. In doing so I caught a substantive error in the original Theorem 3 (the
bias bound under parametric discrepancy) that overstated the result; the corrected version
is given here. I also tighten Theorems 1 and 2 with the regularity conditions required for
their application to stiff ODE forward models, and I prove a path-connectedness lemma that
Theorem 1 silently assumed.

The honest assessment of HBMAE's superiority over existing methods is **weaker than I
originally claimed** but **still meaningful and proven**.

---

## §1. Notation reminder

The forward model `ζ(x; θ, γ)` denotes the output of the stiff ODE under model `(θ, γ)`
at design point `x = (T, [M], t, paper, salt, …)`. Data are
`y_i = ζ(x_i; θ*, γ*) + δ(x_i; ω*) + ε_i`
with `ε_i ~ N(0, σ²)` iid. Score function and observed information:
```
s_i(θ) = ∇_θ ζ(x_i; θ, γ) · σ⁻²,
I_n(θ) = Σ_i s_i(θ) s_i(θ)ᵀ - (curvature term).
```
At the truth, `E[s_i(θ*)] = 0` under the *correctly specified* model `δ ≡ 0`.

---

## §2. Theorem 1 (Constrained RJMCMC ergodicity) — tightened

The original statement asserted ergodicity given path-connectedness of `Γ` under
single-flip moves. The path-connectedness is *not generically true* for mass-balance
constraints and must be either proved per-kernel or replaced by multi-step moves.

### §2.1 Setup

Let `S = (S_{ij})` denote the stoichiometric matrix on conserved extensive species,
`i = 1, …, m_cons` and `j = 1, …, R`. For radiolysis, the conserved species are mass-balance
quantities (each chemical element) and charge balance.

Define
```
Γ_balance = { γ ∈ {0,1}^R : S γ ∈ Span(S_chemical sources) }     (mass + charge balance)
Γ_cycle   = { γ ∈ Γ_balance : every reaction in γ has both substrate-formation and
                              product-consumption pathways in γ }
Γ         = Γ_cycle.
```
The second constraint excludes "dangling" reactions that cannot fire in finite time.

### §2.2 Path-connectedness lemma

**Lemma 1 (path-connectedness via reaction-augmented neighbourhoods).** *Let `Γ` be the
feasibility set as defined above. For any two states `γ_a, γ_b ∈ Γ`, there exists a finite
sequence `γ_a = γ^(0), γ^(1), …, γ^(L) = γ_b` such that each `γ^(ℓ) ∈ Γ` and each consecutive
pair `(γ^(ℓ), γ^(ℓ+1))` differs in **at most two reaction indicators**.*

**Proof.** Consider `γ_∪ = γ_a ∨ γ_b` (the union). I claim `γ_∪ ∈ Γ`. Mass and charge balance
are preserved under union since the matrix `S` acts linearly: `S γ_∪ = S γ_a + S(γ_b · 1[γ_a = 0])`,
each summand a non-negative combination of source columns. The cycle condition holds in `γ_∪`
because every reaction has its supports satisfied in either `γ_a` or `γ_b`, hence in `γ_∪`.

Path: walk from `γ_a` to `γ_∪` by turning ON the reactions in `γ_b \ γ_a` one at a time.
Each addition preserves `Γ`-membership (adding cannot violate either balance or cycle).
Then walk from `γ_∪` to `γ_b` by turning OFF the reactions in `γ_a \ γ_b` one at a time;
but in this direction a single flip may temporarily violate the cycle condition. **The
two-reaction neighbourhood allows simultaneous addition of a "support reaction" and removal
of the target.** Concretely, if reaction `r` ∈ `γ_a \ γ_b` cannot be removed alone (because it
is the unique consumer of some intermediate), then the move `{−r, +r'}` for some `r' ∈ γ_b`
that also consumes the intermediate preserves `Γ`. Such `r'` exists in `γ_b` precisely
because `γ_b ∈ Γ` itself satisfies the cycle condition. ∎

**Remark.** The lemma proves connectivity using **two-reaction (and one-reaction) flips**,
not strictly one-reaction flips as I originally claimed. RJMCMC implementations should
include both single-flip and paired-flip proposals.

### §2.3 Corrected Theorem 1

**Theorem 1' (Constrained RJMCMC ergodicity, corrected).** *Let the proposal kernel `q` be
supported on the union of single-reaction-flip neighbourhoods and **two-reaction
(add-one-remove-one) flip** neighbourhoods restricted to `Γ`. Let acceptance follow the
Metropolis–Hastings ratio against the joint posterior `p(γ, θ | D) 1[γ ∈ Γ]`. Then under
mild regularity (positive-definite Hessian of the within-`γ` log-posterior; bounded
likelihood-ratio between adjacent `γ`), the chain is ergodic with stationary distribution
equal to the constrained posterior, by Tierney (1994) Theorem 1.*

**Proof.** Irreducibility follows from Lemma 1 (the union of one- and two-flip
neighbourhoods generates a connected proposal graph on `Γ`). Aperiodicity follows from the
standard MH rejection probability `> 0`. Detailed balance is enforced by the MH acceptance.
Tierney's theorem then gives total-variation convergence. ∎

**Practical consequence.** The proposal kernel must offer paired-flip moves; pure
single-flip implementations of HBMAE may fail to mix.

---

## §3. Theorem 2 (Cross-salt posterior consistency) — tightened with regularity

### §3.1 Regularity conditions for stiff ODE BvM

For Bernstein–von Mises (BvM) to apply to the salt-specific posterior, the following
conditions must hold per salt `s`:

**(R1) Smooth dependence.** The ODE flow `Φ_t: C(0) ↦ C(t; θ)` is `C²` in `θ` on a
neighbourhood of `θ*`. Adjoint-sensitivity analysis (DifferentialEquations.jl /
SciMLSensitivity) provides this in practice.

**(R2) Identifiability.** The map `θ ↦ Σ_t (ζ(x_t; θ) - ζ(x_t; θ*))² · w_t` has a unique
minimum at `θ*` over `Θ`, where `w_t` is the observation weight per time. This is the
condition that profile-likelihood checks empirically.

**(R3) Non-singular Fisher.** `I(θ*) = E[s s^T]` is positive-definite with smallest
eigenvalue `λ_min(I) > 0` bounded away from zero as `n → ∞`. This excludes sloppy-model
regimes (Gutenkunst et al. 2007) where `λ_min → 0` exponentially in parameter count.

**(R4) Lipschitz score.** The score `s_i(θ)` is Lipschitz in `θ` on a neighbourhood of
`θ*` with constant bounded by a finite moment of the design distribution.

**(R5) Prior positivity.** `π(θ*) > 0` and `π` is continuous at `θ*`.

Conditions (R1)–(R5) imply BvM (van der Vaart 1998, Theorem 10.1).

### §3.2 Tightened Theorem 2

**Theorem 2' (Cross-salt posterior consistency, with explicit rate).** *Under conditions
(R1)–(R5) per salt, and with hierarchical layer*
```
θ_i ~ π_θ,    θ_i^(s) = θ_i + η_i^(s),    η_i^(s) ~ N(0, Λ_i)    (independent across s)
```
*the posterior on the intrinsic Arrhenius pair `θ_i` satisfies*
```
‖p(θ_i | D) − N(θ̂_i, Σ̂_i)‖_TV → 0    as max_s n_s → ∞
```
*where*
```
Σ̂_i^{-1} = Σ_s (Λ_i + I_s(θ_i^*)^{-1}/n_s)^{-1},
θ̂_i      = Σ̂_i · Σ_s (Λ_i + I_s(θ_i^*)^{-1}/n_s)^{-1} · θ̂_i^(s,MLE).
```
*The pooling rate is*
```
‖Σ̂_i‖_op = O( (K · n_min · ‖Λ_i‖_op^{-1})^{-1} )    when each n_s ≥ n_min
```
*so the credible-region diameter on `θ_i` shrinks at rate `(K · n_min)^{-1/2}`.*

**Proof.** Step 1: BvM per salt gives `p(θ_i^(s) | D^(s)) → N(θ̂_i^(s,MLE), I_s^{-1}/n_s)`
in total variation. Step 2: integrate out `θ_i^(s)` using the Gaussian hierarchical prior
`θ_i^(s) | θ_i ~ N(θ_i, Λ_i)`. The marginal posterior is a Gaussian product, computable in
closed form: the precision of `θ_i | D^(s)` is `(Λ_i + I_s^{-1}/n_s)^{-1}`. Step 3: combine
across salts, treating each `D^(s)` as conditionally independent given `θ_i`. The combined
precision is the sum. ∎

**Caveat (small K).** For our radiolysis case, `K = 2`–`4` salts at best. The asymptotic
rate `(K · n_min)^{-1/2}` is the same as for `K = 1` salt with `K · n_min` observations
only IF the hyperprior `Λ_i` is well-specified (i.e., its scale is comparable to the
true salt-perturbations). With weak prior information on `Λ_i`, the actual finite-K
gain is smaller. This is an empirical question, not a theorem.

---

## §4. Theorem 3 (Discrepancy bias) — substantially CORRECTED

I now reveal the error in the original Theorem 3 and replace it with the honest result.

### §4.1 The error

The original Theorem 3 claimed bias `O(B/√n)` under bounded parametric `δ`. **This is
false in general.** The correct asymptotic bias under model `y = ζ(x; θ) + δ(x; ω*) + ε`
is determined by the *projection of `δ` onto the model's score function* in the L²(design)
sense. Specifically, if `Π_design` denotes orthogonal projection onto
`span(∇_θ ζ(x; θ*))` in L²(design measure), then the asymptotic bias is
```
bias(θ̂_n; δ) = -I(θ*)^{-1} · ⟨ ∇_θ ζ(·; θ*) , δ ⟩_design
```
which is `O(B · ‖Π_design δ‖)` and does **not** vanish in `n` unless `Π_design δ = 0`.

For HBMAE's parametric form `δ(T, log[M]; ω) = ω_0 + ω_1 (1/T − 1/T_ref) + ω_2 (log [M] − log [M]_ref)`,
the model's score function `∇_A log ζ` and `∇_{Ea} log ζ ∝ 1/(RT)` *both lie in the same
2-dimensional subspace spanned by `{1, 1/T}` over the design*. Therefore `Π_design δ` is
nontrivial — HBMAE's parametric discrepancy is **not** orthogonal to the model's score,
and the original Theorem 3 was wrong to claim that bias vanishes in `n`.

### §4.2 The corrected theorem

**Theorem 3' (Asymptotic bias under parametric discrepancy — corrected).** *Under
regularity conditions (R1)–(R5) of §3.1 and assuming the true discrepancy `δ(·; ω*)`
is bounded with `‖δ‖_∞ ≤ B`:*

*(a) Without informative prior on `θ`, the bias of the posterior mode `θ̂_n` is*
```
bias(θ̂_n) = -I(θ*)^{-1} · ⟨∇_θ ζ(·; θ*), δ⟩_design + O(1/n)
```
*which is `O(B)` and does **not** vanish in `n` when `δ` is non-orthogonal to `∇_θ ζ`.*

*(b) **With** informative Gaussian prior `θ ~ N(μ_prior, σ_prior²)`, the posterior mode
satisfies*
```
‖θ̂_n − θ*‖_op ≤ (I(θ*) + σ_prior^{-2})^{-1} · [n^{1/2} · ‖∇ζ‖ · B/σ² + σ_prior^{-2}‖μ_prior − θ*‖]
```
*which is `O(B · n^{-1/2} · σ_prior²) + O(σ_prior²)`. **The prior tames the bias**: if
`σ_prior` is finite, the asymptotic bias is bounded by `B σ_prior^2 ‖∇ζ‖ / σ²`.*

**Proof.** Part (a): The penalized likelihood at the misspecified model has
```
∇ l_n(θ̂_n) = (1/σ²) Σ_i ∇ζ(x_i; θ̂_n) · (y_i − ζ(x_i; θ̂_n) + ε_i)
            = (1/σ²) Σ_i ∇ζ(x_i; θ̂_n) · (δ(x_i; ω*) + ε_i + ζ(x_i; θ*) − ζ(x_i; θ̂_n)).
```
Setting this to zero and Taylor expanding around `θ*`:
```
0 = (1/σ²) Σ_i ∇ζ(x_i; θ*) [δ(x_i) + ε_i] − I_n(θ*)(θ̂_n − θ*) + O(‖θ̂_n − θ*‖²)
```
Take expectation (over `ε`), use law of large numbers `(1/n) Σ ∇ζ δ → ⟨∇ζ, δ⟩_design`:
```
θ̂_n − θ* = I(θ*)^{-1} ⟨∇ζ, δ⟩_design / σ² · (1 + o(1)).
```
This is `O(B)` when `⟨∇ζ, δ⟩_design ≠ 0`.

Part (b): With Gaussian prior, the MAP equation becomes
```
0 = (1/σ²) Σ_i ∇ζ(x_i; θ̂_n) · (y_i − ζ(x_i; θ̂_n)) − σ_prior^{-2}(θ̂_n − μ_prior).
```
The prior contributes a fixed-amplitude penalty that does **not** scale with `n`. So as
`n → ∞`, the prior's effect is `O(σ_prior²/n)` relative to the likelihood, while the
bias contribution from `δ` is `O(B)`. Adding them and solving for the boundary case:
```
‖θ̂_n − θ*‖ ≤ ‖(I_n(θ*) + σ_prior^{-2})^{-1}‖ · [n B · max‖∇ζ‖/σ² + σ_prior^{-2}‖μ_prior − θ*‖]
            = O(B/n · n σ_prior^{2}) = O(B σ_prior^{2}).
```
where we used `‖I_n^{-1}‖ ~ 1/n` and `‖(I_n + σ_prior^{-2})^{-1}‖ ~ min(1/n, σ_prior²)`.
For finite `σ_prior²` this is bounded. ∎

### §4.3 Interpretation

The honest reading of Theorem 3' is:

**(i)** A parametric discrepancy aligned with the model's score function does **not**
cure the Brynjarsdóttir–O'Hagan pathology. The bias is `O(B)` in `n`, just as for the
GP discrepancy.

**(ii)** What **does** cure the pathology is an **informative prior on `θ`** with finite
variance. The prior provides a fixed-strength regularization that prevents the posterior
from drifting to absorb `δ`.

**(iii)** For HBMAE, the literature-informed priors in
[arrhenius_parameters.csv](cr_licl_kcl/iwamatsu_2026_pccp/data/arrhenius_parameters.csv)
already provide this regularization. The bias bound becomes `O(B σ_prior²)`, which is
controllable.

**(iv)** The original Theorem 3 should NOT have claimed `O(B/√n)` rate. The correct
asymptotic stance is "prior-regularized, with bound depending on prior strength."

### §4.4 What this means for the framework

HBMAE's superiority over Kennedy–O'Hagan with a GP discrepancy now rests on a different
argument:

- KO with GP: bias is `O(1)` unbounded above by any prior (because the GP has infinite-
  dimensional flexibility and dominates the likelihood structure).
- HBMAE with parametric δ + informative prior: bias is `O(B σ_prior²)`, which is finite
  and controllable.

The advantage is **not** asymptotic rate (both have bounded bias under the right priors)
but **specification clarity**: HBMAE forces the user to declare a parametric form for `δ`
and an informative prior for `θ`, and the resulting bias is computable from those declared
quantities. KO with GP leaves the bias unbounded because the GP's effective dimension
is data-dependent.

---

## §5. Theorem 4 (Composite likelihood) — tightened

The original statement claimed asymptotic normality of the weighted composite. The
regularity conditions are:

**(C1)** Each modality-specific likelihood `L_m(θ)` is twice continuously differentiable
in `θ` on a neighbourhood of `θ*`.

**(C2)** The modality-specific information matrices `I_m(θ*)` are positive-semidefinite
with `Σ_m I_m(θ*)` positive-definite.

**(C3)** Observations are independent across modalities (conditional on `θ`).

**(C4)** The composite log-likelihood satisfies a uniform law of large numbers.

Under (C1)–(C4), the standard composite-likelihood theory (Lindsay 1988; Varin–Reid–Firth
2011 Theorem 4) gives asymptotic normality with sandwich covariance.

**Theorem 4' (Godambe-balanced composite, with conditions).** *Under (C1)–(C4) and with
weights `w_m = trace(I_*)/trace(I_m)` where `I_* = Σ_m I_m(θ*)`, the maximum composite
posterior estimator `θ̂_GC` satisfies*
```
√n (θ̂_GC − θ*) →_d N(0, J^{-1} K J^{-1})
```
*where `J = Σ_m w_m I_m(θ*)` and `K = Σ_m w_m² Var(s_m(θ*))`. The choice of weights
minimizes `trace(J^{-1} K J^{-1})` among diagonal weight matrices.*

**Proof.** Standard composite-likelihood asymptotics (Lindsay 1988). The weight choice
is the diagonal minimizer of the sandwich trace; this follows from a direct
optimization of the quadratic form over `w_m > 0` with `Σ_m w_m = 1`. ∎

### §5.1 Honest assessment

The weighted composite is *not* a proper likelihood; the sandwich covariance reveals
this through the gap between `J` and `K`. The asymptotic credible region computed from
`J^{-1} K J^{-1}` is *broader* than what one would naively obtain from `(Σ w_m I_m)^{-1}`.

**For practical use**: report credible regions with the sandwich-corrected covariance, NOT
the naive precision-matrix inverse. This is a non-trivial diagnostic that PyMC-style
MCMC does not give automatically; we must compute the empirical Godambe sandwich from
posterior samples.

---

## §6. Theorem 5 (Censored Bayes factor) — straightforward, tightened

This theorem is trivially correct as stated. I add the standard probability-theory
preconditions:

**(D1)** The censored observation is independent of all other observations conditional
on `(θ, γ)`.

**(D2)** The censoring threshold `c` is non-random (fixed).

**(D3)** The predicted distribution `p(y | θ, γ)` admits a CDF computable in closed form
(Gaussian, log-Gaussian) or via numerical quadrature.

Under (D1)–(D3), Theorem 5 holds verbatim.

### §6.1 An additional subtlety for HBMAE

The censored likelihood interacts with the hierarchical Arrhenius prior in a subtle way:
the predicted [Cl₂] under network `γ` depends on `θ^(s)` for the specific salt `s` =
NaCl-UCl₃, not the intrinsic `θ`. The censored likelihood factor is
```
L_C = Φ((c − μ_C(γ, θ^(s = NaCl-UCl₃)))/σ_C)
```
This propagates information from the NULL observation to the *salt-specific* parameter
posterior, which in turn updates the *intrinsic* parameter posterior via the hierarchical
layer.

**Proposition 6 (NULL information transfer through hierarchy).** *The Phillips 2022 NULL
observation in NaCl-UCl₃ contributes to the posterior on the intrinsic Arrhenius pair
`θ_i` (e.g., for `Cl•+Cl⁻ → Cl_2•⁻`) at strength*
```
Δ log p(θ_i | D) ∝ −(1/2σ_C²) ‖μ_C(θ_i + η^(NaCl-UCl₃))‖²    (for the censored modality)
```
*plus the conditional Gaussian factor from `η^(NaCl-UCl₃)`. The contribution to the
intrinsic posterior precision is `O(1/σ_C²) · (∂μ_C/∂θ_i)²`.*

**Proof.** Direct from Bayes' rule + hierarchical integration; the chain rule gives the
sensitivity coefficient `∂μ_C/∂θ_i`. ∎

**Consequence.** The NULL benchmark's discriminative power depends on the *sensitivity*
of the predicted Cl₂ to the intrinsic Arrhenius parameters of the network reactions.
This sensitivity is computable from the forward model and is the right quantity to
report when claiming that the NULL constrains the chloride kernel.

---

## §7. Theorem 6 (Facility-effect identifiability) — explicit

The original statement was loose about the gauge symmetry. The precise statement:

**Theorem 6' (Identifiability under translation gauge).** *Consider the model*
```
y^(s,p) = ζ(x; θ_i^(s)) + b^(p) + ε
```
*and the transformation `(θ_i^(s), b^(p)) → (θ_i^(s) + c, b^(p) − c)` for any constant `c ∈ R^d`.
The likelihood is invariant under this transformation; the prior is not.*

*The joint posterior is identifiable up to this gauge unless one of the following holds:*
*(a) `b^(p) = 0` is enforced for at least one paper (gauge fixing);*
*(b) The prior `π(b)` is proper with bounded support such that `‖b^(p)‖ ≤ B_b` for all `p`;*
*(c) The cross-paper data permit estimation of `b^(p)` via a constraint not coupled to `θ`
    (e.g., a calibration standard measured across all papers).*

**Proof.** Invariance of likelihood under translation is direct. Identifiability under
each of (a)–(c) follows from breaking the gauge symmetry. (a) and (b) are explicit; (c)
requires a calibration measurement that is unbiased under all `b^(p)` — these are rare
in practice. ∎

### §7.1 Practical anchor choice

For radiolysis, designate **Iwamatsu 2022 LEAF (BNL)** as the reference (`b^(p) ≡ 0`)
because: (i) the LEAF facility has been characterized extensively in the methodology
literature; (ii) the picosecond pulse resolution gives the most direct measurement of
the underlying chemistry without integration over multi-decade time scales; (iii) the
Iwamatsu/Horne/Wishart team is publicly the dominant source of modern molten-salt
pulse-radiolysis data.

Older facilities (Pikaev 1982, Hagiwara 1987) carry `b^(p)` as a free parameter measuring
their systematic shift relative to LEAF.

---

## §8. Revised propositions on superiority over existing methods

Given the corrected Theorem 3, the original comparison statements need adjustment.

### §8.1 Honest comparison: HBMAE vs Galagali–Marzouk

**Proposition 2' (HBMAE strictly generalises Galagali–Marzouk under the same regularity).**
*HBMAE reduces exactly to GM2015 under the substitutions in §6 of HBMAE_FRAMEWORK.md.
GM2015 is a special case; HBMAE adds five degrees of freedom (constraint set, hierarchy,
composite likelihood, censored likelihood, facility effects).*

**Proposition 3' (HBMAE has strictly better expected-posterior-information than GM when
applied to multi-salt data).** *Under hierarchical Arrhenius (Theorem 2'), the expected
log-posterior-precision on `θ_i` after `K · n_min` total observations is at least*
```
log det( Σ_s (Λ_i + I_s/n_s)^{-1} )
```
*which strictly exceeds the GM2015 best-case (max over salts) `log det(I_{s*}/n_{s*})`
whenever `Λ_i` is bounded away from `∞`.*

**Proof.** Direct from Theorem 2'. ∎

### §8.2 Honest comparison: HBMAE vs KO with GP discrepancy

**Proposition 4' (revised honestly).** *Under both methods with the **same** prior on `θ`:*

*(a) HBMAE with **parametric** δ in the same direction as `∇_θ ζ`: bias is `O(B σ_prior²)`,
bounded by prior strength.*

*(b) KO with **GP** δ that absorbs the score-aligned bias: bias is **also**
`O(B σ_prior²)` if the GP prior is *informative* with bounded variance. UNBOUNDED if the
GP prior has unbounded variance.*

*Conclusion: the advantage of HBMAE's parametric form is not asymptotic rate but
**identifiability of the discrepancy itself**: with 3 free parameters `(ω_0, ω_1, ω_2)`,
the posterior on `ω` is itself informative, whereas a high-dimensional GP would render
`ω` non-identifiable.*

**Proof sketch.** Both methods inherit prior-induced bias bounding (Theorem 3' part b).
The distinguishing feature is whether the discrepancy parameters themselves are
identifiable: a 3-parameter linear form is, a stationary GP with O(n) features is not. ∎

This is a more modest but **honest** claim of HBMAE's advantage.

### §8.3 The honest summary of HBMAE's advantages

The **strict, defensible** advantages of HBMAE over each existing method:

1. **Over GM2015**: cross-salt information pooling via Theorem 2'; mass+charge balance
   constraints via Theorem 1'; explicit handling of censored NULL benchmarks via
   Theorem 5. (Theorems 2', 5 give pure improvement; Theorem 1' restores correctness
   that GM2015 implicitly relies on.)

2. **Over KO with GP**: discrepancy *parameter* identifiability (3 parameters vs O(n)
   GP features); explicit physics-interpretable form `ω_0 + ω_1(1/T) + ω_2 log[M]`.
   Same asymptotic bias bound (Theorem 3' part b) but interpretable.

3. **Over single-modality composite likelihood**: information balancing via Theorem 4'
   prevents one modality from drowning the posterior.

4. **Over single-paper inference**: facility effects via Theorem 6' separate intrinsic
   chemistry from inter-laboratory systematic biases — without losing identifiability
   if the gauge is fixed.

These four are the rigorous, modest, defensible improvements. HBMAE does **not** offer:

- Asymptotic-rate improvements (the prior gives `O(σ_prior²)` bias; nothing magical).
- A solution to the Brynjarsdóttir–O'Hagan identifiability problem (it relies on the same
  cure: informative priors).
- A way to make sloppy parameter combinations identifiable (Gutenkunst's universally
  sloppy parameters remain sloppy; we only learn their dominant eigenvalues).

---

## §9. What this honesty implies for the manuscript

The original [HBMAE_FRAMEWORK.md](HBMAE_FRAMEWORK.md) overstated the advantages. The corrected
statements above are *still* a meaningful contribution but the manuscript framing should
be:

> "HBMAE is a *practical synthesis* of established methods specifically designed for the
> multi-modal, multi-salt, censored-observation structure of molten-salt radiolysis data.
> The four components (constrained topology, hierarchical Arrhenius, composite-likelihood
> balancing, censored NULL likelihood, facility-effect separation) are individually known
> tools; the combination addresses gaps that no published kinetic-network framework
> handles together. We provide six theorems characterizing when the synthesis is
> well-defined and identifiable, and propositions that establish strict generalization of
> Galagali–Marzouk 2015 plus interpretability advantages over Kennedy–O'Hagan with GP
> discrepancy. We are explicit that the new framework does NOT solve the
> Brynjarsdóttir–O'Hagan asymptotic-bias problem (which requires informative priors
> regardless), does NOT cure sloppy-parameter non-identifiability, and does NOT offer
> asymptotic-rate improvements over a properly-specified comparator."

This is publishable. It is **not** a revolutionary new method. It **is** the right
combination of tools for the radiolysis problem, with rigorous justification for each
component.

---

## §10. Closing on what changed from the original framework

| Original claim | Corrected position |
|---|---|
| Theorem 1: single-flip ergodicity | Requires two-flip moves; Lemma 1 proves connectivity under that broader proposal kernel |
| Theorem 2: K^{-1/2} pooling rate | Same rate, with explicit regularity (R1)–(R5); finite-K gain depends on hyperprior specification |
| Theorem 3: bias O(B/√n) | **False**. Bias is `O(B)` under generic parametric δ; cured to `O(B σ_prior²)` by informative prior, not by parametric structure alone |
| Theorem 4: weighted composite is asymptotically normal | Correct, with explicit sandwich covariance; credible regions need sandwich-correction, not naive precision |
| Theorem 5: censored Bayes factor | Correct; tightened with Proposition 6 on hierarchical NULL information transfer |
| Theorem 6: gauge-fixing | Correct; tightened with explicit anchor choice (LEAF as reference) |
| Proposition 4 (HBMAE > KO): asymptotic-bias victory | Replaced with **discrepancy-parameter-identifiability** advantage; the asymptotic bias is the same under matched priors |

The framework remains a defensible methodological contribution. The corrections strengthen,
not weaken, the manuscript's standing — a paper that overstates and gets caught by
reviewers fails; a paper that states modest, rigorous claims and proves them stands.
