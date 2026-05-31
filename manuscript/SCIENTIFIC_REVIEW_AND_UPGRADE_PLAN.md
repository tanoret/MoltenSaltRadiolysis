# Scientific review and upgrade plan

Date: 2026-05-29

Scope: critique `article_paper1.tex` and `article_paper2.tex` as a two-paper
series on HBMAE and molten-salt radiolysis, with emphasis on making them
credible top-tier submissions.

## Style diagnosis

The writing style is clear, pedagogical, and unusually direct. That is a real
asset. The strongest passages explain why a mathematical object exists before
giving the equation. Keep that.

For journal submission, the style should become "pedagogical but not tutorial":
retain plain-language bridges and interpretive paragraphs, but remove textbook
detours, exercises, and claims that sound like internal project notes. The
reader should feel guided, not taught from scratch.

## Literature map that should anchor the introductions

The field history should be framed in three layers.

1. Legacy molten-halide radiation chemistry.

   - Pikaev et al. 1982: broad pulse-radiolysis review of solvated electrons
     and halogen molecular radical anions in molten alkali halides.
   - Makarov et al. 1982: oxidizing products of alkali-halide melt radiolysis,
     including X2 radical-anion spectra and rates.
   - Hagiwara et al. 1987: LiCl-KCl pulse radiolysis, assigning transients to
     `e_s-` and `Cl2.-`.
   - Akiyama et al. 1994: only known fluoride-melt pulse-radiolysis transient
     paper for LiF-KF and FLiNaK-like systems.

2. Modern LEAF and MSR-relevant chloride chemistry.

   - Iwamatsu/Horne/Wishart 2022: Zn2+ kinetics in molten LiCl-KCl and the
     key critique of legacy microsecond pulse widths.
   - Conrad et al. 2023: iodide perturbs the transient pool by adding ICl.-.
   - Iwamatsu et al. 2026: Cr(II)/Cr(III) redox kinetics with both `e_s-`
     and `Cl2.-`.
   - Castro Baldivieso et al. 2026: Nd kinetics in LiCl-KCl.
   - Rotermund et al. 2024: aqueous chloride radical/actinide analogue data,
     useful but not molten-salt evidence.

3. UQ and kinetic-network methodology.

   - Galagali and Marzouk 2015/2019: closest Bayesian reaction-network
     inference precedent.
   - Kennedy and O'Hagan 2001 plus Brynjarsdottir and O'Hagan 2014: model
     discrepancy and the identifiability hazard.
   - Frenklach/B2BDC/Oreluk-Hegde-Packard: deterministic consistency and model
     discrepancy as a comparator, not a Bayesian substitute.
   - Varin/Lindsay composite likelihood, Vehtari/Gelman/Gabry PSIS-LOO, Raue
     profile likelihood, and sloppy-model literature should be used as
     supporting method foundations.

## Critical issues to fix first

### 1. Paper 2 abstract overclaims the chemistry-feature layer

`article_paper2.tex` says the meta-hierarchical layer recovers held-out Cf and
Cr within +/-0.3 log units, but the LOMO table reports Cr median error 2.41
log units and Cf error 4.35 log units. This is a direct contradiction.

Fix: rewrite the abstract to say the layer works inside the Pikaev-like
interpolation envelope and fails, usefully and visibly, for LEAF/aqueous
out-of-envelope metals. The failure is scientifically useful. Do not hide it.

### 2. Paper 2 uses the wrong logarithm for the Phillips Bayes factor

The text reports `log10 K12 = 316.9`, then converts that to `~10^138`. These
cannot both be true. The validation reports show `log L = -316.9` in natural
log units, so `ln K ~= 316.9` and `log10 K ~= 137.6`.

Fix: use one convention everywhere. Recommended wording:
`ln K12 = 316.9`, equivalently `log10 K12 = 137.6`.

### 3. The composite-likelihood weights are internally inconsistent

In Paper 1, the framework defines
`w_m = tr(I_*) / tr(I_m)`, but Theorem 4 states the optimal simplex weight as
`w_m = tr(I_m) / sum_m tr(I_m)`. These are opposite directions.

Fix: decide the actual objective. If the goal is to prevent dense transient
traces from dominating, use normalized inverse information or effective-sample
information weights and prove that objective. If the goal is minimum sandwich
variance under a standard simplex constraint, use the theorem's weights and
accept that high-information modalities receive more weight.

### 4. Theorem 2 currently overstates the cross-salt contraction rate

With random host effects `theta_s = theta + eta_s`, fixed nonzero `Lambda`
means that as within-host sample size grows, each host reveals `theta + eta_s`,
not `theta` itself. The intrinsic posterior variance approaches `Lambda / K`
up to constants; it does not continue contracting as `1/(K n)` unless `Lambda`
is negligible, shrinks with `n`, or the theorem is stated in a measurement-error
dominated regime.

Fix: state two regimes:

- measurement-error dominated: `n_s^{-1} I_s^{-1} >> Lambda`, pooling behaves
  approximately like `1/(sum_s n_s)`;
- host-heterogeneity dominated: `n_s^{-1} I_s^{-1} << Lambda`, intrinsic
  precision grows as `K Lambda^{-1}`.

This correction will make the theorem much more credible.

### 5. Theorem 3 should not claim fixed informative priors cure asymptotic bias

A fixed Gaussian prior does not dominate the likelihood as `n -> infinity`.
It regularizes finite-data inference, and it can bound a prior-dominated MAP
regime, but it does not generically remove asymptotic discrepancy bias unless
the discrepancy is constrained, orthogonalized to the score, or the prior
precision scales with data size.

Fix: present this as finite-sample regularization plus domain-restricted
discrepancy, not as an asymptotic cure.

### 6. Theorem 1 should become either kernel-specific or assumption-based

The proof that the feasible network set is connected under two-flip moves is
not generic for arbitrary mass/charge/cycle constraints. The local enumeration
is valuable, but it proves connectivity for the 10-reaction kernel, not for all
future kernels.

Fix: make graph connectivity an explicit assumption for the general theorem,
then present the enumeration as constructive verification for the chloride-Cr
kernel.

### 7. MCMC diagnostics need to be stronger

The papers report `emcee` walker counts and effective samples, but a top-tier
reviewer will ask whether walkers were treated as independent chains. Use
rank-normalized split R-hat on independent chains, autocorrelation estimates,
ESS for every reported parameter, posterior predictive checks per dataset,
and prior-sensitivity runs saved as reproducible outputs.

## Paper 1 upgrade path

Target identity: rigorous methods paper with a molten-salt radiolysis case
study.

Required upgrades:

1. Lead the introduction with the general problem: heterogeneous calibration
   of stiff chemical kinetic networks with topology uncertainty, censored data,
   facility bias, and cross-host transfer. Introduce molten salts as the
   motivating system.
2. Replace "HBMAE is uniquely superior" language with "HBMAE is the first
   demonstrated synthesis for this radiolysis use case."
3. Add a formal model table: data units, likelihoods, priors, gauge choices,
   hyperparameters, and whether each parameter is learned, fixed, or by analogy.
4. Correct the theorem statements before expanding proofs.
5. Expand the comparator section so each comparator has a fair implementation,
   a mathematical scope, and a documented failure mode on your data.
6. Make validation reproducibility reviewer-proof: exact scripts, random seeds,
   chain diagnostics, posterior predictive plots, and sensitivity outputs.

## Paper 2 upgrade path

Target identity: application paper that provides the first uncertainty-quantified
multi-salt radiolysis kernel, while being painfully honest about what is and is
not calibrated.

Required upgrades:

1. Reframe "full body of literature" as "the digitized corpus assembled here"
   unless every cited paper contributes numerical calibration data.
2. Separate chloride and fluoride confidence levels. Chloride has genuine
   transient kinetics and scalar rates. Fluoride is currently a static screening
   kernel anchored by Davis and Toth-Felker, with Akiyama as an unmet transient
   target.
3. Add a calibrated/assumed/analogy flag to every reaction in the chloride and
   fluoride tables.
4. Treat Rotermund Cf as aqueous chloride analogue evidence, not molten-salt
   calibration evidence.
5. Recast the chemistry-feature layer as a domain-of-applicability tool, not a
   reliable extrapolator. Add Mahalanobis-distance or leverage diagnostics in
   feature space.
6. Recast operational predictions as scenario analyses until the engineering
   assumptions have stronger support: dose deposition, redox buffer depletion,
   cover-gas volume, Henry constants, mass transfer, fission-product inventory,
   and corrosion threshold.
7. Replace "notional safety limit" claims with "screening threshold" unless an
   accepted design-basis basis is cited.

## Best next work packages

1. Correct the three hard contradictions: LOMO abstract, Bayes-factor base, and
   composite weights.
2. Rewrite both introductions using the literature map above.
3. Correct Theorems 1-4 before polishing prose.
4. Add a "claim audit" table: every headline claim, evidence artifact, script,
   data source, and uncertainty caveat.
5. Re-run or at least regenerate diagnostics for Tier 2 and Tier 4 after the
   time-sort bug fix, using independent chains.
6. Add reaction provenance tables with calibrated / literature / analogy /
   placeholder labels.
7. Strengthen Paper 2's operational scenario with sensitivity to redox-buffer
   depletion, dose rate, Henry constant, cover-gas volume, and detection limit.

## Bottom line

The core idea is strong enough for a top article: molten-salt radiolysis badly
needs an uncertainty-quantified, multi-source kinetic calibration framework, and
HBMAE is a plausible first synthesis. The current drafts are not yet top-tier
because a few mathematical statements and headline claims outrun the evidence.

The path forward is not to make the papers louder. It is to make them more
strict: narrower claims, cleaner theorem regimes, stronger diagnostics, and a
more honest separation between calibrated facts, analogical assumptions, and
operational scenario extrapolations.
