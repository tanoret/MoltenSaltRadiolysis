# Article plan: pedagogical journal article on HBMAE

This document is the planning artifact for converting the work in
[validation/](../validation/) and [manuscript/methods.tex](methods.tex) into a
publishable, pedagogically-accessible journal article.

---

## §1. Field scan summary — what we are landing among

I surveyed three classes of comparable work to anchor the article's tone and structure.

### §1.1 Pedagogical / tutorial methodology papers (style template)

| Paper | Venue | What we learn from it |
|---|---|---|
| Schnoerr, Sanguinetti, Grima. "Approximation and inference methods for stochastic biochemical kinetics — a tutorial review." *J. Phys. A: Math. Theor.* **2017**, 50, 093001. [arXiv:1608.06582](https://arxiv.org/abs/1608.06582) | Self-contained intro → CRN → master equation → approximations → inference. Builds from first principles. Numerical case study to compare methods. **Best closest template for our target tone.** |
| Higdon, Williams, Gattiker (LANL). Tutorial on Bayesian calibration of computer models. Slides at [bayesint.github.io](http://bayesint.github.io/slides/bayesCal-INT16.pdf) | Pedagogical examples with a simple physics problem (drop height ↔ drop time). Shows that a tutorial Bayesian-calibration paper can keep the formalism light when motivated by a worked example. |
| Vehtari, Gelman, Gabry. PSIS-LOO. *Stat. Comput.* **2017** | Reference for how a methodological paper can land in a stats journal yet be readable by applied users, with explicit "how to use this" boxes. |

### §1.2 Methodological precedents in our specific niche

| Paper | Venue | Relation to ours |
|---|---|---|
| Galagali, Marzouk. "Bayesian inference of chemical kinetic models from proposed reactions." *Chem. Eng. Sci.* **2015**, 123, 170. [doi:10.1016/j.ces.2014.11.030](https://doi.org/10.1016/j.ces.2014.11.030) | Closest existing methodology. Our paper extends this with hierarchy + composite likelihood + censored data + facility effects. CES is one possible target venue. |
| Galagali, Marzouk. *J. R. Soc. Interface* **2019**, 16, 20180766. [doi:10.1098/rsif.2018.0766](https://doi.org/10.1098/rsif.2018.0766) | RJMCMC scaling via network topology. JRSI is a possible cross-disciplinary target. |
| Kennedy, O'Hagan. *J. R. Stat. Soc. B* **2001**, 63, 425. [doi:10.1111/1467-9868.00294](https://doi.org/10.1111/1467-9868.00294) | Foundational discrepancy framework that we extend. Cited heavily. Not a methodology comparator in itself — it is the floor. |
| Brynjarsdóttir, O'Hagan. *Inverse Problems* **2014**, 30, 114007. [doi:10.1088/0266-5611/30/11/114007](https://doi.org/10.1088/0266-5611/30/11/114007) | Counter-example to naive GP discrepancy. Central to our Theorem 3'. |
| Hegde, Li, Oreluk, Packard, Frenklach. B2BDC. *Combust. Flame* **2018**, 196, 509. | Deterministic UQ comparator. We position HBMAE as the Bayesian complement that absorbs many of the same ideas (consistency check, feasibility set) in a different framework. |

### §1.3 Application precedents in molten-salt radiation chemistry

| Paper | Venue | Relation to ours |
|---|---|---|
| Iwamatsu, Horne, …, Wishart. Zn²⁺ kinetics. *PCCP* **2022**, 24, 25088. [doi:10.1039/D2CP01194H](https://doi.org/10.1039/D2CP01194H) | The primary multi-T dataset we calibrate to. PCCP-style: standard "intro / methods / results / discussion" for a physical-chemistry experimental paper. **Not the right format for OUR paper** — they are reporting measurements; we are reporting a UQ framework. |
| Iwamatsu, Horne, …, Wishart. Cr(II)/Cr(III) kinetics. *PCCP* **2026**, 28, 2061. [doi:10.1039/D4CP04190A](https://doi.org/10.1039/D4CP04190A) | Same as above; second data source. |
| Conrad, Cook et al. Iodide impurity. *PCCP* **2023**. [doi:10.1039/D3CP01477K](https://doi.org/10.1039/D3CP01477K) | Same; third data source. |
| Phillips et al. INL/RPT-22-66727 (NaCl-UCl₃ NULL). [OSTI 1874817](https://www.osti.gov/biblio/1874817) | The NULL benchmark we use. INL report format. |
| Davis et al. F₂ G-values. *Nucl. Sci. Eng.* **2022**, 197(4), 633. | An NSE precedent for radiation-chemistry methodology papers (not Bayesian, but published in NSE — useful venue signal). |

### §1.4 Recent UQ-in-MSR papers (venue-scope check)

| Paper | Venue | Relevance |
|---|---|---|
| BEPU 2024 special issue: Bayesian calibration of fiber-optic temperature sensing. *NSE* **2024**. [tandfonline](https://www.tandfonline.com/doi/full/10.1080/00295639.2025.2528506) | Shows NSE actively publishes Bayesian-calibration methodology in nuclear context. **Strongest venue signal.** |
| Nuclear data uncertainty propagation for MSFR design. *NSE* **2023**, 197(12). | Same; methodology articles welcome in NSE. |
| Multi-physics + UQ for MSR design via DAKOTA. *Applied Sciences* **2024**, 14, 7615. | Open-access methodology venue; less rigorous than NSE but high turnaround. |

**Venue signal**: *Nuclear Science and Engineering* (NSE) is the strongest single target. They have a BEPU-special-issue track that explicitly invites Bayesian-calibration methodology applied to nuclear systems. Open-access option available.

---

## §2. Proposed structure

The paper is a **methodology + application + pedagogy** hybrid. Following Schnoerr et al.'s
template, scaled to the size of our framework and with explicit pedagogical features.

### §2.1 Single-paper vs. two-paper series

**Option A — single comprehensive article** (~30-35 pp double-spaced):
- All five HBMAE components in one paper.
- All four validation case studies (Tier 1 → Tier 4) in one paper.
- Full theorem treatment in appendix.
- Strong cross-references but one self-contained narrative.
- **Risk**: long enough that reviewers may push back on scope creep.
- **Reward**: a single high-citation reference for the whole framework.

**Option B — two-paper series**:
- Paper I "Methodology": HBMAE framework + theorems + controlled simulation (Caveat 3
  contraction study). ~20 pp. Target: NSE methodology special issue OR Chemical Engineering Science.
- Paper II "Application to LiCl-KCl radiation chemistry": Tier 2 calibration of Cr,
  Tier 3 multi-paper Zn, Tier 4 Phillips NULL + integrated MCMC. ~20 pp. Target: NSE
  or PCCP.
- **Risk**: review burden 2×; some reviewers may want to see application in the
  methodology paper.
- **Reward**: cleaner pedagogical focus per paper.

**Recommendation**: **Option A** for the first submission. The unified narrative is
stronger for the user's "pedagogical for early masters students" intent — students
benefit from seeing how the methodology and application connect in one place. The
length concern (35 pp) is within NSE's typical scope.

### §2.2 Detailed outline

The following section structure is for Option A. Page counts are double-spaced
estimates (≈ 350 words/page); 30-35 pages = 10500-12250 words.

**Title (working)**: "Hierarchical Bayesian Mechanism-Adequacy Estimation for Molten
Salt Radiolysis Networks: A Practical Framework with Validation Against the
Pulse-Radiolysis Literature"

**Abstract** (~200 words; structured): Problem → Method → Theorems → Validation →
Headline result.

**§1. Introduction** (~2-3 pp)
- 1.1 Molten salt reactors and the radiolysis bottleneck.
- 1.2 The state of the kinetic-network UQ problem: what's been done, what's missing.
- 1.3 Why molten-salt radiolysis exposes the gaps: time-scale heterogeneity, multi-
  paper data, censored benchmarks, cross-salt mechanism transfer.
- 1.4 Contribution and roadmap.

**§2. Background — building the reader up from first principles** (~6-7 pp)
- 2.1 A worked example of a chloride radiolysis network (3-4 reactions; what the ODE
  looks like; what an Arrhenius parameter is). *Tutorial box: "Reading a Stoichiometric Equation."*
- 2.2 Bayesian inference in 1.5 pages: prior → likelihood → posterior. Single-reaction
  Arrhenius worked example. *Tutorial box: "Choosing a Prior — the Three Choices."*
- 2.3 Model adequacy: definition, why it matters, the Kennedy-O'Hagan formulation.
  *Tutorial box: "When the Model is Wrong."*
- 2.4 Heterogeneous data: what makes radiation chemistry hard. Transient time-series,
  scalar rates, censored bounds, cross-paper drift.

**§3. State of the art and the gaps it leaves open** (~3 pp)
- 3.1 Galagali–Marzouk RJMCMC (with simple schematic diagram of the spike-and-slab prior).
- 3.2 Kennedy–O'Hagan discrepancy and the Brynjarsdóttir caveat.
- 3.3 Frenklach B2BDC: the deterministic alternative.
- 3.4 The five gaps for radiation chemistry (numbered, mapped to the five HBMAE components).

**§4. The HBMAE framework** (~7-8 pp; the central methodological section)
- 4.1 Notation (single table; reused throughout).
- 4.2 The five-layer architecture: schematic figure plus prose.
- 4.3 Constrained topology priors: feasibility set Γ, two-flip connectivity. *Tutorial box: "How do we know which networks are physically possible?"*
- 4.4 Hierarchical Arrhenius: physical motivation (solvent-shift bound from Marcus theory),
  hyperprior, K-scaling intuition.
- 4.5 Composite likelihood with Godambe weights: heuristic from "balancing modalities
  with different signal-to-noise."
- 4.6 Parametric discrepancy: the bounded form and why it avoids the Brynjarsdóttir trap.
- 4.7 Facility-effect terms and the gauge-fixing requirement.
- 4.8 The unified posterior: equation displayed once; product structure highlighted.
- 4.9 Algorithm pseudocode (Algorithm 1, with annotations for the reader).
- 4.10 Software stack used (PyMC vs Stan vs emcee/Julia DiffEqBayes — concrete recommendation).

**§5. Theorems and their plain-language interpretations** (~4-5 pp)
- For each of the 6 theorems + 4 propositions: a one-paragraph statement, an
  intuition box, and a single representative-line proof sketch. **Full proofs deferred
  to Appendix A.**
- Order of presentation: ergodicity → consistency → bias → composite → censored → identifiability.

**§6. Worked Example I: a one-reaction, one-salt calibration** (~3 pp)
- Cr²⁺ + e_s⁻ in LiCl-KCl with the Iwamatsu transients only.
- Step-by-step: write the ODE, set the priors, run NUTS, plot the posterior, interpret
  the credible interval. **The simplest meaningful HBMAE application.**

**§7. Worked Example II: cross-paper resolution** (~3 pp)
- Pikaev (1982) vs Iwamatsu (2022) for Zn²⁺ + e_s⁻.
- The χ² consistency test fails (Tier 1 finding).
- Adding the facility-effect term b^(Pikaev) resolves the inconsistency.
- Pedagogical point: this is where the framework earns its keep.

**§8. Validation case studies** (~6 pp)
- 8.1 Tier 2 Cr calibration (9 transients, 14 parameters, full posterior).
- 8.2 Tier 3 multi-paper Zn (5-method comparison vs SOTA).
- 8.3 Tier 3 Phillips NULL (censored Bayes factor 10¹³⁸).
- 8.4 Tier 4 Theorem 2 K⁺¹ contraction (controlled simulation; figure).
- 8.5 Tier 4 integrated end-to-end (24-parameter joint posterior; table).
- 8.6 (optional) Tier 4 slow-manifold rigour and U-redox sensitivity.

**§9. Comparison with state-of-the-art methods** (~3 pp)
- Table: HBMAE vs M₀–M₄ on the multi-paper Zn case.
- Discussion: when each method wins; when each method fails.
- The honest statement of HBMAE's strict advantages (small, defensible) and what it
  does NOT provide.

**§10. Discussion** (~2-3 pp)
- 10.1 What the radiation chemistry community gains.
- 10.2 What other communities (combustion, biochemistry, environmental) can learn.
- 10.3 Honest limitations:
  - Requires informative priors — not a free lunch.
  - Computational cost scales with ODE solve × parameter dimension × MCMC steps.
  - Discrepancy-class specification still depends on physics knowledge.
  - Small-K identifiability of the facility hierarchy.
- 10.4 Future work: extension to multi-modality combustion data; FLiBe + actinides;
  Bayesian optimal experimental design.

**§11. Conclusion** (~1 pp)
- One paragraph each on:
  - The framework
  - The validation
  - The comparison
  - The pedagogical takeaway

**Appendices** (~10 pp; not in page count above)
- **A. Full proofs of Theorems 1–6 and Propositions 2–5.** With cross-references back to the main text.
- **B. Computational details.** Software versions, ODE solver tolerances, MCMC convergence diagnostics, computational cost per case.
- **C. Code listing pointers.** Reference to GitHub repository with full reproducible scripts. Each Tier (1–4) is a notebook plus a script.
- **D. Glossary of terms.** ~30 terms defined briefly; cross-referenced from main text.
- **E. Suggested exercises for students.** ~10 exercises ranging in difficulty from
  "compute the posterior of a 2-parameter Arrhenius given 5 data points" to "extend
  HBMAE to a new metal in a new salt".

**Supplementary information**
- All raw digitized CSVs (already in [validation/](../validation/)).
- All Python scripts.
- The four narrative reports ([TIER1_NARRATIVE.md](../validation/TIER1_NARRATIVE.md),
  [HBMAE_THEOREMS_TIGHTENED.md](../validation/HBMAE_THEOREMS_TIGHTENED.md),
  [TIER3_EXTENSIONS_REPORT.md](../validation/TIER3_EXTENSIONS_REPORT.md),
  [TIER4_DEEP_CAVEATS.md](../validation/TIER4_DEEP_CAVEATS.md)).
- Posterior chains (`.npy` files; with a script to reproduce visualizations).

### §2.3 Pedagogical features to include throughout

These distinguish a "pedagogical" methodology paper from a standard one:

| Feature | Implementation |
|---|---|
| **Plain-language summary** at start of each section | One paragraph in italics or a styled box, before the technical content. |
| **Tutorial boxes** | ~8-12 boxes throughout, each ~150 words, explaining one concept (e.g. "What is the posterior predictive?"). |
| **Numbered footnotes** for asides early MS students might skip | Standard footnote treatment. |
| **Glossary** in Appendix D | ~30 terms with one-line definitions; cross-referenced from main text. |
| **"Common pitfalls" callouts** | Boxed asides at points where novices typically go wrong (e.g. "Why setting σ_prior → ∞ is not 'uninformative' — it's improper"). |
| **Numerical examples are reproducible** | All data + code are SI; reader can recompute every number in every table. |
| **Figure annotations** | Figure captions tell the story without requiring the reader to find the body text. |
| **Citation density** | Higher for pedagogical references (textbooks: van der Vaart, Gelman BDA, Bishop ML); lower for cutting-edge methodology in the body, deferred to references list. |

### §2.4 Figures plan (target: ~10-12 main figures)

1. Schematic of HBMAE architecture (the 5 layers).
2. Worked Example I: posterior corner on (log A, Eₐ) for Cr²⁺ + e_s⁻.
3. Tier 1 profile-likelihood plot showing identifiable / ridge-degenerate / unidentifiable
   classification.
4. Tier 1 consistency check: standardized residuals across Pikaev + Iwamatsu.
5. Schematic illustrating how the facility-effect term b^(p) resolves the inconsistency.
6. Tier 2 posterior predictive overlay on the 9 Iwamatsu Cr transients.
7. Tier 3 method comparison bars (in-sample, held-out, bias).
8. Tier 3 Phillips NULL censored Bayes factor visualisation.
9. Tier 4 contraction simulation plot.
10. Tier 4 integrated MCMC corner plot.
11. (optional) Tier 4 slow-manifold trajectory.
12. (optional) HBMAE workflow flowchart for the reader.

### §2.5 Tables plan (target: ~6-8 main tables)

1. Notation summary (§2).
2. Comparison of methodologies (§3): SOTA capabilities matrix.
3. Identifiability classification of reactions (§6 or §8.1).
4. Tier 2 posterior summary (Cr).
5. Tier 3 5-method comparison (Zn).
6. Tier 3 Phillips NULL Bayes-factor table.
7. Tier 4 K-scaling exponents.
8. Tier 4 integrated posterior summary (24 parameters).

---

## §3. Venue selection — recommended ordering

| Choice | Journal | Target audience | Fit | Risk |
|---|---|---|---|---|
| **1st choice** | *Nuclear Science and Engineering* (NSE) | Nuclear engineers, mixed MS/PhD. | Strong: BEPU special-issue track explicitly invites methodology; Davis 2022 F₂ paper sets a precedent for radiation-chemistry papers in NSE. | Acceptance rate ~30-40% for methodology submissions; review cycle 6-9 months. |
| 2nd choice | *Reliability Engineering & System Safety* | Broader UQ community, mixed MS/PhD. | Strong: explicitly publishes Bayesian-calibration methodology. Pedagogical articles common. | Less nuclear-domain visibility. |
| 3rd choice | *Chemical Engineering Science* | Chem-eng grads. | Strong: where Galagali-Marzouk 2015 appeared. Methodology + worked example fits. | Less nuclear-domain visibility. |
| 4th choice | *Journal of Computational Chemistry* / *AIChE Journal* | Comp-chem / chem-eng. | Mixed: more code/algorithm orientation. | Less nuclear-domain visibility. |
| Backup | arXiv preprint + journal | Open access | n/a | Always do this. |

**Recommendation**: arXiv preprint immediately on submission; primary target NSE,
secondary target Reliability Engineering & System Safety.

---

## §4. Timeline (suggested)

Assuming the user has ~1 month/week part-time on this:

| Week | Deliverable | Owner |
|---|---|---|
| 1 | Lock target journal + scope (single vs. two-paper). | User decision. |
| 2-3 | Write §1 Introduction + §2 Background from existing materials. ~10-12 pp. | Drafting work. |
| 4 | Write §3-5 Methodology + Theorems (mostly pulling from HBMAE_THEOREMS_TIGHTENED). | Drafting. |
| 5-6 | Write §6-7 Worked Examples I+II. New figures needed. | Drafting + figures. |
| 7 | Write §8 Validation case studies. Pull tables from validation/ directly. | Drafting. |
| 8 | Write §9-11 Comparison + Discussion + Conclusion. | Drafting. |
| 9 | Write appendices A-E. Add glossary + exercises. | Drafting. |
| 10 | Internal review pass. Check pedagogical features are in place. | Self-review or co-author. |
| 11 | External review by 2-3 trusted colleagues. | INL + BNL contacts? |
| 12 | Revise; submit to arXiv + NSE simultaneously. | Final pass. |

Total: ~12 weeks at 8-10 hours/week. Plausible if the user can commit ~80 hours.

---

## §4b. Locked decisions (user-selected 2026-05-28)

| Question | Decision |
|---|---|
| Q1 Scope | **Single comprehensive article** (~40 pp double-spaced after pedagogical depth bump from Q5). |
| Q2 Target journal | **Nuclear Science and Engineering (NSE)** — primary. Reliability Engineering & System Safety as backup. arXiv preprint on submission day. |
| Q3 Co-authors | **Solo with citation** — INL-only authorship. Cite Iwamatsu/Horne/Wishart cluster heavily; do not engage as co-authors. |
| Q5 Pedagogical depth | **No Bayesian background assumed** — full primer on Bayesian thinking required. ~3-4 extra pages in §2; tutorial boxes for prior, likelihood, posterior, MCMC. |

These decisions change the page budget and the §2 structure as follows.

**Revised page budget** (double-spaced, ~350 words/page):

| Section | Pages | Notes |
|---|---|---|
| §1 Introduction | 2-3 | unchanged |
| §2 Background (expanded for Bayesian primer) | **10-11** | up from 6-7 |
| §3 SOTA + gaps | 3 | unchanged |
| §4 HBMAE framework | 7-8 | unchanged |
| §5 Theorems + plain-language interpretations | 4-5 | unchanged |
| §6 Worked Example I | 3 | unchanged |
| §7 Worked Example II | 3 | unchanged |
| §8 Validation case studies | 6 | unchanged |
| §9 SOTA comparison | 3 | unchanged |
| §10 Discussion | 2-3 | unchanged |
| §11 Conclusion | 1 | unchanged |
| **Total main text** | **~44 pp** | up from ~32 pp |
| Appendices A-E | ~12 pp | unchanged |

**Implication**: 44 pp is on the upper edge of NSE's "regular article" length budget
(typical 30-40 pp double-spaced; some longer methodology articles do appear). The
specific BEPU/Bayesian-calibration track papers I surveyed all sit around 30-40 pp.
Strategies if reviewers push back on length:
- Move §5 theorem proofs entirely to Appendix A (saves ~2 pp).
- Move Worked Example II into Validation Case Studies (saves ~1 pp).
- Trim §8 from 6 pp to 4 pp by combining Tier 3 sub-cases.
- Trim §2.1 from worked-example-style to definition-style if Bayesian primer dominates.

I'll write to 44 pp in the first pass and budget-cut if reviewers request.

**Revised §2 structure** (to accommodate no-Bayesian-background readers):

- 2.1 The radiolysis problem in molten chloride salts (1-1.5 pp; physical picture, what
  is observed, why it matters for MSRs).
- 2.2 Kinetic networks as ODEs (1.5 pp; worked example: 3-reaction chloride kernel,
  Arrhenius temperature dependence, what the trajectory looks like). *Tutorial box:
  "Reading a Stoichiometric Equation."*
- 2.3 **Probabilistic thinking 101** (2 pp; new): what a probability distribution is in the
  context of "I don't know the rate constant, but I know it's around 1.7×10¹³ with some
  spread." Discrete vs continuous; PDF and CDF; basic operations. *Tutorial box: "Why
  Bayesian and not Frequentist?"*
- 2.4 **Bayes' rule and its three pieces** (2 pp; new): prior, likelihood, posterior.
  Worked Arrhenius example with two data points. *Tutorial box: "Choosing a Prior — the
  Three Choices."* *Tutorial box: "Common Pitfall — improper 'uninformative' priors."*
- 2.5 **Posterior computation** (1.5 pp; new): why MCMC; what walkers and chains are;
  what "convergence" means. *Tutorial box: "Reading a Corner Plot."*
- 2.6 Model adequacy: what it means; the Kennedy-O'Hagan formulation; why a "perfect-
  model" assumption fails. *Tutorial box: "When the Model is Wrong."*
- 2.7 Heterogeneous radiation-chemistry data: time-resolved transients, scalar rates,
  censored bounds, cross-paper drift. (~1.5 pp).

This brings §2 to ~10-11 pages with 6 tutorial boxes integrated.

## §5. Open strategic questions for the user

Before drafting begins, four decisions to lock down:

**Q1. Single-paper vs. two-paper series?**
- My recommendation: single paper (Option A) for the pedagogical intent.
- But two-paper has the cleaner methodology / application separation.

**Q2. Target journal?**
- My recommendation: NSE primary, RESS secondary.
- Alternative: PCCP for the chemistry-focused readership (but PCCP is more
  "data paper" than "methodology paper" in style).

**Q3. Co-authors?**
- The framework was developed here; the data is from Iwamatsu/Horne/Wishart cluster.
  Should we approach them for co-authorship?
- Pros: stronger acceptance odds, additional radiation-chemistry credibility,
  potential to use their unpublished data.
- Cons: longer timeline; potential scope disagreements.

**Q4. Code / data release strategy?**
- Recommendation: GitHub repo with all scripts + raw CSVs + posterior chains, made
  public on submission day. Cite the repo with a Zenodo DOI.
- Risk: limited; INL has clearance procedures for code release that need ~2 weeks lead time.

**Q5. Pedagogical depth — how deep into the basics do we go?**
- Option (a) "Self-contained for an MS student who has had one course in Bayesian
  statistics": defines posterior, likelihood, MCMC, but assumes Arrhenius/ODE
  background.
- Option (b) "Self-contained for an MS student who has had no Bayesian course":
  more textbook-style introduction to Bayesian thinking; longer §2.
- My recommendation: option (a) — depth (b) would push the article into
  textbook territory and unbalance the methodology/application focus.

---

## §6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Length exceeds journal limits (~30 pp for NSE) | Move detailed proofs and computational details to appendices. Possibly drop §8.6 if not needed. |
| Reviewers push back on "methodology + applied combo" being out of journal scope | Pre-submit query to NSE editorial; backup plan to split into Option B if rejected. |
| Iwamatsu/Horne/Wishart team objects to using their digitized data | Reach out early; they will likely co-author rather than object. |
| INL clearance for code release is slow | Start clearance process Week 1; release on submission day. |
| Pedagogical level slips toward textbook OR toward expert-only | Have one MS student + one PhD reviewer read drafts. |
| The framework itself fails review for "yet another Bayesian hierarchical model" | Strong emphasis on the FIVE-COMPONENT NOVELTY: no other framework simultaneously has constraint + hierarchy + composite + parametric discrepancy + facility effects + the rigorous theorems for each. |

---

## §7. What I need from you to begin drafting

In order of priority:

1. **Decision on Q1 (single vs. two papers)**. I lean single; happy to do whichever.
2. **Decision on Q2 (target journal)**. NSE is my recommendation.
3. **Decision on Q3 (co-authors)**. INL/BNL relationship: do we reach out?
4. **Decision on Q5 (pedagogical depth)**. Option (a) is my recommendation.

Code release strategy (Q4) can be deferred; not blocking.

Once these are locked, the drafting plan above becomes the concrete work-breakdown.
