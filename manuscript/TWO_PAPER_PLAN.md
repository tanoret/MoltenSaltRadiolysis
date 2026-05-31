# Two-paper restructure — execution plan

Date: 2026-05-28

Authoring decision: split the current 54-page combined draft into two independent journal articles.
Both papers must stand alone and pass peer review on their own merits.

---

## Tone and rigor calibration

Confirmed from user feedback:
- **Keep the writing style** (tone, accessibility of prose, clarity of figures).
- **Remove** tutorial boxes, exercises, "Reading a stoichiometric equation"-style
  pedagogical asides. These belong in a textbook, not a journal article.
- **Add** rigor: formal definitions, full proofs, all assumptions stated explicitly,
  numerical diagnostics for every MCMC, sensitivity analyses for every prior,
  cross-validation for every claim.

Each paper must read like a rigorous methods/applications paper, not a tutorial.

---

## Phase A — Close the remaining data gap (priority: critical, ~1 hour)

Two new papers landed in `validation/papers/`:

| Paper | DOI | Role |
|---|---|---|
| Makarov, Zhukova, Pikaev, Spitsyn 1982 (`BF00949993.pdf`) | [10.1007/BF00949993](https://doi.org/10.1007/BF00949993) | **Source of the ε(e_s⁻) ≈ 8000 M⁻¹ cm⁻¹** that Iwamatsu 2022 cites. Closes Tier A1 gap. Also reports F₂⁻/Br₂⁻/I₂⁻/Cl₂⁻ molar absorptivities and disproportionation rates. |
| Iwamatsu, Horne et al. 2022 (`d2cp01194h.pdf`) | [10.1039/D2CP01194H](https://doi.org/10.1039/D2CP01194H) | Full text of the Zn/Cl₂⁻ paper. Cross-check the existing extracted Arrhenius parameters; verify ε values used; extract any missed numerical entries. |

Tasks:
A1. Digitize Makarov 1982: extract every reported molar absorptivity ε(λ_max, T), every k value, every G-value. Build `validation/oxidants_halide_melts/makarov_1982_bull/data/*.csv`.
A2. Re-verify Iwamatsu 2022 PCCP digitization against the full PDF. Cross-check Arrhenius parameters; extract ε values; capture all transient-spectrum and pseudo-first-order tabulated entries.
A3. Update `CHLORIDE_DATA_AUDIT.md` to mark Tier A1 closed; flag any Tier B gaps that now matter.
A4. Convert the scale-free likelihood mode (currently used in Tier 2/4 because ε was missing) to an absolute-concentration likelihood using the calibrated ε(λ, T) prior. Re-run the integrated HBMAE MCMC.

---

## Phase B — Paper 1: methodology

**Working title** (tentative): *"Hierarchical Bayesian Mechanism-Adequacy Estimation: A Multi-Modality, Multi-Host Framework for Calibrating Chemical Kinetic Networks Against Heterogeneous Experimental Data"*

**Candidate journals** (decision required from user):
- *SIAM/ASA Journal on Uncertainty Quantification* — preferred for the rigorous theorem treatment
- *Journal of Computational Physics* — broader audience, applied UQ focus
- *Annals of Applied Statistics* — most theoretical, would lean further on the theorems

**Structure**:

| § | Title | Notes |
|---|---|---|
| 1 | Introduction | Drop the textbook framing. Lead with the calibration problem in chemical kinetics; describe the gap; give the contributions list (5 layers + 6 theorems + the multi-modality composite likelihood). Add the comparison table of methods *up front* (currently in §3.6); state the unique advantages of HBMAE. |
| 2 | Background and notation | Compress current §2.1–2.6. Drop the "tutorial" boxes. Keep the stiff-ODE setup, the inference problem, and the Bayes/MCMC formalism in 4–5 pages. Add a *Definitions* subsection (with formal definitions of every object used later: $\mathcal{G}$, $\boldsymbol\theta$, $\boldsymbol\gamma$, $\eta^{(s)}$, $b^{(p)}$, the likelihood modalities). |
| 3 | State-of-the-art comparator methods | Expand. For each comparator (Galagali–Marzouk RJMCMC, Kennedy–O'Hagan + GP, B2BDC, PSIS-LOO, Frenklach-style data-collaboration) give: (a) the mathematical formulation, (b) the asymptotic guarantees they provide, (c) the failure modes they exhibit when applied to multi-host molten-salt data. Conclude with a formal comparison table. |
| 4 | The HBMAE framework | The core methods section. Expand each of the five layers from the current ~2 pages to ~4 pages with: (i) the formal model statement, (ii) the role and motivation, (iii) the prior specification and its data-availability dependence, (iv) the identifiability properties. Add an explicit *algorithm box* showing the full sampling step with all conditionals. |
| 5 | Theorems and proofs | Full statements with all assumptions inline; full proofs in the body (not appendix). For each theorem: state the assumption set $A_k$; state the conclusion $C_k$; prove via a stated argument structure (often: contraction + a Bayes-factor or BvM step); state the consequence in plain language at the end. Add numerical illustrations (small simulation studies where the rate predicted by the theorem is measured empirically). |
| 6 | Worked Example I — single-reaction calibration on Cr | Keep the Cr Tier 2 example but expand with: (a) prior-sensitivity analysis (vary the prior on each Arrhenius parameter by ±1σ; show the posterior is robust); (b) full convergence diagnostics ($\hat{R}$, $N_{\mathrm{eff}}$, traces, autocorrelations); (c) posterior-predictive checks per trace. |
| 7 | Worked Example II — multi-paper Zn cross-validation | Expand similarly with prior-sensitivity, diagnostics, and posterior-predictive. Add formal verification that the facility-effect resolution is identifiable (not just plausible). |
| 8 | Validation | Tiers 1–4 expanded. Add (a) Tier 1: full identifiability analysis with profile likelihood for every parameter; (b) Tier 2: full diagnostics + prior-sensitivity; (c) Tier 3: Godambe sandwich variance verification + cross-paper consistency χ² tests; (d) Tier 4: full 24-parameter MCMC with posterior decomposition, ELPD-LOO contributions per modality, and posterior-predictive checks. |
| 9 | Comparison with comparator methods | Expand. Re-run the M₀–M₄ comparison with full diagnostics. Add a *synthetic* benchmark with known ground truth where each method's bias is measurable. Add a formal comparison criterion table (ELPD, posterior-predictive RMSE, parameter recovery RMSE, prior-sensitivity score, computational cost). |
| 10 | Discussion | Where HBMAE wins, where it loses. The Brynjarsdóttir–O'Hagan tension. Comparison to recent advances. |
| 11 | Conclusion | Short. Recap, point to Paper 2 for the application. |
| App. A | Detailed proofs of any technical lemmas not in the main text | |
| App. B | Computational implementation details (algorithm pseudocode, MCMC tuning, software stack) | |
| App. C | Synthetic-benchmark results, including the M₂ vs M₃ contraction-rate simulation | |

**Estimated length**: 35–45 pages including appendices.

**Sections to drop entirely** (relocate to Paper 2):
- Existing §10 (predictive-for-any-element)
- The fluoride F₂-kernel material
- The meta-hierarchical chemistry-feature regression
- The PSIS-LOO network selection (it stays only as a *comparator* in §9, not as an example)
- The train/val split figure (move to Paper 2 as a tier-2 validation panel)

---

## Phase C — Paper 2: application

**Working title** (tentative): *"A Calibrated Multi-Salt Radiolysis Model for Molten Chloride and Fluoride Reactor Fuels"*

**Candidate journals** (decision required from user):
- *Journal of Nuclear Materials* — strong fit for radiolysis + reactor materials
- *Nuclear Science and Engineering* — close fit, slightly more reactor-engineering
- *Annals of Nuclear Energy* — broad-scope nuclear
- *Journal of Nuclear Engineering and Radiation Science* — application focus

**Structure**:

| § | Title | Notes |
|---|---|---|
| 1 | Introduction | The MSR radiolysis problem. Chloride and fluoride salts as primary candidates. Practical questions (Cl₂ buildup in MCFR, F₂ buildup in FLiBe reactors, redox state evolution, corrosion coupling). Why a unified calibrated model is needed. Brief reference to Paper 1 for the framework, but the paper must stand alone. |
| 2 | The experimental database | Detailed inventory: every digitized paper, every CSV, every metal and host. Composition tables. Temperature ranges. Modality classification (transient time-series, scalar rate, censored detection, Arrhenius pair, G-value). Provide a *one-figure* summary of the data landscape. |
| 3 | Chloride radiolysis network | The mechanistic kernel: R1–R8 (radical formation), R9–R12 (Cl₂⁻ chemistry), R13–R20 (metal redox: Cr, Zn, Nd, Cf, U). Justify each reaction from the literature. Give the full ODE system. |
| 4 | Fluoride radiolysis network | The static F₂-production kernel (Davis G + Toth/Felker Ea). Comparison with the limited transient data (Akiyama 1994 referenced but undigitized — flag as future work). |
| 5 | Calibration via HBMAE | Brief recap of HBMAE pointing to Paper 1. Per-metal calibration results in a single comprehensive table (Arrhenius parameters with 90% CIs, host effects $\eta^{(s)}$, facility offsets $b^{(p)}$). Convergence diagnostics summary. |
| 6 | Meta-hierarchical chemistry-feature layer | The "predictive-for-any-element" extension. Chemistry features (z, r_ion, χ, log(z/r)). Posterior regression coefficients. LOMO cross-validation table. The interpretation of when the layer extrapolates well vs poorly. |
| 7 | **Validation against all experimental data** (the centerpiece) | One figure per paper, model prediction overlaid on every data point. Tabulated residuals. χ² statistics per dataset. A meta-figure showing the model's residuals across all 33+ observations in a single panel (residual vs predicted log_10 k, separated by metal/host). This is the section the user explicitly asked to "plot ALL the experimental data". |
| 8 | Predictive application to operational MSR conditions | (a) NaCl-UCl₃ MCFR (Phillips conditions, 31 MGy dose, 75–600 °C): predicted Cl₂ inventory over lifetime with full uncertainty propagation. (b) FLiBe-UF₄ MSR: predicted F₂ inventory under operational conditions. (c) Comparison to design-basis safety limits. (d) Identification of the dominant model-uncertainty contributors and recommended experiments to tighten them. |
| 9 | Discussion | Model strengths, limitations, recommended experiments. Honest assessment of where the model is and is not predictive. |
| 10 | Conclusion | |
| App. A | Full network species/reaction tables | |
| App. B | All extracted data CSVs cross-referenced to source papers | |
| App. C | MCMC diagnostics for the production calibration | |
| App. D | Operational-MSR prediction details and sensitivity analyses | |

**Estimated length**: 30–40 pages including appendices.

---

## Phase D — Cross-references and consistency

- Paper 2 cites Paper 1 in §1 and §5 for the HBMAE framework.
- Both papers share the literature bibliography but have non-overlapping content.
- Worked Example II (Zn cross-paper) stays in Paper 1 because the *method* is the point; the Zn results are reproduced briefly in Paper 2's §5 calibration table.
- The Phillips NULL benchmark appears in *both* papers but for different purposes:
  in Paper 1 as a validation case (Tier 3) demonstrating the censored Bayes factor;
  in Paper 2 as the calibration anchor for the U-redox sub-kernel and as the dominant
  driver of the operational MCFR Cl₂ prediction.

---

## Execution sequence

I recommend:

1. **Phase A** (1 hr): close the data gap with Makarov + Iwamatsu re-verification. This is critical for *both* papers since absolute concentrations now become available.
2. **Phase B** (most work): execute Paper 1 first. Strip pedagogical material; expand theory.
3. **Phase C**: execute Paper 2 using Paper 1's calibrated outputs.
4. **Phase D**: final cross-reference pass; recompile both PDFs.

Alternative sequence: do A first, then B and C in parallel (would let me alternate
between methods/application without keeping both contexts loaded). Slightly higher
risk of inconsistency between the two papers.

---

## Resolved decisions (2026-05-28)

1. **Paper 1 venue**: SIAM/ASA Journal on Uncertainty Quantification.
2. **Paper 2 venue**: Journal of Nuclear Materials.
3. **Paper 2 §8 scope**: Full operational-MSR application — NaCl-UCl₃ MCFR Cl₂ inventory
   over reactor lifetime + FLiBe-UF₄ F₂ inventory + comparison to design-basis safety
   limits.
4. **Execution sequence**: Paper 1 first, then Paper 2 (lower inconsistency risk).
5. **Authorship**: solo (M. Tanore-Tamales, INL) unless the user adds co-authors later.
