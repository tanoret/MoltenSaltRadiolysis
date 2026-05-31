# Referee Report — Paper 2 (Multi-salt radiolysis model)

**Manuscript:** *A Calibrated Multi-Salt Radiolysis Model for Molten Chloride and Fluoride Reactor Fuels: Application of Hierarchical Bayesian Mechanism-Adequacy Estimation to Operational MSR Predictions* (M. Tanore-Tamales), `article_paper2.tex`, dated 2026-05.

**Target journal:** *Journal of Nuclear Materials*.

## Recommendation
**Major revisions.**

The manuscript is potentially a strong fit for *JNM*, but the headline operational prediction (cover-gas P(Cl₂) ≈ 2.2 × 10⁻⁸ Pa, ten orders of magnitude below the engineering screening threshold) rests on a single censored observation and on an unmeasured U(III)/U(IV) rate constant whose prior dominates the posterior. The chemistry of §3 is incomplete in ways that matter for the MCFR claim, the corpus inventory has internal counting inconsistencies (14 sources in Table 1, 23 claimed in the body, 11+ in the prompt), and the abstract over-states what the data can support. None of these problems are fatal — they require honest re-framing of the MCFR result as a posterior *upper bound* propagation rather than a "prediction", a more complete chloride network, and tightening of the bookkeeping. With those changes the paper is publishable.

## Summary (≤200 words)

The paper assembles a digitized multi-modality corpus (pulse-radiolysis transients, scalar rate constants, derived Arrhenius pairs, censored detection limits, and steady-state G-values) and feeds it through the HBMAE framework of the companion paper (`Tanore2026PaperI`) to produce a 24-parameter integrated posterior for chloride and fluoride radiolysis kernels in molten salt reactor fuels. The novel application is two operational forward predictions: (i) cover-gas Cl₂ inventory in a 500-MW(th) NaCl-UCl₃ fast reactor over 60 years, predicted ten orders of magnitude below a 100-Pa screening threshold; and (ii) steady-state cover-gas F₂ in a FLiBe-UF₄ thermal loop at 31 [19, 44] Pa, within a factor of three of the same threshold. The headline contributions — the Pikaev-Iwamatsu inter-facility offset b ≈ −3.35, the Phillips-2022 censored Bayes factor ln K = 316.9 endorsing the U(III)/U(IV) sink, the calibrated Cr/Zn Arrhenius pairs, and the LOMO-validated chemistry-feature layer — are real and well-executed. The MCFR ten-orders-of-magnitude margin, however, is essentially a propagation of the Phillips upper bound and should be reported as such.

## Major Concerns

### 1. The MCFR ten-orders-of-magnitude margin is a prior-shaped artifact, not a prediction. (§8.1, lines 1205–1238)

The posterior on η ≡ chain-propagation efficiency is set by Phillips 2022 alone (a single 31-MGy gamma irradiation with no detection above ≈ 1000 ppm); the kernel R12 (Cl₂•⁻ + U(III) → U(IV)) has no measurement — its prior is "diffusion-limited centered, two-decade lognormal spread" (§3.4 R12 paragraph at line 553–563). The Phillips data set is one *one-sided* observation; it constrains η from above but does not measure it. Reporting a posterior median of 2.2 × 10⁻⁸ Pa and a 90% CI of [3.7 × 10⁻⁹, 8.5 × 10⁻⁸] Pa is therefore misleading — what the data actually support is "η ≤ 10⁻¹⁸ at 90% credibility under the stated kernel; the unobserved tail can be many orders of magnitude tighter or looser without changing the data likelihood." The "ten orders of magnitude below the limit" claim in the abstract and §1 (lines 87–90) propagates this artifact.

**Fix.** Reframe §8.1 as "Phillips-bounded screening" with the posterior interpreted as *the prior pushed forward under a censored constraint*. Show a sensitivity test where R12 is two decades slower than the prior median (the lower edge of the lognormal): the MCFR P(Cl₂) margin will likely collapse by ≥ 10×. State explicitly that the predicted median is a *Bayesian propagation of the U-buffer prior*, not an empirical extrapolation.

### 2. The chloride radiolysis network (§3, Table 3) is missing reactions that materially affect the MCFR claim.

Specifically:
- **Oxide-impurity scavenging.** Sims and Forsberg both emphasize that ppm-level O²⁻/OH⁻ contamination is the dominant Cl₂ sink in NaCl-UCl₃ at operational conditions. The current network has no oxide channel — UO₂Cl₂ formation, hydrolysis at the cover-gas interface (H₂O ingress + Cl₂ → HCl + ½ O₂), or the U(IV) + O²⁻ → UO₂²⁺ + 2 e⁻ pseudo-redox couple. If oxide scavenging is even partially active, η would be smaller still; if oxide scavenging *fails* (well-dried salt), η rises by orders of magnitude.
- **Redox-active fission products.** Cs⁺, Cs⁰ (volatile), I⁻, I⁰, Te, Sr²⁺/Sr⁰ are all redox-active on minute-to-hour timescales and would compete with U(III)/U(IV) for the Cl₂•⁻ pool. The paper invokes the chemistry-feature layer for these but does not include any of them in the operational ODE.
- **Cl•/Cl₂•⁻ disproportionation to ClO•** in the presence of residual oxide.
- **R10** (`Cl₂•⁻ + Cl⁻ → Cl₃⁻ + e⁻`) is marked "assumed slow" with no rate constant. Yet the reverse of R10 is essentially R4 + ionization; the asymmetry deserves justification.

**Fix.** Add a §3.5 "Reactions omitted from the kernel and their possible impact" that lists each (oxide, fission products, ClO•) with an order-of-magnitude estimate of the perturbation each would induce on the MCFR P(Cl₂) prediction. This is exactly the place where the U-buffer assumption needs to be challenged on the page.

### 3. Source-counting and corpus claims do not internally reconcile.

- Abstract line 64: "33 rate-constant observations across nine metal scavengers and five hosts, four fluoride G-value compositions" — so ≈ 38 observations.
- §2 line 248: "23 digitized literature sources".
- Table 1 (lines 296–309): 14 sources listed.
- §9.1 line 1428: "73 observations drawn from 13 sources within the larger 23-source digitized database."
- Reviewer prompt: "11+ digitized papers."

Either Table 1 is incomplete (the other 9 sources of the "23" total are not displayed) or the 23 figure is wrong. Either is correctable, but in the present state a careful reviewer cannot verify "every digitized observation" against the corpus. The discrepancy between 33 (abstract) and 73 (discussion) is presumably a per-rate-constant vs. per-(rate-constant + Arrhenius + transient-trace + G-value + censored) count, but this needs to be stated explicitly with a single audit table.

**Fix.** Either expand Table 1 to all 23 sources, or correct "23" to the actual table count (≈ 14). Add a single audit row in Appendix B that adds up obs by modality to the 73 cited in §9.1.

### 4. Internal inconsistencies in the operational scenarios.

- **NaCl-UCl₃ composition.** Table 2 (line 366) lists the Phillips host as "67-33 mol%" but the operational MCFR plant model (line 1180) uses "60-40 mol%". Real MCFR designs are typically nearer 60-40 (UCl₃-rich for breeding), so the 60-40 is plausible — but Table 2 and §8.1 disagree. The 40 mol% UCl₃ choice is also at the high end of solubility for UCl₃ in NaCl at 873 K and deserves a citation.
- **U(III) inventory.** Line 1184 states "Initial U(III) inventory: 40 mol% UCl₃ in the salt loop (∼ 2.5 × 10⁵ mol total U(III))." But UCl₃ contains U(III) — so saying U(III) = 40 mol% is dimensionally fine, but the calculation "2.5 × 10⁵ mol" needs unit checking against 135 t / 50 m³ / molar mass. A back-of-envelope: 135 t × 0.4 (mol-fraction UCl₃) × (1 / M_UCl₃,avg molar mass of salt ~ 100 g/mol) yields ~ 5 × 10⁵ mol of *salt* per mol fraction — the bookkeeping needs to show the work.
- **Dose rate.** Line 1187: "10 kGy/h … moderately lower than the Phillips 2022 13 kGy/h research-irradiation rate." This justification is hand-waving. A 500-MW(th) salt-fueled fast reactor with a 50 m³ salt loop has fission-power volumetric heating dominated by *gamma deposition in the salt*; the appropriate scaling argument (kerma factors, photon penetration depth in NaCl-UCl₃) would give a number in the 50–500 kGy/h range, not 10. Justify with a citation or a back-of-envelope.
- **Henry constant.** Line 1227 cites Sohal et al. 2010 for K_H(Cl₂) = 2 × 10⁻⁵ mol/(m³·Pa) in NaCl-UCl₃ at 600 °C. The Sohal report is a general property compilation; the value for Cl₂ in NaCl-UCl₃ specifically is to my knowledge *not* measured. State whether this is an inferred value (and from what reference reaction) or measured.

### 5. The chemistry-feature LOMO behavior (§6, Table 6) under-predicts the LEAF data by 2–4 orders of magnitude for held-out Cr/Nd/Cf, which is essentially the *entire LEAF dataset for those metals*.

A 2–4 log-unit miss with 0 of 5 coverage (Cr³⁺, Nd³⁺, Cf³⁺ rows in Table 6) is not a "domain of applicability" issue — it is a failure mode of the regression. The reported explanation ("Cr, Nd, Cf were measured at LEAF picosecond; the regression cannot resolve the facility offset for the held-out metal") is exactly correct, but the implication is bigger than the paper acknowledges: it means the chemistry-feature layer cannot be used to predict the kinetics of *any* fission-product element that would be measured at LEAF if such data existed. Practically, every relevant operational element (Cs, Sr, Te, Eu, Pu, Am, Cm) falls into the LEAF-extrapolation regime.

**Fix.** Move §6.2 (LOMO) limitation into a top-level "the chemistry-feature layer is currently descriptive, not predictive, for operational fission-product elements" caveat. Either fit a hierarchical facility-effect into the chemistry-feature regression (analog of b^Pikaev) or restrict the deployment claim to hosts and metals already in the corpus.

### 6. The Toth-Felker recombination kinetics may not transfer to operational FLiBe-UF₄.

Toth & Felker (1990) measured F₂ recombination on *MSRE-composition* salt (≈ LiF-BeF₂-ZrF₄-UF₄) at temperatures *below* 150 °C with a recombination mechanism dominated by surface chemistry at the metal walls. The paper applies the resulting (A_rec, Ea_rec) to a FLiBe-UF₄ loop at 600 °C — extrapolating Ea = 39 kJ/mol from 423 K to 873 K is a 6×-T extrapolation over an Arrhenius factor of ≈ 200. If the mechanism changes (different rate-limiting step at 600 °C, e.g. bulk-phase radical-radical vs surface-trap-mediated), the extrapolation may be wrong by orders of magnitude. The paper's §4.1 acknowledges this but the operational §8.2 result does not propagate any structural uncertainty.

**Fix.** Either inflate the Ea posterior to reflect the 6×-T extrapolation (e.g. assign ±10 kJ/mol structural uncertainty rather than the ±2 kJ/mol experimental), or report a second result with a mechanistically alternative Arrhenius (constant-Ea vs. T-dependent Ea) and bracket the FLiBe-UF₄ P(F₂).

### 7. The Iwamatsu-2022 cross-paper density correction (§9.2, line 1500) is treated as silently absorbed, but it shifted Ea(Zn) by ~ 5 kJ/mol and is a structural correction the reader should see.

The paper states "if a future audit reveals additional density-correction inconsistencies, the calibration would need to be redone." That is a fair statement but the *current* calibration uses the corrected values without showing what changes. A side-by-side comparison (Table) of Iwamatsu-2022 pre-correction, Iwamatsu-2026 density-corrected, and the present posterior would let the reader see how much each layer of correction moves the answer.

### 8. The companion-paper citation (`Tanore2026PaperI`) is leaned on heavily but the methodology is not minimally summarized.

The reader is told (line 187, line 235, line 661, line 1133, line 1512) to consult the companion paper for definitions, theorems, and figure reproductions. Two figures (Fig 2 fluoride kernel; Fig 3 meta-hierarchical layer) are even "reproduced from Tanore2026PaperI". This is acceptable for a true companion submission only if the companion paper is published *simultaneously*. If the companion paper is still under review, two operational-relevance items (Theorem 5 slow-manifold reduction; the four state-of-the-art benchmarks) must be stated as one-paragraph self-contained summaries here. As written, Paper 2 is not readable on its own.

## Specific Comments

- **Abstract, lines 64–72.** "33 rate-constant observations across nine metal scavengers and five hosts, four fluoride G-value compositions, and the Phillips 2022 NaCl-UCl₃ non-detection benchmark; contextual legacy and theoretical sources are retained as validation or provenance anchors." This sentence is the only abstract-level inventory and it conflicts with §2 ("23 sources"), Table 1 (14 rows), and §9.1 ("73 observations from 13 sources within the larger 23-source database"). Rewrite the abstract to give one unambiguous count.

- **Abstract lines 87–90.** "predicts a posterior-median cover-gas P(Cl₂) of 2.2 × 10⁻⁸ Pa at 60 years of operation … about ten orders of magnitude below a notional 100 Pa cover-gas screening threshold." Replace "predicts" with "yields, under the Phillips-bounded U(III)/U(IV) kernel," and add "this margin is dominated by the assumed U-redox sink and would collapse if the U-buffer is incomplete (see §8.1.4)." Then *add* §8.1.4.

- **§1 line 122.** "TerraPower2021, Mausolf2024" — TerraPower is a tech brief (not peer reviewed). Mausolf et al. 2024 in *JNM* is the right primary citation for MCFR radiolysis context; lead with it.

- **§2 Table 1 line 296.** Pikaev N=21 with modality K — but the row also lists "\eS + 6 metals". The 21 observations are *eS + metal* pairs, but the paper says elsewhere there are 9 metal scavengers. Reconcile.

- **§2 Table 1 line 297.** Makarov N=32, modality "K, S". Modality "S" is not defined in Table M (only T, K, A, C, G in the modality table). Define "S" or fix.

- **§2.2 Table 2.** Add the NaCl-UCl₃ and FLiBe-UF₄ density and dose-rate ranges used in the operational predictions, so Tables 2 and 8/9 are cross-checkable.

- **§3 Table 3 R5.** "$\eS \to (\text{bg. impurity})$" is a pure nuisance reaction that has no chemical interpretation. State explicitly that R5 is a first-order loss to *unmodeled* impurity scavenging and that its calibrated rate $\log_{10} k_{bg} \approx 16.5$ (Table 5) means the e_s⁻ lifetime against impurities is ≈ 10⁻¹⁶.⁵ s — physically nonsensical unless k_bg is in non-standard units. Check units.

- **§3 Table 3 R13–R15.** Reported as "$\log_{10} A = 30.51$" etc. with A in M⁻¹ s⁻¹. The corresponding rate at 600 °C is $\log_{10} A − E_a / (2.303 R T) = 30.5 − 33500/(2.303 × 8.314 × 873) = 30.5 − 2.0 ≈ 28.5$, i.e. $k \approx 3 × 10^{28}$ M⁻¹ s⁻¹. That is twenty orders of magnitude above the diffusion limit (≈ $10^{10}$–$10^{11}$ M⁻¹ s⁻¹ in molten salts). The $\log_{10} A = 30.5$ must be in some other base unit — perhaps it has been re-parameterized in M⁻¹ s⁻¹ × dimensionless reference, or expressed as $A / (e^{-T_\mathrm{ref}/T})$. **This needs immediate clarification**; as stated the table is unphysical. The Appendix table (Table 18, line 1654) gives R13 as $A = 3.2 × 10^{13}$ M⁻¹ s⁻¹ — which is at most 1000× above diffusion limit, still suspicious but not catastrophically wrong. So the body table's $\log_{10} A = 30.5$ does not equal $\log_{10}(3.2 × 10^{13}) = 13.5$. Resolve.

- **§5.1 line 712.** "80 walkers, 1500 burn-in samples and 1500 production samples each (total 22,400 production samples after thinning)." With 80 walkers × 1500 production = 120,000 raw samples; thinning to 22,400 is a factor of ~ 5, which is fine but should be stated as such. Also state the thinning ratio.

- **§5.3 line 826.** "Acceptance fraction 0.41 (within 0.2–0.5)" — for emcee-style stretch-move, 0.2–0.5 is the conventional target. Cite Goodman-Weare 2010 or Foreman-Mackey 2013 (already in bib).

- **§7 line 1015.** Master validation figure: 16 panels including panel (m) "Phillips 2022 NaCl-UCl₃ NULL benchmark and censored Bayes-factor cone" — but a censored data point cannot be "shown" in a meaningful overlay panel. Either clarify what is plotted (likely a 1D posterior on $G$ given the NULL) or move to a separate figure.

- **§7.1 Table 7.** $\chi^2_\nu$ = 0.00 for Toth-Felker (line 1072). Reporting a chi-square of identically zero for a 1-observation fit suggests the fit is exact by construction (i.e., $A_\mathrm{rec}$ was calibrated to satisfy the Toth-Felker balance condition exactly, so the residual is zero). Make this point explicit: this row is not validation, it is fitting.

- **§7.3 line 1145.** "$\ln K_{12} = 316.9$, equivalently $\log_{10} K_{12} = 137.6$, in favour of the kernel with the U(III)/U(IV) sink against the kernel without it." A Bayes factor of $10^{137}$ is astronomical; Kass-Raftery's scale tops out at $\ln K > 5$ as "decisive". Such an extreme number signals that the *no-U-redox* kernel predicts orders-of-magnitude excess Cl₂ far above 1000 ppm, so the data are infinitely impossible under that kernel. This is fine as a binary plausibility check but is a *category error* if presented as a quantitative measure of evidence. Reframe: "the data exclude the un-buffered kernel at any plausible significance; the quantitative Bayes-factor is therefore reported only to highlight kernel-A's qualitative consistency with the null."

- **§8.1 line 1218.** "$\eta_{\max} \approx 1 × 10^{-18}$." Derive this number step by step. The reader cannot reproduce $10^{-18}$ from $\dot D$, $G$, $V_\mathrm{salt}$, and the censored bound without seeing the algebra.

- **§8.1 line 1227.** Citation of Sohal 2010 for $K_H(\text{Cl}_2)$ in NaCl-UCl₃ — flag if measured vs. inferred (see Major Concern 4).

- **§8.1 line 1247.** "even if the chain efficiency were only $1 × 10^{-12}$ (six orders of magnitude above the Phillips posterior bound), the predicted 60-year Cl₂ partial pressure would be 20 Pa, still below the 100 Pa limit by a factor of five." Good — this is the right kind of bounding argument. Make this number the *headline* result instead of $10^{-8}$ Pa, because it is the operationally meaningful margin under a plausible relaxation of the U-buffer.

- **§8.2 line 1322.** "$P_{\mathrm{F}_2}^{ss}(\SI{873}{\kelvin}) = 31\,[19,\,44]\,\mathrm{Pa}$" — within a factor of three of 100 Pa. Note explicitly that *one* +1σ shift in either Ea or G(F₂) puts the prediction *above* 100 Pa. The current text downplays this.

- **§8.2 line 1355.** "designs operating below 500 °C would require additional F₂-management systems (a nickel-fluoride scrubber stage)" — provide a citation for nickel-fluoride scrubbers in MSR off-gas designs (Williams 2006 or Holcomb 2012 would work).

- **§8.3 Recommended experiments (lines 1377–1406).** Add a fourth recommended experiment: *direct measurement of U(IV)+e_s⁻ and U(III)+Cl_2^•⁻ rate constants*, since these are the only U-redox kinetics that are completely unconstrained by data.

- **§9.1 line 1428.** "73 observations drawn from 13 sources" — see Major Concern 3.

- **§9.2 line 1467.** "Single-source MCFR constraint" — this is the right honest statement; move it earlier in the paper so the reader sees it before §8.

- **§9.3 line 1510.** "Relation to companion methodology paper" — this section is too short. Expand to clarify which results are taken from the companion (and at what level — definitional, theorem, numerical) vs. computed fresh here.

- **§10 line 1556.** "To our knowledge, the model is the first calibrated multi-salt radiolysis kernel of this breadth" — Mausolf et al. 2024 *JNM* does a related (if smaller) calibration. State the differences explicitly.

- **Appendix A Table 18 line 1662–1667.** Rate constant pre-exponentials $A = 3.2 × 10^{13}$, $4.5 × 10^{13}$ M⁻¹ s⁻¹ are 100× above the diffusion limit in LiCl-KCl. Even with the Marcus/non-adiabatic acceleration, these are large; either the pre-exponential should be lower (in which case Ea is lower) or there is a unit issue. Cross-check against the published Iwamatsu values directly.

- **Appendix B (Cross-referenced data inventory, line 1685–1726).** The directory listing has 13 source directories. This matches "13 sources" in §9.1 but not "23 sources" in §2 — see Major Concern 3.

- **Appendix C Table 23 line 1755.** The autocorrelation time τ ~ 35 on the slowest parameter, with 1500 production steps × 80 walkers = 1.2 × 10⁵ raw samples. Number of independent samples per chain ≈ 1500 / 35 ≈ 43, times 80 walkers = 3,440 — much smaller than the claimed ESS of 510–5290. Reconcile (likely the 22,400 thinned samples are not actually independent).

- **Appendix D Table 24 line 1814.** $P(\mathrm{Cl}_2)$ varies by exactly the same factor (3.5×) for $G$ at ±1σ; this is suspicious for a posterior that is supposedly skewed. Check that the table is showing the right perturbation (1σ in $\log G$ vs. $\sigma$ in $G$).

- **Conventions, throughout.** Mix of 90% CIs (operational) and 95% CIs (Table 5 posterior). Pick one and use throughout, or state both as default and exception.

## Minor Comments

- **Line 84.** "Cr, Nd, and Cf analogues" — italicize element symbols or use chemical formatting consistently. The bib has `\citep{CastroBaldivieso2026}` while the prose has both "Castro-Baldivieso" and "Castro Baldivieso" — pick one.

- **Lines 92–98.** Abstract's last clause runs on. Split into two sentences.

- **Line 188.** "we refer to that paper for the five-layer …" — comma after "paper".

- **Line 232.** "(B) a CSV-level cross-reference of every digitized observation" — Appendix B is a directory listing, not an obs-level cross-reference. Either rename or add a per-obs table.

- **Line 442.** "R8 & $\mathrm{Cl}^- + h^+ \to \mathrm{Cl}^{\bullet}$" — the symbol $h^+$ for "hole" is fine but should be defined where introduced.

- **Line 605.** Eq. (3) "static" kernel for fluoride — describe the integration to steady state explicitly; the result $P^{ss} = G \dot D c_0 / k_\mathrm{rec}$ in Eq. (4) is the steady state, but Eq. (3) is the rate equation. Numbering of equations is inconsistent (Eq. (eq:fkernel) does not appear in (3)).

- **Line 668.** "Akiyama 1994 reports transients but no usable Arrhenius parameters" — this is also true of Hagiwara 1987's e_s⁻ kinetics; the prose has been more generous to Hagiwara than to Akiyama. Reconcile.

- **Line 740.** "consistent with zero at the 95% level for LiCl-KCl, with modest negative shifts for NaCl ($\eta_{\log A} = -0.56$)" — Table 5 line 773 gives $\eta_{\log A}^\mathrm{NaCl} = -0.56\,[-1.31, +0.19]$; that interval *includes* zero too. Restate as "modest, statistically marginal".

- **Line 869.** "where tildes denote z-score-standardized features" — this is already stated in the preceding sentence (line 859). De-duplicate.

- **Line 905.** $b^\mathrm{NaCl} = -2.13\,[-2.52, -1.70]$ in the meta-hier table; compare with $\eta_{\log A}^\mathrm{NaCl} = -0.56$ in Table 5. These are conceptually distinct (one is a *facility offset by host*, the other is a *host correction to the Arrhenius layer*) but readers will conflate them — add a clarifying sentence.

- **Line 1057.** Table 6 — add a column showing per-metal number of LOMO holdouts and reference to which paper the data come from.

- **Line 1167.** "These thresholds are not regulatory design-basis limits" — good, but also state they are not from any specific licensing authority. Add a sentence on what 100 Pa actually represents (corrosion onset? Off-gas system safety? Hot-cell criterion?).

- **Line 1351.** "At a representative operating temperature of 500 °C the predicted P(F₂)_ss is ~ 60 Pa" — this is *above* the 31 Pa median at 600 °C. The text should call out that the operational temperature is itself a uncertainty contributor.

- **Line 1569.** Acknowledgments thank Brookhaven LEAF "for providing primary-data access" — is this *digitized* data or *true primary* data? Clarify. If it is digitized from published figures, no permission was needed (it is OCR'd from published papers).

- **References.** Typos: "Holcomb2012" entry has `year = {2010}`. Mausolf et al. 2024 entry has `pages = {154974}` but volume is 593 — verify against the journal. `Castro Baldivieso` vs `CastroBaldivieso` author key inconsistency. `Tanore2026PaperII` cites itself (Paper 2) inside Paper 2 — remove the self-citation or convert to footnote.

## Strengths

1. **Honest single-paragraph statement of limitations (§9.2)** is unusual for an application paper and should be retained. Particularly the "single-source MCFR constraint" and "chemistry-feature LOMO bias" items are exactly what a referee wants to see.
2. **The Pikaev-Iwamatsu inter-facility offset is the right kind of calibration result** — a structural ambiguity that has been hand-waved in the LEAF literature is here resolved as a single quantitative parameter with a credible interval that excludes zero by ≈ 5σ.
3. **The censored modality treatment for Phillips 2022** is methodologically sound (Tobit likelihood with the slow-manifold reduction), and the Bayes-factor stability test in Table 25 is a nice touch (provided the BF value itself is reinterpreted — see Specific Comment on §7.3).
4. **The actionable experimental roadmap (§8.3)** with three specific experiments, each with a quantified posterior-variance reduction, is the kind of thing that *JNM* readers (designers) want from a calibration paper.
5. **Per-dataset $\chi^2$ table (Table 7)** showing the model is mildly conservative ($\chi^2_\nu = 0.74$, coverage 88%) honestly reports what the calibration achieves.
6. **The MCMC diagnostics in Appendix C** are at the standard expected for a 24-dim Bayesian inference paper; min ESS 510 and $\hat R < 1.05$ are fine.
7. **The 60-year extrapolation is rate-limited by a single steady-state algebra** (Eq. 8.1 / Eq. 8.5), not by any "extrapolating an ODE for $10^{17}$ time steps" — so the *kinetic* extrapolation is well-posed even if the *physical* extrapolation needs more justification.

## References that should be added

1. **Sims, H. E. et al.** Multiple papers on Ce/Eu/Pm radiolysis in chloride and aqueous-Cl⁻ systems (Sims is a co-author on the Castro-Baldivieso 2026 paper but his lanthanide/actinide work in chloride media is independently cited in the LEAF literature). Without these, the chemistry-feature regression for f-element scavengers is under-anchored.
2. **Forsberg, C. W.** "Safety and licensing aspects of the molten salt reactor" or the 2019 *Nuclear Technology* MSR safety review — for the cover-gas Cl₂ tolerance argument and the 100-Pa threshold justification.
3. **Williams, D. F. and Toth, L. M.** "Chemistry of the molten Li_2BeF_4 / LiF-BeF_2 system in the AHTR" or the ORNL TM series on FLiBe corrosion chemistry. The Davis 2022 paper is the *only* FLiBe G(F₂) data set cited; Williams's measurements of redox potential and HF/H₂ equilibria provide structural-uncertainty bracketing.
4. **Andresen, P. L. or Macdonald, D.** on radiolysis-driven IGSCC in LWR primary coolant — not direct, but useful for the methodological argument that chronic-irradiation radiolysis matters at low partial pressures.
5. **Janz, G. J.** *Molten Salts Handbook* / NIST molten-salt property compilations — for the Henry constant $K_H$ value and NaCl-UCl₃ density.
6. **Romatoski & Hu (2017)** is already cited but underused — it has the appropriate FLiBe density and viscosity for §8.2 mass-transfer estimates.
7. **Bersuker (2004) or Banks & Hutchinson (2010)** on Marcus theory in molten salts — for the chemistry-feature regression's "electronegativity" predictor justification.
8. **Mausolf et al. (2024, *JNM*)** is cited in §1 but not used as a benchmark. Compare directly: their G(Cl₂) estimate vs. the posterior here, in §5 or §7.
9. **Joseph, B. and Suresh, A.** on actinide redox in molten chloride — provides published $E°$ for U(III)/U(IV) in NaCl-CsCl and LiCl-KCl that anchors the *thermodynamic* feasibility of R12.
10. **Kennedy & O'Hagan (2001)** is in the bib but not actually cited — either cite or remove.

## Suggested rewrites

### Rewrite 1 — Abstract, lines 87–98 (the MCFR claim)

**Current:**
> "We deploy the calibrated model to two operational MSR design problems: a 500 MW(th) NaCl-UCl₃ fast-spectrum reactor predicts a posterior-median cover-gas $P_{\mathrm{Cl}_2}$ of $2.2 \times 10^{-8}$ Pa at 60 years of operation (90% credible interval $[3.7 \times 10^{-9}, 8.5 \times 10^{-8}]$ Pa), about ten orders of magnitude below a notional 100 Pa cover-gas screening threshold; …"

**Proposed:**
> "We deploy the calibrated model to two operational MSR design problems. A 500 MW(th) NaCl-UCl₃ fast-spectrum reactor, under the assumption that the U(III)/U(IV) redox couple is the dominant Cl₂ sink and that no significant non-uranium scavenger is present, has its cover-gas Cl₂ partial pressure constrained by the Phillips 2022 censored benchmark to lie below approximately $10^{-7}$ Pa over 60 years of operation — far below any plausible engineering screening threshold but driven entirely by the U-buffer prior rather than by direct measurement. A FLiBe-UF₄ salt loop predicts a steady-state cover-gas $P_{\mathrm{F}_2}$ of 31 Pa (90% CI [19, 44] Pa), within a factor of three of the 100 Pa screening anchor and dominated by the Davis-2022 G(F₂) and Toth-Felker recombination Arrhenius uncertainty."

### Rewrite 2 — §3.4 R12 paragraph, lines 553–563 (U(III) sink justification)

**Current:**
> "This is the central reaction for predicting the Phillips 2022 NULL benchmark. The U(III)/U(IV) redox couple acts as a buffered sink that suppresses Cl₂ release into the cover gas. No direct rate constant has been measured because (i) U(III) is too oxygen-sensitive for routine pulse radiolysis and (ii) the LiCl-KCl chemistry is incompatible with U(III) at moderate temperatures. We assign a prior centered at the diffusion-limited value with two-decade lognormal spread; …"

**Proposed:**
> "R12 is the dominant proposed Cl₂ sink in the operational MCFR prediction and is the single most influential unmeasured reaction in the model. No direct rate constant has been measured: U(III) is too oxygen-sensitive for routine pulse radiolysis, and the LiCl-KCl chemistry where the LEAF facility operates is incompatible with U(III) at moderate temperatures. We assign a wide-band prior centered at the diffusion-limited value ($\log_{10} k = 10$, two-decade lognormal spread) and rely on the Phillips 2022 censored upper bound on cover-gas Cl₂ to constrain the *product* $G(\mathrm{Cl}^\bullet) \cdot \eta$, where $\eta$ collapses the U-buffer efficiency. We emphasize that the resulting posterior on $\eta$ is structurally an *upper bound*: the Phillips data are silent on the *magnitude* of $k_{R12}$ as long as the U(III)/U(IV) cycle is at least as efficient as the prior median. A direct LEAF-facility measurement of $k_{R12}$ in a NaCl-UCl₃ matrix (recommended in §8.3) would convert this upper bound to a posterior point estimate and is the highest-priority gap in the current model."

### Rewrite 3 — §8.1 paragraph on the result, lines 1230–1238

**Current:**
> "Figure 8 shows the predicted cover-gas Cl₂ partial pressure over the 60-year lifetime with full posterior uncertainty propagation. The posterior median at year-60 is $P_{\mathrm{Cl}_2}^{(50)} = 2.2 \times 10^{-8}$ Pa, with 90% credible interval $[3.7 \times 10^{-9}, 8.5 \times 10^{-8}]$ Pa. This is about ten orders of magnitude below the notional 100 Pa screening threshold; the posterior-predictive probability of exceeding the threshold at 60 years is $< 10^{-12}$."

**Proposed:**
> "Figure 8 shows the predicted cover-gas Cl₂ partial pressure over the 60-year lifetime under three nested assumption sets: (A) the U(III)/U(IV) sink kernel calibrated against Phillips 2022 as the only Cl₂ scavenger; (B) the same kernel with a conservative 1000× weakening of the U-buffer prior; and (C) a no-U-buffer kernel with literature-prior $G(\mathrm{Cl}^\bullet)$. Under assumption (A), the posterior median at year 60 is $2.2 \times 10^{-8}$ Pa, with 90% credible interval $[3.7 \times 10^{-9}, 8.5 \times 10^{-8}]$ Pa — about ten orders of magnitude below the 100-Pa engineering threshold, but the *width* of this interval reflects only the parametric posterior on the censored bound and *does not include structural uncertainty* in the U-buffer assumption. Under (B), the 60-year median rises to ~ 20 Pa, still below threshold but by less than a decade. Under (C), the 60-year median is ~ $10^{10}$ Pa (above atmospheric pressure within months), confirming that *some* radiolytic sink is required and that the Phillips data exclude the no-sink scenario decisively. We adopt (A) as the headline result, noting that the actual margin against the threshold is bounded above by the (A) prediction and below by approximately the (B) prediction."

### Rewrite 4 — §3 introduction, lines 386–409 (network completeness)

**Current:**
> "The chloride radiolysis kernel comprises 20 elementary reactions partitioned into three functional groups: …"

**Proposed insert at end (new paragraph):**
> "The 20-reaction kernel above is the minimum closure needed to reproduce the pulse-radiolysis transient data in §2; for chronic-irradiation predictions in §8, three additional reaction classes are present in any real MCFR salt but are *not* included in the kernel:
> (i) oxide-impurity scavenging (the reaction of dissolved O²⁻ with Cl_2 or Cl_2^•⁻ to form UO₂Cl₂ or analogous oxychloride species), which is plausibly the dominant Cl₂ sink in a non-perfectly-anhydrous salt and would *tighten* the predicted margin;
> (ii) redox-active fission products (Cs⁺/Cs⁰, I⁻/I•, Te²⁻, Sr²⁺/Sr⁰, …) which would compete with U(III) for the Cl_2^•⁻ pool over the reactor lifetime;
> (iii) hydrolytic chemistry at the cover-gas interface (Cl_2 + H_2O → HCl + ½ O_2 catalyzed by reduced metals), which would shift Cl_2 to HCl with consequences for both off-gas and corrosion.
> Excluding these is appropriate for the calibration set, where the experimental data do not constrain them; it is a *significant assumption* for the operational predictions, where they could perturb the U-buffer balance by orders of magnitude. We return to this in §9.2 (Limitations)."

### Rewrite 5 — §7.3 Bayes factor interpretation, lines 1145–1158

**Current:**
> "The numerical result is $\ln K_{12} = 316.9$, equivalently $\log_{10} K_{12} = 137.6$, in favour of the kernel with the U(III)/U(IV) sink against the kernel without it. This is the strongest formal evidence available from the Phillips data: it constitutes a quantitative endorsement of the U-buffered redox network rather than the un-buffered alternative."

**Proposed:**
> "The numerical result is $\ln K_{12} = 316.9$. This is far beyond the Kass-Raftery 'decisive' threshold ($\ln K > 5$) and signals that the un-buffered chloride kernel predicts a 60-year cover-gas $P(\mathrm{Cl}_2)$ that is astronomical relative to the Phillips 1000-ppm detection limit. The numerical value of $\ln K$ is therefore *not a quantitative measure of evidence* in the sense that, say, $\ln K = 10$ would be — under the un-buffered kernel the data are essentially impossible. We report it for completeness but emphasize that what the data establish is the *binary* fact that some efficient Cl_2 sink must be present; the *quantitative* magnitude of the U-buffer rate constant remains constrained only by the prior."

### Rewrite 6 — §6.2 LOMO interpretation, lines 962–984 (chemistry-feature reach)

**Current:**
> "The LOMO performance is bimodal: metals that lie within the chemistry-feature interpolation envelope (Cd, Tl, Ag, Ca, Ba, Sr, Zn) are recovered within ±1 log unit with 100% coverage of the 90% credible band. Metals that lie outside that envelope (Cr, Nd, Cf) are systematically under-predicted by 2–4 log units."

**Proposed:**
> "The LOMO performance reveals a structural limitation rather than a true 'envelope' effect. The metals that LOMO recovers well (Cd, Tl, Ag, Ca, Ba, Sr, Zn) are all in the *Pikaev microsecond corpus*, and the metals that LOMO under-predicts by 2–4 log units (Cr, Nd, Cf) are all in the *Brookhaven LEAF picosecond corpus*. The chemistry-feature regression has no facility-effect parameter analogous to $b^{(\mathrm{Pikaev})}$ in the integrated chain; consequently, when a LEAF-measured metal is held out, the regression must predict its rate constants without the systematic LEAF-vs-Pikaev correction and necessarily under-predicts.

> The operational implication is that the chemistry-feature layer is reliable for predicting rates of fission-product elements that would (counterfactually) be measured by Pikaev-class apparatus, but unreliable for elements whose actual operational chemistry is closer to LEAF picosecond conditions. Practically, every operationally-relevant fission product (Cs, Sr, Te, I, Eu, Pu, Am, Cm) would, if measured today, be measured at LEAF. The layer is therefore *descriptive* — it identifies the dominant chemistry features ($\beta_\chi$ statistically significant; $\beta_z$ marginal) — but is not yet *predictive* for unmeasured operational elements. A direct extension of the layer with a facility-effect intercept (or a hierarchical layer over $\{\mathrm{Pikaev}, \mathrm{LEAF}\}$ where the LEAF-shift is identified from the in-corpus Cr/Zn data) would resolve this gap; we recommend it as a follow-on to Tanore2026PaperI."

### Rewrite 7 — §1 contributions list, lines 201–227

**Current:**
> Item 4: "A comprehensive validation of model predictions against every digitized observation, including a multi-panel figure and a tabulated per-dataset $\chi^2$ and posterior coverage analysis (§7)."

**Proposed:**
> "A comprehensive validation of model predictions against the 73 likelihood-bearing observations digitized from the corpus (out of a total 14 source papers; see §2 and Appendix B), including a 16-panel master overlay figure with shaded 90% credible bands per dataset, a per-dataset reduced chi-square table that honestly reports a mild conservativeness ($\chi^2_\nu = 0.74$, coverage 88%), and an explicit LOMO-vs-direct partition (§7)."

(And similarly for item 5: clarify the *operational MSR* claims are bounded by the U(III)/U(IV) prior, not by direct measurement.)

### Rewrite 8 — §9.2 Limitations, lines 1467–1507 (re-ordering and strengthening)

**Current ordering:** Single-source MCFR → static fluoride kernel → LOMO bias → host corrections → composition dependence → density correction.

**Proposed reorder (impact-weighted):**
> 1. *Single-source MCFR constraint and U-buffer prior dependence.* (As currently written, but acknowledge that the 10⁻⁸ Pa median is a prior-shaped propagation, not an empirical extrapolation.)
> 2. *Chloride network completeness.* The kernel omits oxide-impurity scavenging, fission-product redox competitors, and cover-gas hydrolytic chemistry; each could perturb the operational P(Cl_2) prediction by orders of magnitude in either direction.
> 3. *Static fluoride kernel.* (As currently written.)
> 4. *Toth-Felker extrapolation to operational T.* The recombination Arrhenius is extrapolated from $T < 423$ K to 873 K, a 6× temperature ratio; the resulting Arrhenius factor of $\approx 200$ is structurally uncertain.
> 5. *Chemistry-feature LOMO bias.* (As currently written but with the facility-effect framing.)
> 6. *Single-host calibration of host corrections.* (As currently written.)
> 7. *G-value composition vs. temperature dependence.* (As currently written.)
> 8. *Cross-paper density correction inheritance.* (As currently written.)

### Rewrite 9 — §10 Conclusion, lines 1533–1564

**Current opening:**
> "We have constructed and validated a calibrated multi-salt radiolysis model …"

**Proposed:**
> "We have constructed a calibrated multi-salt radiolysis model that simultaneously honors transient pulse-radiolysis kinetics, scalar Arrhenius literature values, censored chronic-irradiation upper bounds, and gamma-induced G-values. The 24-parameter HBMAE posterior recovers published Cr and Zn Arrhenius parameters within $1\sigma$ and identifies the long-standing Pikaev-Iwamatsu inter-facility offset as $b^{(\mathrm{Pikaev})} = -3.35\, [-4.53, -2.44]$. Applied to a 500-MW(th) NaCl-UCl_3 MCFR, the model — *conditional on the assumed U(III)/U(IV) buffer kernel and on Phillips 2022 as the only chronic data anchor* — predicts cover-gas Cl_2 partial pressure approximately ten orders of magnitude below a 100-Pa screening threshold over 60 years. This margin is a Bayesian propagation of the U-buffer prior under a censored upper bound, not an empirical measurement; the headline is therefore that *no plausible relaxation of the U-buffer brings the prediction within an order of magnitude of the threshold*, but the *quantitative* margin is the propagated prior width.

> Applied to a FLiBe-UF_4 thermal loop at 600 °C, the model predicts steady-state cover-gas $P(\mathrm{F}_2) = 31 \,[19, 44]$ Pa, within a factor of three of the 100-Pa screening threshold; the dominant uncertainty contributors are $G(\mathrm{F}_2)$ from Davis 2022 (71%) and the Toth-Felker recombination $E_a$ extrapolation to 873 K (21%).

> Three direct experiments — extended-duration NaCl-UCl_3 chronic irradiation with order-of-magnitude tighter detection, LEAF-equivalent FLiBe pulse-radiolysis, and direct measurement of U(III)+Cl_2^•⁻ and U(IV)+e_s⁻ rate constants — would together tighten the operational predictive variance by approximately three orders of magnitude (MCFR Cl_2) and one order of magnitude (FLiBe F_2). The companion methodology paper Tanore2026PaperI develops the HBMAE formalism in full."

## Final assessment

The paper makes a real contribution: a multi-modality, multi-host Bayesian calibration of molten-salt radiolysis is genuinely new, the Pikaev-Iwamatsu offset is the kind of clean structural result that justifies the framework, and the censored Phillips treatment is methodologically clean. The honest reporting of LOMO failure modes and per-dataset $\chi^2$ is encouraging.

However, the MCFR P(Cl_2) headline (ten orders of magnitude below threshold) needs to be reframed as a Bayesian upper-bound propagation rather than a "prediction"; the chloride network must acknowledge what it omits (oxide, fission-product, hydrolytic chemistry); the corpus inventory must be made self-consistent (the 14/23/33/73 mismatch); the §3 Table 3 pre-exponentials need a units sanity check (the body table gives $\log_{10} A = 30$ while the appendix gives $A = 10^{13}$ M⁻¹ s⁻¹ — those are the same parameter and they disagree); and the companion-paper dependency must be made minimally self-contained. The FLiBe F_2 prediction at 31 [19, 44] Pa is honest and well-bounded; that section needs only minor tightening on the Toth-Felker extrapolation.

With these changes — none of which require new science, all of which require clearer presentation and more honest framing — the paper is in scope for *JNM* and should be published. I recommend **Major Revisions** with a target re-review on the rewritten Sections 3, 8, and 9.

— Reviewer (anonymous), 2026-05.
