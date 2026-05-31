# Style Diff: Current Drafts vs. Mauricio Tano's Published Voice

This document highlights specific passages in `article_paper1.tex` and
`article_paper2.tex` that read as LLM-generated or formulaic relative to
Tano's published style (see `MAURICIO_STYLE_GUIDE.md`), and proposes
in-voice rewrites.

## 1. Global pattern findings

### 1.1 Over-use of "catastrophically", "the unique method", "uniquely"

`article_paper1.tex` uses *"catastrophically fails"* / *"catastrophically"*
three times (lines 2080, 2157, 2160) and *"the unique method that
simultaneously..."* (line 117) in the abstract. Tano writes
quantitatively about failure modes ("M_1 exhibits ~35 kJ/mol E_a bias",
"$\chi^2 = 96.5$ on 5 dof, $p<10^{-17}$"); the word "catastrophically" is
not in his vocabulary. The "unique method that simultaneously [verb 1],
[verb 2], [verb 3], [verb 4], and [verb 5]" construction is a hallmark
LLM-abstract clause.

**Fix:** Replace "catastrophically fails" with the numerical failure
metric. Replace "the unique method that simultaneously [list]" with
"the only framework that combines [feature], [feature], [feature]" or
better, drop the unique claim and let the comparison table carry it:
*"HBMAE is the only entry in Table 1 with full coverage of the five
capability axes; each comparator addresses a strict subset."*

### 1.2 "Defensible" appears 5 times

The word *"defensible"* (papers 1 and 2 combined: 6 occurrences) is doing
heavy semantic work as a placeholder for "validated", "auditable",
"licensing-grade", or "with quantified uncertainty". Tano does use
"defensible" occasionally but not as a recurring leitmotif. Most of these
should be replaced with the specific property being claimed.

| Line | Current | Suggested |
|------|---------|-----------|
| p1:146 | "defensible, validated predictions" | "validated predictions" (drop one) |
| p2:58 | "defensible predictions of how" | "quantitative predictions of how" |
| p2:123 | "a defensible understanding" | "a quantitative understanding" |
| p2:176 | "a defensible cover-gas Cl_2 ... prediction with quantified uncertainty" | "an auditable, uncertainty-quantified cover-gas Cl_2 prediction" |
| p2:842 | "within a defensible uncertainty" | "within a documented uncertainty domain" |
| p2:1200 | "require defensible justification" | "require explicit justification" |

### 1.3 "Five concrete ways" / "five distinct methodological problems" / "five-layer framework" / "six theorems"

`article_paper1.tex` lines 84-100, 178, 211-217, 226-228 stack the numeric
inventory: *"five distinct methodological problems"*, *"five concrete
ways"*, *"five-layer framework"*, *"six theorems and four propositions"*.
This rhythm of *"N distinct X..., N-layer Y..., M theorems..."* across the
abstract and intro is a classic LLM scaffolding pattern.

**Fix:** Keep the numbered enumeration *once* (the introduction bullet
list of five gaps). In the abstract, name two or three of the gaps prose-
style rather than restating the count. Tano-style abstract phrasing is
descriptive: *"The framework combines a constrained-topology prior with a
hierarchical Arrhenius layer that couples shared reactions across host
salts, a composite likelihood across observation modalities, and an
explicit facility-effect hierarchy."* (No numeric scaffolding.)

### 1.4 "We prove six theorems characterizing..."

`article_paper1.tex` line 99 leads with the claim. Tano does not advertise
*"we prove"* in abstracts; even when he proves things, he writes more
conservatively. Replace with *"Six theorems characterize..."* or just
remove the meta-statement and let Section 5 carry the proofs.

### 1.5 "To our knowledge, ... is the first..."

`article_paper2.tex` lines 1418 and 1556 both open with *"To our knowledge,
the integrated model presented in this paper is the first..."* This is one
of the most reliable LLM tells. Tano writes about first-of-kind work without
the introductory disclaimer; he says it as a flat fact and lets the
literature review back it up: *"No published compilation simultaneously
calibrates against the chloride and fluoride corpora; the present model
does."* or *"The compilation in Section 2 is the first to assemble
chloride and fluoride pulse-radiolysis kinetics in a single calibrated
framework."*

---

## 2. Per-passage rewrites

### 2.1 Abstract, paper 1 (lines 82–123)

**Current opener (lines 84–95):**

> Bayesian calibration of chemical-kinetic networks against heterogeneous,
> multi-facility, multi-host experimental data raises five distinct
> methodological problems that no single existing framework simultaneously
> resolves: (i) the data combine continuous time-series, scalar rate
> constants, one-sided detection limits, and derived Arrhenius pairs from
> independent laboratories; (ii) the same elementary reactions are reported
> in chemically distinct host environments; (iii) documented inter-laboratory
> systematic biases require explicit gauge separation; (iv) the candidate
> reaction-network topology is itself uncertain; and (v) the parametric
> model is misspecified relative to the true physical chemistry, requiring
> an identifiability-safe discrepancy term.

This is 80+ words in one sentence, with five Roman-numeral sub-clauses,
itself a hallmark LLM-abstract device. Tano's published abstracts open
with a 20-word problem-description sentence and then layer.

**Tano-voice rewrite:**

> Bayesian calibration of chemical-kinetic networks for molten-salt
> radiolysis must handle heterogeneous data: dense pulse-radiolysis
> transients, scalar Arrhenius pairs from multiple laboratories, and
> one-sided detection limits from chronic-irradiation experiments. The
> same elementary reactions appear in chemically distinct host salts, the
> reaction-network topology itself is uncertain, and inter-laboratory
> systematic biases are documented but not formally separated from the
> intrinsic chemistry. No published framework addresses these features
> together. We present Hierarchical Bayesian Mechanism-Adequacy Estimation
> (HBMAE), which combines a constrained-topology prior, a hierarchical
> Arrhenius layer that pools shared reactions across host salts, a
> composite likelihood across observation modalities, a parametric
> Kennedy--O'Hagan discrepancy with bounded flexibility, and an explicit
> facility-effect hierarchy.

### 2.2 Abstract close, paper 1 (lines 116–123)

**Current:**

> Synthetic-data benchmarks against reversible-jump MCMC, Kennedy--O'Hagan
> with Gaussian-process discrepancy, B2BDC, and PSIS-LOO demonstrate that
> HBMAE is the unique method that simultaneously recovers ground-truth
> parameters, handles censored modalities, identifies facility offsets,
> and produces calibrated uncertainty across all five problem features.

**Issue:** "the unique method that simultaneously [4 verbs] across all five
problem features" — LLM-stack of verbs and over-claim. The numerical
results elsewhere in the abstract are strong; this sentence weakens them.

**Tano-voice rewrite:**

> On a synthetic benchmark with known ground truth, HBMAE recovers the
> Arrhenius parameters within 1% bias and attains nominal 90% credible-
> interval coverage; reversible-jump MCMC with weakly informative priors
> exhibits a $\sim 35$~kJ/mol $E_a$ bias of the same sign and magnitude
> as the published-data result. Direct comparison against Kennedy--O'Hagan
> with Gaussian-process discrepancy, B2BDC, and PSIS-LOO is given in
> Section~\ref{sec:comparison}.

### 2.3 Paper 1 line 211–217 (contribution paragraph)

**Current:**

> This paper introduces HBMAE, a Bayesian framework that composes five
> ingredients to address these five gaps. Each ingredient has well-
> established precedent in the broader uncertainty-quantification
> literature; the contribution of the paper is the *synthesis* that
> targets molten-salt radiation chemistry specifically, together with the
> theorems characterizing when the synthesis is well-defined and
> identifiable, and the empirical validation against four published data
> sources.

**Issue:** "composes five ingredients to address these five gaps" — the 5+5
numeric parallelism is too tidy to be human-written. Also "ingredients" is
an unusual choice for "components".

**Tano-voice rewrite:**

> This paper introduces HBMAE, a Bayesian framework whose components are
> individually well-established in the uncertainty-quantification
> literature. The contribution is the synthesis specific to molten-salt
> radiation chemistry, the identifiability and contraction theorems of
> Section~\ref{sec:theorems}, and the validation against four published
> data sources.

### 2.4 Paper 2 abstract (lines 58–98)

**Current opener (line 58):** *"The molten salt reactor (MSR) concept
requires defensible predictions of how the salt's chemical inventory
evolves under sustained radiation..."*

This is closer to Tano's voice — but the abstract then goes on to deliver
five distinct numeric results in 30 lines, which is unusually dense. Compare
to Tano's actual abstracts (Walker/Tano 2023 abstract is 270 words and
delivers one main result with two sub-claims).

**Suggested compression:** Move the chemistry-feature-layer LOMO regression
detail out of the abstract and into Section 5. Keep the abstract focused on
(a) the construction, (b) the operational predictions (MCFR and FLiBe), (c)
the dominant uncertainty contributors.

### 2.5 Paper 2 Discussion opener (line 1418)

**Current:**

> To our knowledge, the integrated model presented in this paper is the
> first calibrated multi-salt radiolysis kernel against which an MSR
> designer can run a screening-level cover-gas Cl_2 or F_2 prediction
> with quantified uncertainty.

**Tano-voice rewrite:**

> The model presented in Sections~\ref{sec:chloride}--\ref{sec:operational}
> assembles the first integrated chloride+fluoride radiolysis kernel
> calibrated against the combined pulse-radiolysis and chronic-irradiation
> literature. An MSR designer can use the kernel for screening-level
> cover-gas Cl_2 or F_2 predictions with quantified uncertainty.

### 2.6 Paper 2 Conclusion opener (line 1533)

**Current:** *"We have constructed and validated a calibrated multi-salt
radiolysis model for molten chloride and fluoride MSR fuels by combining a
digitized pulse-radiolysis and chronic-irradiation corpus with the
Hierarchical Bayesian Mechanism-Adequacy Estimation framework."*

This is in Tano's voice already — leave it. A typical Tano conclusion
recapitulates the construction in past tense ("we have constructed and
validated"), states the principal numeric findings, and points to follow-on
work or open problems. The current draft matches this pattern.

---

## 3. Patterns to systematically purge

Run a global edit pass for these phrases:

| Phrase to remove | Replacement strategy |
|------------------|----------------------|
| *"the unique method that simultaneously..."* | name the feature combination prose-style or refer to the comparison table |
| *"catastrophically fails"* | give the numerical failure metric (bias, coverage, RMSE) |
| *"To our knowledge,"* | drop the disclaimer; state as fact |
| *"five distinct ... five concrete ... five-layer ... six theorems"* (in close proximity) | keep the numeric scaffolding once, not three times |
| *"defensible"* (>2 occurrences) | substitute with the specific property (validated, auditable, quantitative) |
| *"comprehensive validation"* | "validation against every digitized observation" or similar concrete phrasing |
| *"the strict, defensible advantages"* | "the advantages of HBMAE over each comparator" |
| *"composes N ingredients to address N gaps"* | "addresses the gaps of Section X through Y" |
| *"simultaneously resolves"*, *"simultaneously handles"* | "jointly addresses" or just "addresses" |

---

## 4. Patterns to *add* (under-represented in current drafts)

These Tano signatures are absent or thin in the current drafts:

1. **Variance-attribution paragraph.** Paper 2 has one (Section 6's
   "Uncertainty attribution" enumeration) but paper 1 does not. Adding a
   2-sentence posterior-variance attribution near the end of the validation
   section would feel more Tano-like.
2. **"This work represents a fundamental step toward..."** This phrase or
   its variants should appear once, typically in the Conclusion or the
   abstract closing. Neither draft currently has it.
3. **Named-tool sentences.** Tano typically lists the specific software
   stack (MOOSE, Griffin, Pronghorn, Thermochimica, OpenMC). Paper 1
   doesn't name a stack because the work is statistical/Bayesian; that is
   fine. Paper 2 does mention Thermochimica indirectly. Consider adding a
   sentence in Paper 2's Methods naming the computational stack used to
   integrate the ODEs and run the MCMC (emcee, scipy, the specific BDF
   solver), with citations.
4. **Italicized short item-leads in itemize lists.** Paper 2's strengths/
   limitations subsections already use this. Paper 1's Discussion section
   could be edited similarly: each limitation item should open with an
   italicized 2–5 word phrase summarizing it ("*Brynjarsdottir--O'Hagan
   asymptotic limitation.*", "*Parametric discrepancy specification.*",
   "*Small-K facility identifiability.*").

---

## 5. The five most impactful changes (priority order)

Apply these first; they will swing the perceived authorship most.

1. **Remove "the unique method that simultaneously..." and similar
   superlative-stacking** from the paper 1 abstract (line 117) and any
   other location. Replace with descriptive prose or a table reference.
2. **Replace "To our knowledge, ... is the first..."** (paper 2 lines 1418
   and 1556) with a flat-fact construction. This single change visibly
   de-LLMs the Discussion and Conclusion sections.
3. **Replace "catastrophically fails"** (paper 1, three instances) with
   the numerical failure metric the sentence is already carrying anyway.
4. **De-scaffold the numeric inventory** in paper 1's abstract (lines 84–95
   and 116–123). The current "five distinct problems ... five-layer
   framework ... six theorems ... four comparators" cadence is the single
   most LLM-flavored feature of the draft. One Roman-numeral enumeration
   per abstract is the upper bound for Tano voice.
5. **Add a Limitations subsection to paper 1's Discussion section that
   matches the format of paper 2's**: numbered list, italicized short
   item-leads, 2–4 sentences per item. Paper 1's current limitations
   (lines 2188–) are written as a `\paragraph{Limitations.}` block of
   prose with three numbered points; convert these to the same itemized
   format as paper 2 lines 1466–1508.

