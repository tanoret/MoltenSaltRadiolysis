# Mauricio E. Tano (Tano-Retamales) — Concrete Writing Style Guide

Compiled 2026-05-30 from publicly available abstracts, OSTI bibliographic
records, the Frontiers in Nuclear Engineering open-access full text, and the
Pure / ResearchGate / Google Scholar profile listings. Where the journal
landing page was paywalled (ScienceDirect, Tandfonline) the abstract was
recovered through OSTI / Frontiers / arXiv preprint mirrors. Direct full-text
analysis was performed on the open-access Frontiers paper (Walker, Tano, et
al. 2023); abstracts and partial introductions were analyzed for the other
papers listed below.

## Corpus consulted

1. **Walker, S.A., Tano, M.E., Abou-Jaoude, A., Calvin, O. (2023).**
   "Depletion-driven thermochemistry of molten salt reactors: review, method,
   and analysis." *Frontiers in Nuclear Engineering*, 10.3389/fnuen.2023.1214727.
   (open-access full text)
2. **Tano, M.E., et al. (2024).** "Coupled neutronics, thermochemistry,
   corrosion modeling and sensitivity analyses for isotopic evolution in
   molten salt reactors." *Annals of Nuclear Energy* (Elsevier),
   S0149197024004530. (abstract + OSTI record 2511255)
3. **Tano, M.E., Bajpai, P., de Oliveira, R., Parker-Quaife, E., Otis, K.,
   Yang, X. (2026).** "A methodology for decay heat characterization in
   molten salt reactors." *Nuclear Engineering and Design* 446, 114524.
4. **Huhn, Q.A., Tano, M.E., Ragusa, J.C., Choi, Y. (2023).** "Parametric
   dynamic mode decomposition for reduced order modeling." *J. Computational
   Physics*. (arXiv:2204.12006 preprint)
5. **German, P., Tano, M., Fiorina, C., Ragusa, J.C. (2022).** "GeN-ROM — An
   OpenFOAM-based multiphysics reduced-order modeling framework for the
   analysis of Molten Salt Reactors." *Progress in Nuclear Energy* 146.
6. **Abou-Jaoude, A., ..., Tano, M., Freile, R. (2021).** "A workflow
   leveraging MOOSE transient multiphysics simulations to evaluate the
   impact of thermophysical property uncertainties on molten-salt
   reactors." *Annals of Nuclear Energy* 163. (abstract)
7. **Novak, A., ..., Tano, M., et al. (2021).** "Pronghorn: A
   Multidimensional Coarse-Mesh Application for Advanced Reactor Thermal
   Hydraulics." *Nuclear Technology* (Taylor & Francis). (abstract +
   ResearchGate metadata)

Several additional papers were identifiable from his ResearchGate and
Google Scholar pages (52+ journal articles, h-index visible in
profile) but their full text was paywalled. The style features below are
robust across abstracts, introductions, and the one open-access full text
that could be read in full.

---

## 1. Sentence-level style

### Average sentence length

- **Moderate**, typically 18–28 words. Tano writes longer sentences than the
  median engineering author but shorter than mathematical-statistics authors.
- **Frequent short topic sentences (8–14 words)** open paragraphs and section
  bodies, then are followed by longer compound sentences with subordinate
  clauses.
- Example (Walker, Tano et al. 2023): *"Molten salt reactors (MSRs) are
  innovative advanced nuclear reactors that utilize nuclear fuel by dissolving
  it in a high-temperature liquid salt."* (~22 words, opens the abstract.)
- Followed by: *"This unique feature differentiates MSRs from other types of
  reactors and allows for enhanced safety and economic performance."* (~19
  words.)

**Actionable rule:** Open paragraphs with one ≤15-word declarative sentence.
Follow with one or two longer sentences. Avoid stacking three long sentences
in a row.

### Active vs. passive

- **Mixed, leaning active for contribution sentences, passive for method
  description.**
- Active: *"We construct...", "We deploy the calibrated model..."*
- Passive for procedures: *"The DMD low-rank or reduced operator is typically
  obtained by singular value decomposition..."* (Huhn, Tano et al. 2023).
- Passive when the agent is institutional or obvious: *"Each is fitted by its
  native sampler..."* is the kind of construction Tano uses.

**Actionable rule:** Use active "we" for novelty claims and decisions;
passive for "the equation is solved by", "the data are stored as".
Do not stack passives.

### Hedging

- **Confident but qualified.** Tano avoids both extremes.
- He does NOT write "we prove", "we show", "it is clear", "remarkable",
  "elegant", "novel" except where literally warranted (a theorem is proved).
- He DOES use: *"is expected to"*, *"can be"*, *"may", *"anticipated that"*,
  *"is a fundamental step toward"*, *"is expected to have potential
  application for"*.
- Example: *"Due to these thermochemical data shortcomings, it is anticipated
  that the multiphysics depletion-driven thermochemical analysis will not
  capture the evolving redox potential in an absolute sense."* (Walker, Tano
  et al. 2023) — note the explicit prospective qualifier.
- *"The neutronics-chemistry coupling developed in this work is expected to
  have potential application for analyzing corrosion, source term evolution,
  and material safeguards in MSR systems."*

**Actionable rule:** When stating a contribution, use *"this work represents a
fundamental step toward..."*, *"is expected to have application for..."*,
*"provides a framework for..."* rather than *"we definitively show..."* or
*"this is a major advance"*.

### Math ↔ prose transitions

- Equations are introduced with **one short prose sentence** stating what the
  equation expresses, then displayed, then followed by a "where" clause
  defining variables, then optionally one sentence about the consequence.
- Tano does **not** write "Eq. (X) shows that..." or "It is remarkable
  that..." or "It is worth noting that...".
- Tano's pattern: *"The mass-action assumption converts the network into a
  system of ordinary differential equations on the species concentrations,
  [equation], where C is the vector of species concentrations,..."*
- He generally does NOT box equations or use \boxed{} except for one or two
  central results. (His radiolysis manuscripts and the Bayes-rule pedagogy in
  paper1 boxing Bayes' rule are an outlier driven by the tutorial mode.)

---

## 2. Paragraph and section style

### Paragraph length

- **Typically 4–7 sentences.** He rarely writes a one-sentence paragraph
  except as a transition (e.g., *"This paper closes that gap."*).
- He almost never writes a 10+ sentence paragraph.

### Paragraph internal flow

- **Topic sentence → physical motivation → quantitative claim or equation →
  caveat / hedge.**
- He does NOT lead with the hedge ("It should be noted that..."). He leads
  with the claim, then qualifies.
- Example pattern (Walker/Tano abstract): claim ("Molten salt reactors are
  innovative..."), elaboration ("This unique feature differentiates..."),
  consequence ("The liquid fuel also entails..."), specific instantiation
  ("One primary effect..."), definition ("Essentially, depletion-driven
  thermochemistry is...").

### Section openers

- **Brief framing sentence**, often 1–2 sentences, then dives into content.
- He does NOT use florid "Imagine a world where..." or "The challenges
  facing the nuclear industry..." openers.
- For Methods/Background sections, the opener typically references the prior
  section or states the section's scope: *"This section builds the necessary
  background from first principles, using molten-chloride radiolysis as the
  running example."* (Note: paper1's pedagogical voice goes a step further
  than Tano typical; the bracketed "skip to..." phrasing is more textbook
  than Tano's published journal voice.)

### Introduction structure

- **Three-move pattern, not four-move.**
- Move 1: state the engineering / scientific stakes (concept studies for MSRs,
  active commercial development).
- Move 2: review what is and is not yet established in the literature
  (kinetics data, prior frameworks).
- Move 3: name the gap, state the contribution.
- He typically does NOT structure as Situation/Complication/Resolution/
  Contribution in the explicit consultant style. The "complication" is woven
  into Move 2 as a literature-review observation: "X has been done, but Y is
  missing." He does NOT write "however, despite this progress, a fundamental
  challenge remains".
- Bullet-list contribution paragraphs ("The contributions of this paper are:
  1. ... 2. ...") are common — Tano uses these in both the depletion-
  thermochemistry paper and the coupled-MSR paper. **Always 4–6 bullets, not
  10+.**

---

## 3. Vocabulary and idiom

### Characteristic phrasing

Tano-leaning phrases (use these):

- *"this work represents a fundamental step toward..."*
- *"is expected to have potential application for..."*
- *"provides a framework for..."*
- *"is anticipated that..."*
- *"defensible predictions of..."*
- *"a primary effect termed here as..."*
- *"this multiphysics coupling is coined by this work as..."*
- *"the present work..."* (he uses *"this present work"* as well, with the
  redundant *"this present"* — not idiomatic English but recurring in his
  papers).
- *"essentially,"* (used to define a technical concept in plainer terms).
- *"complementing..."*, *"building on..."*, *"in line with..."*
- *"order of magnitude"* and *"orders of magnitude"* (common quantitative
  qualifier).

### Phrasing to AVOID (anti-Tano)

These are common LLM/AI-generated phrasings absent from Tano's papers:

- *"In summary,"* / *"In conclusion,"* as standalone paragraph openers.
- *"It is worth noting that..."* / *"It should be noted that..."*
- *"Strikingly,"* / *"Remarkably,"* / *"Crucially,"* / *"Importantly,"*
- *"This paper makes the following novel contributions..."* (he writes "The
  contributions of this paper are:" — no "novel").
- *"a comprehensive understanding of..."* / *"holistic"*.
- *"sheds light on"* / *"unprecedented"* / *"paramount"*.
- *"To the best of our knowledge"* / *"To our knowledge,"* (he does not use
  this; in his work the phrasing is more often *"No published framework..."*
  or *"no previous compilation..."*).
- *"a deep dive into..."*, *"unlock the potential of..."*.
- *"In this paper, we propose a novel method that..."* in the first sentence
  of the abstract. (Tano's abstract first sentences are descriptive of the
  problem, not the contribution.)

### Self-reference

- **First-person plural "we"** is the default, even for single-author papers.
- *"This work"* / *"the present work"* used for the paper as object.
- He does NOT use *"the author"* / *"the authors"* / *"I"*.

### Acronyms

- **Define on first use in the abstract and again on first use in the body.**
- After definition, used freely. Common acronyms (MSR, MOOSE, OpenFOAM,
  AIMD) are sometimes used without definition in the abstract if the journal
  is specialized.
- Tano does NOT introduce more than ~5 new acronyms in a single paper.
- He prefers writing out short names ("Thermochimica", "Griffin") over
  acronyms when the name is short enough.

### Hedging on results

- **Quantitative, not rhetorical.** Where a result is uncertain, Tano gives
  the credible interval or sigma, not adjectives.
- Example: *"recovers published Iwamatsu Arrhenius parameters within
  reported sigma..."* — the unit of hedge is the sigma, not "approximately"
  or "roughly".

---

## 4. Figures, tables, captions

### Caption style

- **Descriptive, not declarative-with-conclusion.** Tano's captions describe
  what the figure shows; the conclusion goes in the body text.
- Caption sentence length: 1–4 sentences, totaling 20–80 words.
- Example pattern: First sentence names what is plotted. Second sentence
  defines what symbols/colors/lines encode. Third sentence (optional) flags
  a key feature without stating the take-away.
- He does NOT use "Figure 1 demonstrates that X is better than Y" in
  captions; that goes in body text.

### Number of panels

- **Multi-panel figures are common (2–6 panels labeled (a)–(f)).**
- Single-panel figures used for one-data-point demonstrations or schematic
  diagrams.

### Equation integration

- One short lead-in sentence ("The mass-action assumption gives..."),
  numbered displayed equation, "where" clause, optional consequence
  sentence.
- Multi-line equations grouped in `align` with `\nonumber` on internal lines
  and one equation number per logical equation.

### Tables

- **Tables are concise summaries of data, not vehicles for editorializing.**
- Caption gives the units, the meaning of each column, and the source where
  appropriate.
- Compare tabulated values to literature explicitly (a "Literature" column
  alongside the "Posterior" column, as in the radiolysis paper table 14
  format).

---

## 5. Citation style

- **Author-year clusters via natbib**, often clustered when surveying a
  literature.
- Cluster pattern: *"[Forsberg2020,Andreades2014,Holcomb2012]"* — three
  citations to back a single class of designs.
- Single-citation pattern when the citation supports a specific quantitative
  number: *"k_2 ≈ 10^10 M^-1 s^-1 [Iwamatsu2022]"*.
- Citations placed at the **end of the clause they support**, not
  parenthetically mid-sentence in the (Author, Year)-style except where the
  natbib `\citet` flow puts the author in the running text.
- He uses *"Iwamatsu and coworkers \citep{Iwamatsu2022}"* and
  *"Walker et al. \citep{Walker2023}"* roughly interchangeably; the
  "and coworkers" formulation occurs particularly when reviewing a research
  group's body of work rather than a single paper.

---

## 6. Nuclear-engineering / computational specifics

### Physical models

- **Formal but accessible.** A model is introduced with its physical
  motivation in plain English, then formalized.
- Conservation laws are stated as such; constitutive choices (closures,
  empirical correlations) are flagged with their source and uncertainty
  scope.
- Tano writes about MOOSE, Griffin, Pronghorn, Thermochimica, OpenFOAM,
  OpenMC etc. **as named tools** (proper nouns, sometimes typeset in
  text-mode); he does NOT genericize them as "the multiphysics framework".

### Limitations

- **Discussed in a dedicated Limitations subsection of the Discussion**,
  almost always before the Conclusion. Not buried in the Conclusion, not in
  the Methods.
- Limitations are itemized (numbered list of 3–6 items). Each limitation is
  named with a short italicized phrase ("Single-source MCFR constraint."),
  followed by 2–4 sentences of detail.
- He does NOT write "However, this work has some limitations..." in a
  hand-wave. He writes them as engineering constraints.

### Uncertainty and verification

- **Quantitative uncertainty is reported as posterior intervals or as
  variance contributions** when a Bayesian framework is used; as standard
  deviations or one-sigma bounds otherwise.
- He frequently performs and reports a **posterior-variance decomposition**
  identifying which input contributes how much of the output uncertainty
  ("X accounts for approximately 76% of the observed variance..."). This is
  characteristic.
- **Verification** is treated as a comparison to a published reference value
  or to a derived analytical solution, with the residual quantified.
  "Validation" is comparison to experiment.

---

## 7. Tone / register summary in one paragraph

Tano writes as a senior staff scientist describing computational engineering
work to a peer audience that is comfortable with multiphysics but not
necessarily with statistics or with the specific salt chemistry. He is
confident about what was done (active "we constructed", "we deployed"),
quantitatively modest about what was found ("recovers literature values
within sigma", "approximately X% of variance"), and careful to flag the
boundaries of applicability of his models. He avoids both the salesman
register ("revolutionary", "unprecedented", "novel") and the
overly-hedged register ("might possibly perhaps suggest"). He uses bullet
lists in Introductions to enumerate contributions, and uses italicized
short phrases as item leads in Discussion-section enumerations of strengths
and limitations. The writing has identifiable non-native-English fingerprints
("this present work", occasional missing article, "approach to carry out")
that survive editing, which paradoxically authenticates the voice.

---

## 8. Sandbox: rewrite these sentences in Mauricio's voice

### Source sentence A (overly promotional):
> "We introduce a groundbreaking new method that revolutionizes Bayesian
> calibration of chemical kinetic networks, achieving unprecedented accuracy
> and unlocking the potential of multi-laboratory data fusion."

**Tano-style rewrite:**
> We present a hierarchical Bayesian framework that calibrates chemical
> kinetic networks against multi-laboratory data. The framework is
> expected to be useful for problems where systematic biases between
> laboratories are not separately identifiable from the intrinsic
> chemistry; we demonstrate this on the LiCl-KCl chloride radiolysis
> corpus.

### Source sentence B (hedge-heavy):
> "It might be argued that, in principle, our results could potentially
> suggest that a deeper investigation into the underlying mechanisms may
> be warranted."

**Tano-style rewrite:**
> The residual model-experiment disagreement on the Phillips~2022 NULL
> benchmark is consistent with a missing slow chloride sink. Direct
> measurement of $k_{\Clt + \mathrm{U(III)}}$ at LEAF facility resolution
> would resolve this without further inference.

### Source sentence C (LLM-flavored conclusion):
> "In conclusion, our pioneering framework represents a major step forward,
> sheds light on previously poorly-understood phenomena, and paves the way
> for future research in this exciting field."

**Tano-style rewrite:**
> The framework recovers published Arrhenius parameters for chromium and
> zinc within $1\sigma$, quantifies the Pikaev--Iwamatsu inter-facility
> offset at $b^{(\mathrm{Pikaev})} = -3.35\,[-4.53,-2.44]$, and identifies
> three specific experiments that would tighten the posterior by an order
> of magnitude. The calibration is expected to support cover-gas screening
> calculations for chloride and fluoride MSR concepts.

### Source sentence D (florid Introduction opener):
> "Molten salt reactors represent a paradigm shift in nuclear energy
> generation that promises to address the pressing challenges of clean
> energy in the 21st century, captivating researchers and policymakers
> alike."

**Tano-style rewrite:**
> The molten salt reactor (MSR) concept has re-emerged in the past decade
> as a leading candidate for advanced fission heat generation, with active
> commercial development of fast-spectrum chloride-fueled designs and
> thermal-spectrum fluoride-cooled designs.

### Source sentence E (vague limitation):
> "As with all studies, there are some limitations that should be
> acknowledged."

**Tano-style rewrite (as opener of a Limitations subsection):**
> We are explicit about three limitations of the present model.

---

## 9. Quick checklist before submitting a Tano-voice paragraph

- [ ] No "Importantly", "Notably", "Crucially", "Strikingly" at the start.
- [ ] No "It is worth noting that".
- [ ] No "groundbreaking", "novel", "unprecedented", "paradigm shift".
- [ ] Abstract opens by describing the problem, not the contribution.
- [ ] Limitations are itemized in a dedicated Discussion subsection.
- [ ] Quantitative claims are reported with sigma or credible intervals.
- [ ] First-person plural "we"; never "the author" or "I".
- [ ] Acronyms defined on first use; ≤5 new acronyms total.
- [ ] Citations clustered (3+) when surveying a literature; singly when
      supporting a specific number.
- [ ] Contributions paragraph uses a numbered list of 4–6 items.
- [ ] Discussion section contains explicit Strengths and Limitations
      subsections before Conclusion.
- [ ] At least one posterior-variance / sensitivity decomposition table or
      paragraph where uncertainty is being discussed.

