---
title: "Predicting Radiolytic Gas Generation in Disposed Molten Salt Reactor Fuels"
subtitle: "A Calibrated Bayesian Framework for Long-Timescale Container-Safety Assessment"
author:
  - "Mauricio E. Tanore-Tamales"
  - "Idaho National Laboratory, Idaho Falls, ID 83415, USA"
date: "May 2026"
geometry: margin=1in
fontsize: 11pt
mainfont: "Calibri"
---

# Predicting Radiolytic Gas Generation in Disposed Molten Salt Reactor Fuels

**Mauricio E. Tanore-Tamales — Idaho National Laboratory** — *Working summary, May 2026.*

## 1. The disposal problem

Molten salt reactors (MSRs) introduce a new spent-fuel form: the salt itself
is the fuel matrix. When a chloride or fluoride salt loop is taken out of
service, the salt inventory must be conditioned, packaged, and either stored
on site for the medium term or moved into a geological repository. Two
features distinguish this waste form from oxide pellets in zircaloy
cladding.

The first is that the radiolytic chemistry continues after discharge. The
β/γ activity from cesium-137 and strontium-90 deposits energy into the salt
for centuries; the α activity from plutonium-239, americium-241, curium-244,
and the higher actinides deposits energy into the salt for tens of
thousands of years. In a sealed disposal container neither energy source
can escape, and the radiolytic radicals it produces — solvated electrons,
chlorine atoms, fluorine atoms, halogen radical anions — drive
chain-propagation reactions whose end products are molecular Cl$_2$ or
F$_2$ gas. These gases pressurize the container, oxidize the wall, and
in the limiting case can rupture the package.

The second feature is that the U(III)/U(IV) redox buffer, which dominates
the active-reactor radiolytic mass balance and keeps the cover-gas Cl$_2$
inventory at the ppb level during normal operation, is not available
during disposal. The buffer is either consumed during operation,
deliberately extracted before disposal, or rendered unreactive by
solidification. The chain-propagation efficiency $\eta$ in the disposal
regime is therefore set by the residual chemistry of the host salt and
the long-lived fission-product / actinide cocktail, not by the engineered
U cycle.

A defensible prediction of Cl$_2$ or F$_2$ inventory in a sealed disposal
container, as a function of time after discharge, salt composition,
container temperature, and accumulated dose, is the central scientific
question this summary addresses. The present work is a transfer of the
calibrated radiolysis framework developed for active-reactor predictions
[1, 2] into the disposal regime.

![Specific dose rate as a function of time after discharge from a
representative MSR fuel salt. The β/γ contribution follows the
Cs-137/Sr-90 decay; the α contribution is dominated by Cm-244 in the
first decades and by Am-241 ingrowth and Pu-239 longer-lived activity
thereafter. The cross-over near ~150 years marks the transition from
γ-dominant to α-dominant radiolysis chemistry — a transition the
current pulse-radiolysis literature does not cover.](figures/fig_disposal_dose_trajectory.pdf){#fig:dose width=100%}

## 2. Approach: HBMAE-calibrated kernel extrapolated to disposal conditions

The framework underlying the predictions is the hierarchical Bayesian
mechanism-adequacy estimation (HBMAE) calibration of the chloride and
fluoride radiolysis kernels against the digitized pulse-radiolysis and
chronic-irradiation literature [1]. The calibrated kernels reproduce the
published Arrhenius parameters of the dominant Cr$^{3+}$, Zn$^{2+}$,
Nd$^{3+}$, Cf$^{3+}$ reduction reactions and capture the cross-laboratory
systematic offset of the Pikaev 1982 Soviet measurements relative to
the modern Brookhaven LEAF data [2]. With the kernels in hand the
operational MCFR cover-gas Cl$_2$ inventory prediction was carried out
in [2]; the present work extends the same posterior into the disposal
regime.

The extrapolation is not trivial. Pulse-radiolysis experiments are
conducted at salt-loop temperatures (400-600 $^\circ$C) and at dose rates
on the order of 1 kGy per pulse (instantaneous) or 1-10 kGy/h (chronic
gamma irradiation). Sealed disposal containers operate at 100-300
$^\circ$C and, after the first decade, at integrated dose rates that
decline through 1 $\mu$Gy/h and below. The kernel must therefore
extrapolate over four to five decades in dose rate, over more than 300
$^\circ$C in temperature, and from a γ-dominant to an α-dominant
energy-deposition regime. The first two extrapolations are amenable to
the Arrhenius and dose-rate scaling already encoded in the kernel; the
third is not, because the secondary-electron / radical inventory
produced by α tracks differs from the inventory produced by γ tracks
through linear-energy-transfer (LET) effects.

Figure&nbsp;[@fig:regime] situates the operational, laboratory, and
disposal regimes on a (temperature, dose-rate) plane, with the
extrapolation gap explicit. The figure makes precise the calibration
question: which posterior parameters are constrained by the existing
data within the operational box, and which require LET-corrected
α-radiolysis data or low-temperature low-dose-rate confirmation
experiments to be transferable to disposal conditions.

![Operational, laboratory, and disposal radiolysis regimes on a
(temperature, dose-rate) plane. The blue contour is a schematic of the
calibration-coverage posterior under the HBMAE framework: the central
region (high dose rate, $T \sim 700\,\mathrm{K}$) is anchored by Iwamatsu
Cr/Zn data, with the Pikaev cross-laboratory anchor at higher
temperatures. The operational MSR regime sits inside the calibration
envelope; the disposal regime (200 °C, $\sim 1$ $\mu$Gy/h after one
millennium) is four decades below the experimental floor in dose rate
and α-dominant rather than γ-dominant. The arrow indicates the
extrapolation that the calibrated kernel must cover for a defensible
disposal prediction.](figures/fig_disposal_regime_map.pdf){#fig:regime width=100%}

The Bayesian framework propagates this extrapolation explicitly. Each
posterior sample on the chloride or fluoride Arrhenius pair produces a
corresponding sample of the disposal-regime gas-generation rate; the
predictive band quoted below is the 5--95th percentile of that ensemble,
not a single nominal trajectory.

## 3. Findings

Two illustrative cases are presented in
Figure&nbsp;[@fig:inventory].

For a chloride fuel without active U(III)/U(IV) buffering — that is,
either a buffer-exhausted salt or a chloride waste-form in which U has
been separated — the predicted Cl$_2$ partial pressure in a 1 m$^3$
sealed disposal container holding two metric tonnes of salt and 100 L
of head-space reaches the 10$^5$ Pa container-design pressure at
approximately one to ten years after discharge in the posterior median,
and as early as approximately three months in the upper-90th-percentile
tail. The chain-propagation efficiency $\eta$ without U buffering is
inferred at the level of $10^{-1}$ from the unbuffered terms of the
chloride kernel; a 10 ppm-level scavenger (Te, S, O$^{2-}$ from
impurities or fission products) reduces this by a factor of $10^3$ or
better, but in any case the absence of the U cycle removes the buffer
that kept the operational MCFR prediction ten orders of magnitude below
the limit. The implication for disposal package design is that the
chloride waste form requires either an engineered redox buffer in the
solidified matrix, a venting design that releases Cl$_2$ to a getter,
or a hot-cell handling step that reduces the inventory by a controlled
fluorination.

For a fluoride fuel held at 200 $^\circ$C in the same sealed-container
geometry, the F$_2$ steady-state partial pressure is approximately
$10^3$ Pa median, with a 90% posterior band spanning $3\times 10^2$ Pa
to $3\times 10^3$ Pa, and is essentially time-independent over the
disposal interval because F$_2$ recombination at metal walls (Toth and
Felker 1990 activation energy 39 kJ/mol, calibrated pre-exponential
$A_\mathrm{rec} \approx 250$ h$^{-1}$) is fast on the disposal
timescale. The fluoride margin against the $10^5$ Pa container limit is
therefore approximately two orders of magnitude, comfortably below the
operational design pressure in the median but not in the upper-tail
posterior. The dominant uncertainty contributor is the
$6\times$-temperature extrapolation of the Toth-Felker recombination
Arrhenius from 423 K (where it is calibrated) to 473 K (the disposal
target). A confirmation experiment at $T \sim 200\,^\circ$C in a fluoride
salt analog would compress this uncertainty.

![Posterior-predictive gas inventory in a sealed 1 m$^3$ disposal
container holding 2 metric tonnes of (a) chloride fuel salt without
U(III)/U(IV) buffer and (b) fluoride fuel salt at 200 $^\circ$C, as a
function of time after discharge. Median predictions are heavy lines;
the 90% band reflects the joint posterior on $G_\mathrm{Cl^\bullet}$,
chain-propagation efficiency $\eta$, $G_\mathrm{F_2}$, and the
calibrated recombination kinetics. The chloride case crosses the
$10^5$ Pa container-design pressure at $\sim 1-10$ years in the
posterior median; the fluoride case sits at $\sim 10^3$ Pa indefinitely.
The figure assumes no engineered buffer and an integrated dose-rate
trajectory consistent with
Fig.&nbsp;[@fig:dose].](figures/fig_disposal_container_inventory.pdf){#fig:inventory width=100%}

The dose-rate cross-over at $\sim$ 150 years (Fig.&nbsp;[@fig:dose])
introduces a regime change in the chloride prediction: between one
year and 150 years the β/γ flux dominates and the calibrated kernel
applies essentially unmodified, while beyond 150 years the α activity
takes over and the linear-energy-transfer correction enters. The
correction is between a factor of one and a factor of five depending
on the radical-yield ratio assumed for high-LET tracks. The fluoride
prediction is less sensitive to the cross-over because the steady-state
balance is set by the recombination Arrhenius, not by the integrated
dose.

## 4. Implications for the disposal program

The principal conclusions for a fuel-salt disposal program are:

- **Chloride and fluoride disposal forms face qualitatively different
radiolytic-gas constraints.** Chloride salt without an engineered
buffer accumulates Cl$_2$ approximately linearly with integrated dose;
the standard $10^5$ Pa container design pressure is exceeded within
years in posterior median. Fluoride salt reaches a much lower
steady-state F$_2$ pressure set by metal-wall recombination
and remains close to two orders of magnitude below container limit.
A defensible disposal program for chloride fuel therefore requires
*either* engineered post-discharge reduction of U(IV) back to U(III)
to restore the redox buffer in the disposal matrix, *or* a venting
design coupled to a Cl$_2$ getter, *or* a solidification process that
converts the chloride waste form to a stable chemical species
(e.g.\ U-glass or apatite incorporation) before sealing the
container.

- **The dose-rate cross-over near 150 years matters.** Disposal-program
calculations that use only the prompt β/γ dose rate and assume
exponential decay to negligible levels will under-predict the
α-dominant tail. A geological-repository assessment must propagate
the actinide α-radiolysis explicitly for $10^4$-$10^5$ years; the
α/γ G-value ratio under track-structure radiolysis is the key
parameter that current data do not constrain.

- **The calibrated kernel is ready for screening-level disposal
predictions.** The HBMAE-calibrated chloride and fluoride kernels
recover published Arrhenius parameters within published $\sigma$
across the operational regime [1, 2], and the Bayesian posterior
propagates the disposal-regime predictions with explicit uncertainty.
The dominant uncertainty contributors are identified
(chain-propagation efficiency $\eta$ in chloride, the
6$\times$-temperature recombination-Arrhenius extrapolation in
fluoride, the LET-correction across the α/γ cross-over) and a
prioritized list of confirmation experiments is provided below.

- **Three confirmation experiments would tighten the disposal
predictions by more than an order of magnitude:** (i) a chloride-salt
gamma-irradiation at 200 $^\circ$C without a U buffer, to directly
measure $\eta$ in the absence of the redox cycle; (ii) a fluoride-salt
recombination-kinetics experiment in the 150-250 $^\circ$C window, to
confirm the Toth-Felker Arrhenius outside its calibration band; and
(iii) an α-emitter doped salt (Cm-244 spike) at the 100 $\mu$Gy/h
level, to provide the first LET-corrected radiolytic G-value in a
disposal-relevant matrix. Each can be performed within existing
INL hot-cell capability over a 12- to 24-month program.

## 5. Open work

Several extensions are outside the present screening framework and
constitute the next stages of the work. The first is the coupling of
the gas-phase radiolysis prediction to the container-corrosion
chemistry — Cl$_2$ does not merely pressurize the container, it also
attacks the inner wall through Cl + Fe or Cl + Ni surface reactions
whose rates set the cumulative loss of confinement over the disposal
timescale. The second is the explicit inclusion of fission-product
chemistry in the disposal-phase ODE: Cs, I, Te, Sr, and the
lanthanides are redox-active on the disposal-time scale, and their
participation in the Cl$_2$/F$_2$ mass balance is not currently
modeled. The third is the validation of the framework against an
independent dataset: the present HBMAE calibration uses the entire
digitized pulse-radiolysis and chronic-irradiation corpus, so an
honest hold-out for disposal-regime prediction requires either a
new experiment (per the confirmation experiments listed above) or
the digitization of an additional independent legacy dataset, of
which the most promising is the Akiyama 1994 LiF-NaF-KF
pulse-radiolysis data.

## References

[1] Tanore-Tamales, M. (2026). *Hierarchical Bayesian Mechanism-Adequacy
Estimation: A Multi-Modality, Multi-Host Framework for Calibrating
Chemical Kinetic Networks against Heterogeneous Experimental Data.*
Manuscript in preparation; target *SIAM/ASA Journal on Uncertainty
Quantification*.

[2] Tanore-Tamales, M. (2026). *A Calibrated Multi-Salt Radiolysis Model
for Molten Chloride and Fluoride Reactor Fuels.* Manuscript in
preparation; target *Journal of Nuclear Materials*.
