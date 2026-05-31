# Cover letter — Paper 1

**To:** Editorial office, *SIAM/ASA Journal on Uncertainty Quantification*

**Re:** *Hierarchical Bayesian Mechanism-Adequacy Estimation: A Multi-Modality, Multi-Host Framework for Calibrating Chemical Kinetic Networks against Heterogeneous Experimental Data.*

Dear Editor,

Please consider the attached manuscript for publication in *SIAM/ASA Journal on Uncertainty Quantification*.

The paper introduces Hierarchical Bayesian Mechanism-Adequacy Estimation (HBMAE), a Bayesian framework that calibrates chemical-kinetic networks against the heterogeneous datasets typical of radiation chemistry: pulse-radiolysis transients, scalar rate constants, derived Arrhenius pairs, and one-sided detection limits, each carrying its own facility-dependent systematic offset. The framework combines a constrained-topology prior, a hierarchical Arrhenius layer that pools shared reactions across host salts, a composite likelihood across observation modalities, a parametric Kennedy–O'Hagan discrepancy with bounded flexibility, and an explicit facility-effect hierarchy. Six theorems establish ergodicity of the constrained-topology sampler, cross-host Bernstein–von Mises with explicit host-heterogeneity floor, finite-prior discrepancy-bias control, Godambe sandwich asymptotics for the composite estimator, a censored Bayes-factor lower bound for one-sided benchmarks, and gauge identifiability of the facility hierarchy.

The novelty of the contribution is the synthesis: each layer has a precedent in the broader uncertainty-quantification literature, but the combined posterior with explicit identifiability and contraction theorems, and an empirical validation against the BNL LEAF pulse-radiolysis corpus plus the Phillips 2022 NaCl-UCl₃ chronic-irradiation NULL benchmark, is, to our reading, not in the existing literature.

The empirical results recover published Arrhenius parameters within the reported $\sigma$ on four chromium and zinc reactions, quantify a previously qualitative Pikaev–Iwamatsu inter-laboratory offset at $b^{(\mathrm{Pikaev})} = -3.35\,[-4.53,-2.44]$, and exclude the no-U-redox chloride kernel against the Phillips 2022 non-detection benchmark.

A companion application paper (in preparation for *Journal of Nuclear Materials*) uses the calibrated framework to predict operational cover-gas Cl₂ and F₂ inventories in MSR fuels.

We have no conflicts of interest to declare. The manuscript has not been submitted elsewhere.

Sincerely,

M. Tanore-Tamales\
Idaho National Laboratory\
mauricio.tanoretamales@inl.gov
