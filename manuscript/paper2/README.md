# Paper 2 — Multi-salt radiolysis model

A calibrated multi-salt radiolysis model for molten chloride and fluoride reactor fuels: HBMAE calibration of a chloride and a fluoride radiolysis kernel against the digitized pulse-radiolysis and chronic-irradiation corpus, with operational predictions for cover-gas Cl₂ in a NaCl-UCl₃ MCFR and F₂ in a FLiBe-UF₄ loop.

**Target journal:** Journal of Nuclear Materials.

## Files

| File | Purpose |
|---|---|
| `article.tex` | Manuscript source |
| `references.bib` | Cited references (40 entries) |
| `figures/fig_data_landscape.{pdf,png}` | Calibration-corpus map (1000/T vs log k) |
| `figures/fig_fluoride_kernel.{pdf,png}` | Davis G + Toth/Felker recombination, 4 panels |
| `figures/fig_master_validation.{pdf,png}` | 16-panel model-vs-data overlay |
| `figures/fig_meta_hier.{pdf,png}` | Chemistry-feature regression + LOMO |
| `figures/fig_residuals_all.{pdf,png}` | LOMO residuals scatter |
| `figures/fig_mcfr_cl2_lifetime.{pdf,png}` | NaCl-UCl₃ MCFR 60-y prediction |
| `figures/fig_flibe_f2_lifetime.{pdf,png}` | FLiBe-UF₄ F₂ inventory vs T |
| `cover_letter.md` | Editor cover letter |
| `article.pdf` | Compiled manuscript (37 pages) |

## Build

```bash
tectonic article.tex
```

## Companion paper

This application paper has a companion methodology paper at `../paper1/article.tex`, *Hierarchical Bayesian Mechanism-Adequacy Estimation: A Multi-Modality, Multi-Host Framework for Calibrating Chemical Kinetic Networks against Heterogeneous Experimental Data* (target *SIAM/ASA Journal on Uncertainty Quantification*). The framework's six theorems and the algorithmic details are in the companion paper; the present manuscript invokes the framework but stands alone for *Journal of Nuclear Materials* readers, with a one-page methodology summary in §1.4.

## Reproducibility

Calibrated parameters, the integrated MCMC chain, and the operational-MSR-prediction figures are all reproducible from `../../scripts/paper2_*.py` and `../../validation/`. The Zenodo-DOI release will accompany the journal submission.
