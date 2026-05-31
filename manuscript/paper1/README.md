# Paper 1 — HBMAE methodology

Hierarchical Bayesian Mechanism-Adequacy Estimation for chemical-kinetic networks: a multi-modality, multi-host framework calibrating against heterogeneous experimental data, with six theorems and a validation case on molten-salt pulse radiolysis.

**Target journal:** SIAM/ASA Journal on Uncertainty Quantification.

## Files

| File | Purpose |
|---|---|
| `article.tex` | Manuscript source |
| `references.bib` | Cited references (53 entries) |
| `figures/fig_method_comparison.{pdf,png}` | Comparator-method benchmark figure |
| `figures/fig_worked1_corner.{pdf,png}` | Cr Tier 2 posterior corner plot |
| `cover_letter.md` | Editor cover letter |
| `article.pdf` | Compiled manuscript (66 pages) |

## Build

```bash
tectonic article.tex
```

The build needs only TeX Live (or tectonic) and the included `references.bib`. No external scripts or data are required to produce the PDF.

## Companion paper

This methodology paper has a companion application paper: `../paper2/article.tex`, *A calibrated multi-salt radiolysis model for molten chloride and fluoride reactor fuels* (target *Journal of Nuclear Materials*). The two papers are independent and stand alone, with the application paper citing the present methodology paper (`Tanore2026PaperI`) where the HBMAE framework is invoked.

## Reproducibility

All numerical artifacts in this paper (MCMC chains, posterior summaries, validation tables, benchmark synthetic-data results) are derived from scripts and data in the parent project repository at `../../` (`scripts/tier*.py`, `validation/`). The Zenodo-DOI release will accompany the journal submission; the manuscript references the placeholder `[repository URL to be provided at submission]` which must be replaced before mailing.
