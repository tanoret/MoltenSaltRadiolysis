# Submission packages

This `manuscript/` directory contains two independent journal articles, each
self-contained in its own subdirectory with sources, figures, bibliography,
cover letter, and Makefile.

| Folder | Title | Target journal | Pages |
|---|---|---|---|
| [`paper1/`](paper1/) | Hierarchical Bayesian Mechanism-Adequacy Estimation: A Multi-Modality, Multi-Host Framework for Calibrating Chemical Kinetic Networks against Heterogeneous Experimental Data | *SIAM/ASA Journal on Uncertainty Quantification* | 66 |
| [`paper2/`](paper2/) | A Calibrated Multi-Salt Radiolysis Model for Molten Chloride and Fluoride Reactor Fuels | *Journal of Nuclear Materials* | 37 |

The two papers are companions, citing each other where the methodology /
application division benefits the reader. Each can be reviewed independently.

## Build instructions

In either folder:

```bash
cd paper1     # or paper2
make          # produces article.pdf via tectonic
```

The Makefile cleans up auxiliary files (`make clean`) and rebuilds the PDF
from `article.tex` + `references.bib` + `figures/`. No external scripts or
data are needed at build time; the figure PDFs are pre-generated and shipped
in `figures/`.

Build dependency: `tectonic` (recommended) or any modern TeX Live distribution
with `pdflatex` + `biber`/`bibtex`. The Makefile defaults to `tectonic`; set
`ENGINE=pdflatex` to switch.

## Supporting materials

Pre-generated supporting artifacts kept at the project root level:

- `../validation/` — digitized experimental data CSVs, posterior chains,
  validation tables (parent project repository).
- `../scripts/` — Python scripts that produced the calibrated parameters,
  validation tables, and the operational-MSR forward predictions.
- `../REVIEW_PAPER1.md`, `../REVIEW_PAPER2.md` — self-administered referee
  reports used to revise the manuscripts; kept for review-trail traceability.
- `../MAURICIO_STYLE_GUIDE.md`, `../STYLE_DIFF.md` — writing-style guides
  used to align the manuscripts with the author's published voice.

At submission, the supporting materials should be uploaded to Zenodo (or the
journal's data-availability repository) and the manuscripts updated with the
DOI in place of the `[repository URL to be provided at submission]`
placeholder. Both manuscripts already contain a "Code and data availability"
statement that links to this Zenodo DOI.

## Author

Mauricio E. Tanore-Tamales, Idaho National Laboratory.
`mauricio.tanoretamales@inl.gov`.
