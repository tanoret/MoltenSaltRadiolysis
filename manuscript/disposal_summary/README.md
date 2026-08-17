# Disposal-oriented summary

A 4-page working summary of the radiolysis program oriented toward **molten
salt fuel disposal** (sealed-container gas generation under long-timescale
α-dominant dose-rate regimes, no engineered U(III)/U(IV) buffer).

This summary is independent of the two journal manuscripts in
`../paper1/` (HBMAE methodology) and `../paper2/` (operational-reactor
predictions). It targets a different audience — disposal program managers
and waste-form designers — and uses three original figures generated
specifically for this document.

## Files

| File | Purpose |
|---|---|
| `disposal_summary.md` | Source markdown |
| `disposal_summary.docx` | Word document (compiled via pandoc with embedded PNG figures) |
| `figures/fig_disposal_dose_trajectory.{pdf,png}` | β/γ vs α dose-rate trajectory over $10^{-1}$–$10^6$ years; cross-over at ~150 y |
| `figures/fig_disposal_container_inventory.{pdf,png}` | Cl₂ and F₂ inventory in a sealed disposal container vs. time after discharge |
| `figures/fig_disposal_regime_map.{pdf,png}` | Operational vs. disposal radiolysis regimes on (T, dose-rate) plane; extrapolation gap |
| `scripts/build_disposal_figures.py` | Reproduces the three figures from the calibrated kernel |

## Build

```bash
# Regenerate figures from the calibrated kernel
python scripts/build_disposal_figures.py

# Regenerate the Word document
pandoc disposal_summary.md -o disposal_summary.docx
```

## Status

Standalone short summary. The three figures are original and do not appear
in either of the journal manuscripts. Posterior numerical estimates are
drawn from the integrated HBMAE chain in `../../validation/`.
