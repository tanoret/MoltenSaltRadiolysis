# Validation literature inventory

Status:
- **digitized** — quantitative data CSVs in `data/` and a manifest pointing to them
- **partial** — text values + computed Arrhenius extracted; figure-trace digitization still useful
- **mechanistic** — computational/theory paper; informs model structure but holds no experimental data
- **stubbed** — directory exists with citation; awaiting paper access or figure digitization
- **candidate** — found in literature, no folder yet

For each digitized/partial case, the manifest's `data/` directory holds some combination of:
- `Source.txt` — full citation + experimental conditions
- `reported_rate_constants.csv` — single-T rate constants reported in paper text
- `arrhenius_parameters.csv` — A, Ea, ΔH‡, ΔS‡ values
- `k_vs_T_from_arrhenius.csv` — per-T 2nd-order k from the published Arrhenius with propagated σ
- `vision_fig*.csv` — vision-extracted figure points (~5-10% precision; for sanity check)
- `compositions_table1.csv`, `experimental_conditions.csv`, etc. — paper-specific tables
- raw transient CSVs (where available)

---

## Chloride-melt pulse radiolysis (LiCl–KCl eutectic)

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **digitized** | Iwamatsu, Horne, Ramos-Ballesteros, …, Wishart. "Kinetics of radiation-induced Cr(II)/Cr(III) redox chemistry in molten LiCl–KCl eutectic." *PCCP* **2026**, 28, 2061. | [10.1039/D4CP04190A](https://doi.org/10.1039/D4CP04190A) | [cr_licl_kcl/iwamatsu_2026_pccp/](cr_licl_kcl/iwamatsu_2026_pccp/) — 9 transient CSVs + reported_rate_constants + arrhenius + k_vs_T + vision_fig2B/3B |
| **partial** | Iwamatsu, Horne, …, Wishart. "Radiation-induced reaction kinetics of Zn²⁺ with eₛ⁻ and Cl₂•⁻ in molten LiCl–KCl eutectic at 400–600 °C." *PCCP* **2022**, 24, 25088. | [10.1039/D2CP01194H](https://doi.org/10.1039/D2CP01194H) | [zn_licl_kcl/horne_2022_pccp/](zn_licl_kcl/horne_2022_pccp/) — reported + arrhenius + k_vs_T + vision_fig4B (no transients yet) |
| **partial** | Conrad, Cook, et al. "Impact of iodide ions on the speciation of radiolytic transients in molten LiCl–KCl eutectic." *PCCP* **2023**. | [10.1039/D3CP01477K](https://doi.org/10.1039/D3CP01477K) | [i_licl_kcl/horne_2023_pccp/](i_licl_kcl/horne_2023_pccp/) — compositions_table1 + reported_kinetics + vision_fig5 (SI Table S1 still needed) |
| **stubbed** | Iwamatsu, Horne, et al. "Influence of Nd(II)/Nd(III) on radiolytic transients in LiCl–KCl." *Inorg. Chem.* **2026**. | [10.1021/acs.inorgchem.5c04788](https://doi.org/10.1021/acs.inorgchem.5c04788) | [nd_licl_kcl/horne_2026_ic/](nd_licl_kcl/horne_2026_ic/) — k(eₛ⁻+Nd³⁺) = 4.54e10, k(Nd²⁺+Cl₂•⁻) = 1.72e10 at 673 K (text only; no OSTI preprint accessible) |
| **stubbed** | Hagiwara, Sawamura, Sumiyoshi, Katayama. "Pulse radiolysis study of transient species in LiCl–KCl melt." *Rad. Phys. Chem.* **1987**, 30(2), 143. | — | [licl_kcl_baseline/hagiwara_1987_rpc/](licl_kcl_baseline/hagiwara_1987_rpc/) — historical baseline values already embedded in Zn manifest |

## Chloride-melt pulse radiolysis (other / neat)

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **mechanistic** | Nguyen, Gibson, Emerson, et al. "Chlorine gas and anion radical reactivity in molten salts and the link to chlorobasicity." *PCCP* **2025**, 27, 4290. **(Computational MD — no rate data.)** | [10.1039/D4CP03285C](https://doi.org/10.1039/D4CP03285C) | [cl3_chlorobasicity/horne_2025_pccp/](cl3_chlorobasicity/horne_2025_pccp/) — qualitative Cl2/Cl3⁻ branching only |
| **mechanistic** | Nguyen, Bryantsev, Margulis. "Are high-temperature molten salts reactive with excess electrons? Case of ZnCl₂." *JPC B* **2023**, 127(42). **(Computational MD — no rate data.)** Author attribution corrected (was mis-listed as Gibson). | [10.1021/acs.jpcb.3c04210](https://doi.org/10.1021/acs.jpcb.3c04210) | [zncl2_neat/gibson_2023_jpcb/](zncl2_neat/gibson_2023_jpcb/) — three e⁻ states identified |
| **stubbed** | Horne, et al. "Chemistry of solvated electrons in molten alkali chloride salts." *JPC C* **2018**. | [10.1021/acs.jpcc.8b05716](https://doi.org/10.1021/acs.jpcc.8b05716) | [ealkali_chlorides/horne_2018_jpcc/](ealkali_chlorides/horne_2018_jpcc/) — **needed for e_s⁻ molar absorptivity ε(λ,T) in LiCl-KCl** |
| **stubbed** | Pikaev, Makarov, Zhukova. "Solvated electron in irradiated melts of alkaline halides." *Rad. Phys. Chem.* **1982**, 19, 377. | [10.1016/0146-5724(82)90005-X](https://doi.org/10.1016/0146-5724(82)90005-X) | [alkali_halide_review/pikaev_1982_rpc/](alkali_halide_review/pikaev_1982_rpc/) — companion Russian papers listed |
| **stubbed** | Makarov, Zhukova, Pikaev, Spitzyn. "Oxidizing agents produced by radiolysis of alkali-metal halide melts." *Bull. Acad. Sci. USSR Div. Chem. Sci.* **1982**, 31(4), 662. | [10.1007/BF00949993](https://doi.org/10.1007/BF00949993) | [oxidants_halide_melts/makarov_1982_rcb/](oxidants_halide_melts/makarov_1982_rcb/) |

## Bromide and iodide neat melts

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **stubbed** | Sawamura, Gebicki, Mayer, Kroh. "Pulse radiolysis of LiBr–KBr melts." *Rad. Phys. Chem.* **1990**, 46(2), 433. | — | [libr_kbr_pulse_rad/sawamura_1990_rpc/](libr_kbr_pulse_rad/sawamura_1990_rpc/) — needs bromide kernel addition |

## Fluoride melts (FLiBe, FLiNaK, MSRE fuel)

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **stubbed** | (Anon.) "Fluorine generation by gamma radiolysis of a fluoride salt mixture." *Rad. Eff. Def. Solids* **1990**, 112(4). LiF-BeF₂-ZrF₄-UF₄ (MSRE composition); **Ea = 39 kJ/mol** for F-recombination. | [10.1080/10420159008213046](https://doi.org/10.1080/10420159008213046) | [flibe_msre_f2/heron_1990_redds/](flibe_msre_f2/heron_1990_redds/) — **HIGH PRIORITY**: real MSRE-relevant F₂ generation kinetics |
| **stubbed** | Davis, Hania, Boomstra, …, Riedel. "Radiolytic Production of Fluorine Gas from MSR Relevant Fluoride Salts." *NSE* **2022**, 197(4), 633. G-values: LiF~0.004, BeF₂~0.009, ThF₄~0.021, FLiBe-UF₄~0.005. | [10.1080/00295639.2022.2129951](https://doi.org/10.1080/00295639.2022.2129951) | [eflibe_f2_yield/davis_2022_nse/](eflibe_f2_yield/davis_2022_nse/) |
| **stubbed** | Akiyama, Kitaichi, Fujiwara, Sawamura. "Short-lived species in pulse-irradiated LiF-KF and LiF-NaF-KF melts." *J. Nucl. Sci. Technol.* **1994**, 31(3), 250. | — | [flinak_pulse_rad/akiyama_1994_jnst/](flinak_pulse_rad/akiyama_1994_jnst/) — **THE fluoride pulse-rad reference** |
| **stubbed** | Haubenreich, Williams, Icenhour. ORNL-MSR-69-46 series. MSRE fuel-salt storage F₂ release < 150 °C. | — | [flibe_msre/haubenreich_msr_69_46/](flibe_msre/haubenreich_msr_69_46/) |

## Chloride fast-spectrum / MCFR

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **partial** | **Phillips**, Cao, Warmann, Mohr, Lovel, Core. "Gamma Irradiation of NaCl-UCl₃ Salt for the Molten Chloride Fast Reactor." INL/RPT-22-66727 (2022). **(Author attribution corrected from earlier "Karlsson".)** Conclusion: Cl₂ < detection limit (~1000 ppm) at all 4 capsules (75–600 °C, 31 MGy total dose). | [OSTI 1874817](https://www.osti.gov/biblio/1874817) | [mcfr_uc13_irradiation/karlsson_2022_inl/](mcfr_uc13_irradiation/karlsson_2022_inl/) — **NULL benchmark** with experimental_conditions.csv |
| **stubbed** | Ramos-Ballesteros et al. "Radiation-induced long-lived transients and metal-particle formation in solid KCl–MgCl₂." *JPC C* **2022**, 126, 9820. | [10.1021/acs.jpcc.2c01725](https://doi.org/10.1021/acs.jpcc.2c01725) | [kcl_mgcl2_solid/ramos_ballesteros_2022_jpcc/](kcl_mgcl2_solid/ramos_ballesteros_2022_jpcc/) — solid-phase EPR characterization |

## Spectroscopic calibration (enables absorbance → concentration)

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **digitized** | Moon, J.; Chidambaram, D. "Near-infrared spectra and molar absorption coefficients of trivalent lanthanides in molten LiCl-KCl eutectic." *Prog. Nucl. Energy* **2022**, 152, 104375. **(Authorship corrected from earlier "Lazaridis".)** | [10.1016/j.pnucene.2022.104375](https://doi.org/10.1016/j.pnucene.2022.104375) | [lanthanide_epsilon_licl_kcl/moon_chidambaram_2022_pnse/](lanthanide_epsilon_licl_kcl/moon_chidambaram_2022_pnse/) — 20 ε values for Sm³⁺, Nd³⁺, Dy³⁺ at 500 °C with σ |

## Cross-cutting / EFRC

| Status | Paper | DOI | Local |
|--------|-------|-----|-------|
| **candidate** | Horne, Dias, et al. "Radiation-Assisted Formation of Metal Nanoparticles in Molten Salts." *JPC Lett.* **2020**, 12, 157. | [10.1021/acs.jpclett.0c03231](https://doi.org/10.1021/acs.jpclett.0c03231) | — |
| **candidate** | Phillips et al. "High-temperature furnace and cell holder for in situ spectroscopic, electrochemical, and radiolytic investigations of molten salts." *Rev. Sci. Instrum.* **2020**, 91, 083105. | [10.1063/1.5140463](https://doi.org/10.1063/1.5140463) | — instrumentation |
| **candidate** | Roy et al. "X-ray scattering reveals ion clustering of dilute chromium species in molten chloride medium." *Chem. Sci.* **2021**, 12, 8026. | [10.1039/d1sc01224j](https://doi.org/10.1039/d1sc01224j) | — |
| **candidate** | Gill et al. "Speciation and solubility of Ni(II) and Co(II) in molten ZnCl₂." *JPC B* **2020**, 124, 1253. | [10.1021/acs.jpcb.0c00195](https://doi.org/10.1021/acs.jpcb.0c00195) | — |
| **candidate** | Wishart (PI). Molten Salts in Extreme Environments EFRC technical summary, 2024. | [BES summary](https://science.osti.gov/-/media/bes/efrc/pdf/technical-summaries/2025/MSEE-Wishart-BNL-Technical-Summary-202410.pdf) | — |

---

## Inventory summary

| Bucket | Count | Status |
|--------|-------|--------|
| Fully digitized (quantitative CSVs) | 3 | Cr (Iwamatsu 2026), MCFR (Phillips 2022), lanthanide ε (Moon 2022) |
| Partial digitized | 2 | Zn (Iwamatsu 2022), iodide (Conrad 2023) |
| Mechanistic / theory-only | 2 | Both Nguyen MD papers; no experimental data |
| Stubbed (citation captured, data not yet extracted) | 10 | Nd, Hagiwara, Pikaev, Makarov, Sawamura LiBr, Davis F₂, Akiyama FLiNaK, Heron 1990 F₂, Haubenreich MSRE, Ramos-Ballesteros |
| Candidate (no folder yet) | 5 | Horne nano, Phillips instr, Roy XRD, Gill ZnCl₂, Wishart EFRC |
| **Total cases tracked** | **22** | |

## Most impactful next digitization targets (ranked)

1. **Horne 2018 JPC C** — gives e_s⁻ molar absorptivity ε(λ,T) in LiCl-KCl. Unlocks absolute concentration overlays for every Iwamatsu/Conrad case. Currently the harness uses `absorbance_scale_free` because ε is unknown.
2. **Akiyama 1994 JNST** — only available fluoride-melt pulse-radiolysis paper. The fluoride kernel in `database.yaml` currently has placeholder G-values and 2 toy reactions.
3. **Heron 1990 RedDS** — MSRE-composition F₂ generation; Ea = 39 kJ/mol already extracted from abstract. Full PDF would give A and yield rates.
4. **Davis 2022 NSE** — concrete G-values for F₂ production from solid fluoride salts (LiF, BeF₂, UF₄, ThF₄, FLiBe-UF₄).
5. **Conrad 2023 SI Table S1** — pseudo-1st-order eS⁻ decay rate vs [KI] across 400-700 °C (validates the iodine extension once added to the chloride kernel).

## Bayesian-readiness checklist

For each "digitized" or "partial" case, the data structures support a Bayesian model-adequacy workflow:

- **Likelihood**: pair `harness.run_case()` output with each trace's experimental observable (Gaussian σ ~ experimental noise).
- **Priors**: use Arrhenius A and Ea from `arrhenius_parameters.csv` as informative Gaussian priors with reported σ.
- **Posterior**: marginalize over A, Ea, and any nuisance parameters (pulse dose, impurity scavenging baseline).
- **Posterior predictive vs `k_vs_T_from_arrhenius.csv`** — directly comparable per-T 2nd-order rate constants.
- **Posterior predictive vs `vision_fig*.csv`** — qualitative overlay sanity check.

The MCFR Phillips 2022 case is a posterior-predictive NULL test: the model must not predict
detectable Cl₂ under the documented dose/T conditions.

## Caveats / known limitations

- **Vision-extracted CSVs** carry ~5-10 % precision vs ~few percent experimental σ. Prefer
  `k_vs_T_from_arrhenius.csv` for inference; use `vision_fig*.csv` as overlay sanity checks.
- **SI tables** for the modern Iwamatsu/Conrad papers contain per-T rate constant values that are NOT
  in the OSTI preprints. For full-precision per-T data, fetch the publisher SI.
- **Paywalled papers**: Akiyama 1994 JNST, Heron 1990 RedDS, Davis 2022 NSE, Iwamatsu 2026 Inorg Chem.
  Library access at INL should unlock all four.
- **Author attributions** previously had errors: Phillips et al. ≠ Karlsson; Moon & Chidambaram ≠ Lazaridis;
  Nguyen et al. ≠ Gibson for the ZnCl₂ MD paper. All corrected above.

## Still to search

- Hu & Steinberg–era 1960s-70s NaCl/KCl gamma irradiation literature.
- Other 1990s-era Sawamura/Akiyama follow-ups not in the Conrad 2023 reference list.
- Actinide-bearing molten chlorides (UCl₃, PuCl₃) under irradiation beyond Phillips 2022.
- Modern (2024-2026) MSEE EFRC outputs as they're released — Wishart group prolific.
