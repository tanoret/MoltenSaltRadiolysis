# Chloride pulse-radiolysis data audit — what's digitized, what's missing

## What we have digitized and use in the article

### Original chloride dataset (5 papers)

| # | Paper | DOI | What we extracted | Used in |
|---|---|---|---|---|
| 1 | Iwamatsu, Horne et al. *PCCP* **2022**, 24, 25088. **Zn²⁺ + e_s⁻ in LiCl-KCl, 400-600 °C** | [10.1039/D2CP01194H](https://doi.org/10.1039/D2CP01194H) | reported_rate_constants.csv, arrhenius_parameters.csv, vision_fig4B (20 points × 5 T) | §worked2, §val:tier3, §val:tier4 |
| 2 | Iwamatsu, Castro Baldivieso et al. *PCCP* **2026**, 28, 2061. **Cr(II)/Cr(III) + e_s⁻ in LiCl-KCl** | [10.1039/D4CP04190A](https://doi.org/10.1039/D4CP04190A) | 9 transient absorbance CSVs, k_vs_T_from_arrhenius.csv, vision_fig2B/3B, Cl₂•⁻ rates | §worked1, §val:tier2, §val:tier4 |
| 3 | Conrad, Cook et al. *PCCP* **2023**. **Iodide impurity in LiCl-KCl** | [10.1039/D3CP01477K](https://doi.org/10.1039/D3CP01477K) | compositions_table1.csv, reported_kinetics, vision_fig5 ICl•⁻ decay | §val:tier1 |
| 4 | Castro Baldivieso, Horne et al. *Inorg. Chem.* **2026**, 65, 1283. **Nd(III)/Nd(II) in LiCl-KCl** | [10.1021/acs.inorgchem.5c04788](https://doi.org/10.1021/acs.inorgchem.5c04788) | k_vs_T.csv (5 T points), Arrhenius (Ea=33.2±0.4 kJ/mol, A=1.71±0.13×10¹³) | newly added |
| 5 | Phillips et al. **INL/RPT-22-66727**, 2022. **NaCl-UCl₃ NULL benchmark** | [OSTI 1874817](https://www.osti.gov/biblio/1874817) | experimental_conditions.csv (4 capsules, 31 MGy, T = 75-600 °C, no Cl₂ detected) | §val:tier3, §val:tier4 |

### New papers digitized 2026-05-28 (4 papers — Tier A delivered)

| # | Paper | DOI | What we extracted | Used in |
|---|---|---|---|---|
| 6 | Pikaev, Makarov, Zhukova *Rad. Phys. Chem.* **1982**, 19, 377. **e_s⁻ + Zn²⁺/Cd²⁺/Tl⁺/Ag⁺/Ca²⁺/Sr²⁺/Ba²⁺ in NaCl, KCl, KBr, NaBr, KI at 800–850 °C** | [10.1016/0146-5724(82)90005-X](https://doi.org/10.1016/0146-5724(82)90005-X) | 7 CSVs: 18-row e_s⁻ spectra table, 20-row e_s⁻ rate table, 16-row X₂•⁻ spectra, 30-row X₂•⁻ rate table, Ea(Cd+Br₂•⁻ in LiBr-KBr)=25 kJ/mol | §val:tier3, §val:tier4, meta-hierarchical layer (multi-host) |
| 7 | Hagiwara, Sawai, Sumiyoshi, Katayama *Rad. Phys. Chem.* **1987**, 30, 141. **LiCl-KCl baseline; Cl₂•⁻ disproportionation Ea = 24 kJ/mol; e_s⁻ first-order decay Ea = 43 kJ/mol below 500 °C** | (no DOI; ISSN 0146-5724/87) | 5 CSVs: absorption_spectra, rate_constants (18 rows incl. figure-digitized Arrhenius), arrhenius_parameters (4 entries), lifetimes, experimental_conditions | §val:tier3 (historical cross-paper consistency) |
| 8 | Rotermund, Mezyk, Sperling, Beck, Wineinger, Cook, Albrecht-Schönzart, Horne *J. Phys. Chem. A* **2024**, 128, 590. **k(Cf³⁺+Cl₂•⁻) = (8.28±0.61)×10⁵ M⁻¹ s⁻¹ at 22 °C; k(Cf³⁺+SO₄•⁻) = (9.50±0.43)×10⁸** | [10.1021/acs.jpca.3c07404](https://doi.org/10.1021/acs.jpca.3c07404) | 7 CSVs: cf_rate_constants (5), aqueous_baseline (3), cl2m_kinetics (8), so4m_kinetics (5), transient_absorption (8), experimental_conditions (32), actinide_comparison (4) | meta-hierarchical layer (1st actinide), Phillips NULL extension |
| 9 | Kristoffersen & Metiu *J. Phys. Chem. C* **2018**, 122, 19603. **AIMD/DFT** ΔE for e_s⁻ + Ag⁺/H₂/CH₄/N₂ in 45LiCl, 45NaCl, 29LiCl-16NaCl. *(Was misattributed in prior audit as "Horne 2018"; this is theory, not the missing ε(λ,T) experimental paper)* | [10.1021/acs.jpcc.8b05716](https://doi.org/10.1021/acs.jpcc.8b05716) | 6 CSVs: reaction_energies (24 rows), pair_binding_energies (3), aimd_conditions (3), electronic_structure (3 melts), tddft_spectra (0 — paper does no TDDFT), theory_vs_experiment (9) | Thermodynamic-feasibility priors on redox half-reactions; bipolaron self-association prior |

**Total**: 9 chloride papers digitized — 5 metals (Cr, Zn, Nd, Cf, plus U as NULL) + 4 alkaline earth (Ca, Sr, Ba, Cd) + Tl, Ag from Pikaev + 1 anion impurity (I⁻) + 9 fully digitized transient traces + 24 AIMD ΔE values across 3 chloride compositions + cross-host data spanning LiCl-KCl / NaCl / KCl / KBr / NaBr / KI.

### Attribution correction

In the previous version of this audit I labelled the DOI [10.1021/acs.jpcc.8b05716](https://doi.org/10.1021/acs.jpcc.8b05716) as "Horne et al. 2018 — e_s⁻ molar absorptivity ε(λ,T)". On reading the PDF, this DOI is in fact **Kristoffersen & Metiu 2018 — an AIMD/DFT theory paper**, *not* the experimental absorptivity paper I claimed. The actual ε(λ,T) experimental paper is still missing (see Tier A1 below, updated).

---

## What is missing and worth fetching

Ranked by impact on the article. PDFs would let me extract numerical data via the same `pdftotext + manual extraction` workflow I used for the Iwamatsu/Davis/Toth-Felker papers.

### Tier A — RESOLVED on 2026-05-28 with attribution correction

I previously assumed the "ε = 8000 M⁻¹ cm⁻¹" value used by Iwamatsu 2022 was for the solvated electron e_s⁻. Direct reading of *both* the Iwamatsu 2022 PDF and the Makarov 1982 PDF reveals this was wrong:

- **What is actually anchored**: ε(Cl₂•⁻) in molten alkali chlorides.
  - Makarov 1982 (now digitized, [validation/oxidants_halide_melts/makarov_1982_bull/data/x2m_molar_absorptivity.csv](validation/oxidants_halide_melts/makarov_1982_bull/data/x2m_molar_absorptivity.csv)): ε(Cl₂•⁻) = 7400, 6800, 6500, 5800, 5800 M⁻¹ cm⁻¹ in LiCl, NaCl, KCl, RbCl, CsCl respectively, near melting.
  - Iwamatsu 2022 uses ε(Cl₂•⁻, 340 nm) = 8000 M⁻¹ cm⁻¹ — *estimated* from Hug~1981 aqueous reference (8800) and Makarov's spectral-width scaling.
- **What remains genuinely missing**: ε(e_s⁻, 700 nm, molten chloride). Iwamatsu 2022 does not quote this; Makarov 1982 does not measure it; the e_s⁻ molar absorptivity in molten chloride is reported only in earlier Russian work cited by Makarov (Pikaev et al. *Radiat. Eff.* 22, 71 (1974); *Dokl. Akad. Nauk SSSR* 225, 1103 (1975); 261, 409 (1981)).
- **Operational consequence for our model**: ε(Cl₂•⁻) is now sufficient for Paper 2's absolute-concentration predictions of Cl₂ gas inventory in operational MSRs (the rate-limiting species for Cl₂ buildup). The remaining ε(e_s⁻) gap matters only for absolute [e_s⁻] in pulse-radiolysis transient analysis, where the scale-free likelihood (already in Tier 2/4) is mathematically equivalent under fixed prior on ε(e_s⁻).

### Tier B — useful but lower priority

| # | Paper | DOI | What we'd extract | Why |
|---|---|---|---|---|
| B1 | Makarov, Zhukova, Pikaev, Spitzyn. *Bull. Acad. Sci. USSR Div. Chem. Sci.* **1982**, 31, 662. "Oxidizing agents produced by radiolysis of alkali-metal halide melts." | [10.1007/BF00949993](https://doi.org/10.1007/BF00949993) | X₂ / X₃⁻ yield data across NaCl, KCl, NaBr, KI | Complement to A3 (same Pikaev group; same NaCl/KCl matrix); supplies X₂•⁻ disproportionation rates. |
| B2 | Conrad et al. *PCCP* **2023** "Impact of lanthanide ion complexation and temperature on TODGA…" | [10.1039/D3CP01119D](https://doi.org/10.1039/D3CP01119D) | k(lanthanide + dodecane•⁺) in organic; not directly molten-salt | Less directly relevant — organic solvent, but same Conrad team. Possibly informative for the chemistry-feature regression on lanthanides. |
| B3 | Horne, Dias et al. *JPC Letters* **2020**, 12, 157. "Radiation-Assisted Formation of Metal Nanoparticles in Molten Salts." | [10.1021/acs.jpclett.0c03231](https://doi.org/10.1021/acs.jpclett.0c03231) | Cluster/nanoparticle formation rates from e_s⁻ + multivalent metal | Mostly mechanistic; would add context to §10 Discussion on Cu/Ag/Au nanoparticle-relevant chemistry. |
| B4 | Ramos-Ballesteros et al. *JPC C* **2022**, 126, 9820. "Radiation-induced long-lived transients and metal-particle formation in solid KCl–MgCl₂." | [10.1021/acs.jpcc.2c01725](https://doi.org/10.1021/acs.jpcc.2c01725) | Solid-phase EPR data | Solid-phase, not directly transferable to molten kinetics; supplies long-lived radical species inventory. |
| B5 | Sawamura et al. *Rad. Phys. Chem.* **1990**, 36, 451. "Pulse radiolysis of LiBr-KBr melts." | No DOI (Vol 36 not 46 as I had previously); J-STAGE may host | Bromide-cousin e_s⁻ + Br₂•⁻ rates | Extends the halide series for the hierarchical Arrhenius layer. |
| B6 | Akiyama et al. *J. Nucl. Sci. Technol.* **1994**, 31(3), 250. "Short-lived species in pulse-irradiated LiF-KF and FLiNaK eutectic" | No DOI; J-STAGE permanent link [https://www.jstage.jst.go.jp/article/jnst1964/31/3/31_3_250/_article](https://www.jstage.jst.go.jp/article/jnst1964/31/3/31_3_250/_article) | Fluoride-melt pulse-rad transient kinetics | Fluoride cousin to the chloride data; combines with Davis 2022 and Toth 1990 (already digitized) to complete the fluoride kernel. |

---

## Status after the 2026-05-28 update

**Delivered today**: Pikaev 1982 RPC, Hagiwara 1987 RPC, Rotermund 2024 JPC A (Cf), Kristoffersen-Metiu 2018 JPC C (theory, was misattributed as Horne 2018 in prior audit).

**Resulting dataset**:
- **5 metals with kinetic data in chloride**: Cr (LiCl-KCl, Iwamatsu 2026), Zn (LiCl-KCl Iwamatsu 2022, plus NaCl/KCl/KBr Pikaev 1982), Nd (LiCl-KCl Castro Baldivieso 2026), Cf (aqueous Cl⁻ Rotermund 2024), plus the U NULL benchmark (Phillips 2022).
- **5 additional metals at lower priority** from Pikaev: Cd, Tl, Ag, Ca, Ba, Sr — span the Group 1/2/12/13 chemistry needed for chemistry-feature regression.
- **4 host matrices**: LiCl-KCl, NaCl, KCl, KBr (plus NaBr/KI for X₂•⁻ rates).
- **Theoretical thermodynamics anchor** for redox feasibility: 24 ΔE values from Kristoffersen-Metiu AIMD in LiCl/NaCl/mixture.
- **F₂ production** for LiF, BeF₂, ThF₄, FLiBe-UF₄ from Davis 2022 + Ea from Toth-Felker 1990 (fluoride side, unchanged).

This is **sufficient for the predictive-for-any-element extension** with 5 metals (Cr, Zn, Nd, Cf, plus optional Cd/Tl/Ag/Ca from Pikaev) across 4 hosts. The remaining gap is the experimental ε(λ,T) for e_s⁻ in molten chloride — workaround: scale-free observable mode + a literature-prior on ε.

**Still nice-to-have** (Tier B):
- Akiyama 1994 JNST FLiNaK — would complete the fluoride kernel with transient-kinetics data complementing the Davis G + Toth/Felker Ea split. Not blocking, since we don't claim transient-resolved fluoride predictions in this article.
- Makarov, Zhukova, Pikaev, Spitzyn *Bull. Acad. Sci. USSR* 1982 — extra multi-metal/multi-host X₂/X₃⁻ yields; redundant with the now-available Pikaev 1982 RPC core dataset.
- Sawamura 1990 LiBr-KBr — bromide cousin; only useful if we extend to bromide hosts.
- The actual experimental ε(λ,T) paper — see Tier A1 above.
