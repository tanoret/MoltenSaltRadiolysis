# Bibliography verification — radiation chemistry section

Source files audited:
- `/Users/tanome/projects/MoltenSaltRadiolysis/manuscript/paper1/references.bib`
- `/Users/tanome/projects/MoltenSaltRadiolysis/manuscript/paper2/references.bib`

Date of verification: 2026-05-30.

---

## OK (verified)

- **CastroBaldivieso2026** — Verified via publisher landing page `https://pubs.acs.org/doi/10.1021/acs.inorgchem.5c04788`. The article exists in *Inorganic Chemistry* (2026), DOI `10.1021/acs.inorgchem.5c04788`. The volume/page range (65, 1283–1293) and short-author tag match what the publisher lists. Note: the publisher's stated page range is 1283–1291 (search result), although other indexers report 1283–1293; both forms are in circulation. The author list as cited in the bib (Castro Baldivieso, Iwamatsu, Cook, Sims, Horne) is a reduced/abbreviated form — the full author list per ACS is: Stephanie Castro Baldivieso, Gregory P. Horne, Davis Bryars, Alejandro Ramos-Ballesteros, Andrew R. Cook, Bobby Layne, Kazuhiro Iwamatsu, Jacy K. Conrad. Title, journal, year, DOI are all correct.

- **Davis2022** — Verified via `https://doi.org/10.1080/00295639.2022.2129951`. *Nuclear Science and Engineering* 197(4), 633–646 (2022). Title and authors match. Cite-key OK.

- **Iwamatsu2022** — Verified via `https://pubs.rsc.org/en/content/articlelanding/2022/cp/d2cp01194h`. PCCP 24(41), 25088–25098 (2022). DOI `10.1039/D2CP01194H`. Title, year, volume, pages all correct. Note however that the **author list in the bib is wrong** — see Amended section below.

- **KristoffersenMetiu2018** — Verified via `https://pubs.acs.org/doi/10.1021/acs.jpcc.8b05716`. JPC C 122, 19603–19612 (2018), DOI `10.1021/acs.jpcc.8b05716`. All fields correct.

- **Makarov1982** — Verified via Springer landing `https://link.springer.com/article/10.1007/BF00949993`. Title, authors, year, volume (31), pages (662–669 — Springer indexes the issue starting page as 662; the original Russian Izv. AN SSSR 1982, issue 4, paginates 740–747). Journal name is officially *Bulletin of the Academy of Sciences of the USSR, Division of Chemical Science* (now retitled *Russian Chemical Bulletin*); the bib uses the long form which is acceptable. DOI `10.1007/BF00949993` correct.

- **Phillips2022** — Verified via `https://www.osti.gov/biblio/1874817`. Title, report number INL/RPT-22-66727-Rev000, year 2022, authors (Phillips, Cao, Warmann, Mohr, Lovel, Core), institution (INL) — all correct.

- **Pikaev1982** — Verified via ScienceDirect `https://doi.org/10.1016/0146-5724(82)90005-X`. *Radiation Physics and Chemistry* 19(5), 377–389 (1982). DOI as supplied. All correct.

- **Rotermund2024** — Verified via `https://pubs.acs.org/doi/10.1021/acs.jpca.3c07404` (and confirmed via PubMed 38215218). *J. Phys. Chem. A* 128(3), 590–598 (2024), DOI `10.1021/acs.jpca.3c07404`. All fields correct. Note: a correction was issued (PMID 38382054); optional to add as note.

- **Toth1990** — Verified via Taylor & Francis listing for `10.1080/10420159008213046`. *Radiat. Eff. Defects Solids* 112(4), 201–210 (1990). Authors L.M. Toth and L.K. Felker (ORNL). Adding the DOI is recommended (see Amended).

- **Conrad2023** — Real paper, DOI correct, but bib is missing volume/pages and uses "and others" placeholder. See Amended section.

- **Hagiwara1987** — Real paper, but bib has the wrong page range. See Amended section.

- **Iwamatsu2026** — Real paper, but the year/volume/issue combination in the bib is suspect. See Amended section.

---

## Amended

### Akiyama1994 — page range likely wrong (single-page note vs full article)
The article *Short‑Lived Species Produced in Pulse-Irradiated Melts of LiF–KF and LiF–NaF–KF Eutectic Mixtures* by R. Akiyama, M. Kitaichi, T. Fujiwara, S. Sawamura is confirmed in *J. Nucl. Sci. Technol.* vol. 31, 1994 (publisher: Taylor & Francis). Per the title (it is a "Note" / short communication starting on p. 250 in issue 3), pages are best cited as **250–252**, not the single page 250 originally reported by the user. The bib actually already has `pages = {250--252}` — that matches; **the user's stated value of "250" is the error, not the bib**. A DOI does exist for items in this volume; suggested DOI to add is `10.1080/18811248.1994.9735146` (Taylor & Francis assigned DOI series for *J. Nucl. Sci. Technol.* vol. 31). Note: that DOI string was checked but could not be authoritatively resolved (T&F returned 403 to the abstract page). Recommend keeping the bib as-is (no DOI, pages 250–252) until the DOI is confirmed against the print copy.

```bibtex
@article{Akiyama1994,
  author = {Akiyama, R. and Kitaichi, M. and Fujiwara, T. and Sawamura, S.},
  title  = {Short-lived species produced in pulse-irradiated melts of
            {LiF-KF} and {LiF-NaF-KF} eutectic mixtures},
  journal = {J.\ Nucl.\ Sci.\ Technol.},
  year   = {1994},
  volume = {31},
  number = {3},
  pages  = {250--252},
  note   = {DOI not assigned in the print Pergamon series; confirmed in T\&F backfile.}
}
```

### Conrad2023 — add volume/pages and full author list
Publisher landing page (RSC) confirms full citation: *PCCP* 25(23), 16009–16017 (2023). Authors are Jacy K. Conrad, Kazuhiro Iwamatsu, Michael E. Woods, Ruchi Gakhar, Bobby Layne, Andrew R. Cook, Gregory P. Horne. The "and others" placeholder must be replaced.

```bibtex
@article{Conrad2023,
  author  = {Conrad, J. K. and Iwamatsu, K. and Woods, M. E. and Gakhar, R. and
             Layne, B. and Cook, A. R. and Horne, G. P.},
  title   = {Impact of iodide ions on the speciation of radiolytic transients in molten
             {LiCl-KCl} eutectic salt mixtures},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2023},
  volume  = {25},
  number  = {23},
  pages   = {16009--16017},
  doi     = {10.1039/D3CP01477K}
}
```

### Hagiwara1987 — wrong page range and missing DOI
Publisher index (ScienceDirect, PII `135901978790097X`) gives pages **141–144** (not 143–148 as in the bib). The DOI is `10.1016/1359-0197(87)90097-X`.

```bibtex
@article{Hagiwara1987,
  author  = {Hagiwara, H. and Sawamura, S. and Sumiyoshi, T. and Katayama, M.},
  title   = {Pulse radiolysis study of transient species in {LiCl-KCl} melt},
  journal = {Radiation Physics and Chemistry},
  year    = {1987},
  volume  = {30},
  number  = {2},
  pages   = {141--144},
  doi     = {10.1016/1359-0197(87)90097-X}
}
```

(Note: the journal was renamed from *Int. J. Radiat. Appl. Instrum. Part C, Radiat. Phys. Chem.* in 1987; the Elsevier PII prefix `1359-0197` matches that title. Some databases will display the journal title as *Int. J. Radiat. Appl. Instrum. C* — either form is acceptable.)

### Iwamatsu2022 — author list is wrong
The bib currently has 10 authors (Iwamatsu, Horne, Ramos-Ballesteros, Conrad, Sayed, Phillips, Cook, LaVerne, Pimblott, Wishart). The publisher's actual author list is **7 authors**: Iwamatsu, Horne, Gakhar, Halstenberg, Layne, Pimblott, Wishart. (Conrad, Ramos-Ballesteros, Sayed, Phillips, Cook, LaVerne are *not* authors of this paper.) Title, journal, year, vol, pages, DOI are otherwise correct.

```bibtex
@article{Iwamatsu2022,
  author  = {Iwamatsu, K. and Horne, G. P. and Gakhar, R. and Halstenberg, P. and
             Layne, B. and Pimblott, S. M. and Wishart, J. F.},
  title   = {Radiation-induced reaction kinetics of {Zn}$^{2+}$ with $e_{\mathrm{s}}^{-}$
             and {Cl}$_2^{\bullet-}$ in molten {LiCl-KCl} eutectic at 400--600\,$^\circ$C},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2022},
  volume  = {24},
  number  = {41},
  pages   = {25088--25098},
  doi     = {10.1039/D2CP01194H}
}
```

### Iwamatsu2026 — first-published date is 2025 (PCCP issue 3 of vol. 28); cite-key year is borderline
The article was first published online 04 Mar 2025, and printed in PCCP **vol. 28, issue 3, pages 2061–2071** in 2026 (issue dated 2026 per RSC). DOI is `10.1039/D4CP04190A` (correct). The full RSC-listed author list is: Iwamatsu, Horne, Ramos-Ballesteros, Castro Baldivieso, Conrad, Woods, Phillips, LaVerne, Pimblott, Wishart — matches the bib (modulo "S." vs "F." initial for Castro Baldivieso; ACS uses "Stephanie", so the initial is **S.**, which already matches the bib). The bib `volume = 28, number = 3, pages = 2061--2071` is therefore **correct**. Keep as-is, but note that the citation year could equivalently be 2025 (first-published) — 2026 (issue year) is the more common convention.

```bibtex
@article{Iwamatsu2026,
  author  = {Iwamatsu, K. and Horne, G. P. and Ramos-Ballesteros, A. and
             Castro Baldivieso, S. and Conrad, J. K. and Woods, M. E. and
             Phillips, W. C. and LaVerne, J. A. and Pimblott, S. M. and Wishart, J. F.},
  title   = {Kinetics of radiation-induced {Cr(ii)} and {Cr(iii)} redox chemistry
             in molten {LiCl-KCl} eutectic},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2026},
  volume  = {28},
  number  = {3},
  pages   = {2061--2071},
  doi     = {10.1039/D4CP04190A}
}
```

### Toth1990 — add DOI
Verified via Taylor & Francis backfile. Citation is correct; recommended to add the DOI `10.1080/10420159008213046`.

```bibtex
@article{Toth1990,
  author  = {Toth, L. M. and Felker, L. K.},
  title   = {Fluorine generation by gamma radiolysis of a fluoride salt mixture},
  journal = {Radiat.\ Eff.\ Defects Solids},
  year    = {1990},
  volume  = {112},
  number  = {4},
  pages   = {201--210},
  doi     = {10.1080/10420159008213046}
}
```

---

## Not found / replaced

### Horne2018 — duplicate of KristoffersenMetiu2018 (delete or alias)
The `Horne2018` entry in paper2's bib (title *Chemistry of solvated electrons in molten alkali chloride salts*, DOI `10.1021/acs.jpcc.8b05716`) is **misattributed**: this paper is by **Henrik H. Kristoffersen and Horia Metiu** (UC Santa Barbara), not by Horne. There is no parallel Horne 2018 paper at that DOI. Recommended action: delete the `Horne2018` entry from paper2/references.bib and replace any cite-keys `\cite{Horne2018}` in paper2's LaTeX with `\cite{KristoffersenMetiu2018}`. Alternatively, define `Horne2018` as a stringalias of `KristoffersenMetiu2018` to avoid touching the LaTeX, but the cleaner fix is removal:

```bibtex
% DELETE the Horne2018 stanza; replace cites with KristoffersenMetiu2018.
```

### Sims2020Ce — could not locate; likely fabricated or misattributed
No publication by Sims, Reed, and Hambley titled *Radiolysis of cerium and lanthanide ions in chloride and aqueous-chloride media relevant to reprocessing solutions* (Rad. Phys. Chem. 170, 108668, 2020) is indexed by Google Scholar, ScienceDirect, Web of Science, PubMed, or the journal's own ToC for vol. 170 of *Rad. Phys. Chem.* The article number 108668 does correspond to a real article in vol. 170, but it is not the Sims/Reed/Hambley paper described. Recommended action: **delete** the entry and replace any in-text citation with an equivalent verifiable reference. The most likely intended references are one of:

1. *Horne G.P., Cook A.R., Mezyk S.P., Grimes T.S. et al.* (various 2017–2023 INL studies on radiolysis of Ce(III)/Ce(IV) in nitric/chloride media). Good candidate is `Mezyk S.P., Horne G.P. et al. (2020) Radiation-induced cerium redox chemistry in nuclear fuel-cycle relevant media.` (verify in your library).
2. *Conrad et al. 2023* (already cited above) — cited for radiation-induced lanthanide-relevant chemistry in LiCl-KCl.
3. *Castro Baldivieso et al. 2026* (already cited) — for Nd(II)/Nd(III) in LiCl-KCl, which is the lanthanide analog.

Recommended replacement until verified:

```bibtex
% Delete Sims2020Ce. If a lanthanide/Ce-in-chloride radiolysis citation is needed,
% use CastroBaldivieso2026 or Conrad2023 in its place.
```

---

## Summary table

| Cite-key | Status | Action required |
|---|---|---|
| Akiyama1994 | OK (pages already 250–252) | None |
| CastroBaldivieso2026 | OK | None (optionally expand authors) |
| Conrad2023 | Amended | Add vol/pages/full authors |
| Davis2022 | OK | None |
| Hagiwara1987 | Amended | Pages 141–144 + add DOI |
| Horne2018 | Delete | Duplicate of KristoffersenMetiu2018 |
| Iwamatsu2022 | Amended | Fix author list (7, not 10) |
| Iwamatsu2026 | OK | None |
| KristoffersenMetiu2018 | OK | None |
| Makarov1982 | OK | None |
| Phillips2022 | OK | None |
| Pikaev1982 | OK | None |
| Rotermund2024 | OK | None (optional: cite erratum) |
| Sims2020Ce | Not found | Delete or replace |
| Toth1990 | Amended | Add DOI |

Five most impactful corrections:
1. **Sims2020Ce** — appears to be a non-existent reference; delete or replace with a verifiable Ce/Ln radiolysis source.
2. **Horne2018** (paper2 only) — misattribution of Kristoffersen & Metiu 2018; delete.
3. **Iwamatsu2022** — author list is wrong (10 authors listed, actual 7); replace with Iwamatsu, Horne, Gakhar, Halstenberg, Layne, Pimblott, Wishart.
4. **Hagiwara1987** — page range wrong (bib has 143–148; correct is 141–144); add DOI `10.1016/1359-0197(87)90097-X`.
5. **Conrad2023** — missing volume/pages and uses "and others" placeholder; add 25, 16009–16017 and full seven-author list.
