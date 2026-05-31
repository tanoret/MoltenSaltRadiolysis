# Bibliography verification — MSR / nuclear engineering section

Verification of MSR design, nuclear engineering, molten-salt properties, and safety review entries from
`paper1/references.bib` and `paper2/references.bib`.

---

## OK (verified)

### Briggs1964
- ORNL-3708, R. B. Briggs, "Molten-Salt Reactor Program Semiannual Progress Report for Period
  Ending July 31, 1964", November 1964. Confirmed via OSTI.GOV biblio 4676587 and UNT Digital Library.
  Source: https://www.osti.gov/biblio/4676587/  — entry is correct as cited.

### Forsberg2020
- C. W. Forsberg, "Market Basis for Salt-Cooled Reactors: Dispatchable Heat, Hydrogen, and Electricity
  with Assured Peak Power Capacity", Nuclear Technology 206 (11), pp. 1659–1685, 2020.
  DOI: 10.1080/00295450.2020.1743628 (Taylor & Francis). Entry correct; DOI could be added.

### Olander2002
- D. Olander, "Redox condition in molten fluoride salts: Definition and control", J. Nucl. Mater.
  300 (2002) 270–272. DOI: 10.1016/S0022-3115(01)00742-5. Entry correct.

### Romatoski2017
- R. R. Romatoski and L.-W. Hu, "Fluoride salt coolant properties for nuclear reactor applications:
  A review", Annals of Nuclear Energy 109 (2017) 635–647. DOI: 10.1016/j.anucene.2017.05.036.
  Entry correct (DOI already present).

### Williams2006
- D. F. Williams, L. M. Toth, K. T. Clarno, "Assessment of Candidate Molten Salt Coolants for
  the Advanced High-Temperature Reactor (AHTR)", ORNL/TM-2006/12, March 2006.
  Confirmed via ORNL publications portal and OSTI. Entry correct.

### Sohal2010
- M. S. Sohal, M. A. Ebner, P. Sabharwall, P. Sharpe, "Engineering Database of Liquid Salt
  Thermophysical and Thermochemical Properties", INL/EXT-10-18297, March 2010 (Rev. 1 issued
  June 2013). Confirmed via INL Digital Library and OSTI biblio 980801. Entry correct.

### Holcomb2012 (note: misnamed key — it is a 2010 report)
- D. E. Holcomb and S. M. Cetiner, "An Overview of Liquid-Fluoride-Salt Heat Transport Systems",
  ORNL/TM-2010/156, September 2010. Confirmed via ORNL Pub25407 and OSTI biblio 990239.
  The bib entry has year=2010, which is correct; only the BibTeX key (`Holcomb2012`) is
  cosmetically wrong. Either rename the key to `Holcomb2010` or leave it (functionally fine).

### IAEA2013MSR
- IAEA, "Challenges Related to the Use of Liquid Metal and Molten Salt Coolants in Advanced
  Reactors", IAEA-TECDOC-1696, IAEA, Vienna, 2013. Confirmed via iaea.org/publications/8942.
  Entry correct.

### TerraPower2021 (gray literature — caveat)
- TerraPower's "Molten Chloride Fast Reactor (MCFR)" materials exist as company tech briefs
  (terrapower.com/our-work/molten-chloride-fast-reactor-technology/) and an NRC pre-application
  presentation (ML21228A222, 2021). Citation as "Technical brief, TerraPower LLC, 2021" is
  acceptable but is gray literature. **Recommendation:** if the citation supports a technical
  claim, replace with peer-reviewed substitute such as Andreades2014 (PB-FHR design) or
  Latkowski/Cisneros TerraPower MCFR ICAPP/PHYSOR proceedings. Otherwise retain as misc with
  an accessed URL.

---

## Amended

### Andreades2014 — year and pagination need fixing
**Current:** year = 2014, volume 195, pages 223–238.
**Correct:** *Nuclear Technology* 195 (3), 223–238, **September 2016** (received 2014, published
2016). DOI: 10.13182/NT16-2.

Corrected BibTeX:
```bibtex
@article{Andreades2016,
  author  = {Andreades, C. and Cisneros, A. T. and Choi, J. K. and Chong, A. Y. K. and
             Fratoni, M. and Hong, S. and Huddar, L. R. and Huff, K. D. and Kendrick, J. and
             Krumwiede, D. L. and Laufer, M. R. and Munk, M. and Scarlat, R. O. and
             Zweibaum, N. and Greenspan, E. and Peterson, P. F.},
  title   = {Design Summary of the {Mark-I} Pebble-Bed, Fluoride Salt--Cooled,
             High-Temperature Reactor Commercial Power Plant},
  journal = {Nuclear Technology},
  year    = {2016},
  volume  = {195},
  number  = {3},
  pages   = {223--238},
  doi     = {10.13182/NT16-2}
}
```
Update citation key `Andreades2014` → `Andreades2016` (or keep the key but fix `year=2016`).

### Carotti2017 — title and venue look fabricated
**Current entry:** Carotti, Liu, Scarlat, "A review of corrosion mechanisms and chemistry control
strategies for fluoride salt-cooled high-temperature reactors", *Annals of Nuclear Energy* 110
(2017) 1051–1057.
**Findings:** No paper with this exact title/venue exists. Carotti+Scarlat 2017 publications are
electrochemistry in JES (Carotti, Wu, Scarlat, "Characterization of a Thermodynamic Reference
Electrode for Molten LiF-BeF2 (FLiBe)", *J. Electrochem. Soc.* 164(12), H854–H861, 2017).
The likely intended citation is the closely related multi-author review:
Zhang J., Forsberg C. W., Simpson M. F., Guo S., Lam S. T., Scarlat R. O., Carotti F.,
Chan K. J., Singh P. M., Doniger W., Sridharan K., Keiser J. R., "Redox potential control in
molten salt systems for corrosion mitigation", *Corrosion Science* 144 (2018) 44–53.
DOI: 10.1016/j.corsci.2018.08.035.

Corrected BibTeX (replace Carotti2017 with):
```bibtex
@article{Zhang2018Redox,
  author  = {Zhang, J. and Forsberg, C. W. and Simpson, M. F. and Guo, S. and Lam, S. T. and
             Scarlat, R. O. and Carotti, F. and Chan, K. J. and Singh, P. M. and Doniger, W. and
             Sridharan, K. and Keiser, J. R.},
  title   = {Redox potential control in molten salt systems for corrosion mitigation},
  journal = {Corrosion Science},
  year    = {2018},
  volume  = {144},
  pages   = {44--53},
  doi     = {10.1016/j.corsci.2018.08.035}
}
```

### Janz1988 — venue is JPCRD Supplement, not NSRDS-NBS 61
**Current entry:** lists "series = NSRDS-NBS, NSRDS-NBS 61, Part II".
**Correct identification:** G. J. Janz, "Thermodynamic and Transport Properties for Molten Salts:
Correlation Equations for Critically Evaluated Density, Surface Tension, Electrical Conductance,
and Viscosity Data", *Journal of Physical and Chemical Reference Data* **17**, Supplement No. 2,
1988 (ACS / AIP for NBS). NSRDS-NBS 61 (Parts I, II, IV) is a different series published 1978/1979/1981.

Corrected BibTeX:
```bibtex
@article{Janz1988,
  author  = {Janz, G. J.},
  title   = {Thermodynamic and Transport Properties for Molten Salts: Correlation Equations
             for Critically Evaluated Density, Surface Tension, Electrical Conductance, and
             Viscosity Data},
  journal = {Journal of Physical and Chemical Reference Data},
  year    = {1988},
  volume  = {17},
  number  = {Supplement 2},
  publisher = {American Chemical Society / American Institute of Physics (for the National
               Bureau of Standards)}
}
```

---

## Not found / replaced

### ForsbergSafety2019 — could not verify
**Current entry:** Forsberg, "Safety and licensing aspects of the molten salt reactor",
*Nuclear Technology* 206 (2019) 1685–1709.
**Findings:** No Forsberg paper with that exact title exists in *Nuclear Technology*. The page
range overlaps suspiciously with Forsberg2020 (1659–1685). Real adjacent Forsberg articles in
*NT* vol. 206 (2019/2020) include:
- "Heat-Pipe Heat Exchangers for Salt-Cooled Fission and Fusion Reactors..." NT 206 (2020),
  DOI 10.1080/00295450.2019.1681222.
- "Fusion Blankets and FHRs with Flibe Salt Coolant..." NT 206 (2020) 1778,
  DOI 10.1080/00295450.2019.1691400.

Recommendation: replace with one of the verified Forsberg NT 206 articles, or with the
peer-reviewed safety/licensing review by Andreades/Forsberg/Greenspan, **or** with the ORNL
report "Molten Salt Reactor Fundamental Safety Function PIRT" (Holcomb et al., ORNL Pub165504).
If the citation context is generic "MSR safety", an alternative is:
- Elsheikh, B. M., "Safety assessment of molten salt reactors in comparison with light water
  reactors", *J. Radiation Research and Applied Sciences* 6 (2013) 63–70.
If the citation is non-load-bearing, **delete it** rather than fabricate.

### Mausolf2024 — could not locate
**Current entry:** Mausolf, Cao, Phillips et al., "Radiolytic effects in molten chloride salts
for advanced reactor applications", *J. Nucl. Mater.* 593 (2024) 154974.
**Findings:** *J. Nucl. Mater.* 593 (May 2024) exists, but no Mausolf-authored article on
molten chloride radiolysis was returned by Scholar, OSTI, or ScienceDirect search.
Article number 154974 does not resolve to a Mausolf paper. The closest *extant* work by the
named co-authors is:
- Phillips, W. C., Cao, G., Warmann, S. A., Mohr, B. C., et al., **"Gamma Irradiation of NaCl-UCl3
  Salt for the Molten Chloride Fast Reactor"**, INL/RPT-22-66727, INL, 2022.
  OSTI biblio 1874817; https://www.osti.gov/servlets/purl/1874817.

Recommendation: **replace** Mausolf2024 with the Phillips2022 INL report (which is in fact
already cited elsewhere in the bib as `Phillips2022`). Likely the Mausolf2024 entry is a
hallucinated duplicate of Phillips2022. **Delete** Mausolf2024 and reuse Phillips2022.

Suggested replacement BibTeX (if not already present):
```bibtex
@techreport{Phillips2022,
  author      = {Phillips, W. C. and Cao, G. and Warmann, S. A. and Mohr, B. C. and others},
  title       = {Gamma Irradiation of {NaCl-UCl3} Salt for the Molten Chloride Fast Reactor},
  institution = {Idaho National Laboratory},
  number      = {INL/RPT-22-66727},
  year        = {2022},
  url         = {https://www.osti.gov/biblio/1874817}
}
```

### Tomeck2024MCFR — fabricated; recommend deletion
**Current entry:** Tomeck, J., et al., "Survey of molten chloride fast reactor designs and
NaCl-UCl3 composition windows", *Progress in Nuclear Energy* 171 (2024) 105180.
**Findings:** No author "Tomeck" appears in Scopus, Scholar, or ScienceDirect with any MCFR
publication. *Progress in Nuclear Energy* vol. 171 exists but does not contain this paper.
This is almost certainly a hallucinated reference produced by an earlier agent.

Closest extant works:
- Wang, J., Yan, R., Zuo, X., Zhou, B., Cai, X.-Z., "Design space analysis for a NaCl-UCl3 based
  breed and burn reactor system", *Annals of Nuclear Energy* (2025).
- Yu, C., et al., "Design and assessment of a molten chloride fast reactor",
  *Nuclear Engineering and Design* (2021), DOI 10.1016/j.nucengdes.2021.111180.
- Latkowski, J. F., et al. (TerraPower), "TerraPower's Molten Chloride Fast Reactor Technology"
  (industry talk / proceedings).

Recommendation: **delete** `Tomeck2024MCFR`. Replace with Wang2025 or Yu2021 if a peer-reviewed
MCFR design citation is needed.

### Yoshikawa2018 — could not verify; likely fabricated
**Current entry:** Yoshikawa, A. and Tropf, D. R. et al., "Equilibrium considerations of molten
fluoride salt redox chemistry relevant to molten salt reactor design", *Nuclear Engineering and
Design* 340 (2018) 284–292.
**Findings:** No Yoshikawa+Tropf paper on this topic appears in any database. *Nuclear
Engineering and Design* vol. 340 (2018) does not contain this article. The author "Tropf, D. R."
returns no results.

Closest extant works on fluoride salt redox 2017-2018:
- Guo, S., Zhang, J., Wu, W., Zhou, W., "Corrosion in the molten fluoride and chloride salts and
  materials development for nuclear applications", *Progress in Materials Science* 97 (2018) 448–487.
- Wang, Y.-L. et al., "Measurement of europium(III)/europium(II) couple in fluoride molten salt
  for redox control in a molten salt reactor concept", *J. Nucl. Mater.* 496 (2017) 197–206.
  DOI: 10.1016/j.jnucmat.2017.09.038.
- Olander2002 (already in bib).

Recommendation: **delete** Yoshikawa2018 and replace with Wang2017 (the Eu(III)/Eu(II) JNM
paper) or Olander2002 if the cited context is redox potential / equilibrium in fluoride salts.

---

## Summary

- **Verified as cited:** Briggs1964, Forsberg2020, Olander2002, Romatoski2017, Williams2006,
  Sohal2010, Holcomb2012/2010 (cosmetic key mismatch only), IAEA2013MSR, TerraPower2021
  (gray-lit caveat).
- **Need fixing:** Andreades2014 (year 2016, full author list, add DOI), Janz1988 (cite JPCRD
  Suppl. 17(2) not NSRDS-NBS 61), Carotti2017 (replace with Zhang et al. 2018 *Corros. Sci.* 144).
- **Likely fabricated / not found:** ForsbergSafety2019, Mausolf2024, Tomeck2024MCFR,
  Yoshikawa2018.

### Five highest-impact corrections
1. **Tomeck2024MCFR** — delete; this citation is fabricated.
2. **Mausolf2024** — delete; merge into already-existing Phillips2022 INL report.
3. **Yoshikawa2018** — delete; replace with Wang et al. 2017 *J. Nucl. Mater.* 496 if needed.
4. **Andreades2014** — change year to 2016, add full author list and DOI (10.13182/NT16-2).
5. **Carotti2017** — replace with Zhang et al., "Redox potential control in molten salt systems
   for corrosion mitigation", *Corros. Sci.* 144 (2018) 44–53,
   DOI 10.1016/j.corsci.2018.08.035.

Additional cosmetic fix: rename key `Holcomb2012` → `Holcomb2010` (the report is from 2010).
Additional cosmetic fix: Janz1988 should be a JPCRD article, not an NSRDS-NBS series book.
