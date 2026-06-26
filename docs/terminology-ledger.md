# Paper 2 terminology ledger

This ledger records the terminology used consistently across the Figure 1-3 Results sections and should guide subsequent manuscript writing.

| Canonical term | Use | Avoid or reserve |
|---|---|---|
| human faecal metabolome atlas | Complete integrated Paper 2 annotation resource | human gut metabolome atlas, unless referring more broadly to the gut |
| faecal | Preferred spelling throughout the manuscript | fecal |
| multi-dimensional fractionation | Overall extraction and fractionation strategy | multidimensional fractionation |
| primary fractionation | Initial polar or non-polar flash-chromatography fractionation | fractionation collection when the processing stage must be explicit |
| subfractionation | Further separation of selected primary fractions | secondary fractionation |
| biological pool | One of the four pooled sample groups | donor pool when discussing pooled analytical material |
| crude profile | Integrated profile from crude processing levels | crude dataset when referring to the combined profile |
| fraction-resolved profile | Profile defined by biological pool, fraction type and fraction sequence | fractionated donor profile |
| processing level | Methanol crude, DCM phase, aqueous phase, non-polar fractions or polar fractions | extraction level |
| non-repeating annotation | Distinct `annotation_id`, defined by compound name where available and otherwise SMILES | unique compound, unless chemical identity is independently established |
| processing-level-specific annotation | Annotation detected at one processing level and absent from all others | unique annotation without specifying the comparison scope |
| platform-unique annotation | Annotation detected by one chromatographic and ionization platform only | unique annotation without specifying the platform scope |
| pool-unique annotation | Annotation detected in one biological pool only | donor-unique annotation |
| annotation coverage | Number or breadth of retained annotations | metabolite coverage when referring specifically to annotation counts |
| platform-resolved record | Retained annotation record within a chromatographic and ionization platform | unique annotation |
| analysis feature | Quantitative feature retained separately by platform and processing condition | annotation when discussing PCA or differential-feature counts |
| NPClassifier superclass/subclass | Chemical classification used in current Figures 1-3 | ClassyFire class or superclass |
| Phe-Hex | Phenyl-hexyl chromatography platform | PheHex |
| HILIC | Hydrophilic interaction liquid chromatography platform | Hilic |
| SAX | Strong anion-exchange chromatography platform | strong anion exchange when used as the platform label |
| positive-ion mode / negative-ion mode | Ionization polarity in prose | positive mode / negative mode |
| confidence level 1-3 | Annotation confidence categories | level one / level two / level three |
| MSNovelist | Excluded annotation source | MSNOVELIST, MS Novelist |
| Metabolite Annotation Propagation and Synthesis (MAPS) | Full name at first use, then MAPS | Metabolite Annotation and Propagation |
| FASST public repositories | Spectrum-level public-repository comparison | public FASST repo |
| atlas-derived in-house spectral library | Platform- and ion-mode-specific MSP library generated from the human faecal metabolome atlas | in-house lib, atlas library when the derivation or analytical scope is unclear |
| MZmine reprocessing | Reanalysis of an LC-MS/MS dataset in MZmine using its corresponding atlas-derived in-house spectral library | rerun, reprocess at mzmine |
| resolved comparison identity | Figure 5 identity after standardised-name matching and SMILES fallback | annotation ID when referring to cross-database comparison |
| atlas-only identity | Non-lipid resolved identity present in the atlas and absent from healthy HMDB and MiMeDB | novel metabolite, unless novelty is independently established |

## Statistical terminology

| Canonical term | Use |
|---|---|
| principal component analysis (PCA) | Unsupervised ordination used for crude and fraction-resolved quantitative profiles |
| PERMANOVA | Test of association between pool membership and multivariate composition |
| differential feature | Analysis feature meeting both the effect-size and BH-adjusted FDR thresholds |
| BH-adjusted FDR | Multiple-testing-adjusted significance measure used for differential analysis |

## Interpretation boundaries

- Crude-profile analyses retain individual-donor information.
- Fraction-resolved analyses describe pooled material and do not provide independent donor-level replication.
- Explained-variance percentages from the crude and fraction-resolved PCAs should not be compared directly because the analyses use different samples and feature matrices.
- Use `annotation` rather than `metabolite` when the statement depends specifically on computational annotation identity.
