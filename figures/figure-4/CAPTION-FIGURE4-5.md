# Merged Figure 4/5. Atlas-derived libraries enable public faecal dataset remining and atlas-specific network interpretation

**a,** Platform-matched MZmine reprocessing of crude atlas datasets using atlas-derived in-house MSP libraries, compared with the corresponding MAPS-processed outputs. MAPS annotations were restricted to confidence levels 1-3 and MSNovelist annotations were excluded. Stacked bars count distinct identities within each annotation source.

**b,** Public faecal dataset remining with atlas-derived in-house MSP libraries. Bars show non-repeating atlas-library annotations recovered from the public BeefDiet positive-ion dataset (`HGMD_0357`), public BeefDiet negative-ion dataset (`HGMD_0358`) and NIST human faecal dataset (`HGMD_0359`). Labels show non-repeating annotations and, in parentheses, the annotation rate calculated as non-repeating atlas-library annotations divided by the total number of resolved MZmine features in the corresponding `ms2 or ion identity` feature list. These public searches were performed without retention-time matching and should therefore be interpreted as spectral-library remining rather than retention-time-validated identification.

**c,** Venn diagram showing identity overlap among the human faecal atlas, healthy-faecal HMDB records and MiMeDB when lipid and steroid-associated identities are retained. Standardised `compound.name` was used as the primary identity key, with SMILES used as fallback when standardised names were unavailable.

**d,** Chemical-class breakdown of atlas-only identities from the lipid-included comparison in **c**. Bars show the most abundant NPClassifier superclasses among identities found in the atlas but absent from both healthy-faecal HMDB and MiMeDB. Unmatched structures are shown as `Unclassified`.

**e,** Full NIST faecal molecular network annotated with the atlas-derived MSP library. Node colour denotes the best atlas MSP match confidence level; node size denotes the number of NIST subjects in which the feature was detected after collapsing technical replicates; gold outer marks indicate atlas-matched identities absent from both HMDB and MiMeDB; and edge width encodes MS/MS cosine similarity. Text labels are shown for the 30 level 1 or level 2 atlas MSP anchors with the highest network degree.

For all panels, non-repeating annotations were counted after exclusion of MSNovelist annotations where relevant. Identity-level overlap analyses use resolved comparison identities, not the full atlas-wide `annotation_id` count.
