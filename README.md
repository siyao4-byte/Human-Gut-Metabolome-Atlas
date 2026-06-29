# Human Faecal Metabolome Atlas

This repository contains the public code, figures, captions and atlas-derived spectral libraries associated with a deep fractionation study of the human faecal metabolome. The resource was designed to expose faecal chemical space that is compressed in conventional crude LC-MS/MS profiles by linking retained metabolite annotations to their extraction, fractionation, chromatographic platform and ion-mode provenance.

The repository is intended as a reproducible companion to the manuscript and as a practical entry point for reuse of the atlas-derived in-house MSP libraries. It does not include raw LC-MS/MS data, private sample-management exports, donor workflow notes, local QC files or manuscript Results/Discussion drafts.

## Repository Contents

The public release includes:

- `README.md`: overview of the atlas resource and repository layout.
- `docs/terminology-ledger.md`: controlled terminology used across the manuscript, figures and scripts.
- `config/paths.example.yml`: template for local paths required to rerun analyses.
- `scripts/`: analysis and figure-generation scripts.
- `figures/`: final assembled figures, panel-level figure exports and concise captions.
- `libraries/`: six atlas-derived MSP libraries for platform- and ion-mode-specific spectral matching.
- `supplementary_material.xlsx`: curated supplementary tables prepared for manuscript release.

The release intentionally excludes private/local files through `.gitignore`, including raw outputs, local metadata exports, working Word documents, manuscript prose drafts, donor workflow summaries and local path configuration.

## Atlas Scope

The atlas integrates annotations from crude extracts, biphasic phases, primary fractions, selected subfractions, complementary chromatographic platforms and positive- and negative-ion modes. Retained annotations are reported using confidence levels 1-3 and counted as non-repeating annotation identities after exclusion of MSNovelist-only annotations.

The resource should be interpreted as a lineage-resolved faecal metabolomics atlas rather than a population prevalence database. Some analyses use pooled biological material, and confidence categories should remain visible when interpreting or reusing the annotations.

## Data Products

### Figures

Final assembled figures are provided in `figures/`:

- `figures/figure-1/figure1-assembled.png`
- `figures/figure-2/figure2-assembled.png`
- `figures/figure-3/figure3-assembled.png`
- `figures/figure-4/figure4-5-atlas-utility-assembled.png`

Each figure folder also contains its corresponding `CAPTION.md` and selected panel-level exports used during assembly.

### MSP Libraries

Atlas-derived in-house spectral libraries are provided in `libraries/`:

- `PheHex-pos-paper2-rt075.msp`
- `PheHex-neg-paper2-rt075.msp`
- `HILIC-pos-paper2-rt075.msp`
- `HILIC-neg-paper2-rt075.msp`
- `SAX-pos-paper2-rt075.msp`
- `SAX-neg-paper2-rt075.msp`

These libraries were generated from atlas annotations and spectra after platform-specific library input preparation. They are intended for faecal LC-MS/MS reannotation and public dataset remining. Retention-time use depends on the compatibility of the target chromatographic method; cross-study matching should be interpreted primarily as spectral-library transfer unless retention time is experimentally validated.

## Repository Structure

```text
.
|-- config/
|   `-- paths.example.yml
|-- docs/
|   `-- terminology-ledger.md
|-- figures/
|   |-- figure-1/
|   |-- figure-2/
|   |-- figure-3/
|   `-- figure-4/
|-- libraries/
|-- scripts/
|   |-- 00_setup/
|   |-- 01_metadata/
|   |-- 02_merge_annotations/
|   |-- 03_classification/
|   |-- 04_figure1/
|   |-- 05_figure2/
|   |-- 06_figure3/
|   |-- 07_figure4/
|   `-- 08_figure5/
`-- supplementary_material.xlsx
```

## Running the Scripts

The scripts are written primarily in R, with Python used for MSP library generation. To rerun analyses locally, copy the example path file and edit it for your own system:

```powershell
Copy-Item config\paths.example.yml config\paths.yml
```

Then update `config/paths.yml` with local paths to annotation folders, metadata files and external comparison resources. Private institutional paths are not included in this public repository.

Example R execution:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\01_metadata\run_metadata.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\02_merge_annotations\run_merge.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\03_classification\run_classification.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\04_figure1\run_figure1.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\05_figure2\run_figure2.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\06_figure3\run_figure3.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\08_figure5\run_figure5.R'
```

Figure 4 library generation and atlas-utility analyses are organized under `scripts/07_figure4/`. Some steps require external MZmine, GNPS/MASSIVE or local annotation exports that are not redistributed here.

## Software Requirements

The analysis scripts use R and Python packages for data handling, statistics and plotting. Core R packages include:

- `tidyverse`
- `readxl`
- `yaml`
- `ggplot2`
- `patchwork`
- `ComplexHeatmap`
- `circlize`
- `igraph`
- `ggraph`
- `ggrepel`

Python is used for MSP-library construction from prepared annotation inputs and feature USIs. Install package versions compatible with your local LC-MS/MS processing workflow.

## Interpretation Notes

- `HGMH_xxxx`, `HGME_xxxx`, `HGMF_xxxx`, `HGMA_xxxx` and `HGMD_xxxx` identifiers are deidentified internal sample, fraction, annotation and dataset identifiers.
- Diet, gender and ethnicity metadata are retained where required for figure interpretation.
- Real donor names, private book codes, raw sample-management exports and local file paths are excluded from the public release.
- Atlas-only identities should not be interpreted as confirmed novel metabolites unless independently validated.
- Network-neighbouring features are prioritization candidates, not structural identifications.
- Public-dataset remining with the MSP libraries should be interpreted according to spectral match quality, confidence level and chromatographic compatibility.

## Citation

Please cite the associated manuscript when using this repository or the MSP libraries. Citation details will be added after publication.

## Contact

For questions about the atlas, scripts or MSP libraries, please contact the corresponding manuscript authors.

Run scripts with:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\01_metadata\run_metadata.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\02_merge_annotations\run_merge.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\03_classification\run_classification.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\04_figure1\run_figure1.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\05_figure2\run_figure2.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\08_figure5\run_figure5.R'
```
