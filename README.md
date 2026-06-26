# Paper 2: Human Fecal Fractionation Atlas

This folder is the working workspace for the fractionation paper.

The original Word documents remain at the project root:

- `Paper-2_Fecal-Fractionation_v1-4.docx`: current manuscript draft
- `Data set relevant to paper2 and example plots.docx`: analysis instructions

## Workspace layout

| Folder | Purpose |
|---|---|
| `config/` | External paths and analysis settings |
| `data/raw/` | Read-only source data, or links/manifests pointing to source data |
| `data/metadata/` | Paper-specific dataset and sample metadata |
| `data/intermediate/` | Rebuildable intermediate datasets |
| `data/processed/` | Analysis-ready datasets, including `paper2_total_list.csv` |
| `scripts/` | Numbered analysis scripts in execution order |
| `scripts/reference/` | Copies of useful existing scripts; do not run as the final pipeline |
| `outputs/tables/` | Final CSV and supplementary tables |
| `outputs/qc/` | Data checks, exclusions, and merge audit reports |
| `figures/` | Final and draft figures grouped by figure number |
| `docs/` | Analysis plan, decisions, and data requirements |

## Planned execution order

1. Build the Paper 2 dataset manifest from the HGM workbook.
2. Build `HGMA -> parent -> extract/fraction -> HGMH/pool` metadata.
3. Merge selected annotation datasets into `data/processed/paper2_total_list.csv` using Annotation Merger logic for Figure 1.
4. Build the four-pool Figure 2 dataset and Figure 3 individual-crude/pool-fraction dataset using combining-datasets logic and `area > 0` presence.
5. Append NPClassifier superclass and subclass.
6. Generate figure-ready tables before generating plots.

Figure 4 library-generation inputs and Figure 5 cross-collection benchmarking are implemented.

See `docs/data-requirements.md` and `docs/analysis-plan.md` before starting analysis.

## Current pipeline status

Updated on 2026-06-10:

- Figure 1 candidate manifest: 140 datasets after explicitly excluding small-scale and lipidomic datasets, including small-scale bile-acid dataset `HGMD_0294`
- `data/processed/paper2_total_list.csv`: 16,926 merged level 1-3 platform-resolved annotation rows after MSNovelist exclusion
- `data/processed/paper2_total_list_classified.csv`: current NPClassifier-appended total list
- Figure 1 review tables and first-draft plots
- Four-pool sample lineage for all 24 Figure 2/3 datasets
- `data/intermediate/figure2-four-pool-presence.csv`: 16,355 distinct pool/processing-level annotation records using `area > 0`
- Figure 2 annotation-gain, processing-level unique annotation, and class-enrichment tables
- Demo-guided complete Figure 1 and Figure 2 assemblies, with individual subfigures and PDF exports
- Nature-style Figure 3 analysis and assembly, including crude-profile clustering, crude/fractionated PCA, matched differential analysis, prevalence and chemical-class summaries
- Figure 5 identity-resolved comparison of the non-lipid human faecal atlas, healthy-faecal HMDB, and MiMeDB collections; FASST public-repository evidence is reported separately because an identity-level FASST result list is unavailable

See `docs/figure1-2-specifications.md` for the implemented logic of every subfigure.
Draft publication captions are in `docs/figure-captions.md` and copied alongside each assembled figure as `CAPTION.md`.

Run scripts with:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\01_metadata\run_metadata.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\02_merge_annotations\run_merge.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\03_classification\run_classification.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\04_figure1\run_figure1.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\05_figure2\run_figure2.R'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' 'scripts\08_figure5\run_figure5.R'
```
