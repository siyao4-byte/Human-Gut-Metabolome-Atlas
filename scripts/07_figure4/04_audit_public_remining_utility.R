source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(scales)
})

# Figure 4 public remining utility audit.
# Purpose: summarize atlas-library MZmine reprocessing outputs for internal
# crude datasets and public faecal datasets, recover sample identifiers from
# MZmine USI strings, and scan the HGM workbook for matching public dataset
# metadata. This script does not change assembled figures.

reprocessed_root <- file.path(project_root, "outputs", "mzmine reprocessing")
table_dir <- file.path(project_root, "outputs", "tables")
qc_dir <- file.path(project_root, "outputs", "qc")
figure_dir <- file.path(project_root, "figures", "figure-4")
ensure_dirs(table_dir, qc_dir, figure_dir)
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

internal_ids <- c("HGMD_0314", "HGMD_0315", "HGMD_0318", "HGMD_0319")
public_ids <- c("HGMD_0357", "HGMD_0358", "HGMD_0359")
dataset_ids <- c(internal_ids, public_ids)

manual_metadata <- tibble::tribble(
  ~dataset.ID, ~utility_role, ~collection_label, ~public_source_hint,
  "HGMD_0314", "Mine crude atlas datasets", "Phe-Hex positive crude", "",
  "HGMD_0315", "Mine crude atlas datasets", "Phe-Hex negative crude", "",
  "HGMD_0318", "Mine crude atlas datasets", "HILIC positive crude", "",
  "HGMD_0319", "Mine crude atlas datasets", "HILIC negative crude", "",
  "HGMD_0357", "Remine public datasets without RT matching", "Public BeefDiet positive-ion dataset", "MASSIVE/MSV public repository",
  "HGMD_0358", "Remine public datasets without RT matching", "Public BeefDiet negative-ion dataset", "MASSIVE/MSV public repository",
  "HGMD_0359", "Remine public datasets without RT matching", "Human faecal NIST dataset", "NIST public human faecal dataset"
)

find_reprocessed_csv <- function(dataset_id) {
  files <- list.files(
    file.path(reprocessed_root, dataset_id),
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(files)) return(NA_character_)
  files[[which.max(file.info(files)$size)]]
}

make_annotation_id <- function(compound_name, smiles) {
  compound_name <- na_if(clean_text(compound_name), "")
  smiles <- na_if(clean_text(smiles), "")
  coalesce(compound_name, smiles)
}

read_reprocessed <- function(dataset_id) {
  path <- find_reprocessed_csv(dataset_id)
  if (is.na(path)) return(tibble())
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    transmute(
      dataset.ID = dataset_id,
      source_file = basename(path),
      source_row = row_number(),
      annotation_id = make_annotation_id(compound_name, smiles),
      compound.name = na_if(clean_text(compound_name), ""),
      smiles = na_if(clean_text(smiles), ""),
      score = suppressWarnings(as.numeric(score)),
      precursor_mz = suppressWarnings(as.numeric(precursor_mz)),
      rt = suppressWarnings(as.numeric(rt)),
      query_spectrum_usi = clean_text(query_spectrum_usi)
    ) %>%
    filter(!is.na(annotation_id))
}

parse_usi_samples <- function(annotation_rows) {
  annotation_rows %>%
    filter(query_spectrum_usi != "") %>%
    tidyr::separate_rows(query_spectrum_usi, sep = ";") %>%
    mutate(
      query_spectrum_usi = clean_text(query_spectrum_usi),
      sample_id = str_match(query_spectrum_usi, "^mzspec:[^:]+:([^:]+):(.+)$")[, 2],
      scan_id = str_match(query_spectrum_usi, "^mzspec:[^:]+:([^:]+):(.+)$")[, 3],
      sample_id = na_if(clean_text(sample_id), ""),
      scan_id = na_if(clean_text(scan_id), "")
    ) %>%
    filter(!is.na(sample_id)) %>%
    mutate(
      public_collection = case_when(
        str_detect(sample_id, regex("BeefDiet", ignore_case = TRUE)) ~ "BeefDiet MASSIVE dataset",
        str_detect(sample_id, regex("^NIST", ignore_case = TRUE)) ~ "NIST human faecal dataset",
        str_detect(dataset.ID, "HGMD_03(14|15|18|19)") ~ "Internal crude atlas dataset",
        TRUE ~ "Unknown public dataset"
      ),
      inferred_polarity = case_when(
        str_detect(sample_id, regex("MS2pos|_POS_", ignore_case = TRUE)) ~ "positive",
        str_detect(sample_id, regex("MS2neg|_NEG_", ignore_case = TRUE)) ~ "negative",
        TRUE ~ NA_character_
      ),
      nist_subject = str_match(sample_id, "NIST_[A-Z]+_Samp_(\\d+)-\\d+")[, 2],
      nist_replicate = str_match(sample_id, "NIST_[A-Z]+_Samp_\\d+-(\\d+)")[, 2],
      beefdiet_sample_number = str_match(sample_id, "MS2(?:pos|neg)_(\\d+)_")[, 2],
      beefdiet_fraction = str_match(sample_id, "_(L\\d+-\\d+)$")[, 2]
    )
}

scan_hgm_workbook <- function(terms) {
  if (!file.exists(paths$hgm_workbook)) return(tibble())
  sheets <- readxl::excel_sheets(paths$hgm_workbook)
  purrr::map_dfr(sheets, function(sheet) {
    x <- suppressMessages(readxl::read_excel(paths$hgm_workbook, sheet = sheet, col_types = "text"))
    if (!nrow(x)) return(tibble())
    names(x) <- make.names(names(x), unique = TRUE)
    row_text <- apply(as.data.frame(x, stringsAsFactors = FALSE), 1, function(row) {
      paste(clean_text(row), collapse = " | ")
    })
    hit <- str_detect(row_text, regex(paste(terms, collapse = "|"), ignore_case = TRUE))
    if (!any(hit)) return(tibble())
    tibble(
      sheet = sheet,
      row_number = which(hit) + 1L,
      matched_text = row_text[hit]
    )
  })
}

annotation_rows <- bind_rows(lapply(dataset_ids, read_reprocessed)) %>%
  left_join(manual_metadata, by = "dataset.ID")

usi_samples <- parse_usi_samples(annotation_rows)

dataset_summary <- annotation_rows %>%
  group_by(dataset.ID, utility_role, collection_label, public_source_hint, source_file) %>%
  summarise(
    annotation_rows = n(),
    non_repeating_annotations = n_distinct(annotation_id),
    unique_compound_names = n_distinct(compound.name[!is.na(compound.name)]),
    unique_smiles = n_distinct(smiles[!is.na(smiles)]),
    median_score = median(score, na.rm = TRUE),
    rt_available_rows = sum(!is.na(rt)),
    .groups = "drop"
  ) %>%
  left_join(
    usi_samples %>%
      group_by(dataset.ID) %>%
      summarise(
        usi_links = n(),
        unique_sample_ids = n_distinct(sample_id),
        public_collection = paste(sort(unique(public_collection)), collapse = "; "),
        inferred_polarity = paste(sort(unique(na.omit(inferred_polarity))), collapse = "; "),
        .groups = "drop"
      ),
    by = "dataset.ID"
  ) %>%
  arrange(dataset.ID)

sample_summary <- usi_samples %>%
  distinct(dataset.ID, sample_id, public_collection, inferred_polarity, nist_subject, nist_replicate, beefdiet_sample_number, beefdiet_fraction) %>%
  arrange(dataset.ID, sample_id)

public_sample_counts <- sample_summary %>%
  filter(dataset.ID %in% public_ids) %>%
  group_by(dataset.ID, public_collection, inferred_polarity) %>%
  summarise(
    unique_sample_ids = n_distinct(sample_id),
    nist_subjects = n_distinct(nist_subject[!is.na(nist_subject)]),
    beefdiet_sample_numbers = n_distinct(beefdiet_sample_number[!is.na(beefdiet_sample_number)]),
    beefdiet_fractions = n_distinct(beefdiet_fraction[!is.na(beefdiet_fraction)]),
    .groups = "drop"
  )

workbook_hits <- scan_hgm_workbook(c(dataset_ids, "MSV", "MassIVE", "MASSIVE", "NIST", "BeefDiet"))

analysis_suggestions <- tibble::tribble(
  ~priority, ~analysis, ~rationale, ~required_inputs,
  1L, "Public-data remine yield", "Show that the atlas-derived library can annotate public faecal LC-MS/MS data without RT matching.", "Reprocessed annotation CSVs for HGMD_0357, HGMD_0358 and HGMD_0359; count non-repeating annotation_id/compound_name and score distributions.",
  2L, "Cross-resource overlap", "Test whether public-remine hits recover metabolites already present in the atlas and identify public-only candidates.", "Standardized compound.name or SMILES from public reprocessing outputs plus paper2_total_list_classified.csv.",
  3L, "Sample-level prevalence", "Use the USI-derived sample IDs to estimate which atlas-library metabolites recur across public samples/subjects.", "query_spectrum_usi sample IDs; for NIST, subject and replicate parsed from NIST_POS_Samp_##-##.",
  4L, "Diet or cohort contrast if metadata are available", "The BeefDiet filenames suggest a diet study; if group metadata can be recovered, compare compound classes or named bile-acid/fatty-acid metabolites by diet.", "Original MassIVE metadata table mapping sample IDs to diet/group and acquisition polarity.",
  5L, "Technical reproducibility in NIST", "NIST sample IDs appear to encode subject and replicate; replicate concordance can be used as a conservative validation of library utility.", "NIST sample IDs, replicate labels, and feature-level abundance/export table if available.",
  6L, "Class-level public remine profile", "Summarize whether public-remine hits are dominated by expected faecal metabolite classes rather than random spectral matches.", "NPClassifier/classified names via existing paper2_total_list_classified mapping or external classifier cache."
)

write_csv_stable(dataset_summary, file.path(table_dir, "figure4-remining-utility-dataset-summary.csv"))
write_csv_stable(sample_summary, file.path(table_dir, "figure4-remining-utility-sample-ids.csv"))
write_csv_stable(public_sample_counts, file.path(table_dir, "figure4-remining-utility-public-sample-counts.csv"))
write_csv_stable(workbook_hits, file.path(qc_dir, "figure4-hgm-workbook-public-dataset-hits.csv"))
write_csv_stable(analysis_suggestions, file.path(table_dir, "figure4-public-dataset-analysis-suggestions.csv"))

plot_data <- dataset_summary %>%
  mutate(
    dataset_label = paste0(dataset.ID, "\n", collection_label),
    utility_role = factor(
      utility_role,
      levels = c("Mine crude atlas datasets", "Remine public datasets without RT matching")
    )
  )

p_utility <- ggplot(plot_data, aes(dataset_label, non_repeating_annotations, fill = utility_role)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = fmt_int(non_repeating_annotations)),
    vjust = -0.35,
    size = 4.6,
    family = "Arial",
    fontface = "bold"
  ) +
  facet_wrap(~utility_role, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = c(
      "Mine crude atlas datasets" = "#404B74",
      "Remine public datasets without RT matching" = "#9A4049"
    ),
    guide = "none"
  ) +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Atlas-derived in-house libraries support crude mining and public-data remining",
    x = NULL,
    y = "Non-repeating annotations"
  ) +
  theme_classic(base_size = 17, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(size = 21, face = "bold"),
    strip.text = element_text(size = 14.5, face = "bold"),
    axis.text.x = element_text(size = 10.5, angle = 25, hjust = 1),
    plot.margin = margin(10, 12, 10, 10)
  )

ggsave(
  file.path(figure_dir, "figure4-inhouse-library-utility-overview.png"),
  p_utility,
  width = 9,
  height = 3.9,
  dpi = 300,
  bg = "white"
)
