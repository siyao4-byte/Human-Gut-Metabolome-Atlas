# Figure 3 data build.
# Purpose: assemble level 1-3, non-MSNovelist feature annotations for the four
# manually defined pools, then create normalized quantitative profiles and
# pool-level presence tables used by the downstream PCA, volcano, and prevalence
# plots.
source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(limma)
  library(vegan)
})

dataset_manifest <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure23-dataset-manifest.csv"),
  show_col_types = FALSE
)
sample_map <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure23-sample-lineage.csv"),
  col_types = readr::cols(fraction_sequence = readr::col_character()),
  show_col_types = FALSE
)
pool_defs <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure2-combined-pool-definitions.csv"),
  show_col_types = FALSE
)
pool_members <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure2-combined-pool-members.csv"),
  show_col_types = FALSE
)

pool_lookup <- setNames(pool_members$combined_pool_id, pool_members$HGMH.ID)
canonical_lookup <- setNames(pool_defs$combined_pool_id, pool_defs$hgmh_ids)

annotation_type_rank <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, stringr::regex("authentic|mzmine", ignore_case = TRUE)) ~ 1L,
    stringr::str_detect(x, stringr::regex("gnps", ignore_case = TRUE)) ~ 2L,
    stringr::str_detect(x, stringr::regex("ms2query", ignore_case = TRUE)) ~ 3L,
    TRUE ~ 99L
  )
}

read_best_annotations <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      feature.ID = as.character(feature.ID),
      confidence.level = suppressWarnings(as.numeric(confidence.level)),
      confidence.score = suppressWarnings(as.numeric(confidence.score)),
      id.prob = suppressWarnings(as.numeric(id.prob)),
      compound.name = na_if(clean_text(compound.name), ""),
      smiles = na_if(clean_text(smiles), ""),
      annotation_id = coalesce(compound.name, smiles),
      annotation_rank = annotation_type_rank(annotation.type)
    ) %>%
    # Count and model non-repeating annotations by annotation_id. This is
    # compound.name when available, otherwise smiles; MSNovelist-only calls and
    # vague analogue/candidate labels are excluded consistently across figures.
    filter(
      confidence.level %in% 1:3,
      !is.na(annotation_id),
      !str_detect(coalesce(annotation.type, ""), regex("^MSNovelist$", ignore_case = TRUE)),
      !str_detect(compound.name, regex("analogue|candidate|PUBCHEM", ignore_case = TRUE))
    ) %>%
    arrange(feature.ID, confidence.level, annotation_rank, desc(id.prob), desc(confidence.score)) %>%
    group_by(feature.ID) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(feature.ID, annotation_id, compound.name, smiles, confidence.level, annotation.type)
}

assign_sample_metadata <- function(dataset_id) {
  # Crude profiles are kept at individual HGMH level; fraction profiles are
  # collapsed to the canonical four-pool definition and blocked by fraction
  # sequence so matching fractions are compared like-for-like.
  sample_map %>%
    filter(dataset.ID == dataset_id) %>%
    mutate(
      hgmh_single = if_else(hgmh_count == 1, hgmh_ids, NA_character_),
      crude_pool = unname(pool_lookup[hgmh_single]),
      canonical_pool = unname(canonical_lookup[hgmh_ids]),
      state = case_when(
        processing_level %in% c("methanol_crude", "dcm_crude", "water_crude") ~ "Crude",
        processing_level %in% c("polar_fraction", "nonpolar_fraction") ~ "Fractionated",
        TRUE ~ NA_character_
      ),
      combined_pool_id = if_else(state == "Crude", crude_pool, canonical_pool),
      profile_id = case_when(
        state == "Crude" ~ hgmh_single,
        state == "Fractionated" ~ paste(combined_pool_id, processing_level, fraction_sequence, sep = "::"),
        TRUE ~ NA_character_
      ),
      block_id = if_else(state == "Fractionated", paste(processing_level, fraction_sequence, sep = "::"), NA_character_)
    ) %>%
    filter(
      !is.na(state), !is.na(combined_pool_id), !is.na(profile_id),
      hgmh_single != "HGMH_0064" | is.na(hgmh_single),
      state != "Fractionated" | str_detect(fraction_sequence, "^\\d+$")
    ) %>%
    select(dataset.ID, HGMA.ID, state, profile_id, combined_pool_id, block_id, processing_level, fraction_sequence)
}

normalise_dataset <- function(dataset_id, annotation_path, abundance_path, column_type, ion_mode) {
  message("Processing Figure 3 dataset ", dataset_id)
  ann <- read_best_annotations(annotation_path)
  meta <- assign_sample_metadata(dataset_id)
  abundance <- readr::read_csv(abundance_path, show_col_types = FALSE, progress = FALSE) %>%
    transmute(
      feature.ID = as.character(feature.ID),
      HGMA.ID = str_remove(clean_text(samples), "\\.raw\\.area$"),
      area = pmax(suppressWarnings(as.numeric(area)), 0, na.rm = TRUE)
    ) %>%
    filter(str_detect(HGMA.ID, "^HGMA_\\d+$"), !is.na(area)) %>%
    group_by(feature.ID, HGMA.ID) %>%
    summarise(area = max(area), .groups = "drop") %>%
    inner_join(ann, by = "feature.ID") %>%
    inner_join(meta, by = "HGMA.ID")

  if (!nrow(abundance)) return(tibble())

  abundance %>%
    group_by(state) %>%
    group_modify(~ {
      samples <- sort(unique(.x$profile_id))
      features <- sort(unique(.x$feature.ID))
      matrix_long <- tidyr::expand_grid(feature.ID = features, profile_id = samples) %>%
        left_join(
          .x %>% group_by(feature.ID, profile_id) %>% summarise(area = max(area), .groups = "drop"),
          by = c("feature.ID", "profile_id")
        ) %>%
        mutate(area = coalesce(area, 0))

      keep_features <- matrix_long %>%
        group_by(feature.ID) %>%
        summarise(zero_prop = mean(area <= 0), .groups = "drop") %>%
        filter(zero_prop <= 0.5) %>%
        pull(feature.ID)
      matrix_long <- matrix_long %>% filter(feature.ID %in% keep_features)
      if (!nrow(matrix_long)) return(tibble())

      feature_min <- matrix_long %>%
        filter(area > 0) %>%
        group_by(feature.ID) %>%
        summarise(min_positive = min(area), .groups = "drop")
      matrix_long <- matrix_long %>%
        left_join(feature_min, by = "feature.ID") %>%
        mutate(area_imputed = if_else(area > 0, area, 0.5 * min_positive))
      sample_medians <- matrix_long %>%
        group_by(profile_id) %>%
        summarise(sample_median = median(area_imputed, na.rm = TRUE), .groups = "drop")
      global_median <- median(sample_medians$sample_median, na.rm = TRUE)

      # Quantitative profiles use half-minimum imputation, sample median
      # normalization to the global median, then log10(x + 1) transformation.
      matrix_long %>%
        left_join(sample_medians, by = "profile_id") %>%
        mutate(value = log10(area_imputed / (sample_median / global_median) + 1)) %>%
        select(feature.ID, profile_id, value)
    }) %>%
    ungroup() %>%
    mutate(dataset.ID = dataset_id, feature_key = paste(dataset_id, feature.ID, sep = "::")) %>%
    left_join(ann, by = "feature.ID") %>%
    left_join(meta %>% distinct(state, profile_id, combined_pool_id, block_id, processing_level), by = c("state", "profile_id")) %>%
    mutate(column.type = column_type, ion.mode = ion_mode) %>%
    distinct()
}

quant_long <- bind_rows(Map(
  normalise_dataset,
  dataset_manifest$HGMD.ID,
  dataset_manifest$annotation_path,
  dataset_manifest$abundance_path,
  dataset_manifest$column.type,
  dataset_manifest$ion.mode
))

presence <- readr::read_csv(
  file.path(project_root, "data", "intermediate", "figure2-four-pool-presence.csv"),
  col_types = readr::cols(fraction_sequence = readr::col_character()),
  show_col_types = FALSE
) %>%
  # Presence summaries use the same four-pool lineage table as Figure 2, but
  # only crude and first-level fractionation states are retained for Figure 3.
  filter(figure2_processing_level %in% c("methanol_crude", "dcm_crude", "water_crude", "polar_fraction", "nonpolar_fraction")) %>%
  mutate(state = if_else(str_detect(figure2_processing_level, "fraction"), "Fractionated", "Crude")) %>%
  distinct(annotation_id, combined_pool_id, state, compound.name, smiles, confidence.level)

write_csv_stable(quant_long, file.path(project_root, "data", "intermediate", "figure3-abundance-long.csv"))
write_csv_stable(presence, file.path(project_root, "data", "intermediate", "figure3-pool-presence.csv"))

qc <- tibble(
  metric = c("quantitative_rows", "quantitative_features", "crude_profiles", "fractionated_profiles", "presence_annotations"),
  value = c(
    nrow(quant_long), n_distinct(quant_long$feature_key),
    n_distinct(quant_long$profile_id[quant_long$state == "Crude"]),
    n_distinct(quant_long$profile_id[quant_long$state == "Fractionated"]),
    n_distinct(presence$annotation_id)
  )
)
write_csv_stable(qc, file.path(project_root, "outputs", "qc", "figure3-data-summary.csv"))
print(qc)
