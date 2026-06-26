source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

nist_dir <- file.path(project_root, "outputs", "mzmine reprocessing", "HGMD_0359")
metadata_path <- file.path(nist_dir, "gnps_metadata.tsv")
annotation_path <- file.path(nist_dir, "data_annotations_NIST.csv")
quant_path <- file.path(nist_dir, "data_iimn_gnps_NIST_quant.csv")
out_dir <- file.path(project_root, "outputs", "tables")

metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    filename = clean_text(filename),
    diet = clean_text(ATTRIBUTE_diet_info),
    sample_type = clean_text(SampleType),
    subject = str_match(filename, "Samp_([0-9]{2})-")[, 2],
    replicate = str_match(filename, "Samp_[0-9]{2}-([0-9]{2})")[, 2],
    quant_col = paste0(str_remove(filename, "\\.mzXML$"), ".mzXML Peak area")
  )

nist_samples <- metadata %>%
  filter(sample_type == "animal", !is.na(subject), diet %in% c("omnivore", "vegetarian")) %>%
  distinct(filename, quant_col, subject, replicate, diet, Sample_ID_NIST, UniqueSubjectID)

metadata_summary <- bind_rows(
  metadata %>%
    count(sample_type, diet, name = "n_files") %>%
    mutate(summary_level = "files"),
  nist_samples %>%
    distinct(subject, diet) %>%
    count(sample_type = "animal_subject", diet, name = "n_files") %>%
    mutate(summary_level = "subjects")
) %>%
  select(summary_level, sample_type, diet, n_files)

annotations <- readr::read_csv(annotation_path, show_col_types = FALSE) %>%
  transmute(
    feature_id = as.character(id),
    compound_name = clean_text(compound_name),
    smiles = clean_text(smiles),
    atlas_score = suppressWarnings(as.numeric(score)),
    rt = suppressWarnings(as.numeric(rt)),
    precursor_mz = suppressWarnings(as.numeric(precursor_mz))
  ) %>%
  filter(compound_name != "") %>%
  group_by(feature_id, compound_name) %>%
  summarise(
    smiles = first(smiles[smiles != ""]),
    atlas_score = max(atlas_score, na.rm = TRUE),
    rt = first(rt),
    precursor_mz = first(precursor_mz),
    .groups = "drop"
  ) %>%
  mutate(atlas_score = if_else(is.infinite(atlas_score), NA_real_, atlas_score))

quant <- readr::read_csv(quant_path, show_col_types = FALSE) %>%
  rename(feature_id = `row ID`) %>%
  mutate(feature_id = as.character(feature_id))

sample_cols <- intersect(nist_samples$quant_col, names(quant))
if (length(sample_cols) == 0) {
  stop("No NIST sample columns in quant table matched gnps_metadata.tsv filenames.")
}

quant_long <- quant %>%
  select(feature_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "quant_col",
    values_to = "peak_area"
  ) %>%
  left_join(nist_samples, by = "quant_col") %>%
  mutate(
    peak_area = suppressWarnings(as.numeric(peak_area)),
    present = !is.na(peak_area) & peak_area > 0,
    log10_area = log10(coalesce(peak_area, 0) + 1)
  )

sample_normalisation <- quant_long %>%
  group_by(quant_col) %>%
  summarise(
    sample_median_log10_area = median(log10_area[peak_area > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_median_log10_area = if_else(
      is.finite(sample_median_log10_area),
      sample_median_log10_area,
      median(sample_median_log10_area[is.finite(sample_median_log10_area)], na.rm = TRUE)
    ),
    global_median_log10_area = median(sample_median_log10_area, na.rm = TRUE)
  )

quant_long <- quant_long %>%
  left_join(sample_normalisation, by = "quant_col") %>%
  mutate(
    normalised_log10_area = log10_area - sample_median_log10_area + global_median_log10_area
  )

feature_subject <- quant_long %>%
  group_by(feature_id, subject, diet) %>%
  summarise(
    replicate_n = n(),
    present_replicates = sum(present, na.rm = TRUE),
    subject_present = present_replicates > 0,
    mean_log10_area = mean(log10_area, na.rm = TRUE),
    mean_normalised_log10_area = mean(normalised_log10_area, na.rm = TRUE),
    median_log10_area = median(log10_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(annotations, by = "feature_id", relationship = "many-to-many")

safe_t <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  if (sum(ok) < 4 || n_distinct(group[ok]) != 2) {
    return(tibble(p_value = NA_real_, statistic = NA_real_))
  }
  group_sizes <- table(group[ok])
  if (any(group_sizes < 2)) {
    return(tibble(p_value = NA_real_, statistic = NA_real_))
  }
  out <- tryCatch(t.test(value[ok] ~ group[ok]), error = function(e) NULL)
  if (is.null(out)) {
    tibble(p_value = NA_real_, statistic = NA_real_)
  } else {
    tibble(p_value = out$p.value, statistic = unname(out$statistic))
  }
}

safe_fisher <- function(present, group) {
  ok <- !is.na(present) & !is.na(group)
  if (sum(ok) == 0 || n_distinct(group[ok]) != 2) {
    return(NA_real_)
  }
  tab <- table(group[ok], present[ok])
  if (nrow(tab) != 2 || ncol(tab) != 2) {
    return(NA_real_)
  }
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

feature_stats <- feature_subject %>%
  group_by(feature_id, compound_name, smiles, atlas_score, rt, precursor_mz) %>%
  summarise(
    n_subjects = n_distinct(subject),
    omnivore_subjects = n_distinct(subject[diet == "omnivore"]),
    vegetarian_subjects = n_distinct(subject[diet == "vegetarian"]),
    omnivore_prevalence = sum(subject_present[diet == "omnivore"], na.rm = TRUE),
    vegetarian_prevalence = sum(subject_present[diet == "vegetarian"], na.rm = TRUE),
    omnivore_mean_normalised_log10_area = mean(mean_normalised_log10_area[diet == "omnivore"], na.rm = TRUE),
    vegetarian_mean_normalised_log10_area = mean(mean_normalised_log10_area[diet == "vegetarian"], na.rm = TRUE),
    normalised_log10_area_difference_omnivore_minus_vegetarian =
      omnivore_mean_normalised_log10_area - vegetarian_mean_normalised_log10_area,
    abundance_test = list(safe_t(mean_normalised_log10_area, diet)),
    prevalence_p_value = safe_fisher(subject_present, diet),
    .groups = "drop"
  ) %>%
  unnest(abundance_test) %>%
  mutate(
    abundance_fdr = p.adjust(p_value, method = "BH"),
    prevalence_fdr = p.adjust(prevalence_p_value, method = "BH")
  ) %>%
  arrange(abundance_fdr, desc(abs(normalised_log10_area_difference_omnivore_minus_vegetarian)))

compound_subject <- feature_subject %>%
  group_by(compound_name, smiles, subject, diet) %>%
  summarise(
    subject_present = any(subject_present),
    mean_normalised_log10_area = max(mean_normalised_log10_area, na.rm = TRUE),
    best_atlas_score = max(atlas_score, na.rm = TRUE),
    feature_n = n_distinct(feature_id),
    .groups = "drop"
  ) %>%
  mutate(
    mean_normalised_log10_area = if_else(is.infinite(mean_normalised_log10_area), NA_real_, mean_normalised_log10_area),
    best_atlas_score = if_else(is.infinite(best_atlas_score), NA_real_, best_atlas_score)
  )

compound_stats <- compound_subject %>%
  group_by(compound_name, smiles) %>%
  summarise(
    feature_n = max(feature_n, na.rm = TRUE),
    best_atlas_score = max(best_atlas_score, na.rm = TRUE),
    n_subjects = n_distinct(subject),
    omnivore_prevalence = sum(subject_present[diet == "omnivore"], na.rm = TRUE),
    vegetarian_prevalence = sum(subject_present[diet == "vegetarian"], na.rm = TRUE),
    omnivore_mean_normalised_log10_area = mean(mean_normalised_log10_area[diet == "omnivore"], na.rm = TRUE),
    vegetarian_mean_normalised_log10_area = mean(mean_normalised_log10_area[diet == "vegetarian"], na.rm = TRUE),
    normalised_log10_area_difference_omnivore_minus_vegetarian =
      omnivore_mean_normalised_log10_area - vegetarian_mean_normalised_log10_area,
    abundance_test = list(safe_t(mean_normalised_log10_area, diet)),
    prevalence_p_value = safe_fisher(subject_present, diet),
    .groups = "drop"
  ) %>%
  unnest(abundance_test) %>%
  mutate(
    abundance_fdr = p.adjust(p_value, method = "BH"),
    prevalence_fdr = p.adjust(prevalence_p_value, method = "BH")
  ) %>%
  arrange(abundance_fdr, desc(abs(normalised_log10_area_difference_omnivore_minus_vegetarian)))

diet_overview <- tibble(
  metric = c(
    "metadata_animal_files",
    "metadata_subjects",
    "omnivore_subjects",
    "vegetarian_subjects",
    "technical_replicates_per_subject",
    "quant_sample_columns_matched",
    "atlas_matched_feature_ids_tested",
    "feature_compound_pairs_tested",
    "non_repeating_compounds_tested",
    "feature_abundance_fdr_lt_0_05",
    "compound_abundance_fdr_lt_0_05",
    "feature_prevalence_fdr_lt_0_05",
    "compound_prevalence_fdr_lt_0_05"
  ),
  value = c(
    nrow(nist_samples),
    n_distinct(nist_samples$subject),
    n_distinct(nist_samples$subject[nist_samples$diet == "omnivore"]),
    n_distinct(nist_samples$subject[nist_samples$diet == "vegetarian"]),
    paste(sort(unique(nist_samples %>% count(subject) %>% pull(n))), collapse = ";"),
    length(sample_cols),
    n_distinct(feature_stats$feature_id),
    nrow(feature_stats),
    n_distinct(compound_stats$compound_name),
    sum(feature_stats$abundance_fdr < 0.05, na.rm = TRUE),
    sum(compound_stats$abundance_fdr < 0.05, na.rm = TRUE),
    sum(feature_stats$prevalence_fdr < 0.05, na.rm = TRUE),
    sum(compound_stats$prevalence_fdr < 0.05, na.rm = TRUE)
  )
)

top_compound_stats <- compound_stats %>%
  filter(is.finite(normalised_log10_area_difference_omnivore_minus_vegetarian)) %>%
  arrange(abundance_fdr, desc(abs(normalised_log10_area_difference_omnivore_minus_vegetarian))) %>%
  slice_head(n = 50)

write_csv_stable(metadata_summary, file.path(out_dir, "figure4-nist-metadata-summary.csv"))
write_csv_stable(diet_overview, file.path(out_dir, "figure4-nist-diet-atlas-match-overview.csv"))
write_csv_stable(feature_stats, file.path(out_dir, "figure4-nist-diet-feature-stats.csv"))
write_csv_stable(compound_stats, file.path(out_dir, "figure4-nist-diet-compound-stats.csv"))
write_csv_stable(top_compound_stats, file.path(out_dir, "figure4-nist-diet-top-compound-stats.csv"))

message("NIST metadata diet analysis complete.")
message("Subjects: ", n_distinct(nist_samples$subject), " (",
        n_distinct(nist_samples$subject[nist_samples$diet == "omnivore"]), " omnivore, ",
        n_distinct(nist_samples$subject[nist_samples$diet == "vegetarian"]), " vegetarian).")
message("Atlas-matched features tested: ", n_distinct(feature_stats$feature_id))
message("Feature-compound pairs tested: ", nrow(feature_stats))
message("Non-repeating compounds tested: ", n_distinct(compound_stats$compound_name))
message("Compound abundance FDR < 0.05: ", sum(compound_stats$abundance_fdr < 0.05, na.rm = TRUE))
