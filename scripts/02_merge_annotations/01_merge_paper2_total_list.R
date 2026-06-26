source(file.path("scripts", "00_setup", "helpers.R"))

manifest_path <- file.path(project_root, "data", "metadata", "paper2-dataset-manifest.csv")
manifest_fallback_path <- file.path(
  project_root, "data", "metadata", "paper2-dataset-manifest.current.csv"
)
manifest_candidates <- c(manifest_path, manifest_fallback_path)
manifest_candidates <- manifest_candidates[file.exists(manifest_candidates)]
if (length(manifest_candidates)) {
  manifest_path <- manifest_candidates[which.max(file.info(manifest_candidates)$mtime)]
}
if (!file.exists(manifest_path)) stop("Run scripts/01_metadata/01_build_dataset_manifests.R first.")
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

required_cols <- c(
  "feature.ID", "rt", "mz", "compound.name", "smiles", "annotation.type",
  "confidence.level", "confidence.score", "id.prob", "CID", "HMDB.ID",
  "Formula", "IUPAC", "Monoisotopic.Mass", "feature.usi", "Samples"
)

read_annotation <- function(dataset_id, annotation_path, column_type, ion_mode, targeted_flag) {
  message("Reading ", dataset_id)
  x <- readr::read_csv(annotation_path, show_col_types = FALSE, progress = FALSE)
  missing <- setdiff(required_cols, names(x))
  for (col in missing) x[[col]] <- NA
  x %>%
    mutate(
      across(any_of(required_cols), as.character),
      dataset.ID = dataset_id,
      column.type = clean_text(column_type),
      ion.mode = clean_text(ion_mode),
      targeted = as.logical(targeted_flag),
      source_row = row_number(),
      confidence.level.numeric = suppressWarnings(as.numeric(confidence.level)),
      confidence.score.numeric = suppressWarnings(as.numeric(confidence.score)),
      id.prob.numeric = suppressWarnings(as.numeric(id.prob)),
      compound.name = na_if(clean_text(compound.name), ""),
      smiles = na_if(clean_text(smiles), ""),
      annotation_id = coalesce(compound.name, smiles)
    )
}

raw <- bind_rows(Map(
  read_annotation,
  manifest$HGMD.ID,
  manifest$annotation_path,
  manifest$column.type,
  manifest$ion.mode,
  manifest$targeted_flag
))

stage_counts <- tibble(stage = "loaded", rows = nrow(raw))

eligible <- raw %>%
  mutate(
    exclusion_reason = case_when(
      is.na(confidence.level.numeric) ~ "missing confidence level",
      confidence.level.numeric > 3 ~ "confidence level above 3",
      is.na(annotation_id) ~ "missing compound name and SMILES",
      str_detect(coalesce(annotation.type, ""), regex("^MSNovelist$", ignore_case = TRUE)) ~
        "MSNovelist annotation source",
      str_detect(compound.name, regex("analogue|candidate", ignore_case = TRUE)) ~
        "analogue/candidate compound name",
      str_detect(compound.name, regex("PUBCHEM", ignore_case = TRUE)) ~
        "PUBCHEM placeholder compound name",
      TRUE ~ NA_character_
    )
  )

excluded_initial <- eligible %>% filter(!is.na(exclusion_reason))
eligible <- eligible %>% filter(is.na(exclusion_reason))
stage_counts <- bind_rows(stage_counts, tibble(stage = "eligible_level_1_to_3", rows = nrow(eligible)))

provenance <- eligible %>%
  group_by(annotation_id, column.type, ion.mode) %>%
  summarise(
    source_dataset_ids = paste(sort(unique(dataset.ID)), collapse = "; "),
    source_dataset_count = n_distinct(dataset.ID),
    source_feature_count = n(),
    .groups = "drop"
  )

ranked <- eligible %>%
  mutate(
    platform_priority = case_when(
      str_detect(column.type, regex("^C18", ignore_case = TRUE)) ~ 1L,
      column.type %in% c("Phe-Hex", "HILIC", "SAX") ~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  arrange(
    annotation_id, column.type, ion.mode,
    confidence.level.numeric,
    desc(id.prob.numeric), desc(confidence.score.numeric),
    dataset.ID, source_row
  ) %>%
  group_by(annotation_id, column.type, ion.mode) %>%
  mutate(best_within_platform = row_number() == 1L) %>%
  ungroup()

excluded_within <- ranked %>%
  filter(!best_within_platform) %>%
  mutate(exclusion_reason = "lower-ranked duplicate within column type and ion mode")

best_platform <- ranked %>% filter(best_within_platform)
stage_counts <- bind_rows(stage_counts, tibble(stage = "best_within_platform", rows = nrow(best_platform)))

with_priority <- best_platform %>%
  group_by(annotation_id, ion.mode) %>%
  mutate(
    c18_available = any(str_detect(column.type, regex("^C18", ignore_case = TRUE))),
    keep_platform_priority = !c18_available |
      str_detect(column.type, regex("^C18", ignore_case = TRUE))
  ) %>%
  ungroup()

excluded_priority <- with_priority %>%
  filter(!keep_platform_priority) %>%
  mutate(exclusion_reason = "C18-priority duplicate across column types")

paper2_total <- with_priority %>%
  filter(keep_platform_priority) %>%
  left_join(provenance, by = c("annotation_id", "column.type", "ion.mode")) %>%
  mutate(
    confidence.level = confidence.level.numeric,
    confidence.score = confidence.score.numeric,
    id.prob = id.prob.numeric
  ) %>%
  select(
    annotation_id, compound.name, smiles, confidence.level,
    confidence.score, id.prob,
    annotation.type, feature.ID, feature.usi, rt, mz, Formula, IUPAC,
    Monoisotopic.Mass, CID, HMDB.ID, dataset.ID, targeted, column.type, ion.mode,
    source_dataset_ids, source_dataset_count, source_feature_count, Samples,
    everything(), -exclusion_reason, -confidence.level.numeric,
    -confidence.score.numeric, -id.prob.numeric
  ) %>%
  arrange(confidence.level, annotation_id, column.type, ion.mode)

stage_counts <- bind_rows(stage_counts, tibble(stage = "paper2_total_list", rows = nrow(paper2_total)))
exclusions <- bind_rows(excluded_initial, excluded_within, excluded_priority) %>%
  transmute(
    dataset.ID, source_row, feature.ID, annotation_id, compound.name, smiles,
    confidence.level = confidence.level.numeric, column.type, ion.mode, exclusion_reason
  )

write_csv_stable(paper2_total, file.path(project_root, "data", "processed", "paper2_total_list.csv"))
write_csv_stable(exclusions, file.path(project_root, "outputs", "qc", "paper2-total-list-exclusions.csv"))
write_csv_stable(stage_counts, file.path(project_root, "outputs", "qc", "paper2-total-list-stage-counts.csv"))
print(stage_counts)
