source(file.path("scripts", "00_setup", "helpers.R"))

dataset_manifest <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure23-dataset-manifest.csv"),
  show_col_types = FALSE
)
sample_map <- readr::read_csv(
  file.path(project_root, "data", "metadata", "figure23-sample-lineage.csv"),
  show_col_types = FALSE
)

split_ids_local <- function(x) {
  x <- clean_text(x)
  x <- x[x != ""]
  if (!length(x)) return(character())
  clean_text(unlist(str_split(x, ";\\s*"), use.names = FALSE))
}

pool_defs <- sample_map %>%
  filter(str_detect(processing_level, "fraction"), hgmh_ids != "") %>%
  count(hgmh_ids, sort = TRUE, name = "lineage_records") %>%
  slice_head(n = 4) %>%
  arrange(hgmh_ids) %>%
  mutate(combined_pool_id = paste0("Pool_", row_number()))

pool_members <- pool_defs %>%
  separate_rows(hgmh_ids, sep = ";\\s*") %>%
  transmute(combined_pool_id, HGMH.ID = clean_text(hgmh_ids)) %>%
  distinct()

assign_pool <- function(hgmh_ids) {
  ids <- split_ids_local(hgmh_ids)
  hit <- pool_members$combined_pool_id[pool_members$HGMH.ID %in% ids]
  hit <- sort(unique(hit))
  if (length(hit) == 1) hit else NA_character_
}

sample_map <- sample_map %>%
  mutate(
    combined_pool_id = vapply(hgmh_ids, assign_pool, character(1)),
    figure2_processing_level = case_when(
      processing_level == "polar_fraction" & str_starts(Parent.Code, "HGME_") ~ "polar_batch_crude",
      processing_level == "nonpolar_fraction" & str_starts(Parent.Code, "HGME_") ~ "nonpolar_batch_crude",
      TRUE ~ processing_level
    )
  )

annotation_type_rank <- function(x) {
  case_when(
    str_detect(x, regex("authentic|mzmine", ignore_case = TRUE)) ~ 1L,
    str_detect(x, regex("gnps", ignore_case = TRUE)) ~ 2L,
    str_detect(x, regex("ms2query", ignore_case = TRUE)) ~ 3L,
    TRUE ~ 99L
  )
}

process_dataset <- function(dataset_id, annotation_path, abundance_path, column_type, ion_mode) {
  message("Combining ", dataset_id)
  ann <- readr::read_csv(
    annotation_path,
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
    filter(
      !is.na(confidence.level), confidence.level <= 3,
      !is.na(annotation_id),
      !str_detect(coalesce(annotation.type, ""), regex("^MSNovelist$", ignore_case = TRUE)),
      !str_detect(compound.name, regex("analogue|candidate|PUBCHEM", ignore_case = TRUE))
    ) %>%
    arrange(annotation_id, confidence.level, annotation_rank, desc(id.prob), desc(confidence.score)) %>%
    group_by(annotation_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(
      feature.ID, annotation_id, compound.name, smiles, confidence.level,
      confidence.score, annotation.type
    )

  abundance <- readr::read_csv(abundance_path, show_col_types = FALSE, progress = FALSE) %>%
    transmute(
      feature.ID = as.character(feature.ID),
      HGMA.ID = str_remove(clean_text(samples), "\\.raw\\.area$"),
      area = suppressWarnings(as.numeric(area))
    ) %>%
    filter(str_detect(HGMA.ID, "^HGMA_\\d+$"), !is.na(area), area > 0) %>%
    group_by(feature.ID, HGMA.ID) %>%
    summarise(area = max(area), .groups = "drop")

  abundance %>%
    inner_join(ann, by = "feature.ID") %>%
    mutate(dataset.ID = dataset_id) %>%
    left_join(
      sample_map %>%
        filter(dataset.ID == dataset_id) %>%
        select(
          dataset.ID, HGMA.ID, combined_pool_id, figure2_processing_level,
          processing_level, pool_id, hgmh_ids, donor_ids, Parent.Code,
          fraction_sequence
        ),
      by = c("dataset.ID", "HGMA.ID")
    ) %>%
    mutate(column.type = column_type, ion.mode = ion_mode) %>%
    filter(!is.na(combined_pool_id))
}

presence <- bind_rows(Map(
  process_dataset,
  dataset_manifest$HGMD.ID,
  dataset_manifest$annotation_path,
  dataset_manifest$abundance_path,
  dataset_manifest$column.type,
  dataset_manifest$ion.mode
))

presence_distinct <- presence %>%
  arrange(annotation_id, combined_pool_id, figure2_processing_level, confidence.level) %>%
  group_by(annotation_id, combined_pool_id, figure2_processing_level) %>%
  slice_head(n = 1) %>%
  ungroup()

qc <- tibble(
  metric = c(
    "four_pool_definitions", "presence_rows", "distinct_pool_processing_annotations",
    "unmapped_pool_rows", "processing_levels"
  ),
  value = c(
    nrow(pool_defs), nrow(presence), nrow(presence_distinct),
    sum(is.na(presence$combined_pool_id)), n_distinct(presence_distinct$figure2_processing_level)
  )
)

write_csv_stable(pool_defs, file.path(project_root, "data", "metadata", "figure2-combined-pool-definitions.csv"))
write_csv_stable(pool_members, file.path(project_root, "data", "metadata", "figure2-combined-pool-members.csv"))
write_csv_stable(presence_distinct, file.path(project_root, "data", "intermediate", "figure2-four-pool-presence.csv"))
write_csv_stable(qc, file.path(project_root, "outputs", "qc", "figure2-presence-summary.csv"))
print(pool_defs)
print(qc)
