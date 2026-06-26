source(file.path("scripts", "00_setup", "helpers.R"))

dataset_ids <- c("HGMD_0222", "HGMD_0223", "HGMD_0230", "HGMD_0231")
manifest <- readr::read_csv(
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.current.csv"),
  show_col_types = FALSE
) %>%
  filter(HGMD.ID %in% dataset_ids)

if (nrow(manifest) != length(dataset_ids)) {
  stop("Expected all four BH01-31/BH01-40 datasets in the current manifest.")
}

annotation_type_rank <- function(x) {
  case_when(
    str_detect(x, regex("authentic|mzmine", ignore_case = TRUE)) ~ 1L,
    str_detect(x, regex("gnps", ignore_case = TRUE)) ~ 2L,
    str_detect(x, regex("ms2query", ignore_case = TRUE)) ~ 3L,
    TRUE ~ 99L
  )
}

read_collection_dataset <- function(dataset_id, annotation_path, abundance_path) {
  annotations <- readr::read_csv(
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
    select(feature.ID, annotation_id, compound.name, smiles, confidence.level)

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
    inner_join(annotations, by = "feature.ID") %>%
    mutate(dataset.ID = dataset_id)
}

if (all(file.exists(manifest$annotation_path)) && all(file.exists(manifest$abundance_path))) {
  combined <- bind_rows(Map(
    read_collection_dataset,
    manifest$HGMD.ID,
    manifest$annotation_path,
    manifest$abundance_path
  ))
} else {
  existing_membership_path <- file.path(
    project_root, "outputs", "tables", "figure2e-bh01-31-collection-annotation-membership.csv"
  )
  if (!file.exists(existing_membership_path)) {
    stop("BH01 collection source files and the audited local membership fallback are both unavailable.")
  }
  message("External BH01 source files unavailable; rebuilding from the audited local collection membership.")
  combined <- readr::read_csv(existing_membership_path, show_col_types = FALSE) %>%
    filter(collection %in% c("Complete BH01-31 collection", "Complete BH01-40 collection")) %>%
    select(dataset.ID, HGMA.ID, annotation_id, compound.name, smiles, confidence.level) %>%
    distinct()
}

collection_definitions <- tibble(
  collection = c(
    "BH01-31 crude reference",
    "Complete BH01-31 collection",
    "BH01-31 fractions 18 and 19",
    "Complete BH01-40 collection"
  ),
  display_level = c(
    "Crude/reference\nsample",
    "Complete primary\nfraction collection",
    "Selected primary\nfractions 18 + 19",
    "Complete C18\nsubfraction collection"
  ),
  dataset_ids = c(
    "HGMD_0222; HGMD_0223",
    "HGMD_0222; HGMD_0223",
    "HGMD_0222; HGMD_0223",
    "HGMD_0230; HGMD_0231"
  ),
  sample_scope = c("HGMA_4255", "all HGMA samples", "parents HGMF_3451 and HGMF_3452", "all HGMA samples")
)

analysis_map <- readxl::read_excel(paths$hgm_workbook, sheet = "A - Analysis") %>%
  mutate(across(everything(), clean_text)) %>%
  separate_rows(dataset.ID, sep = ";\\s*")
selected_fraction_hgma <- analysis_map %>%
  filter(
    dataset.ID %in% c("HGMD_0222", "HGMD_0223"),
    str_detect(Parent.Code, "HGMF_3451|HGMF_3452")
  ) %>%
  pull(HGMA.ID) %>%
  unique()

collection_membership <- bind_rows(
  combined %>%
    filter(dataset.ID %in% c("HGMD_0222", "HGMD_0223"), HGMA.ID == "HGMA_4255") %>%
    mutate(collection = "BH01-31 crude reference"),
  combined %>%
    filter(dataset.ID %in% c("HGMD_0222", "HGMD_0223")) %>%
    mutate(collection = "Complete BH01-31 collection"),
  combined %>%
    filter(
      dataset.ID %in% c("HGMD_0222", "HGMD_0223"),
      HGMA.ID %in% selected_fraction_hgma
    ) %>%
    mutate(collection = "BH01-31 fractions 18 and 19"),
  combined %>%
    filter(dataset.ID %in% c("HGMD_0230", "HGMD_0231")) %>%
    mutate(collection = "Complete BH01-40 collection")
)

collection_counts <- collection_membership %>%
  group_by(collection) %>%
  summarise(
    collection_annotations = n_distinct(annotation_id),
    analysis_samples = n_distinct(HGMA.ID),
    .groups = "drop"
  ) %>%
  left_join(collection_definitions, by = "collection") %>%
  mutate(
    collection = factor(collection, levels = collection_definitions$collection),
    relative_to_crude_percent = 100 * (
      collection_annotations / collection_annotations[collection == "BH01-31 crude reference"] - 1
    )
  ) %>%
  arrange(collection)

write_csv_stable(
  collection_counts,
  file.path(project_root, "outputs", "tables", "figure2e-bh01-31-collection-counts.csv")
)
write_csv_stable(
  collection_membership %>%
    distinct(collection, dataset.ID, HGMA.ID, annotation_id, compound.name, smiles, confidence.level) %>%
    arrange(collection, dataset.ID, HGMA.ID, annotation_id),
  file.path(project_root, "outputs", "tables", "figure2e-bh01-31-collection-annotation-membership.csv")
)

print(collection_counts)
