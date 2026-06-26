source(file.path("scripts", "00_setup", "helpers.R"))

x <- readr::read_csv(
  file.path(project_root, "data", "intermediate", "figure2-four-pool-presence.csv"),
  col_types = readr::cols(fraction_sequence = readr::col_character()),
  show_col_types = FALSE
)

paper_cache_path <- file.path(project_root, "data", "processed", "paper2-npclassifier-cache.csv")
cache_path <- if (file.exists(paper_cache_path)) paper_cache_path else paths$npclassifier_cache
npc <- readr::read_csv(cache_path, show_col_types = FALSE, progress = FALSE)
npc_small <- npc %>%
  transmute(
    smiles = clean_text(if ("SMILES" %in% names(npc)) SMILES else smiles),
    npc_superclass = clean_text(if ("np_superclass" %in% names(npc)) np_superclass else superclass),
    npc_class = clean_text(if ("np_class" %in% names(npc)) np_class else subclass)
  ) %>%
  filter(smiles != "") %>%
  distinct(smiles, .keep_all = TRUE)

x <- x %>%
  mutate(smiles = clean_text(smiles)) %>%
  left_join(npc_small, by = "smiles") %>%
  mutate(
    npc_superclass = if_else(is.na(npc_superclass) | npc_superclass == "", "Unclassified", npc_superclass),
    npc_class = if_else(is.na(npc_class) | npc_class == "", "Unclassified", npc_class)
  )

annotation_gain <- x %>%
  distinct(annotation_id, combined_pool_id, figure2_processing_level, confidence.level) %>%
  count(combined_pool_id, figure2_processing_level, confidence.level, name = "unique_annotations")

processing_presence <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  group_by(annotation_id) %>%
  mutate(processing_level_count = n()) %>%
  ungroup()

unique_annotations <- processing_presence %>%
  filter(processing_level_count == 1) %>%
  select(-processing_level_count) %>%
  left_join(
    x %>% distinct(annotation_id, compound.name, smiles, npc_superclass, npc_class),
    by = "annotation_id"
  ) %>%
  distinct()

class_enrichment <- x %>%
  filter(npc_superclass != "Unclassified") %>%
  distinct(annotation_id, figure2_processing_level, npc_superclass) %>%
  count(figure2_processing_level, npc_superclass, name = "unique_annotations") %>%
  group_by(npc_superclass) %>%
  mutate(total = sum(unique_annotations)) %>%
  ungroup() %>%
  slice_max(total, n = 30, with_ties = TRUE) %>%
  select(-total)

unique_class <- unique_annotations %>%
  filter(npc_superclass != "Unclassified") %>%
  count(figure2_processing_level, npc_superclass, name = "unique_annotations")

write_csv_stable(annotation_gain, file.path(project_root, "outputs", "tables", "figure2-annotation-gain.csv"))
write_csv_stable(unique_annotations, file.path(project_root, "outputs", "tables", "figure2-processing-level-unique-annotations.csv"))
write_csv_stable(class_enrichment, file.path(project_root, "outputs", "tables", "figure2-class-enrichment.csv"))
write_csv_stable(unique_class, file.path(project_root, "outputs", "tables", "figure2-unique-class-identity.csv"))
