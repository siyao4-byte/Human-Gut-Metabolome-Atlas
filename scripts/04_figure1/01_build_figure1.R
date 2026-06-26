source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages(library(ggplot2))

classified_fallback_path <- file.path(
  project_root, "data", "processed", "paper2_total_list_classified.current.csv"
)
classified_path <- file.path(project_root, "data", "processed", "paper2_total_list_classified.csv")
classified_candidates <- c(classified_path, classified_fallback_path)
classified_candidates <- classified_candidates[file.exists(classified_candidates)]
input_path <- classified_candidates[which.max(file.info(classified_candidates)$mtime)]
x <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    confidence.level = as.integer(confidence.level),
    platform = paste(column.type, ion.mode, sep = " / "),
    annotation_source = if_else(
      is.na(annotation.type) | annotation.type == "", "Unspecified", annotation.type
    ),
    superclass = if_else(
      is.na(npc_superclass) | npc_superclass == "", "Unclassified", npc_superclass
    ),
    subclass = if_else(
      is.na(npc_class) | npc_class == "", "Unclassified", npc_class
    )
  )

best_annotation <- x %>%
  arrange(annotation_id, confidence.level, desc(confidence.score)) %>%
  group_by(annotation_id) %>%
  slice_head(n = 1) %>%
  ungroup()

confidence_counts <- best_annotation %>%
  count(confidence.level, name = "unique_annotations")

source_counts <- best_annotation %>%
  count(annotation_source, sort = TRUE, name = "unique_annotations")

composition <- x %>%
  filter(
    confidence.level <= 3,
    smiles != "",
    superclass != "Unclassified",
    subclass != "Unclassified",
    !str_detect(str_to_lower(superclass), "^other"),
    !str_detect(str_to_lower(subclass), "^other")
  ) %>%
  distinct(smiles, superclass, subclass) %>%
  count(superclass, subclass, name = "unique_annotations") %>%
  arrange(desc(unique_annotations))

top_classes <- x %>%
  filter(confidence.level <= 3, subclass != "Unclassified") %>%
  distinct(annotation_id, subclass) %>%
  count(subclass, sort = TRUE) %>%
  slice_head(n = 30) %>%
  pull(subclass)

platform_class <- x %>%
  filter(confidence.level <= 3, subclass %in% top_classes) %>%
  distinct(annotation_id, platform, subclass) %>%
  count(platform, subclass, name = "unique_annotations")

platform_summary <- x %>%
  filter(confidence.level <= 3, subclass != "Unclassified") %>%
  distinct(annotation_id, platform) %>%
  group_by(annotation_id) %>%
  mutate(platform_count = n()) %>%
  ungroup() %>%
  group_by(platform) %>%
  summarise(
    total_unique_annotations = n_distinct(annotation_id),
    platform_only_annotations = n_distinct(annotation_id[platform_count == 1]),
    .groups = "drop"
  ) %>%
  arrange(desc(total_unique_annotations))

write_csv_stable(confidence_counts, file.path(project_root, "outputs", "tables", "figure1b-confidence-level-counts.csv"))
write_csv_stable(source_counts, file.path(project_root, "outputs", "tables", "figure1-annotation-source-counts.csv"))
write_csv_stable(composition, file.path(project_root, "outputs", "tables", "figure1c-chemical-composition.csv"))
write_csv_stable(platform_class, file.path(project_root, "outputs", "tables", "figure1-platform-class-heatmap.csv"))
write_csv_stable(platform_summary, file.path(project_root, "outputs", "tables", "figure1-platform-total-unique.csv"))

p_confidence <- ggplot(confidence_counts, aes(factor(confidence.level), unique_annotations)) +
  geom_col(fill = uom_colours[["heritage"]], width = 0.72) +
  geom_text(aes(label = unique_annotations), vjust = -0.3, size = 5.4) +
  labs(x = "Annotation confidence level", y = "Non-repeating annotations") +
  theme_minimal(base_size = 19) +
  theme(panel.grid.major.x = element_blank())

p_source <- source_counts %>%
  slice_head(n = 20) %>%
  ggplot(aes(unique_annotations, reorder(annotation_source, unique_annotations))) +
  geom_col(fill = uom_colours[["blue_light"]]) +
  geom_text(aes(label = unique_annotations), hjust = -0.1, size = 4.3) +
  labs(x = "Non-repeating annotations", y = "Annotation source") +
  theme_minimal(base_size = 18)

p_heatmap <- ggplot(platform_class, aes(subclass, platform, fill = unique_annotations)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = uom_colours[["heritage"]]) +
  labs(x = "NPClassifier class (top 30)", y = "Platform", fill = "Non-repeating\nannotations") +
  theme_minimal(base_size = 17) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1)
  )

figure_dir <- file.path(project_root, "figures", "figure-1")
ggsave(file.path(figure_dir, "figure1b-confidence-level-counts.png"), p_confidence, width = 5.6, height = 4, dpi = 300)
ggsave(file.path(figure_dir, "figure1-annotation-source-counts.png"), p_source, width = 6.4, height = 5.6, dpi = 300)
ggsave(file.path(figure_dir, "figure1-platform-class-heatmap.png"), p_heatmap, width = 12, height = 7.2, dpi = 300)
