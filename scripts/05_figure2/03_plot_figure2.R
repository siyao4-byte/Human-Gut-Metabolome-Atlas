source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages(library(ggplot2))

gain <- readr::read_csv(
  file.path(project_root, "outputs", "tables", "figure2-annotation-gain.csv"),
  show_col_types = FALSE
)
enrichment <- readr::read_csv(
  file.path(project_root, "outputs", "tables", "figure2-class-enrichment.csv"),
  show_col_types = FALSE
)
unique_class <- readr::read_csv(
  file.path(project_root, "outputs", "tables", "figure2-unique-class-identity.csv"),
  show_col_types = FALSE
)

level_order <- c(
  "methanol_crude", "water_crude", "dcm_crude",
  "polar_fraction", "nonpolar_fraction"
)
gain <- gain %>%
  mutate(figure2_processing_level = factor(figure2_processing_level, levels = level_order))

p_gain <- ggplot(
  gain,
  aes(figure2_processing_level, unique_annotations, fill = factor(confidence.level))
) +
  geom_col() +
  facet_wrap(~ combined_pool_id, ncol = 2) +
  scale_fill_manual(values = uom_confidence, name = "Confidence level") +
  labs(x = "Processing level", y = "Unique annotations") +
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_heatmap <- ggplot(
  enrichment,
  aes(npc_superclass, figure2_processing_level, fill = unique_annotations)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = uom_colours[["heritage"]]) +
  labs(x = "NPClassifier superclass", y = "Processing level", fill = "Unique\nannotations") +
  theme_minimal(base_size = 22.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1)
  )

p_bubble <- ggplot(
  unique_class,
  aes(npc_superclass, figure2_processing_level, size = unique_annotations, color = unique_annotations)
) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = uom_colours[["blue_light"]], high = uom_colours[["heritage"]]) +
  labs(
    x = "NPClassifier superclass", y = "Processing level",
    size = "Unique annotations", color = "Unique annotations"
  ) +
  theme_minimal(base_size = 22.5) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))

figure_dir <- file.path(project_root, "figures", "figure-2")
ggsave(file.path(figure_dir, "figure2-annotation-gain-four-pools.png"), p_gain, width = 12, height = 8, dpi = 300)
ggsave(file.path(figure_dir, "figure2-class-enrichment-heatmap.png"), p_heatmap, width = 15.2, height = 6.4, dpi = 300)
ggsave(file.path(figure_dir, "figure2-unique-class-bubble.png"), p_bubble, width = 15.2, height = 7.2, dpi = 300)
