source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)
  library(png)
  library(grid)
  library(scales)
  library(stringr)
  library(tidyr)
})

# Revised merged Figure 4/5 assembly.
# Panels:
# a, Platform-matched MZmine reprocessing comparison.
# b, Public faecal dataset remining yield with atlas-derived MSP libraries.
# c, Lipid-included atlas/HMDB/MiMeDB identity overlap.
# d, Atlas-only full-chemistry class breakdown beside the overlap plot.
# e, NIST public faecal molecular network mined with the atlas MSP library.

figure_dir <- file.path(project_root, "figures", "figure-4")
ensure_dirs(figure_dir)
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

theme_utility <- function(base_size = 16) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.text = element_text(colour = "#2B2B2B"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

png_panel <- function(path, bg = "white") {
  if (!file.exists(path)) stop("Missing panel image: ", path)
  img <- png::readPNG(path)
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  patchwork::wrap_elements(
    full = grid::grobTree(
      grid::rectGrob(gp = grid::gpar(col = NA, fill = bg)),
      grob
    )
  )
}

public_summary <- readr::read_csv(
  file.path(project_root, "outputs", "tables", "figure4-remining-utility-dataset-summary.csv"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  filter(utility_role == "Remine public datasets without RT matching") %>%
  left_join(
    readr::read_csv(
      file.path(project_root, "outputs", "tables", "figure4-public-remining-annotation-rates.csv"),
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      select(dataset.ID, total_resolved_features, annotation_rate_percent),
    by = "dataset.ID"
  ) %>%
  mutate(
    collection_label = str_replace_all(collection_label, "Public ", ""),
    collection_label = str_replace_all(collection_label, " dataset", ""),
    collection_label = factor(collection_label, levels = collection_label[order(non_repeating_annotations)]),
    rate_label = paste0(fmt_int(non_repeating_annotations), "\n(", round(annotation_rate_percent, 1), "%)")
  )

panel_a <- ggplot(public_summary, aes(collection_label, non_repeating_annotations, fill = inferred_polarity)) +
  geom_col(width = 0.68, colour = "#303030", linewidth = 0.25) +
  geom_text(
    aes(label = rate_label),
    vjust = -0.35,
    family = "Arial",
    fontface = "bold",
    size = 4.6,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(positive = "#617156", negative = "#BE7440", unknown = "#A7A7A7"),
    name = "Ion mode",
    na.value = "#A7A7A7"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.32)), labels = fmt_int) +
  labs(
    title = "Public faecal dataset remining",
    x = NULL,
    y = "Atlas-library annotations"
  ) +
  theme_utility(17) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.ticks.x = element_blank()
  )
panel_reprocessing <- png_panel(file.path(figure_dir, "figure4-mzmine-reprocessing-comparison.png"))

identity_audit <- readr::read_csv(
  file.path(project_root, "outputs", "qc", "figure5-identity-and-lipid-audit.csv"),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(
    collection = clean_text(collection),
    identity_key = clean_text(identity_key),
    npc_superclass = clean_text(npc_superclass)
  ) %>%
  filter(identity_key != "")

atlas_all <- identity_audit %>% filter(collection == "Human faecal atlas") %>% pull(identity_key) %>% unique()
hmdb_all <- identity_audit %>% filter(collection == "Healthy HMDB") %>% pull(identity_key) %>% unique()
mimedb_all <- identity_audit %>% filter(collection == "MiMeDB") %>% pull(identity_key) %>% unique()

panel_b <- ggVennDiagram(
  list(
    "Human faecal atlas" = atlas_all,
    "Healthy HMDB" = hmdb_all,
    "MiMeDB" = mimedb_all
  ),
  label_alpha = 0,
  category.names = c("Atlas", "HMDB", "MiMeDB"),
  set_color = c("#617156", "#BE7440", "#404B74"),
  set_size = 4.4,
  label = "count",
  label_size = 4.2,
  edge_size = 1
) +
  scale_fill_gradient(low = "#F4F4F4", high = "#9A4049") +
  labs(
    title = "Identity overlap including lipids"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    legend.position = "none",
    plot.margin = margin(22, 22, 18, 22)
  )

all_membership <- identity_audit %>%
  distinct(identity_key, collection) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = collection, values_from = present, values_fill = FALSE)

for (nm in c("Human faecal atlas", "Healthy HMDB", "MiMeDB")) {
  if (!nm %in% names(all_membership)) all_membership[[nm]] <- FALSE
}

atlas_unique_all <- all_membership %>%
  filter(.data[["Human faecal atlas"]], !.data[["Healthy HMDB"]], !.data[["MiMeDB"]]) %>%
  pull(identity_key)

atlas_unique_classes_all <- identity_audit %>%
  filter(collection == "Human faecal atlas", identity_key %in% atlas_unique_all) %>%
  mutate(npc_superclass = if_else(npc_superclass == "" | is.na(npc_superclass), "Unclassified", npc_superclass)) %>%
  distinct(identity_key, npc_superclass) %>%
  count(npc_superclass, sort = TRUE, name = "atlas_unique_identities") %>%
  slice_head(n = 18) %>%
  mutate(npc_superclass = factor(npc_superclass, levels = rev(npc_superclass)))

readr::write_csv(
  atlas_unique_classes_all,
  file.path(project_root, "outputs", "tables", "figure5-atlas-unique-all-classes.csv")
)

panel_c <- ggplot(atlas_unique_classes_all, aes(atlas_unique_identities, npc_superclass)) +
  geom_col(fill = "#617156", width = 0.72) +
  geom_text(
    aes(label = fmt_int(atlas_unique_identities)),
    hjust = -0.12,
    family = "Arial",
    fontface = "bold",
    size = 4
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "Atlas-only chemical classes",
    x = "Atlas-only identities",
    y = NULL
  ) +
  theme_utility(14) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = 10.5))

panel_d <- png_panel(file.path(figure_dir, "figure4-nist-atlas-msp-network.png"))

assembled <- (
  (panel_reprocessing | panel_a) + plot_layout(widths = c(4, 1.6))
) /
  (
    (panel_b | panel_c) + plot_layout(widths = c(1.25, 1.75))
  ) /
  panel_d +
  plot_layout(heights = c(1, 1.25, 3.2)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 36, face = "bold", family = "Arial"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

ggsave(
  file.path(figure_dir, "figure4-5-atlas-utility-assembled.png"),
  assembled,
  width = 14.4,
  height = 16.8,
  dpi = 300,
  bg = "white"
)

ggsave(file.path(figure_dir, "figure4-5a-public-remining-yield.png"), panel_a, width = 7.2, height = 3.12, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure4-5b-lipid-included-venn.png"), panel_b, width = 5.4, height = 4.44, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure4-5c-atlas-unique-all-classes.png"), panel_c, width = 5.4, height = 4.44, dpi = 300, bg = "white")

message("Revised merged Figure 4/5 assembly complete.")
