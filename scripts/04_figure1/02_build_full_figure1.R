source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(patchwork)
  library(scales)
  library(colorspace)
})

classified_fallback_path <- file.path(
  project_root, "data", "processed", "paper2_total_list_classified.current.csv"
)
classified_path <- file.path(project_root, "data", "processed", "paper2_total_list_classified.csv")
classified_candidates <- c(classified_path, classified_fallback_path)
classified_candidates <- classified_candidates[file.exists(classified_candidates)]
classified_path <- classified_candidates[which.max(file.info(classified_candidates)$mtime)]
x <- readr::read_csv(classified_path, show_col_types = FALSE) %>%
  mutate(
    confidence.level = as.integer(confidence.level),
    platform = paste(column.type, ion.mode, sep = " / "),
    annotation_source = if_else(is.na(annotation.type) | annotation.type == "", "Unspecified", annotation.type),
    superclass = if_else(is.na(npc_superclass) | npc_superclass == "", "Unclassified", npc_superclass),
    subclass = if_else(is.na(npc_class) | npc_class == "", "Unclassified", npc_class)
  )

best <- x %>%
  arrange(annotation_id, confidence.level, desc(confidence.score)) %>%
  group_by(annotation_id) %>%
  slice_head(n = 1) %>%
  ungroup()

blue <- "#404B74"
red <- "#9A4049"
light_blue <- uom_colours[["blue_light"]]
yellow <- uom_colours[["yellow"]]
orange <- uom_colours[["brown"]]
purple <- uom_colours[["maroon"]]
green <- uom_colours[["green_dark"]]
grey <- uom_colours[["grey_dark"]]
light_grey <- uom_colours[["grey_light"]]

theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#272727"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#272727"),
      axis.text = element_text(colour = "#272727"),
      axis.title = element_text(colour = "#272727"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size - 0.3, face = "bold"),
      legend.text = element_text(size = base_size - 0.8),
      legend.key.size = unit(0.3, "cm"),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

panel_label <- function(label) {
  annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.2, vjust = 1.1,
           fontface = "bold", family = "Arial", size = 5)
}

box <- function(xmin, xmax, ymin, ymax, label, fill, size = 3.4) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = "#606060", linewidth = 0.35),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, fontface = "bold", family = "Arial", size = size)
  )
}

p_a <- ggplot() +
  box(0.2, 2.1, 3.5, 4.25, "Donor and pooled\nfaecal extracts", light_grey, 4.1) +
  box(2.7, 4.6, 3.5, 4.25, "Crude LC-MS/MS", light_blue, 4.1) +
  box(5.2, 7.1, 4.25, 5.0, "Biphasic partition", uom_colours[["yellow"]], 4.1) +
  box(5.2, 7.1, 3.05, 3.8, "Deep fractionation", uom_colours[["yellow"]], 4.1) +
  box(7.7, 9.6, 3.5, 4.25, "Fraction LC-MS/MS", light_blue, 4.1) +
  box(3.0, 7.0, 1.55, 2.3, "MAPS annotation and class curation", uom_colours[["pink"]], 4.1) +
  box(3.0, 7.0, 0.25, 1.0, "Human faecal fraction atlas", uom_colours[["sage_light"]], 4.1) +
  annotate("segment", x = 2.1, xend = 2.7, y = 3.875, yend = 3.875,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 4.6, xend = 5.2, y = 3.875, yend = 4.625,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 6.15, xend = 6.15, y = 4.25, yend = 3.8,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 7.1, xend = 7.7, y = 3.425, yend = 3.875,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 3.65, xend = 4.15, y = 3.5, yend = 2.3,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 8.65, xend = 6.2, y = 3.5, yend = 2.3,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  annotate("segment", x = 5.0, xend = 5.0, y = 1.55, yend = 1.0,
           arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35) +
  coord_cartesian(xlim = c(0, 9.8), ylim = c(0, 5.15), clip = "off") +
  theme_void() +
  labs(title = "Construction of the fraction-resolved atlas") +
  theme(
    plot.title = element_text(size = 17, face = "bold", hjust = 0, family = "Arial"),
    plot.margin = margin(6, 6, 4, 6)
  )

confidence_stack <- best %>%
  count(confidence.level, annotation_source, name = "unique_annotations")
confidence_totals <- confidence_stack %>%
  group_by(confidence.level) %>%
  summarise(unique_annotations = sum(unique_annotations), .groups = "drop")
p_b <- ggplot(confidence_stack, aes(factor(confidence.level), unique_annotations, fill = annotation_source)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  geom_text(
    data = confidence_totals,
    aes(factor(confidence.level), unique_annotations, label = fmt_int(unique_annotations)),
    inherit.aes = FALSE, vjust = -0.25, size = 4.7, family = "Arial"
  ) +
  scale_fill_manual(values = rep(uom_discrete, length.out = n_distinct(confidence_stack$annotation_source)),
                    name = "Annotation evidence") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Annotation confidence", x = "Annotation level", y = "Non-repeating annotations") +
  theme_nature(14) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(12, 8, 8, 8)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(colour = NA)))

sources <- best %>% count(annotation_source, sort = TRUE, name = "unique_annotations")
p_c <- ggplot(sources, aes(unique_annotations, reorder(annotation_source, unique_annotations))) +
  geom_col(fill = grey) +
  geom_text(aes(label = fmt_int(unique_annotations)), hjust = -0.12, size = 5.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Annotation evidence", x = "Non-repeating annotations", y = NULL) +
  theme_classic(base_size = 17) + panel_label("c")

figure1c_smiles_qc <- x %>%
  filter(confidence.level <= 3, smiles != "") %>%
  distinct(smiles, superclass, subclass) %>%
  summarise(
    total_smiles = n_distinct(smiles),
    matched_smiles = n_distinct(smiles[superclass != "Unclassified" & subclass != "Unclassified"]),
    unmatched_smiles = total_smiles - matched_smiles,
    named_class_smiles = n_distinct(smiles[
      superclass != "Unclassified" &
        subclass != "Unclassified" &
        !str_detect(str_to_lower(superclass), "^other") &
        !str_detect(str_to_lower(subclass), "^other")
    ])
  )

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
  count(superclass, subclass, name = "n") %>%
  group_by(superclass) %>%
  mutate(superclass_total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(superclass_total), superclass, desc(n), subclass)

top_super <- composition %>%
  distinct(superclass, superclass_total) %>%
  slice_max(superclass_total, n = 10, with_ties = FALSE) %>%
  pull(superclass)

composition <- composition %>%
  filter(superclass %in% top_super) %>%
  group_by(superclass) %>%
  mutate(subclass_rank = min_rank(desc(n))) %>%
  ungroup() %>%
  filter(subclass_rank <= 8) %>%
  group_by(superclass) %>%
  mutate(superclass_total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(superclass_total), superclass, desc(n), subclass) %>%
  mutate(superclass_plot = superclass, subclass_plot = subclass)

figure1c_smiles_qc <- figure1c_smiles_qc %>%
  mutate(plotted_named_smiles = sum(composition$n))

super_tbl <- composition %>%
  distinct(superclass_plot, superclass_total) %>%
  mutate(
    end = cumsum(superclass_total),
    start = lag(end, default = 0)
  )
outer_tbl <- composition %>%
  left_join(super_tbl %>% select(superclass_plot, super_start = start), by = "superclass_plot") %>%
  arrange(desc(superclass_total), superclass_plot, desc(n), subclass_plot) %>%
  group_by(superclass_plot) %>%
  mutate(end = super_start + cumsum(n), start = lag(end, default = first(super_start))) %>%
  ungroup()

super_tbl <- super_tbl %>% mutate(legend_label = paste0("Superclass: ", superclass_plot))
outer_tbl <- outer_tbl %>% mutate(legend_label = paste0("Subclass: ", subclass_plot))
ring_total <- sum(super_tbl$superclass_total)
super_tbl <- super_tbl %>%
  mutate(
    middle = (start + end) / 2 / ring_total * 2 * pi,
    label_x = sin(middle) * 0.515,
    label_y = cos(middle) * 0.515,
    ring_label = str_wrap(superclass_plot, width = 12)
  )
outer_tbl <- outer_tbl %>%
  mutate(
    middle = (start + end) / 2 / ring_total * 2 * pi,
    label_x = sin(middle) * 0.85,
    label_y = cos(middle) * 0.85,
    ring_label = if_else(n >= 50, str_wrap(subclass_plot, width = 10), "")
  )
pie_family_cols <- c(
  "#404B74", "#9A4049", "#406D80", "#965A78", "#BE7440",
  "#8087A2", "#BB8086", "#809DAA", "#B991A5", "#D4A280"
)
superclass_cols <- setNames(
  rep(pie_family_cols, length.out = nrow(super_tbl)),
  super_tbl$superclass_plot
)
outer_tbl <- outer_tbl %>%
  group_by(superclass_plot) %>%
  mutate(
    subclass_shade = seq(0.2, 0.72, length.out = n()),
    subclass_colour = colorspace::lighten(
      superclass_cols[superclass_plot],
      amount = subclass_shade,
      space = "HCL"
    )
  ) %>%
  ungroup()
classification_cols <- c(
  setNames(superclass_cols[super_tbl$superclass_plot], super_tbl$legend_label),
  setNames(outer_tbl$subclass_colour, outer_tbl$legend_label)
)

p_d <- ggplot() +
  geom_arc_bar(
    data = super_tbl,
    aes(
      x0 = 0, y0 = 0, r0 = 0.35, r = 0.68,
      start = start / sum(superclass_total) * 2 * pi,
      end = end / sum(superclass_total) * 2 * pi,
      fill = legend_label
    ),
    colour = "white", linewidth = 0.35
  ) +
  geom_arc_bar(
    data = outer_tbl,
    aes(
      x0 = 0, y0 = 0, r0 = 0.70, r = 1,
      start = start / sum(n) * 2 * pi,
      end = end / sum(n) * 2 * pi,
      fill = legend_label
    ),
    colour = "white", linewidth = 0.2
  ) +
  geom_text(
    data = super_tbl,
    aes(x = label_x, y = label_y, label = ring_label),
    inherit.aes = FALSE, size = 3.8, fontface = "bold", family = "Arial", lineheight = 0.85
  ) +
  geom_text(
    data = outer_tbl,
    aes(x = label_x, y = label_y, label = ring_label),
    inherit.aes = FALSE, size = 2.7, family = "Arial", lineheight = 0.8
  ) +
  scale_fill_manual(
    values = classification_cols,
    breaks = outer_tbl$legend_label,
    labels = outer_tbl$subclass_plot
  ) +
  coord_fixed() +
  labs(
    title = "Atlas-wide chemical composition",
    fill = "Chemical classification"
  ) +
  theme_void(base_size = 11.5, base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(0.28, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.margin = margin(t = 1, r = 0, b = 0, l = 0),
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_legend(
    title = "Subclass",
    ncol = 2,
    byrow = TRUE,
    override.aes = list(colour = NA)
  ))

platform_long <- x %>%
  filter(confidence.level <= 3, subclass != "Unclassified") %>%
  distinct(annotation_id, platform, subclass)
top_classes <- platform_long %>%
  distinct(annotation_id, subclass) %>%
  count(subclass, sort = TRUE) %>%
  slice_head(n = 30) %>%
  pull(subclass)
heat <- platform_long %>%
  filter(subclass %in% top_classes) %>%
  count(platform, subclass, name = "unique_annotations")
platform_sum <- x %>%
  filter(confidence.level <= 3) %>%
  distinct(annotation_id, platform) %>%
  group_by(annotation_id) %>%
  mutate(platform_n = n()) %>%
  ungroup() %>%
  group_by(platform) %>%
  summarise(total = n_distinct(annotation_id), unique = n_distinct(annotation_id[platform_n == 1]), .groups = "drop") %>%
  pivot_longer(c(total, unique), names_to = "metric", values_to = "n")

p_e_heat <- ggplot(heat, aes(subclass, platform, fill = unique_annotations)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(
    colours = c("#D9E2EA", "#9DB9C7", "#D5B4A3", "#B66F70", "#7F2635"),
    values = c(0, 0.18, 0.43, 0.72, 1),
    trans = "sqrt"
  ) +
  labs(title = "Platform coverage by compound class", x = NULL, y = NULL, fill = "Non-repeating\nannotations") +
  theme_nature(13) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 9.5),
    legend.position = "right"
  )
p_e_side <- ggplot(platform_sum, aes(n, platform, fill = metric)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = fmt_int(n)), position = position_dodge(width = 0.7),
            hjust = -0.05, size = 3.1, family = "Arial") +
  scale_fill_manual(
    values = c(
      total = blue,
      unique = red
    ),
    labels = c(total = "non-repeating", unique = "unique to\nplatform")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = "Annotation count", y = NULL, fill = NULL) +
  theme_nature(13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )
p_e <- p_e_heat + p_e_side + plot_layout(widths = c(3.7, 2.3))
p_e <- wrap_elements(full = p_e)

source(file.path("scripts", "04_figure1", "06_build_circular_atlas_lineage.R"))
p_lineage_assembled <- wrap_elements(full = p_lineage)

p_d_standalone <- p_d +
  guides(fill = guide_legend(
    title = "Subclass",
    ncol = 4,
    byrow = TRUE,
    override.aes = list(colour = NA)
  )) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6.3),
    legend.key.size = unit(0.22, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.margin = margin(t = 1, r = 0, b = 0, l = 0)
  )
p_d_assembled <- p_d +
  guides(fill = guide_legend(
    title = "Subclass",
    nrow = 30,
    byrow = FALSE,
    override.aes = list(colour = NA)
  )) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 6.8),
    legend.key.size = unit(0.26, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.margin = margin(t = 1, r = 0, b = 0, l = 0)
  )
p_d_assembled <- wrap_elements(full = p_d_assembled)
row1 <- p_a
row2 <- (p_lineage_assembled | p_b) +
  plot_layout(widths = c(3, 1))
row3 <- p_d_assembled
row4 <- p_e
full <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(0.55, 1.05, 1.95, 1.15)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 36, face = "bold", family = "Arial")
    )
  ) &
  theme(text = element_text(family = "Arial"))

figure_dir <- file.path(project_root, "figures", "figure-1")
write_csv_stable(confidence_stack, file.path(project_root, "outputs", "tables", "figure1b-confidence-level-source-stacks.csv"))
write_csv_stable(figure1c_smiles_qc, file.path(project_root, "outputs", "tables", "figure1c-smiles-classification-coverage.csv"))
write_csv_stable(heat, file.path(project_root, "outputs", "tables", "figure1d-platform-class-counts.csv"))
write_csv_stable(platform_sum, file.path(project_root, "outputs", "tables", "figure1d-platform-total-and-unique.csv"))
ggsave(file.path(figure_dir, "figure1a-atlas-workflow-flowchart.png"), p_a, width = 10.4, height = 4, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure1c-nested-chemical-composition.png"), p_d_standalone, width = 9.6, height = 9.6, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure1d-platform-contribution-with-summary.png"), p_e, width = 9.6, height = 3.84, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure1-assembled.png"), full, width = 12.8, height = 24.8, dpi = 300, bg = "white")
