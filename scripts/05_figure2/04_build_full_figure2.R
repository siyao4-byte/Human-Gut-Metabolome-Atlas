source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

theme_nature2 <- function(base_size = 21) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "#272727"),
      axis.ticks = element_line(linewidth = 0.4, colour = "#272727"),
      axis.text = element_text(colour = "#272727"),
      plot.title = element_text(size = base_size + 1.5, face = "bold"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank()
    )
}

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
    superclass = clean_text(if ("np_superclass" %in% names(npc)) np_superclass else superclass)
  ) %>%
  filter(smiles != "") %>%
  distinct(smiles, .keep_all = TRUE)
x <- x %>%
  mutate(smiles = clean_text(smiles)) %>%
  left_join(npc_small, by = "smiles") %>%
  mutate(superclass = if_else(is.na(superclass) | superclass == "", "Unclassified", superclass))

blue <- "#404B74"
light_blue <- uom_colours[["blue_light"]]
red <- "#9A4049"
orange <- uom_colours[["brown"]]
teal <- uom_colours[["blue_dark"]]
purple <- uom_colours[["maroon"]]
grey <- uom_colours[["grey_dark"]]
confidence_cols <- c("1" = blue, "2" = "#A0A5B7", "3" = red)
processing_cols <- c(
  methanol_crude = blue,
  dcm_crude = red,
  water_crude = teal,
  nonpolar_fraction = orange,
  polar_fraction = purple
)
processing_cols_light <- c(
  methanol_crude = "#929AB7",
  dcm_crude = "#D3A0A5",
  water_crude = "#A7C8D1",
  nonpolar_fraction = "#E2BD91",
  polar_fraction = "#CDB2C0"
)

layer_order <- c("methanol_crude", "dcm_crude", "water_crude", "nonpolar_fraction", "polar_fraction")
layer_labels <- c(
  methanol_crude = "Methanol\ncrude",
  dcm_crude = "DCM\nphase",
  water_crude = "Aqueous\nphase",
  nonpolar_fraction = "DCM\nfractions",
  polar_fraction = "Aqueous\nfractions"
)
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

layer_counts <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  count(figure2_processing_level, name = "n") %>%
  complete(figure2_processing_level = layer_order, fill = list(n = 0)) %>%
  mutate(label = layer_labels[figure2_processing_level])
atlas_n <- n_distinct(x$annotation_id)
count_for <- function(level) layer_counts$n[match(level, layer_counts$figure2_processing_level)]

box <- function(xmin, xmax, ymin, ymax, label, fill, colour) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = colour, linewidth = 0.9),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, fontface = "bold", family = "Arial", size = 8.7)
  )
}

p_a <- ggplot() +
  box(0.1, 1.7, 1.8, 2.8, paste0("Methanol crude\n", fmt_int(count_for("methanol_crude"))), uom_colours[["blue_light"]], blue) +
  box(2.4, 4.0, 2.8, 3.8, paste0("DCM phase\n", fmt_int(count_for("dcm_crude"))), "#DDBFC2", red) +
  box(2.4, 4.0, 0.8, 1.8, paste0("Aqueous phase\n", fmt_int(count_for("water_crude"))), "#BFCED5", teal) +
  box(4.8, 6.5, 2.8, 3.8, paste0("DCM fractions\n", fmt_int(count_for("nonpolar_fraction"))), uom_colours[["yellow"]], orange) +
  box(4.8, 6.5, 0.8, 1.8, paste0("Aqueous fractions\n", fmt_int(count_for("polar_fraction"))), "#DCC8D2", purple) +
  box(7.3, 8.8, 1.8, 2.8, paste0("Atlas\n", fmt_int(atlas_n)), uom_colours[["grey_light"]], grey) +
  annotate("segment", x = 1.7, xend = 2.4, y = 2.3, yend = 3.3, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 1.7, xend = 2.4, y = 2.3, yend = 1.3, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 4.0, xend = 4.8, y = 3.3, yend = 3.3, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 4.0, xend = 4.8, y = 1.3, yend = 1.3, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 6.5, xend = 7.3, y = 3.3, yend = 2.45, arrow = arrow(length = unit(0.12, "inches"))) +
  annotate("segment", x = 6.5, xend = 7.3, y = 1.3, yend = 2.15, arrow = arrow(length = unit(0.12, "inches"))) +
  coord_cartesian(xlim = c(0, 9), ylim = c(0.4, 4.3), clip = "off") +
  labs(title = "Hierarchical gain from crude extract to phases and fractions") +
  theme_void() +
  theme(plot.title = element_text(size = 28.5, face = "bold", hjust = 0.5, family = "Arial"))

class_counts <- x %>%
  filter(superclass != "Unclassified") %>%
  distinct(annotation_id, figure2_processing_level, superclass) %>%
  count(figure2_processing_level, superclass, name = "n")
top_class <- class_counts %>%
  group_by(superclass) %>% summarise(total = sum(n), .groups = "drop") %>%
  slice_max(total, n = 30, with_ties = FALSE) %>% pull(superclass)
enrich <- class_counts %>%
  filter(superclass %in% top_class) %>%
  complete(figure2_processing_level = layer_order, superclass = top_class, fill = list(n = 0)) %>%
  group_by(figure2_processing_level) %>%
  mutate(level_total = sum(n)) %>%
  ungroup() %>%
  group_by(superclass) %>%
  mutate(class_total = sum(n)) %>%
  ungroup() %>%
  mutate(
    grand_total = sum(n),
    expected_count = level_total * class_total / grand_total,
    enrichment_percent = if_else(expected_count > 0, 100 * (n / expected_count - 1), NA_real_)
  )

confidence_by_level <- x %>%
  distinct(annotation_id, figure2_processing_level, confidence.level) %>%
  count(figure2_processing_level, confidence.level, name = "unique_annotations")
p_b_confidence <- ggplot(
  confidence_by_level,
  aes(figure2_processing_level, unique_annotations, fill = factor(confidence.level))
) +
  geom_col(colour = "white", linewidth = 0.25) +
  scale_x_discrete(labels = layer_labels) +
  scale_fill_manual(
    values = confidence_cols,
    name = "Confidence level",
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(title = "Annotation confidence", x = NULL, y = "Non-repeating annotations") +
  theme_nature2(21) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "bottom",
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

processing_unique_counts <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  group_by(annotation_id) %>%
  mutate(processing_level_count = n_distinct(figure2_processing_level)) %>%
  ungroup() %>%
  group_by(figure2_processing_level) %>%
  summarise(
    `Total detected` = n_distinct(annotation_id),
    `Unique to processing level` = n_distinct(annotation_id[processing_level_count == 1]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(`Total detected`, `Unique to processing level`),
    names_to = "coverage_type",
    values_to = "annotations"
  ) %>%
  mutate(
    coverage_type = factor(
      coverage_type,
      levels = c("Total detected", "Unique to processing level")
    )
  )

p_b_coverage <- ggplot(
  processing_unique_counts,
  aes(figure2_processing_level, annotations, fill = figure2_processing_level)
) +
  geom_col(width = 0.72, colour = "#303030", linewidth = 0.35) +
  geom_text(
    aes(label = fmt_int(annotations)),
    vjust = -0.35, family = "Arial", fontface = "bold", size = 5.4
  ) +
  facet_wrap(~coverage_type, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = layer_labels) +
  scale_fill_manual(values = processing_cols_light, guide = "none") +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Total and processing-level-specific coverage", x = NULL, y = "Annotations") +
  theme_nature2(21) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

p_b <- p_b_coverage +
  plot_annotation(
    title = "Annotation coverage across processing levels",
    theme = theme(
      plot.title = element_text(size = 28.5, face = "bold", hjust = 0.5, family = "Arial")
    )
  )

# UpSet representation of total, unique and overlapping processing-level coverage.
upset_membership <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = figure2_processing_level,
    values_from = present,
    values_fill = FALSE
  )

upset_intersections <- upset_membership %>%
  rowwise() %>%
  mutate(
    signature = paste(layer_order[c_across(all_of(layer_order))], collapse = " + "),
    set_count = sum(c_across(all_of(layer_order)))
  ) %>%
  ungroup() %>%
  count(signature, set_count, name = "intersection_size") %>%
  arrange(desc(intersection_size), set_count, signature) %>%
  mutate(
    intersection_id = row_number(),
    intersection_axis = factor(intersection_id, levels = intersection_id)
  )

upset_matrix <- upset_intersections %>%
  select(intersection_id, intersection_axis, signature) %>%
  crossing(figure2_processing_level = layer_order) %>%
  mutate(
    present = purrr::map2_lgl(
      signature, figure2_processing_level,
      ~ .y %in% str_split(.x, " \\+ ")[[1]]
    ),
    figure2_processing_level = factor(figure2_processing_level, levels = rev(layer_order))
  )

upset_set_sizes <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  count(figure2_processing_level, name = "set_size") %>%
  mutate(figure2_processing_level = factor(figure2_processing_level, levels = rev(layer_order)))

p_b_intersections <- ggplot(
  upset_intersections,
  aes(intersection_axis, intersection_size)
) +
  geom_col(width = 0.72, fill = "#777777", colour = "#303030", linewidth = 0.25) +
  geom_text(
    aes(label = fmt_int(intersection_size)),
    vjust = -0.35, size = 3.7, family = "Arial", fontface = "bold"
  ) +
  scale_x_discrete(drop = FALSE, breaks = NULL, expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.17))) +
  labs(x = NULL, y = "Intersection size") +
  theme_nature2(17) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

p_b_matrix <- ggplot(
  upset_matrix,
  aes(intersection_axis, figure2_processing_level)
) +
  geom_line(
    data = upset_matrix %>% filter(present) %>% group_by(intersection_id),
    aes(group = intersection_id),
    colour = "#9A9A9A", linewidth = 0.65
  ) +
  geom_point(colour = "#D8D8D8", size = 3.1) +
  geom_point(
    data = upset_matrix %>% filter(present),
    aes(colour = figure2_processing_level),
    size = 3.5
  ) +
  scale_colour_manual(values = processing_cols_light, guide = "none") +
  scale_x_discrete(drop = FALSE, breaks = NULL, expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(labels = layer_labels) +
  labs(x = NULL, y = NULL) +
  theme_nature2(17) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank()
  )

p_b_set_sizes <- ggplot(
  upset_set_sizes,
  aes(set_size, figure2_processing_level, fill = figure2_processing_level)
) +
  geom_col(width = 0.68, colour = "#303030", linewidth = 0.25) +
  geom_text(
    aes(label = fmt_int(set_size)),
    hjust = -0.18, size = 4.2, family = "Arial", fontface = "bold"
  ) +
  scale_fill_manual(values = processing_cols_light, guide = "none") +
  scale_x_continuous(
    trans = "reverse", labels = fmt_int,
    expand = expansion(mult = c(0.22, 0.02))
  ) +
  scale_y_discrete(labels = layer_labels) +
  labs(x = "Total annotations", y = NULL) +
  theme_nature2(17) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 14)
  )

p_broad_membership <- x %>%
  distinct(annotation_id, figure2_processing_level) %>%
  mutate(
    collection = case_when(
      figure2_processing_level == "methanol_crude" ~ "Methanol crude",
      figure2_processing_level %in% c("dcm_crude", "water_crude") ~ "Biphasic phases",
      figure2_processing_level %in% c("nonpolar_fraction", "polar_fraction") ~ "Fractions"
    )
  ) %>%
  filter(!is.na(collection)) %>%
  distinct(annotation_id, collection) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = collection, values_from = present, values_fill = FALSE)

broad_set_sizes <- c(
  "Methanol crude" = sum(p_broad_membership$`Methanol crude`),
  "Biphasic phases" = sum(p_broad_membership$`Biphasic phases`),
  "Fractions" = sum(p_broad_membership$Fractions)
)
broad_intersections <- p_broad_membership %>%
  count(
    methanol = `Methanol crude`,
    biphasic = `Biphasic phases`,
    fractions = Fractions,
    name = "n"
  )
intersection_n <- function(m, b, f) {
  broad_intersections %>%
    filter(methanol == m, biphasic == b, fractions == f) %>%
    pull(n) %>%
    sum()
}
venn_circles <- tibble(
  collection = names(broad_set_sizes),
  x0 = c(-0.62, 0.62, 0),
  y0 = c(0.30, 0.30, -0.42),
  r = sqrt(unname(broad_set_sizes) / max(broad_set_sizes)) * 1.22
)
venn_labels <- tibble(
  x = c(-1.08, 1.08, 0, 0, -0.57, 0.57, 0),
  y = c(0.50, 0.50, -1.20, 0.72, -0.37, -0.37, 0.03),
  n = c(
    intersection_n(TRUE, FALSE, FALSE),
    intersection_n(FALSE, TRUE, FALSE),
    intersection_n(FALSE, FALSE, TRUE),
    intersection_n(TRUE, TRUE, FALSE),
    intersection_n(TRUE, FALSE, TRUE),
    intersection_n(FALSE, TRUE, TRUE),
    intersection_n(TRUE, TRUE, TRUE)
  )
)
p_b_venn <- ggplot() +
  ggforce::geom_circle(
    data = venn_circles,
    aes(x0 = x0, y0 = y0, r = r, fill = collection, colour = collection),
    alpha = 0.23, linewidth = 0.8
  ) +
  geom_text(
    data = venn_labels, aes(x, y, label = fmt_int(n)),
    family = "Arial", fontface = "bold", size = 4.2
  ) +
  scale_fill_manual(
    values = c("Methanol crude" = processing_cols_light[["methanol_crude"]],
               "Biphasic phases" = processing_cols_light[["dcm_crude"]],
               "Fractions" = processing_cols_light[["nonpolar_fraction"]])
  ) +
  scale_colour_manual(
    values = c("Methanol crude" = processing_cols[["methanol_crude"]],
               "Biphasic phases" = processing_cols[["dcm_crude"]],
               "Fractions" = processing_cols[["nonpolar_fraction"]])
  ) +
  coord_fixed(xlim = c(-2.0, 2.0), ylim = c(-1.95, 1.75), clip = "off") +
  labs(title = "Broad processing collections", fill = NULL, colour = NULL) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0)
  )

p_b_design <- "
AB
CD
"
p_b <- p_b_venn + p_b_intersections + p_b_set_sizes + p_b_matrix +
  plot_layout(
    design = p_b_design,
    widths = c(0.36, 1),
    heights = c(0.68, 0.32)
  ) +
  plot_annotation(
    title = "Unique and overlapping annotations across processing levels",
    theme = theme(
      plot.title = element_text(size = 28.5, face = "bold", hjust = 0.5, family = "Arial")
    )
  )
p_b <- wrap_elements(full = p_b)

heat_max <- max(enrich$n)
class_order <- enrich %>%
  group_by(superclass) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(total) %>%
  pull(superclass)
enrich <- enrich %>%
  mutate(superclass = factor(superclass, levels = class_order))

p_c <- ggplot(enrich, aes(figure2_processing_level, superclass, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours = c("#EEF1F6", "#B4C0D2", "#D3B5AF", "#B66F70", red),
    values = c(0, 0.18, 0.42, 0.70, 1),
    trans = "sqrt",
    limits = c(0, heat_max)
  ) +
  scale_x_discrete(labels = layer_labels) +
  labs(title = "Class coverage", x = NULL, y = NULL, fill = "Non-repeating\nannotations") +
  theme_nature2(21) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))

presence_by_level <- x %>%
  distinct(annotation_id, figure2_processing_level, superclass) %>%
  group_by(annotation_id) %>%
  mutate(level_n = n_distinct(figure2_processing_level)) %>%
  ungroup() %>%
  filter(level_n == 1, superclass != "Unclassified") %>%
  count(figure2_processing_level, superclass, name = "n")
unique_bubble <- presence_by_level %>%
  filter(superclass %in% class_order) %>%
  mutate(superclass = factor(superclass, levels = class_order))
p_d <- ggplot(unique_bubble, aes(n, superclass, size = n, colour = figure2_processing_level)) +
  geom_point(alpha = 0.9, stroke = 1.1) +
  scale_colour_manual(
    values = processing_cols_light,
    labels = layer_labels,
    name = "Processing level",
    guide = guide_legend(
      override.aes = list(size = 8, alpha = 1),
      keyheight = unit(1.15, "cm")
    )
  ) +
  scale_size_continuous(
    range = c(3.5, 15),
    trans = "sqrt",
    breaks = c(25, 75, 150),
    guide = guide_legend(
      override.aes = list(size = c(8, 14, 20)),
      keyheight = unit(1.5, "cm")
    )
  ) +
  scale_x_sqrt(
    breaks = c(0, 5, 10, 25, 50, 100, 200),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(title = "Unique-to-level annotations", x = "Annotation count", y = NULL, size = "Unique annotations") +
  theme_nature2(21) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.7, "cm"),
    plot.margin = margin(10, 18, 10, 10)
  )

bh_gain <- readr::read_csv(
  file.path(project_root, "outputs", "tables", "figure2e-bh01-31-collection-counts.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    stage = factor(
      display_level,
      levels = display_level
    ),
    stage_colour = case_when(
      collection == "BH01-31 crude reference" ~ "Crude reference",
      collection == "Complete BH01-31 collection" ~ "Primary fractionation",
      collection == "BH01-31 fractions 18 and 19" ~ "Selected fractions",
      TRUE ~ "C18 subfractionation"
    )
  )

p_e_overview_data <- bh_gain %>%
  filter(collection %in% c("BH01-31 crude reference", "Complete BH01-31 collection"))
p_e_zoom_data <- bh_gain %>%
  filter(collection %in% c("BH01-31 fractions 18 and 19", "Complete BH01-40 collection"))

stage_cols <- c(
  "Crude reference" = "#858585",
  "Primary fractionation" = orange,
  "Selected fractions" = "#C8A96B",
  "C18 subfractionation" = "#617156"
)
make_bh_bar <- function(d, title) {
  ggplot(d, aes(stage, collection_annotations, fill = stage_colour)) +
  geom_col(width = 0.68, colour = "#303030", linewidth = 0.4) +
  geom_text(
    aes(label = fmt_int(collection_annotations)),
    vjust = -0.4, family = "Arial", fontface = "bold", size = 6.2
  ) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = title,
    x = NULL, y = "Non-repeating annotations"
  ) +
  theme_nature2(19) +
  theme(
    axis.text.x = element_text(size = 15, face = "bold")
  )
}

p_e_overview <- make_bh_bar(
  p_e_overview_data,
  "Primary fractionation relative to crude baseline"
)
p_e_zoom <- make_bh_bar(
  p_e_zoom_data,
  "Selected primary-fraction branch"
)
p_e_composite <- p_e_overview + p_e_zoom +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "A traceable non-polar branch quantifies primary and deeper fractionation",
    caption = "Positive- and negative-ion datasets are merged within each collection before non-repeating annotations are counted.",
    theme = theme(
      plot.title = element_text(size = 25, face = "bold", family = "Arial", margin = margin(l = 28)),
      plot.caption = element_text(size = 12, colour = "#606060", family = "Arial", hjust = 0, margin = margin(l = 28))
    )
  )
p_e <- wrap_elements(full = p_e_composite)

write_csv_stable(layer_counts, file.path(project_root, "outputs", "tables", "figure2-hierarchy-counts.csv"))
write_csv_stable(enrich, file.path(project_root, "outputs", "tables", "figure2-relative-class-enrichment.csv"))
write_csv_stable(confidence_by_level, file.path(project_root, "outputs", "tables", "figure2-confidence-level-by-processing-level.csv"))
write_csv_stable(processing_unique_counts, file.path(project_root, "outputs", "tables", "figure2-total-and-unique-by-processing-level.csv"))
write_csv_stable(unique_bubble, file.path(project_root, "outputs", "tables", "figure2-unique-class-bubble.csv"))
write_csv_stable(upset_intersections, file.path(project_root, "outputs", "tables", "figure2-processing-level-intersections.csv"))
write_csv_stable(upset_set_sizes, file.path(project_root, "outputs", "tables", "figure2-processing-level-set-sizes.csv"))
for (level in layer_order) {
  unique_level <- x %>%
    distinct(annotation_id, figure2_processing_level, compound.name, smiles) %>%
    group_by(annotation_id) %>% mutate(level_n = n_distinct(figure2_processing_level)) %>% ungroup() %>%
    filter(level_n == 1, figure2_processing_level == level)
  write_csv_stable(unique_level, file.path(project_root, "outputs", "tables", paste0("figure2-", level, "-unique-compounds.csv")))
}

full <- p_a / p_b / p_d / p_e +
  plot_layout(heights = c(0.55, 1.05, 1.15, 0.9)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 36, face = "bold", family = "Arial")
    )
  )
figure_dir <- file.path(project_root, "figures", "figure-2")
ggsave(file.path(figure_dir, "figure2a-processing-workflow-flowchart.png"), p_a, width = 21.6, height = 5.2, dpi = 300)
ggsave(file.path(figure_dir, "figure2b-annotation-overlap-summary.png"), p_b, width = 21.6, height = 8.4, dpi = 300)
ggsave(file.path(figure_dir, "figure2c-unique-class-bubble.png"), p_d, width = 12.8, height = 12.8, dpi = 300)
ggsave(file.path(figure_dir, "figure2d-lineage-subfractionation-comparison.png"), p_e, width = 21.6, height = 8.8, dpi = 300)
ggsave(file.path(figure_dir, "figure2-assembled.png"), full, width = 21.6, height = 28, dpi = 300)
