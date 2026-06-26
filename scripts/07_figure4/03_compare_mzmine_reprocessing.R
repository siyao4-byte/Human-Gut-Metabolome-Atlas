source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

dataset_ids <- c("HGMD_0314", "HGMD_0315", "HGMD_0318", "HGMD_0319")
reprocessed_root <- file.path(project_root, "outputs", "mzmine reprocessing")
maps_root <- path_value("annotation_root")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "figures", "figure-4")
ensure_dirs(table_dir, figure_dir)
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

dataset_metadata <- tibble::tribble(
  ~dataset.ID, ~platform,
  "HGMD_0314", "Phe-Hex positive-ion mode",
  "HGMD_0315", "Phe-Hex negative-ion mode",
  "HGMD_0318", "HILIC positive-ion mode",
  "HGMD_0319", "HILIC negative-ion mode"
) %>%
  mutate(dataset_label = paste0(dataset.ID, "\n", platform))

find_reprocessed_csv <- function(dataset_id) {
  files <- list.files(
    file.path(reprocessed_root, dataset_id),
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(files) != 1) {
    stop("Expected one reprocessed annotation CSV for ", dataset_id, "; found ", length(files))
  }
  files[[1]]
}

make_annotation_id <- function(compound_name, smiles) {
  compound_name <- na_if(clean_text(compound_name), "")
  smiles <- na_if(clean_text(smiles), "")
  coalesce(compound_name, smiles)
}

annotation_type_rank <- function(x) {
  case_when(
    str_detect(x, regex("authentic|mzmine|in-house spectral library", ignore_case = TRUE)) ~ 1L,
    str_detect(x, regex("gnps", ignore_case = TRUE)) ~ 2L,
    str_detect(x, regex("ms2query", ignore_case = TRUE)) ~ 3L,
    TRUE ~ 99L
  )
}

read_reprocessed <- function(dataset_id) {
  path <- find_reprocessed_csv(dataset_id)
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  x %>%
    transmute(
      dataset.ID = dataset_id,
      source = "Atlas-library MZmine reprocessing",
      source_row = row_number(),
      annotation_id = make_annotation_id(compound_name, smiles),
      compound.name = na_if(clean_text(compound_name), ""),
      smiles = na_if(clean_text(smiles), ""),
      score = suppressWarnings(as.numeric(score)),
      annotation_source = "Atlas-derived in-house spectral library",
      confidence.level = 1,
      annotation_rank = 1L,
      id.prob = score,
      confidence.score = score
    ) %>%
    filter(!is.na(annotation_id))
}

read_maps <- function(dataset_id) {
  path <- file.path(maps_root, paste0(dataset_id, ".csv"))
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  x %>%
    transmute(
      dataset.ID = dataset_id,
      source = "MAPS",
      source_row = row_number(),
      annotation_id = make_annotation_id(compound.name, smiles),
      compound.name = na_if(clean_text(compound.name), ""),
      smiles = na_if(clean_text(smiles), ""),
      confidence.level = suppressWarnings(as.numeric(confidence.level)),
      annotation.type = clean_text(annotation.type),
      annotation_source = clean_text(annotation.type),
      annotation_rank = annotation_type_rank(annotation.type),
      id.prob = suppressWarnings(as.numeric(id.prob)),
      confidence.score = suppressWarnings(as.numeric(confidence.score))
    ) %>%
    filter(
      !is.na(annotation_id),
      confidence.level %in% 1:3,
      !str_detect(annotation.type, regex("^MSNovelist$", ignore_case = TRUE))
    )
}

maps_files_available <- all(file.exists(file.path(maps_root, paste0(dataset_ids, ".csv"))))
if (maps_files_available) {
  annotation_rows <- bind_rows(
    lapply(dataset_ids, read_maps),
    lapply(dataset_ids, read_reprocessed)
  )

  summary_table <- annotation_rows %>%
    group_by(dataset.ID, source) %>%
    summarise(
      retained_annotation_rows = n(),
      non_repeating_annotations = n_distinct(annotation_id),
      redundant_rows_collapsed = retained_annotation_rows - non_repeating_annotations,
      .groups = "drop"
    ) %>%
    left_join(dataset_metadata, by = "dataset.ID") %>%
    mutate(
      source = factor(source, levels = c("MAPS", "Atlas-library MZmine reprocessing"))
    ) %>%
    arrange(dataset.ID, source)

  source_summary <- annotation_rows %>%
    arrange(
      dataset.ID, source, annotation_id,
      confidence.level, annotation_rank,
      desc(id.prob), desc(confidence.score), source_row
    ) %>%
    group_by(dataset.ID, source, annotation_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    count(dataset.ID, source, annotation_source, name = "non_repeating_annotations") %>%
    left_join(dataset_metadata, by = "dataset.ID") %>%
    mutate(
      source = factor(source, levels = c("MAPS", "Atlas-library MZmine reprocessing"))
    ) %>%
    arrange(dataset.ID, source, desc(non_repeating_annotations))

  identity_membership <- annotation_rows %>%
    distinct(dataset.ID, source, annotation_id) %>%
    mutate(present = TRUE) %>%
    tidyr::pivot_wider(names_from = source, values_from = present, values_fill = FALSE) %>%
    mutate(
      comparison_category = case_when(
        MAPS & `Atlas-library MZmine reprocessing` ~ "Shared",
        MAPS ~ "MAPS only",
        `Atlas-library MZmine reprocessing` ~ "MZmine reprocessing only"
      )
    )

  comparison_table <- identity_membership %>%
    group_by(dataset.ID) %>%
    summarise(
      maps_annotations = sum(MAPS),
      reprocessed_annotations = sum(`Atlas-library MZmine reprocessing`),
      shared_annotations = sum(MAPS & `Atlas-library MZmine reprocessing`),
      maps_only_annotations = sum(MAPS & !`Atlas-library MZmine reprocessing`),
      reprocessed_only_annotations = sum(!MAPS & `Atlas-library MZmine reprocessing`),
      reprocessed_minus_maps = reprocessed_annotations - maps_annotations,
      percent_change_vs_maps = 100 * reprocessed_minus_maps / maps_annotations,
      .groups = "drop"
    ) %>%
    left_join(dataset_metadata, by = "dataset.ID")

  write_csv_stable(
    summary_table,
    file.path(table_dir, "figure4-mzmine-reprocessing-non-repeating-counts.csv")
  )
  write_csv_stable(
    source_summary,
    file.path(table_dir, "figure4-mzmine-reprocessing-source-counts.csv")
  )
  write_csv_stable(
    comparison_table,
    file.path(table_dir, "figure4-mzmine-reprocessing-overlap-summary.csv")
  )
  write_csv_stable(
    annotation_rows %>% distinct(dataset.ID, source, annotation_id, compound.name, smiles),
    file.path(table_dir, "figure4-mzmine-reprocessing-annotation-identities.csv")
  )
  write_csv_stable(
    identity_membership %>% arrange(dataset.ID, comparison_category, annotation_id),
    file.path(table_dir, "figure4-mzmine-reprocessing-identity-membership.csv")
  )
} else {
  message("External MAPS annotation files unavailable; using cached source-count table.")
  source_summary <- readr::read_csv(
    file.path(table_dir, "figure4-mzmine-reprocessing-source-counts.csv"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      source = factor(source, levels = c("MAPS", "Atlas-library MZmine reprocessing"))
    )
}

source_order <- source_summary %>%
  group_by(annotation_source) %>%
  summarise(n = sum(non_repeating_annotations), .groups = "drop") %>%
  arrange(desc(n)) %>%
  pull(annotation_source)
source_summary <- source_summary %>%
  mutate(annotation_source = factor(annotation_source, levels = rev(source_order)))
source_palette <- setNames(
  c("#9A4049", "#404B74", "#BE7440", "#617156", "#A0A5B7", "#D4A280", "#95A08E", "#BB8086", "#8087A2")[
    seq_along(source_order)
  ],
  source_order
)
dataset_axis_labels <- c(
  HGMD_0314 = "HGMD_0314 Phe-Hex +",
  HGMD_0315 = "HGMD_0315 Phe-Hex -",
  HGMD_0318 = "HGMD_0318 HILIC +",
  HGMD_0319 = "HGMD_0319 HILIC -"
)
source_facet_labels <- c(
  "Atlas-library MZmine reprocessing" = "Atlas-library\nMZmine reprocessing",
  MAPS = "MAPS"
)

plot_totals <- source_summary %>%
  group_by(dataset.ID, source) %>%
  summarise(non_repeating_annotations = sum(non_repeating_annotations), .groups = "drop")

p <- ggplot(
  source_summary,
  aes(dataset.ID, non_repeating_annotations, fill = annotation_source)
) +
  geom_col(width = 0.7) +
  geom_text(
    data = plot_totals,
    aes(dataset.ID, non_repeating_annotations, label = fmt_int(non_repeating_annotations)),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 4.6,
    family = "Arial",
    fontface = "bold"
  ) +
  facet_wrap(~source, nrow = 1, labeller = labeller(source = source_facet_labels)) +
  scale_fill_manual(values = source_palette, name = "Annotation source") +
  scale_x_discrete(labels = dataset_axis_labels) +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Source-resolved MAPS and atlas-library annotations",
    x = NULL,
    y = "Annotations"
  ) +
  theme_classic(base_size = 18, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(size = 17, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10.5),
    legend.title = element_text(size = 11.5, face = "bold"),
    strip.text = element_text(size = 12.5, face = "bold", lineheight = 0.9),
    axis.text.x = element_text(size = 10.5, face = "plain", angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 16, 24, 10)
  )

ggsave(
  file.path(figure_dir, "figure4-mzmine-reprocessing-comparison.png"),
  p,
  width = 9.6,
  height = 4.8,
  dpi = 300
)
