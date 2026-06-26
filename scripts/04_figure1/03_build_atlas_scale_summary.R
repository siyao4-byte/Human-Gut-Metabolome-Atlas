source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

manifest_candidates <- c(
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.current.csv"),
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.csv")
)
manifest_path <- manifest_candidates[file.exists(manifest_candidates)][[1]]
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE) %>%
  filter(include_figure1)
included_datasets <- unique(manifest$HGMD.ID)

analysis <- read_hgm_sheet("A - Analysis") %>%
  mutate(across(everything(), clean_text)) %>%
  separate_rows(dataset.ID, sep = ";\\s*") %>%
  filter(dataset.ID %in% included_datasets)
extracts <- read_hgm_sheet("E - Extracts") %>%
  mutate(across(everything(), clean_text))
fractions <- read_hgm_sheet("F - Fractions") %>%
  mutate(across(everything(), clean_text))
humans <- read_hgm_sheet("H - Human Samples") %>%
  mutate(across(everything(), clean_text))

node_map <- bind_rows(
  extracts %>% transmute(node_id = HGME.ID, parent_id = Parent.Code, node_type = "extract"),
  fractions %>% transmute(node_id = HGMF.ID, parent_id = Parent.Code, node_type = "fraction"),
  humans %>% transmute(node_id = HGMH.ID, parent_id = "", node_type = "human")
) %>%
  filter(node_id != "") %>%
  distinct(node_id, .keep_all = TRUE)

extract_hgm_ids <- function(x) {
  ids <- stringr::str_extract_all(clean_text(x), "HGM[EFH]_\\d+")
  sort(unique(unlist(ids, use.names = FALSE)))
}

trace_ancestors <- function(start_id, max_depth = 20L) {
  visited <- character()
  unresolved <- character()

  walk <- function(ids, depth = 1L) {
    if (depth > max_depth) return(invisible(NULL))
    for (current in extract_hgm_ids(ids)) {
      if (current %in% visited) next
      visited <<- c(visited, current)
      if (stringr::str_starts(current, "HGMH_")) next
      row <- node_map[node_map$node_id == current, , drop = FALSE]
      if (!nrow(row)) {
        unresolved <<- c(unresolved, current)
        next
      }
      walk(row$parent_id[[1]], depth + 1L)
    }
    invisible(NULL)
  }

  walk(start_id)
  list(nodes = sort(unique(visited)), unresolved = sort(unique(unresolved)))
}

analysis_lineage <- lapply(analysis$Parent.Code, trace_ancestors)
all_nodes <- sort(unique(unlist(lapply(analysis_lineage, `[[`, "nodes"), use.names = FALSE)))
unresolved_nodes <- sort(unique(unlist(lapply(analysis_lineage, `[[`, "unresolved"), use.names = FALSE)))

used_fraction_ids <- all_nodes[stringr::str_starts(all_nodes, "HGMF_")]
used_human_ids <- all_nodes[stringr::str_starts(all_nodes, "HGMH_")]
used_donor_ids <- humans %>%
  filter(HGMH.ID %in% used_human_ids, Donor.ID != "") %>%
  distinct(Donor.ID) %>%
  pull(Donor.ID)

has_upstream_fraction <- function(fraction_id) {
  row <- fractions %>% filter(HGMF.ID == fraction_id) %>% slice_head(n = 1)
  if (!nrow(row)) return(NA)
  upstream <- trace_ancestors(row$Parent.Code)$nodes
  any(stringr::str_starts(upstream, "HGMF_"))
}

fraction_classes <- tibble(HGMF.ID = used_fraction_ids) %>%
  mutate(
    is_subfraction = vapply(HGMF.ID, has_upstream_fraction, logical(1)),
    fraction_level = case_when(
      is.na(is_subfraction) ~ "Unresolved fraction hierarchy",
      is_subfraction ~ "Subfractions generated",
      TRUE ~ "Primary fractions generated"
    )
  )

summary_table <- tibble(
  metric = c(
    "Human donors represented",
    "Primary fractions generated",
    "Subfractions generated",
    "LC-MS/MS datasets run"
  ),
  display_label = c(
    "Human donors\nrepresented",
    "Primary fractions\ngenerated",
    "Subfractions\ngenerated",
    "LC-MS/MS datasets\nrun"
  ),
  value = c(
    length(used_donor_ids),
    sum(!fraction_classes$is_subfraction, na.rm = TRUE),
    sum(fraction_classes$is_subfraction, na.rm = TRUE),
    length(included_datasets)
  ),
  definition = c(
    "Unique non-empty Donor.ID values linked through HGMH ancestors to included dataset injections",
    "Unique HGMF records without an upstream HGMF ancestor",
    "Unique HGMF records with an upstream HGMF ancestor",
    "Good-quality untargeted datasets retained in Figure 1"
  ),
  order = seq_len(4)
)

audit <- tibble(
  metric = c(
    "included_datasets",
    "analysis_rows_linked_to_included_datasets",
    "unique_analysis_samples",
    "unique_fraction_records",
    "unresolved_fraction_hierarchy",
    "unique_donor_ids",
    "unique_human_ancestors",
    "unresolved_lineage_nodes"
  ),
  value = c(
    length(included_datasets),
    nrow(analysis),
    n_distinct(analysis$HGMA.ID),
    length(used_fraction_ids),
    sum(is.na(fraction_classes$is_subfraction)),
    length(used_donor_ids),
    length(used_human_ids),
    length(unresolved_nodes)
  )
)

write_csv_stable(summary_table, file.path(project_root, "outputs", "tables", "figure1-atlas-scale-summary.csv"))
write_csv_stable(fraction_classes, file.path(project_root, "outputs", "tables", "figure1-atlas-fraction-hierarchy.csv"))
write_csv_stable(audit, file.path(project_root, "outputs", "qc", "figure1-atlas-scale-summary-audit.csv"))
write_csv_stable(
  tibble(unresolved_node = unresolved_nodes),
  file.path(project_root, "outputs", "qc", "figure1-atlas-scale-summary-unresolved-lineage.csv")
)

card_cols <- c("#404B74", "#617156", "#BE7440", "#9A4049")
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")
summary_table <- summary_table %>%
  mutate(
    x = order,
    label = fmt_int(value),
    metric = factor(metric, levels = metric)
  )

p_summary <- ggplot(summary_table, aes(x, 0)) +
  geom_rect(
    aes(xmin = x - 0.43, xmax = x + 0.43, ymin = -0.72, ymax = 0.72, fill = metric),
    colour = "#303030", linewidth = 0.45
  ) +
  geom_text(
    aes(y = 0.23, label = label),
    family = "Arial", fontface = "bold", size = 12.5, colour = "white"
  ) +
  geom_text(
    aes(y = -0.22, label = display_label),
    family = "Arial", fontface = "bold", size = 5.8, lineheight = 0.92, colour = "white"
  ) +
  annotate(
    "segment", x = 1.45, xend = 1.55, y = 0, yend = 0,
    arrow = arrow(length = unit(0.10, "inches")), linewidth = 0.65, colour = "#606060"
  ) +
  annotate(
    "segment", x = 2.45, xend = 2.55, y = 0, yend = 0,
    arrow = arrow(length = unit(0.10, "inches")), linewidth = 0.65, colour = "#606060"
  ) +
  annotate(
    "segment", x = 3.45, xend = 3.55, y = 0, yend = 0,
    arrow = arrow(length = unit(0.10, "inches")), linewidth = 0.65, colour = "#606060"
  ) +
  scale_fill_manual(values = setNames(card_cols, levels(summary_table$metric))) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-0.9, 1.0), clip = "off") +
  labs(
    title = "Scale of the human faecal fraction atlas",
    caption = "One linked fraction record lacked hierarchy metadata and was excluded from primary/subfraction counts."
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 24, face = "bold", hjust = 0),
    plot.caption = element_text(size = 10.5, colour = "#606060", hjust = 0),
    plot.margin = margin(18, 20, 18, 20),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

figure_dir <- file.path(project_root, "figures", "figure-1")
base <- file.path(figure_dir, "figure1-atlas-scale-summary")
ggsave(paste0(base, ".png"), p_summary, width = 12.8, height = 4.64, dpi = 300, bg = "white")

message("Standalone Figure 1 atlas-scale summary written; assembled Figure 1 was not modified.")
