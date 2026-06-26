source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(graphlayouts)
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
  library(scales)
})

# Figure 4/5 utility candidate: NIST public faecal FBMN network mined with
# the atlas-derived MSP library. The plot is exported as PNG only and is not
# inserted into any assembled manuscript figure.

figure_dir <- file.path(project_root, "figures", "figure-4")
table_dir <- file.path(project_root, "outputs", "tables")
qc_dir <- file.path(project_root, "outputs", "qc")
ensure_dirs(figure_dir, table_dir, qc_dir)

nist_dir <- file.path(project_root, "outputs", "mzmine reprocessing", "HGMD_0359")
graphml_path <- list.files(nist_dir, pattern = "\\.graphml$", full.names = TRUE)
if (length(graphml_path) != 1) {
  stop("Expected one NIST GraphML file in ", nist_dir, "; found ", length(graphml_path))
}
graphml_path <- graphml_path[[1]]
annotation_path <- file.path(nist_dir, "data_annotations_NIST.csv")
quant_path <- file.path(nist_dir, "data_iimn_gnps_NIST_quant.csv")
figure5_audit_path <- file.path(qc_dir, "figure5-identity-and-lipid-audit.csv")
atlas_path <- file.path(project_root, "data", "processed", "paper2_total_list_classified.csv")

for (p in c(annotation_path, quant_path, figure5_audit_path, atlas_path)) {
  if (!file.exists(p)) stop("Missing NIST network input: ", p)
}

norm_name <- function(x) {
  str_replace_all(str_to_lower(clean_text(x)), "[[:space:]]+", " ")
}
norm_smiles <- function(x) {
  str_replace_all(clean_text(x), "[[:space:]]+", "")
}

make_identity_key <- function(name, smiles, smiles_to_name) {
  name_key <- norm_name(name)
  smiles_key <- norm_smiles(smiles)
  mapped_name_key <- smiles_to_name$mapped_name_key[match(smiles_key, smiles_to_name$smiles_key)]
  case_when(
    name_key != "" ~ paste0("name:", name_key),
    !is.na(mapped_name_key) & mapped_name_key != "" ~ paste0("name:", mapped_name_key),
    smiles_key != "" ~ paste0("smiles:", smiles_key),
    TRUE ~ NA_character_
  )
}

figure5_audit <- readr::read_csv(figure5_audit_path, show_col_types = FALSE, progress = FALSE)
smiles_to_name <- figure5_audit %>%
  filter(compound_name != "", smiles != "") %>%
  mutate(smiles_key = norm_smiles(smiles), name_key = norm_name(compound_name)) %>%
  count(smiles_key, name_key, sort = TRUE, name = "mapping_rows") %>%
  group_by(smiles_key) %>%
  arrange(desc(mapping_rows), name_key, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(smiles_key, mapped_name_key = name_key)

reference_membership <- figure5_audit %>%
  mutate(identity_key = make_identity_key(compound_name, smiles, smiles_to_name)) %>%
  filter(!is.na(identity_key)) %>%
  distinct(identity_key, collection) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = collection, values_from = present, values_fill = FALSE)
for (nm in c("Healthy HMDB", "MiMeDB")) {
  if (!nm %in% names(reference_membership)) reference_membership[[nm]] <- FALSE
}

atlas_confidence <- readr::read_csv(atlas_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    identity_key = make_identity_key(compound.name, smiles, smiles_to_name),
    atlas_confidence = suppressWarnings(as.integer(confidence.level))
  ) %>%
  filter(!is.na(identity_key), atlas_confidence %in% 1:3) %>%
  group_by(identity_key) %>%
  summarise(atlas_confidence = min(atlas_confidence, na.rm = TRUE), .groups = "drop")

nist_annotations <- readr::read_csv(annotation_path, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    feature_id = as.character(id),
    identity_key = make_identity_key(compound_name, smiles, smiles_to_name),
    atlas_score = suppressWarnings(as.numeric(score))
  ) %>%
  filter(!is.na(identity_key)) %>%
  left_join(atlas_confidence, by = "identity_key") %>%
  left_join(reference_membership, by = "identity_key") %>%
  mutate(
    `Healthy HMDB` = coalesce(.data[["Healthy HMDB"]], FALSE),
    MiMeDB = coalesce(MiMeDB, FALSE)
  ) %>%
  arrange(feature_id, atlas_confidence, desc(atlas_score)) %>%
  group_by(feature_id) %>%
  summarise(
    compound_name = first(na_if(clean_text(compound_name), "")),
    smiles = first(na_if(clean_text(smiles), "")),
    identity_key = first(identity_key),
    atlas_score = max(atlas_score, na.rm = TRUE),
    atlas_confidence = min(atlas_confidence, na.rm = TRUE),
    in_hmdb = any(.data[["Healthy HMDB"]]),
    in_mimedb = any(MiMeDB),
    .groups = "drop"
  ) %>%
  mutate(
    atlas_confidence = if_else(is.infinite(atlas_confidence), NA_integer_, atlas_confidence),
    atlas_score = if_else(is.infinite(atlas_score), NA_real_, atlas_score),
    atlas_hit = !is.na(identity_key),
    absent_from_hmdb_mimedb = atlas_hit & !in_hmdb & !in_mimedb
  )

quant <- readr::read_csv(quant_path, show_col_types = FALSE, progress = FALSE)
sample_cols <- names(quant)[str_detect(names(quant), "^NIST_POS_Samp_\\d+-\\d+\\.mzXML Peak area$")]
subject_ids <- str_match(sample_cols, "^NIST_POS_Samp_(\\d+)-\\d+\\.mzXML Peak area$")[, 2]
subject_levels <- sort(unique(subject_ids))

prevalence <- quant %>%
  transmute(feature_id = as.character(`row ID`), across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "area") %>%
  mutate(
    subject_id = str_match(sample, "^NIST_POS_Samp_(\\d+)-\\d+\\.mzXML Peak area$")[, 2],
    present = !is.na(area) & area > 0
  ) %>%
  group_by(feature_id, subject_id) %>%
  summarise(subject_present = any(present), .groups = "drop") %>%
  group_by(feature_id) %>%
  summarise(
    n_subjects_present = sum(subject_present),
    subject_prevalence = n_subjects_present / length(subject_levels),
    .groups = "drop"
  )

g <- read_graph(graphml_path, format = "graphml")
vertex_id <- vertex_attr(g, "id")
if (is.null(vertex_id) || !length(vertex_id)) vertex_id <- V(g)$name
V(g)$feature_id <- as.character(vertex_id)

node_attrs <- tibble(feature_id = V(g)$feature_id) %>%
  left_join(nist_annotations, by = "feature_id") %>%
  left_join(prevalence, by = "feature_id") %>%
  mutate(
    network_degree = degree(g),
    n_subjects_present = coalesce(n_subjects_present, 0L),
    subject_prevalence = coalesce(subject_prevalence, 0),
    atlas_hit = coalesce(atlas_hit, FALSE),
    absent_from_hmdb_mimedb = coalesce(absent_from_hmdb_mimedb, FALSE),
    confidence_group = case_when(
      atlas_confidence == 1 ~ "Level 1",
      atlas_confidence == 2 ~ "Level 2",
      atlas_confidence == 3 ~ "Level 3",
      TRUE ~ "No atlas match"
    )
  )

atlas_hit_lookup <- setNames(node_attrs$atlas_hit, node_attrs$feature_id)
atlas_hit_indices <- which(atlas_hit_lookup[V(g)$feature_id])
atlas_neighbour_ids <- unique(unlist(lapply(atlas_hit_indices, function(v_idx) {
  V(g)$feature_id[neighbors(g, v_idx)]
}), use.names = FALSE))
node_attrs <- node_attrs %>%
  mutate(
    atlas_neighbour = feature_id %in% atlas_neighbour_ids & !atlas_hit,
    border_group = case_when(
      atlas_hit ~ "Atlas MSP hit",
      atlas_neighbour ~ "Network neighbour",
      TRUE ~ "Other feature"
    )
  )

label_feature_ids <- node_attrs %>%
  filter(atlas_hit, atlas_confidence %in% 1:2) %>%
  arrange(desc(network_degree), atlas_confidence, desc(subject_prevalence), compound_name) %>%
  slice_head(n = 30) %>%
  pull(feature_id)

node_attrs <- node_attrs %>%
  mutate(
    label_group = feature_id %in% label_feature_ids
  )

for (nm in names(node_attrs)) {
  g <- set_vertex_attr(g, nm, value = node_attrs[[nm]][match(V(g)$feature_id, node_attrs$feature_id)])
}

edge_score <- suppressWarnings(as.numeric(E(g)$score))
if (all(is.na(edge_score))) edge_score <- suppressWarnings(as.numeric(E(g)$EdgeScore))
matched_peaks <- suppressWarnings(as.numeric(E(g)$matched_peaks))
if (all(is.na(matched_peaks))) {
  matched_peaks <- suppressWarnings(as.numeric(str_extract(E(g)$matched_peaks, "\\d+")))
}
E(g)$edge_score_plot <- coalesce(edge_score, 0)
E(g)$matched_peaks_plot <- coalesce(matched_peaks, median(matched_peaks, na.rm = TRUE))
E(g)$matched_peaks_plot[is.na(E(g)$matched_peaks_plot)] <- 1

component_tbl <- tibble(component_id = components(g)$membership) %>%
  count(component_id, name = "component_size")
V(g)$component_size <- component_tbl$component_size[match(components(g)$membership, component_tbl$component_id)]

# Keep the full network visible and encode atlas-library information as node
# metadata so unmatched graph structure remains interpretable context.
keep_component <- tapply(V(g)$atlas_hit, components(g)$membership, any)
V(g)$atlas_component <- keep_component[as.character(components(g)$membership)]

tg <- as_tbl_graph(g)
layout <- create_layout(tg, layout = "stress")

spread_node_coordinates <- function(tbl, min_dist = 1.25, iterations = 120, step = 0.45) {
  x <- tbl$x
  y <- tbl$y
  n <- length(x)
  min_dist_sq <- min_dist^2

  for (iter in seq_len(iterations)) {
    for (i in seq_len(n - 1)) {
      dx <- x[i] - x[(i + 1):n]
      dy <- y[i] - y[(i + 1):n]
      dist_sq <- dx^2 + dy^2
      too_close <- which(dist_sq < min_dist_sq)
      if (!length(too_close)) next

      j_idx <- too_close + i
      dist <- sqrt(pmax(dist_sq[too_close], 1e-8))
      dx_close <- dx[too_close]
      dy_close <- dy[too_close]

      zero_dist <- dist < 1e-4
      if (any(zero_dist)) {
        angle <- (i + j_idx[zero_dist]) * 2.399963
        dx_close[zero_dist] <- cos(angle)
        dy_close[zero_dist] <- sin(angle)
        dist[zero_dist] <- 1
      }

      push <- (min_dist - dist) * step * 0.5
      move_x <- push * dx_close / dist
      move_y <- push * dy_close / dist

      x[i] <- x[i] + sum(move_x)
      y[i] <- y[i] + sum(move_y)
      x[j_idx] <- x[j_idx] - move_x
      y[j_idx] <- y[j_idx] - move_y
    }
  }

  tbl$x <- x
  tbl$y <- y
  tbl
}

expand_large_components <- function(tbl) {
  tbl %>%
    group_by(component) %>%
    mutate(
      component_expand = case_when(
        first(component_size) >= 120 ~ 2.25,
        first(component_size) >= 60 ~ 1.95,
        first(component_size) >= 25 ~ 1.55,
        first(component_size) >= 10 ~ 1.25,
        TRUE ~ 1
      ),
      component_x = mean(x, na.rm = TRUE),
      component_y = mean(y, na.rm = TRUE),
      x = component_x + (x - component_x) * component_expand,
      y = component_y + (y - component_y) * component_expand
    ) %>%
    ungroup() %>%
    select(-component_expand, -component_x, -component_y)
}

layout_tbl <- as_tibble(layout) %>%
  mutate(
    x = x * 6.975,
    y = y * 6.4125
  ) %>%
  expand_large_components() %>%
  spread_node_coordinates(min_dist = 7.5375, iterations = 420, step = 0.90) %>%
  mutate(
    halo_radius = rescale(sqrt(subject_prevalence + 0.03), to = c(0.058, 0.231)),
    label_text = if_else(label_group, coalesce(compound_name, feature_id), NA_character_)
  )

edge_ends_tbl <- as_tibble(ends(g, E(g), names = FALSE), .name_repair = "minimal") %>%
  setNames(c("from_idx", "to_idx")) %>%
  mutate(
    from = V(g)$feature_id[from_idx],
    to = V(g)$feature_id[to_idx],
    edge_score_plot = E(g)$edge_score_plot,
    matched_peaks_plot = E(g)$matched_peaks_plot
  )

plot_node_ids <- layout_tbl$feature_id
layout_plot <- layout_tbl

edge_tbl <- edge_ends_tbl %>%
  left_join(layout_tbl %>% select(feature_id, x, y), by = c("from" = "feature_id")) %>%
  rename(x = x, y = y) %>%
  left_join(layout_tbl %>% select(feature_id, xend = x, yend = y), by = c("to" = "feature_id")) %>%
  mutate(
    edge_alpha = rescale(matched_peaks_plot, to = c(0.18, 0.58)),
    edge_width = rescale(edge_score_plot, to = c(0.16, 0.72))
  )

confidence_cols <- c(
  "Level 1" = "#404B74",
  "Level 2" = "#9A4049",
  "Level 3" = "#BE7440",
  "No atlas match" = "#D6D6D6"
)

p_network <- ggplot() +
  geom_segment(
    data = edge_tbl,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = edge_width, alpha = edge_alpha),
    colour = "#6E7480",
    lineend = "round",
    show.legend = FALSE
  ) +
  ggforce::geom_circle(
    data = layout_plot %>% filter(absent_from_hmdb_mimedb),
    aes(x0 = x, y0 = y, r = halo_radius),
    inherit.aes = FALSE,
    colour = NA,
    fill = "#FFE9A6",
    alpha = 0.55
  ) +
  geom_point(
    data = layout_plot %>% filter(absent_from_hmdb_mimedb),
    aes(x, y, colour = "Absent from HMDB and MiMeDB"),
    shape = 21,
    fill = "#FFE9A6",
    size = 3,
    stroke = 0.9,
    alpha = 0.92,
    show.legend = TRUE
  ) +
  geom_point(
    data = layout_plot,
    aes(x, y, size = subject_prevalence, fill = confidence_group),
    shape = 21,
    colour = "white",
    stroke = 0.16,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = layout_plot %>% filter(!is.na(label_text)),
    aes(x, y, label = label_text),
    family = "Arial",
    size = 2.65,
    colour = "#252525",
    min.segment.length = 0,
    segment.size = 0.18,
    box.padding = 0.62,
    point.padding = 0.50,
    max.overlaps = Inf,
    force = 2.2,
    force_pull = 0.08,
    max.time = 2,
    seed = 1
  ) +
  scale_fill_manual(values = confidence_cols, name = "Atlas MSP match") +
  scale_colour_manual(
    values = c("Absent from HMDB and MiMeDB" = "#D9A441"),
    name = "Reference overlap"
  ) +
  scale_size_continuous(
    range = c(0.42, 3.09),
    breaks = c(0, 6 / 18, 12 / 18, 1),
    labels = c("0/18", "6/18", "12/18", "18/18"),
    name = "NIST subjects\nwith feature"
  ) +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  guides(
    fill = guide_legend(override.aes = list(size = 5), order = 1),
    colour = guide_legend(override.aes = list(fill = "#FFE9A6", linewidth = 1.8, alpha = 0.9), order = 2),
    size = guide_legend(
      override.aes = list(
        fill = "#D6D6D6",
        colour = "#6F6F6F",
        shape = 21,
        stroke = 0.3,
        alpha = 1
      ),
      order = 3
    )
  ) +
  labs(
    title = "Atlas-library annotation of the full NIST faecal molecular network"
  ) +
  coord_equal(clip = "off", expand = FALSE) +
  theme_void(base_size = 15, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 11.5, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(4, 8, 4, 4)
  )

network_summary <- tibble(
  graphml_file = basename(graphml_path),
  nodes = vcount(g),
  edges = ecount(g),
  atlas_msp_hit_nodes = sum(V(g)$atlas_hit),
  no_atlas_match_nodes = sum(!V(g)$atlas_hit),
  atlas_neighbour_nodes = sum(V(g)$atlas_neighbour),
  hmdb_mimedb_absent_atlas_hits = sum(V(g)$absent_from_hmdb_mimedb),
  plotted_nodes = length(plot_node_ids),
  plotted_edges = nrow(edge_tbl),
  level_1_nodes = sum(V(g)$atlas_confidence == 1, na.rm = TRUE),
  level_2_nodes = sum(V(g)$atlas_confidence == 2, na.rm = TRUE),
  level_3_nodes = sum(V(g)$atlas_confidence == 3, na.rm = TRUE),
  nist_subjects = length(subject_levels)
)

write_csv_stable(network_summary, file.path(table_dir, "figure4-nist-network-summary.csv"))
write_csv_stable(
  layout_tbl %>%
    mutate(plotted = feature_id %in% plot_node_ids) %>%
    select(feature_id, plotted, component, mz, rt, compound_name, atlas_score, atlas_confidence, confidence_group,
           atlas_hit, atlas_neighbour, in_hmdb, in_mimedb, absent_from_hmdb_mimedb,
           network_degree, n_subjects_present, subject_prevalence, x, y),
  file.path(table_dir, "figure4-nist-network-node-attributes.csv")
)
write_csv_stable(
  edge_tbl,
  file.path(table_dir, "figure4-nist-network-edge-attributes.csv")
)

ggsave(
  file.path(figure_dir, "figure4-nist-atlas-msp-network.png"),
  p_network,
  width = 10.8,
  height = 7.2,
  dpi = 300,
  bg = "white"
)

message("NIST atlas-MSP network plotted with igraph/tidygraph/ggraph.")
