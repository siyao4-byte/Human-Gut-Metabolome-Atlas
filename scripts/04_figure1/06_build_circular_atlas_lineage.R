source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(ggplot2)
  library(readxl)
})

manifest <- readr::read_csv(
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.current.csv"),
  show_col_types = FALSE
) %>%
  filter(include_figure1)

analysis <- read_hgm_sheet("A - Analysis") %>%
  mutate(across(everything(), clean_text)) %>%
  separate_rows(dataset.ID, sep = ";\\s*") %>%
  filter(dataset.ID %in% manifest$HGMD.ID)
extracts <- read_hgm_sheet("E - Extracts") %>%
  mutate(across(everything(), clean_text))
fractions <- read_hgm_sheet("F - Fractions") %>%
  mutate(across(everything(), clean_text))
humans <- read_hgm_sheet("H - Human Samples") %>%
  mutate(across(everything(), clean_text))

extract_ids <- function(x) {
  unique(unlist(str_extract_all(coalesce(x, ""), "HGM[EFH]_\\d+[A-Z_]*"), use.names = FALSE))
}

raw_edges <- bind_rows(
  extracts %>% transmute(child = HGME.ID, parent_raw = Parent.Code),
  fractions %>% transmute(child = HGMF.ID, parent_raw = Parent.Code)
) %>%
  filter(child != "") %>%
  rowwise() %>%
  mutate(parent = list(extract_ids(parent_raw))) %>%
  ungroup() %>%
  unnest_longer(parent) %>%
  filter(parent != "") %>%
  distinct(parent, child)

start_nodes <- unique(unlist(lapply(analysis$Parent.Code, extract_ids), use.names = FALSE))
used_nodes <- start_nodes
frontier <- start_nodes
for (i in seq_len(30)) {
  parents <- raw_edges %>% filter(child %in% frontier) %>% pull(parent) %>% unique()
  parents <- setdiff(parents, used_nodes)
  if (!length(parents)) break
  used_nodes <- union(used_nodes, parents)
  frontier <- parents
}

edges <- raw_edges %>% filter(parent %in% used_nodes, child %in% used_nodes)
branch_parents <- unique(edges$parent[str_detect(edges$parent, "^HGMF_")])

fraction_meta <- fractions %>%
  transmute(
    id = HGMF.ID,
    parent_key = Parent.Code,
    sequence = suppressWarnings(as.numeric(str_extract(coalesce(
      if ("Fraction Sequence" %in% names(fractions)) .data[["Fraction Sequence"]] else "",
      if ("Full.Book.Code" %in% names(fractions)) .data[["Full.Book.Code"]] else "",
      HGMF.ID
    ), "\\d+(?=\\D*$)"))),
    description = str_to_lower(paste(
      coalesce(if ("Fractionation Method" %in% names(fractions)) .data[["Fractionation Method"]] else "", ""),
      coalesce(if ("Sample.Type" %in% names(fractions)) .data[["Sample.Type"]] else "", ""),
      coalesce(if ("Full.Book.Code" %in% names(fractions)) .data[["Full.Book.Code"]] else "", ""),
      coalesce(if ("Notes" %in% names(fractions)) .data[["Notes"]] else "", "")
    ))
  ) %>%
  filter(id %in% used_nodes) %>%
  group_by(parent_key) %>%
  mutate(
    sibling_rank = row_number(),
    sibling_n = n(),
    keep_sibling = sibling_rank == 1 | sibling_rank == sibling_n | id %in% branch_parents
  ) %>%
  ungroup()

keep_fraction <- fraction_meta %>% filter(keep_sibling) %>% pull(id)
keep_nodes <- used_nodes[!str_detect(used_nodes, "^HGMF_") | used_nodes %in% keep_fraction]
edges <- edges %>%
  filter(parent %in% keep_nodes, child %in% keep_nodes) %>%
  group_by(child) %>%
  slice_head(n = 1) %>%
  ungroup()

used_hgmh <- keep_nodes[str_detect(keep_nodes, "^HGMH_")]
donor_alias <- humans %>%
  filter(HGMH.ID %in% used_hgmh) %>%
  distinct(Donor.ID) %>%
  arrange(Donor.ID) %>%
  mutate(donor_alias = sprintf("Donor_%02d", row_number()))
human_alias <- humans %>%
  filter(HGMH.ID %in% used_hgmh) %>%
  distinct(HGMH.ID, Donor.ID) %>%
  left_join(donor_alias, by = "Donor.ID")
edges <- bind_rows(
  human_alias %>% transmute(parent = donor_alias, child = HGMH.ID),
  edges
) %>%
  distinct(parent, child)

fraction_types <- fraction_meta %>%
  mutate(
    node_group = case_when(
      str_detect(description, "hilic") ~ "HILIC subfraction",
      str_detect(description, "sax") ~ "SAX subfraction",
      str_detect(description, "c18|reverse|reversed|bh01-40") ~ "C18 subfraction",
      str_detect(description, "non.?polar|dcm|phe.?hex|phenyl") ~ "Non-polar fraction",
      str_detect(description, "polar|aqueous|water") ~ "Polar fraction",
      id %in% branch_parents ~ "Further-fractionated node",
      TRUE ~ "Other fraction"
    )
  ) %>%
  select(id, node_group)

nodes <- tibble(id = unique(c(edges$parent, edges$child))) %>%
  left_join(fraction_types, by = "id") %>%
  mutate(
    node_group = case_when(
      str_detect(id, "^Donor_") ~ "Donor",
      str_detect(id, "^HGMH_") ~ "Human sample",
      str_detect(id, "^HGME_") ~ "Extract",
      TRUE ~ coalesce(node_group, "Other fraction")
    ),
    label = if_else(node_group == "Human sample", id, "")
  )

lineage_cols <- c(
  "Donor" = "#303030",
  "Human sample" = "#617156",
  "Extract" = "#A7A7A7",
  "Polar fraction" = "#8DBDCB",
  "Non-polar fraction" = "#D39A72",
  "C18 subfraction" = "#C7A36A",
  "SAX subfraction" = "#A68AB0",
  "HILIC subfraction" = "#7299B5",
  "Further-fractionated node" = "#9A4049",
  "Other fraction" = "#CFCFCF"
)

graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
roots <- V(graph)$name[degree(graph, mode = "in") == 0]
if (length(roots) > 1) {
  graph <- add_vertices(graph, 1, name = "Atlas")
  graph <- add_edges(graph, as.vector(rbind(rep("Atlas", length(roots)), roots)))
  V(graph)$node_group[V(graph)$name == "Atlas"] <- "Donor"
  V(graph)$label[V(graph)$name == "Atlas"] <- "Atlas"
}

p_lineage <- ggraph(graph, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(colour = "#B7B7B7", alpha = 0.62, linewidth = 0.22) +
  geom_node_point(aes(colour = node_group, size = node_group), alpha = 0.95) +
  geom_node_text(
    aes(label = label, angle = -((-node_angle(x, y) + 90) %% 180) + 90),
    size = 1.9, family = "Arial", repel = FALSE
  ) +
  scale_colour_manual(
    values = lineage_cols,
    breaks = names(lineage_cols),
    name = "Lineage node",
    guide = guide_legend(ncol = 1, override.aes = list(size = 3.4, alpha = 1))
  ) +
  scale_size_manual(
    values = c(
      "Donor" = 3.8, "Human sample" = 3.0, "Extract" = 1.4,
      "Polar fraction" = 1.15, "Non-polar fraction" = 1.15,
      "C18 subfraction" = 1.15, "SAX subfraction" = 1.15,
      "HILIC subfraction" = 1.15, "Further-fractionated node" = 2.1,
      "Other fraction" = 0.8
    ),
    guide = "none"
  ) +
  labs(
    title = "Atlas fraction and subfraction lineage"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(l = 8),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 18, 8, 8)
  )

write_csv_stable(nodes, file.path(project_root, "outputs", "tables", "figure1a-circular-lineage-nodes.csv"))
write_csv_stable(edges, file.path(project_root, "outputs", "tables", "figure1a-circular-lineage-edges.csv"))
write_csv_stable(
  fraction_meta,
  file.path(project_root, "outputs", "qc", "figure1a-fraction-type-inference-audit.csv")
)
write_csv_stable(
  human_alias %>% select(HGMH.ID, donor_alias),
  file.path(project_root, "outputs", "qc", "figure1a-deidentified-donor-aliases.csv")
)

ggsave(
  file.path(project_root, "figures", "figure-1", "figure1a-circular-lineage.png"),
  p_lineage, width = 10.4, height = 7.2, dpi = 300, bg = "white"
)
