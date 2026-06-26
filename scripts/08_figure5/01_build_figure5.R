source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(ggalluvial)
  library(ComplexUpset)
  library(patchwork)
  library(scales)
})
fmt_int <- function(x) scales::number(x, accuracy = 1, big.mark = "")

hmdb_path <- path_value("hmdb_healthy_faecal_csv")
mimedb_cid_path <- path_value("mimedb_pubchem_csv")
mimedb_no_cid_path <- path_value("mimedb_no_cid_csv")
atlas_path <- file.path(project_root, "data", "processed", "paper2_total_list_classified.csv")

for (p in c(hmdb_path, mimedb_cid_path, mimedb_no_cid_path, atlas_path)) {
  if (!file.exists(p)) stop("Missing Figure 5 input: ", p)
}

norm_name <- function(x) {
  x <- str_to_lower(clean_text(x))
  x <- str_replace_all(x, "[[:space:]]+", " ")
  x
}

norm_smiles <- function(x) {
  str_replace_all(clean_text(x), "[[:space:]]+", "")
}

lipid_name_regex <- regex(
  paste(
    c(
      "\\blipid\\b", "phosphatid", "lysophosph", "sphing", "ceramide",
      "glycerolipid", "glycerophosph", "triacylglycer", "diacylglycer",
      "monoacylglycer", "cholest", "steroid", "sterol", "bile acid",
      "fatty acid", "fatty acyl", "fatty amide", "fatty ester",
      "acylcarnitine", "eicosanoid", "octadecanoid"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

atlas_lipid_superclasses <- c(
  "Eicosanoids", "Fatty Acids and Conjugates", "Fatty acyls", "Fatty amides",
  "Fatty esters", "Glycerolipids", "Glycerophospholipids", "Octadecanoids",
  "Sphingolipids", "Steroids"
)

mimedb_lipid_classes <- c(
  "Fatty Acyls", "Glycerolipids", "Glycerophospholipids", "Prenol lipids",
  "Sphingolipids", "Steroids and steroid derivatives"
)

atlas <- readr::read_csv(atlas_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    collection = "Human faecal atlas",
    source_id = annotation_id,
    compound_name = clean_text(compound.name),
    smiles = clean_text(smiles),
    name_key = norm_name(compound.name),
    smiles_key = norm_smiles(smiles),
    explicit_lipid = clean_text(npc_superclass) %in% atlas_lipid_superclasses |
      str_detect(clean_text(npc_superclass), regex("terpenoid|prenol", ignore_case = TRUE)),
    npc_superclass = clean_text(npc_superclass)
  ) %>%
  distinct(collection, source_id, .keep_all = TRUE)

hmdb <- readr::read_csv(hmdb_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    collection = "Healthy HMDB",
    source_id = paste0("HMDB:", row_number()),
    compound_name = clean_text(compound.name),
    smiles = clean_text(smiles),
    name_key = norm_name(compound.name),
    smiles_key = norm_smiles(smiles),
    explicit_lipid = FALSE,
    npc_superclass = ""
  )

mimedb_cid <- readr::read_csv(mimedb_cid_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    collection = "MiMeDB",
    source_id = paste0("MiMeDB:", clean_text(id)),
    compound_name = clean_text(compound.name),
    smiles = clean_text(SMILES),
    name_key = norm_name(compound.name),
    smiles_key = norm_smiles(SMILES),
    explicit_lipid = clean_text(classification) %in% mimedb_lipid_classes,
    npc_superclass = ""
  )

mimedb_no_cid <- readr::read_csv(mimedb_no_cid_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    collection = "MiMeDB",
    source_id = paste0("MiMeDB:", clean_text(id)),
    compound_name = "",
    smiles = clean_text(moldb_smiles),
    name_key = "",
    smiles_key = norm_smiles(moldb_smiles),
    explicit_lipid = clean_text(classification) %in% mimedb_lipid_classes,
    npc_superclass = ""
  )

records <- bind_rows(atlas, hmdb, mimedb_cid, mimedb_no_cid) %>%
  filter(name_key != "" | smiles_key != "") %>%
  mutate(name_lipid = str_detect(compound_name, lipid_name_regex))

# Standardised compound names are primary. SMILES only supplies an identity when
# a standardised name is absent; if that SMILES occurs with a named record, the
# missing-name record inherits the named identity.
smiles_to_name <- records %>%
  filter(name_key != "", smiles_key != "") %>%
  count(smiles_key, name_key, sort = TRUE, name = "mapping_rows") %>%
  group_by(smiles_key) %>%
  arrange(desc(mapping_rows), name_key, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(smiles_key, mapped_name_key = name_key)

records <- records %>%
  left_join(smiles_to_name, by = "smiles_key") %>%
  mutate(
    identity_key = case_when(
      name_key != "" ~ paste0("name:", name_key),
      mapped_name_key != "" ~ paste0("name:", mapped_name_key),
      smiles_key != "" ~ paste0("smiles:", smiles_key),
      TRUE ~ ""
    ),
    match_key_type = case_when(
      name_key != "" ~ "standardised compound.name",
      mapped_name_key != "" ~ "SMILES fallback linked to standardised name",
      TRUE ~ "SMILES fallback only"
    )
  ) %>%
  filter(identity_key != "")

identity_lipid <- records %>%
  group_by(identity_key) %>%
  summarise(is_lipid = any(explicit_lipid | name_lipid), .groups = "drop")

records <- records %>%
  left_join(identity_lipid, by = "identity_key")

membership_all <- records %>%
  distinct(identity_key, collection, is_lipid)

membership_nonlipid <- membership_all %>%
  filter(!is_lipid) %>%
  distinct(identity_key, collection)

collection_levels <- c("Human faecal atlas", "Healthy HMDB", "MiMeDB")

comparison_summary <- records %>%
  group_by(collection) %>%
  summarise(
    source_rows = n(),
    distinct_identities_all = n_distinct(identity_key),
    identities_flagged_lipid = n_distinct(identity_key[is_lipid]),
    distinct_identities_nonlipid = n_distinct(identity_key[!is_lipid]),
    standardised_name_rows = sum(match_key_type == "standardised compound.name"),
    smiles_fallback_rows = sum(match_key_type != "standardised compound.name"),
    .groups = "drop"
  )

wide <- membership_nonlipid %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = collection, values_from = present, values_fill = FALSE)

for (nm in collection_levels) if (!nm %in% names(wide)) wide[[nm]] <- FALSE

intersections <- wide %>%
  mutate(
    signature = case_when(
      .data[["Human faecal atlas"]] & .data[["Healthy HMDB"]] & .data[["MiMeDB"]] ~ "Atlas + healthy HMDB + MiMeDB",
      .data[["Human faecal atlas"]] & .data[["Healthy HMDB"]] ~ "Atlas + healthy HMDB",
      .data[["Human faecal atlas"]] & .data[["MiMeDB"]] ~ "Atlas + MiMeDB",
      .data[["Healthy HMDB"]] & .data[["MiMeDB"]] ~ "Healthy HMDB + MiMeDB",
      .data[["Human faecal atlas"]] ~ "Human faecal atlas only",
      .data[["Healthy HMDB"]] ~ "Healthy HMDB only",
      .data[["MiMeDB"]] ~ "MiMeDB only",
      TRUE ~ "Other"
    )
  ) %>%
  count(signature, name = "intersection_size") %>%
  arrange(desc(intersection_size))

atlas_unique_keys <- wide %>%
  filter(.data[["Human faecal atlas"]], !.data[["Healthy HMDB"]], !.data[["MiMeDB"]]) %>%
  pull(identity_key)

atlas_unique_classes <- records %>%
  filter(
    collection == "Human faecal atlas",
    identity_key %in% atlas_unique_keys,
    !is_lipid,
    npc_superclass != "",
    npc_superclass != "Unclassified"
  ) %>%
  distinct(identity_key, npc_superclass) %>%
  count(npc_superclass, sort = TRUE, name = "atlas_unique_annotations")

atlas_unique_keys_all <- membership_all %>%
  tidyr::pivot_wider(names_from = collection, values_from = collection, values_fn = length, values_fill = 0) %>%
  filter(.data[["Human faecal atlas"]] > 0, .data[["Healthy HMDB"]] == 0, .data[["MiMeDB"]] == 0) %>%
  pull(identity_key)

atlas_details <- readr::read_csv(atlas_path, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    annotation_id = clean_text(annotation_id),
    annotation.type = clean_text(annotation.type),
    confidence.level = suppressWarnings(as.numeric(confidence.level)),
    confidence.score = suppressWarnings(as.numeric(confidence.score)),
    npc_superclass = clean_text(npc_superclass),
    canopus.NPC.superclass = clean_text(canopus.NPC.superclass),
    gnps.library.usi = clean_text(gnps.library.usi),
    gnps.in.silico.bile.acid.info = clean_text(gnps.in.silico.bile.acid.info)
  ) %>%
  transmute(
    annotation_id,
    annotation.type,
    confidence.level,
    confidence.score,
    npc_superclass = if_else(npc_superclass != "", npc_superclass, canopus.NPC.superclass),
    gnps.library.usi,
    gnps.in.silico.bile.acid.info
  )

evidence_rank <- c(
  "In silico bile-acid library" = 1L,
  "Authentic standard" = 2L,
  "Spectral library" = 3L,
  "Spectral neighbour" = 4L,
  "In silico structure" = 5L,
  "Other" = 6L
)
evidence_levels <- c(
  "In silico bile-acid library",
  "Authentic standard",
  "Spectral library",
  "Spectral neighbour",
  "In silico structure",
  "Other"
)

atlas_unique_flow_records <- records %>%
  filter(collection == "Human faecal atlas", identity_key %in% atlas_unique_keys_all) %>%
  distinct(identity_key, annotation_id = source_id) %>%
  left_join(atlas_details, by = "annotation_id") %>%
  mutate(
    annotation_source = case_when(
      gnps.in.silico.bile.acid.info != "" ~ "In silico bile-acid library",
      annotation.type == "authentic standard" | confidence.level == 1 ~ "Authentic standard",
      annotation.type %in% c("gnps", "mzmine") | gnps.library.usi != "" ~ "Spectral library",
      annotation.type == "ms2query" ~ "Spectral neighbour",
      annotation.type == "CSI:FingerID" ~ "In silico structure",
      TRUE ~ "Other"
    ),
    evidence_rank = evidence_rank[annotation_source],
    confidence_level = paste0("Level ", confidence.level),
    npc_superclass = if_else(is.na(npc_superclass) | npc_superclass == "", "Unclassified", npc_superclass)
  ) %>%
  group_by(identity_key) %>%
  arrange(evidence_rank, confidence.level, desc(confidence.score), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

top_flow_classes <- atlas_unique_flow_records %>%
  filter(npc_superclass != "Unclassified") %>%
  count(npc_superclass, sort = TRUE, name = "annotations") %>%
  slice_head(n = 10) %>%
  pull(npc_superclass)

atlas_unique_flow <- atlas_unique_flow_records %>%
  mutate(
    superclass = case_when(
      npc_superclass %in% top_flow_classes ~ npc_superclass,
      npc_superclass == "Unclassified" ~ "Unclassified",
      TRUE ~ "Other classified"
    ),
    superclass = factor(superclass, levels = c(top_flow_classes, "Other classified", "Unclassified")),
    annotation_source = factor(annotation_source, levels = evidence_levels),
    confidence_level = factor(confidence_level, levels = paste0("Level ", 1:3))
  ) %>%
  count(superclass, annotation_source, confidence_level, name = "annotations") %>%
  filter(annotations > 0) %>%
  arrange(superclass, annotation_source, confidence_level)

audit <- records %>%
  select(
    collection, source_id, identity_key, match_key_type, compound_name, smiles,
    explicit_lipid, name_lipid, is_lipid, npc_superclass
  )

write_csv_stable(comparison_summary, file.path(project_root, "outputs", "tables", "figure5-collection-summary.csv"))
write_csv_stable(intersections, file.path(project_root, "outputs", "tables", "figure5-nonlipid-intersections.csv"))
write_csv_stable(atlas_unique_classes, file.path(project_root, "outputs", "tables", "figure5-atlas-unique-nonlipid-classes.csv"))
write_csv_stable(atlas_unique_flow, file.path(project_root, "outputs", "tables", "figure5-atlas-unique-sankey-flows.csv"))
write_csv_stable(audit, file.path(project_root, "outputs", "qc", "figure5-identity-and-lipid-audit.csv"))

atlas_set <- membership_nonlipid %>% filter(collection == "Human faecal atlas") %>% pull(identity_key)
hmdb_set <- membership_nonlipid %>% filter(collection == "Healthy HMDB") %>% pull(identity_key)
mimedb_set <- membership_nonlipid %>% filter(collection == "MiMeDB") %>% pull(identity_key)

set_colours <- c(
  "Human faecal atlas" = "#617156",
  "Healthy HMDB" = "#BE7440",
  "MiMeDB" = "#404B74"
)

theme_figure5 <- function(base_size = 18) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.text = element_text(colour = "#252525"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

overview <- tibble::tibble(
  collection = c("Human faecal\natlas", "Healthy\nHMDB", "MiMeDB", "FASST public\nrepositories"),
  x = 1:4,
  y = 1,
  fill = c("#CAD0C6", "#E9D1BF", "#BFC3D1", "#D6D6D6"),
  border = c("#617156", "#BE7440", "#404B74", "#616161"),
  main = c(
    fmt_int(comparison_summary$distinct_identities_nonlipid[comparison_summary$collection == "Human faecal atlas"]),
    fmt_int(comparison_summary$distinct_identities_nonlipid[comparison_summary$collection == "Healthy HMDB"]),
    fmt_int(comparison_summary$distinct_identities_nonlipid[comparison_summary$collection == "MiMeDB"]),
    "Identity list\nunavailable"
  ),
  detail = c(
    "non-lipid identities",
    "non-lipid identities",
    "non-lipid identities",
    "spectrum-level search only"
  )
)

p_a <- ggplot(overview) +
  geom_tile(aes(x, y, fill = fill, colour = border), width = 0.82, height = 0.82, linewidth = 1.2) +
  geom_text(aes(x, y + 0.19, label = collection), family = "Arial", fontface = "bold", size = 6.2) +
  geom_text(aes(x, y - 0.05, label = main), family = "Arial", fontface = "bold", size = 5.4) +
  geom_text(aes(x, y - 0.30, label = detail), family = "Arial", size = 4.0, colour = "#444444") +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(0.45, 1.55), clip = "off") +
  labs(
    title = "Collection comparison and evidence compatibility"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.margin = margin(12, 20, 12, 20)
  )

p_b <- ggVennDiagram(
  list(
    "Human faecal atlas" = atlas_set,
    "Healthy HMDB" = hmdb_set,
    "MiMeDB" = mimedb_set
  ),
  label_alpha = 0,
  category.names = c("Atlas", "HMDB", "MiMeDB"),
  set_color = unname(set_colours),
  set_size = 6,
  label = "count",
  label_size = 5,
  edge_size = 1
) +
  scale_fill_gradient(low = "#F4F4F4", high = "#9A4049") +
  labs(
    title = "Non-lipid metabolite-identity overlap"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    legend.position = "none"
  )

p_c_data <- intersections %>%
  mutate(signature = factor(signature, levels = rev(signature)))
p_c <- ggplot(p_c_data, aes(intersection_size, signature)) +
  geom_col(fill = "#404B74", width = 0.72) +
  geom_text(
    aes(label = fmt_int(intersection_size)),
    hjust = -0.12, family = "Arial", fontface = "bold", size = 4.4
  ) +
  scale_x_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Exact non-lipid intersections",
    x = "Identity count",
    y = NULL
  ) +
  theme_figure5(17) +
  theme(axis.text.y = element_text(size = 12))

p_c_wrapped <- wrap_elements(full = p_c)

p_d <- ggplot(
  atlas_unique_flow,
  aes(y = annotations, axis1 = superclass, axis2 = annotation_source, axis3 = confidence_level)
) +
  ggalluvial::geom_alluvium(aes(fill = superclass), width = 0.10, alpha = 0.76, knot.pos = 0.45) +
  ggalluvial::geom_stratum(width = 0.16, fill = "#F8F8F8", colour = "#333333", linewidth = 0.28) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    family = "Arial",
    size = 4.0,
    lineheight = 0.88
  ) +
  scale_x_discrete(
    limits = c("Superclass", "Annotation source", "Confidence level"),
    expand = c(0.06, 0.03)
  ) +
  scale_fill_manual(
    values = setNames(rep(c("#617156", "#BE7440", "#9A4049", "#404B74", "#406D80", "#965A78", "#8087A2", "#BB8086", "#809DAA", "#B991A5", "#C8C8C8", "#E0E0E0"), length.out = nlevels(atlas_unique_flow$superclass)), levels(atlas_unique_flow$superclass)),
    guide = "none"
  ) +
  labs(
    title = "Atlas-only identities by class, source and confidence",
    x = NULL,
    y = "Atlas-only non-repeating identities"
  ) +
  theme_figure5(18) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

p_a_wrapped <- wrap_elements(full = p_a)
p_b_wrapped <- wrap_elements(full = p_b)
p_d_wrapped <- wrap_elements(full = p_d)

full <- p_a_wrapped / (p_b_wrapped | p_c_wrapped) / p_d_wrapped +
  plot_layout(heights = c(0.62, 1.15, 1.30)) +
  plot_annotation(
    title = "Figure 5. The human faecal atlas expands identity-resolved metabolite coverage",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", family = "Arial"),
      plot.tag = element_text(size = 24, face = "bold", family = "Arial")
    )
  )

figure_dir <- file.path(project_root, "figures", "figure-4")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "figure5a-collection-overview.png"), p_a, width = 10.8, height = 3.48, dpi = 300)
ggsave(file.path(figure_dir, "figure5b-nonlipid-venn.png"), p_b, width = 6.6, height = 5.4, dpi = 300)
ggsave(file.path(figure_dir, "figure5c-nonlipid-upset.png"), p_c, width = 7.8, height = 5.4, dpi = 300)
ggsave(file.path(figure_dir, "figure5d-atlas-unique-sankey.png"), p_d, width = 10.2, height = 7.2, dpi = 300)
ggsave(file.path(figure_dir, "figure5-assembled.png"), full, width = 14.4, height = 16.8, dpi = 300)

message("Figure 5 complete. FASST is shown descriptively because no identity-level membership list was available.")
