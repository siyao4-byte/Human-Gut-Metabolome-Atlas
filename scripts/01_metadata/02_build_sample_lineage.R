source(file.path("scripts", "00_setup", "helpers.R"))

analysis <- readxl::read_excel(paths$hgm_workbook, sheet = "A - Analysis") %>%
  mutate(across(everything(), clean_text)) %>%
  separate_rows(dataset.ID, sep = ";\\s*") %>%
  mutate(dataset.ID = clean_text(dataset.ID))
extracts <- readxl::read_excel(paths$hgm_workbook, sheet = "E - Extracts") %>%
  mutate(across(everything(), clean_text))
fractions <- readxl::read_excel(paths$hgm_workbook, sheet = "F - Fractions") %>%
  mutate(across(everything(), clean_text))
humans <- readxl::read_excel(paths$hgm_workbook, sheet = "H - Human Samples") %>%
  mutate(across(everything(), clean_text))
datasets <- readxl::read_excel(paths$hgm_workbook, sheet = "D - Dataset") %>%
  mutate(across(everything(), clean_text))

figure23_ids <- c(sprintf("HGMD_%04d", 314:321), sprintf("HGMD_%04d", 336:351))
excluded_hgmh_ids <- c("HGMH_0064")
node_map <- bind_rows(
  extracts %>% transmute(node_id = HGME.ID, parent_id = Parent.Code, node_type = "extract"),
  fractions %>% transmute(node_id = HGMF.ID, parent_id = Parent.Code, node_type = "fraction"),
  humans %>% transmute(node_id = HGMH.ID, parent_id = NA_character_, node_type = "human")
) %>%
  filter(!is.na(node_id), node_id != "") %>%
  distinct(node_id, .keep_all = TRUE)

split_ids <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(character())
  out <- unlist(str_split(x, ";\\s*"), use.names = FALSE)
  out <- clean_text(out)
  out[out != ""]
}

resolve_lineage <- function(start_id, max_depth = 12L) {
  visited <- character()
  edges <- character()
  unresolved <- character()

  walk <- function(ids, depth = 1L) {
    if (depth > max_depth) return(invisible(NULL))
    for (current in split_ids(ids)) {
      if (current %in% visited) next
      visited <<- c(visited, current)
      if (str_starts(current, "HGMH_")) next
      row <- node_map[node_map$node_id == current, , drop = FALSE]
      if (!nrow(row)) {
        unresolved <<- c(unresolved, current)
        next
      }
      parents <- split_ids(row$parent_id[[1]])
      if (!length(parents)) next
      edges <<- c(edges, paste0(current, " > ", parents))
      walk(parents, depth + 1L)
    }
    invisible(NULL)
  }

  walk(start_id)
  hgmh <- sort(unique(visited[str_starts(visited, "HGMH_")]))
  tibble(
    lineage = paste(edges, collapse = " | "),
    hgmh_ids = paste(hgmh, collapse = "; "),
    hgmh_count = length(hgmh),
    unresolved_ids = paste(sort(unique(unresolved)), collapse = "; ")
  )
}

analysis_by_hgma <- analysis %>%
  filter(str_starts(HGMA.ID, "HGMA_")) %>%
  group_by(HGMA.ID) %>%
  summarise(
    Parent.Code = first(Parent.Code[Parent.Code != ""]),
    Book.Code = first(Book.Code[Book.Code != ""]),
    source_workbook_dataset_ids = paste(sort(unique(dataset.ID[dataset.ID != ""])), collapse = "; "),
    .groups = "drop"
  )

read_observed_samples <- function(dataset_id) {
  path <- dataset_abundance_path(dataset_id)
  if (is.na(path)) return(tibble(HGMA.ID = character(), dataset.ID = character()))
  readr::read_csv(
    path,
    col_select = any_of(c("samples", "Samples")),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    rename_with(~ "HGMA.ID", matches("^(samples|Samples)$")) %>%
    transmute(
      HGMA.ID = str_remove(clean_text(HGMA.ID), "\\.raw\\.area$"),
      dataset.ID = dataset_id
    ) %>%
    filter(str_detect(HGMA.ID, "^HGMA_\\d+$")) %>%
    distinct()
}

observed_samples <- bind_rows(lapply(figure23_ids, read_observed_samples))

# Positive/negative pairs represent the same injection set. Use their union for
# metadata completeness; actual presence remains defined from each abundance file.
dataset_pairs <- tibble(
  pair_id = rep(seq_len(12), each = 2),
  dataset.ID = figure23_ids
)

observed_samples <- observed_samples %>%
  left_join(dataset_pairs, by = "dataset.ID") %>%
  group_by(pair_id) %>%
  complete(HGMA.ID, dataset.ID = unique(dataset.ID)) %>%
  ungroup() %>%
  select(HGMA.ID, dataset.ID) %>%
  distinct()

base_samples <- observed_samples %>%
  left_join(analysis_by_hgma, by = "HGMA.ID") %>%
  distinct(HGMA.ID, dataset.ID, .keep_all = TRUE)

message("Figure 2/3 observed HGMA + dataset rows found: ", nrow(base_samples))
if (!nrow(base_samples)) {
  stop("No Figure 2/3 analysis rows found in A - Analysis.")
}

lineage_rows <- lapply(base_samples$Parent.Code, resolve_lineage)
lineage_df <- bind_rows(lineage_rows)
if (!"hgmh_ids" %in% names(lineage_df)) {
  stop("Lineage resolver did not produce hgmh_ids. Columns: ", paste(names(lineage_df), collapse = ", "))
}

human_lookup <- humans %>%
  select(HGMH.ID, Donor.ID) %>%
  filter(HGMH.ID != "") %>%
  distinct(HGMH.ID, .keep_all = TRUE)

resolve_donors <- function(hgmh_ids) {
  ids <- split_ids(hgmh_ids)
  donors <- human_lookup$Donor.ID[match(ids, human_lookup$HGMH.ID)]
  paste(sort(unique(donors[!is.na(donors) & donors != ""])), collapse = "; ")
}

sample_map <- bind_cols(base_samples, lineage_df) %>%
  left_join(
    datasets %>% select(HGMD.ID, column.type, ion.mode, quality, dataset_notes = notes),
    by = c("dataset.ID" = "HGMD.ID")
  ) %>%
  left_join(
    extracts %>% select(
      extract_id = HGME.ID, extract_parent = Parent.Code,
      extract_book_code = Full.Book.Code, extract_scale = Scale,
      extract_sample_type = Sample.Type, extraction_solvent = `Extraction Solvent`,
      extract_replicate = Replicate, extract_notes = Notes
    ),
    by = c("Parent.Code" = "extract_id")
  ) %>%
  left_join(
    fractions %>% select(
      fraction_id = HGMF.ID, fraction_parent = Parent.Code,
      fraction_book_code = Full.Book.Code,
      fraction_sequence = Fraction.Subidentifier,
      fraction_sample_type = Sample.Type, fraction_notes = Notes,
      fraction_notes_2 = `Notes 2`
    ),
    by = c("Parent.Code" = "fraction_id")
  ) %>%
  mutate(
    resolved_hgmh_ids = hgmh_ids,
    resolved_hgmh_count = hgmh_count,
    donor_ids = vapply(hgmh_ids, resolve_donors, character(1)),
    processing_level = case_when(
      dataset.ID %in% c("HGMD_0314", "HGMD_0315", "HGMD_0318", "HGMD_0319") ~ "methanol_crude",
      dataset.ID %in% c("HGMD_0316", "HGMD_0317") ~ "dcm_crude",
      dataset.ID %in% c("HGMD_0320", "HGMD_0321") ~ "water_crude",
      dataset.ID %in% sprintf("HGMD_%04d", 336:343) ~ "polar_fraction",
      dataset.ID %in% sprintf("HGMD_%04d", 344:351) ~ "nonpolar_fraction",
      TRUE ~ "review"
    ),
    pool_id = case_when(
      processing_level == "polar_fraction" ~ str_extract(Book.Code, "BH01-\\d+"),
      processing_level == "nonpolar_fraction" ~ str_extract(Book.Code, "SL01-\\d+"),
      TRUE ~ NA_character_
    ),
    lineage_override_reason = case_when(
      processing_level == "nonpolar_fraction" & pool_id == "SL01-96" ~
        "Canonical Pool_1 override: workbook SL01-96 parent graph incorrectly resolves most fractions to HGMH_0053 only",
      TRUE ~ NA_character_
    ),
    hgmh_ids = case_when(
      !is.na(lineage_override_reason) ~ "HGMH_0053; HGMH_0057; HGMH_0068",
      TRUE ~ hgmh_ids
    ),
    hgmh_count = lengths(str_split(hgmh_ids, ";\\s*")),
    donor_ids = vapply(hgmh_ids, resolve_donors, character(1)),
    presence_rule = "area > 0"
  ) %>%
  filter(!vapply(
    hgmh_ids,
    function(ids) any(split_ids(ids) %in% excluded_hgmh_ids),
    logical(1)
  )) %>%
  arrange(dataset.ID, HGMA.ID)

message(
  "Excluded all Figure 2/3 records originating from: ",
  paste(excluded_hgmh_ids, collapse = ", ")
)

qc <- sample_map %>%
  summarise(
    samples = n(),
    datasets = n_distinct(dataset.ID),
    unresolved_hgmh = sum(hgmh_count == 0),
    unresolved_parent = sum(unresolved_ids != ""),
    fraction_samples = sum(str_detect(processing_level, "fraction")),
    fraction_sequence_missing = sum(
      str_detect(processing_level, "fraction") &
        (is.na(fraction_sequence) | fraction_sequence == "")
    ),
    canonical_lineage_overrides = sum(!is.na(lineage_override_reason))
  ) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value")

pool_summary <- sample_map %>%
  filter(str_detect(processing_level, "fraction")) %>%
  count(pool_id, processing_level, dataset.ID, column.type, ion.mode, name = "sample_count")

write_csv_stable(sample_map, file.path(project_root, "data", "metadata", "figure23-sample-lineage.csv"))
write_csv_stable(pool_summary, file.path(project_root, "data", "metadata", "figure2-four-pool-sample-map.csv"))
write_csv_stable(qc, file.path(project_root, "outputs", "qc", "sample-lineage-summary.csv"))
print(qc)
print(pool_summary, n = Inf)
