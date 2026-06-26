source(file.path("scripts", "00_setup", "helpers.R"))

input_path <- file.path(project_root, "data", "processed", "paper2_total_list.csv")
if (!file.exists(input_path)) stop("Run the Paper 2 total-list merger first.")

x <- readr::read_csv(input_path, show_col_types = FALSE)
paper_cache_path <- file.path(project_root, "data", "processed", "paper2-npclassifier-cache.csv")
cache_path <- if (file.exists(paper_cache_path)) paper_cache_path else paths$npclassifier_cache
npc <- readr::read_csv(cache_path, show_col_types = FALSE, progress = FALSE)

pick_col <- function(candidates) {
  hit <- candidates[candidates %in% names(npc)]
  if (length(hit)) hit[[1]] else NA_character_
}

smiles_col <- pick_col(c("SMILES", "smiles"))
pathway_col <- pick_col(c("np_pathway", "pathway", "NP_Pathway"))
superclass_col <- pick_col(c("np_superclass", "superclass", "NP_Superclass"))
class_col <- pick_col(c("np_class", "class", "subclass", "NP_Class"))

if (is.na(smiles_col)) stop("No SMILES column found in NPClassifier cache.")

npc_small <- tibble(
  smiles_key = str_to_lower(clean_text(npc[[smiles_col]])),
  npc_pathway = if (!is.na(pathway_col)) clean_text(npc[[pathway_col]]) else NA_character_,
  npc_superclass = if (!is.na(superclass_col)) clean_text(npc[[superclass_col]]) else NA_character_,
  npc_class = if (!is.na(class_col)) clean_text(npc[[class_col]]) else NA_character_
) %>%
  filter(smiles_key != "") %>%
  distinct(smiles_key, .keep_all = TRUE)

out <- x %>%
  mutate(smiles = clean_text(smiles), smiles_key = str_to_lower(smiles)) %>%
  left_join(npc_small, by = "smiles_key") %>%
  select(-smiles_key) %>%
  mutate(
    npc_pathway = na_if(npc_pathway, ""),
    npc_superclass = na_if(npc_superclass, ""),
    npc_class = na_if(npc_class, "")
  )

qc <- tibble(
  metric = c(
    "rows", "rows_with_smiles", "rows_with_npclassifier_superclass",
    "rows_with_npclassifier_class", "rows_with_npclassifier_pathway"
  ),
  value = c(
    nrow(out), sum(out$smiles != ""), sum(!is.na(out$npc_superclass)),
    sum(!is.na(out$npc_class)), sum(!is.na(out$npc_pathway))
  )
)

classified_path <- file.path(project_root, "data", "processed", "paper2_total_list_classified.csv")
classified_fallback_path <- file.path(
  project_root, "data", "processed", "paper2_total_list_classified.current.csv"
)
classified_written <- tryCatch(
  {
    write_csv_stable(out, classified_path)
    classified_path
  },
  error = function(e) {
    warning(
      "Canonical classified atlas is locked; writing the fresh rebuild to ",
      classified_fallback_path
    )
    write_csv_stable(out, classified_fallback_path)
    classified_fallback_path
  }
)
write_csv_stable(qc, file.path(project_root, "outputs", "qc", "npclassifier-coverage.csv"))
message("Current classified atlas: ", classified_written)
print(qc)
