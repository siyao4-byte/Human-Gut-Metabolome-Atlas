source(file.path("scripts", "00_setup", "helpers.R"))

suppressPackageStartupMessages({
  library(future.apply)
  library(httr)
  library(jsonlite)
})

source_cache_path <- paths$npclassifier_cache
paper_cache_path <- file.path(project_root, "data", "processed", "paper2-npclassifier-cache.csv")
query_log_path <- file.path(project_root, "outputs", "qc", "paper2-npclassifier-query-log.csv")

empty_result <- function(smiles, status, raw = NA_character_) {
  tibble(
    SMILES = smiles,
    np_superclass = NA_character_,
    np_class = NA_character_,
    np_pathway = NA_character_,
    is_glycoside = NA,
    raw = raw,
    cached_at = as.character(Sys.time()),
    query_status = status
  )
}

npclassify_one <- function(smiles, timeout_s = 45, attempts = 3) {
  url <- paste0(
    "https://npclassifier.gnps2.org/classify?smiles=",
    utils::URLencode(smiles, reserved = TRUE)
  )

  for (attempt in seq_len(attempts)) {
    response <- try(httr::GET(url, httr::timeout(timeout_s)), silent = TRUE)
    if (!inherits(response, "try-error") && httr::status_code(response) == 200) {
      raw <- httr::content(response, as = "text", encoding = "UTF-8")
      parsed <- try(jsonlite::fromJSON(raw), silent = TRUE)
      if (!inherits(parsed, "try-error")) {
        return(tibble(
          SMILES = smiles,
          np_superclass = paste(parsed$superclass_results %||% character(), collapse = "; "),
          np_class = paste(parsed$class_results %||% character(), collapse = "; "),
          np_pathway = paste(parsed$pathway_results %||% character(), collapse = "; "),
          is_glycoside = isTRUE(parsed$isglycoside),
          raw = raw,
          cached_at = as.character(Sys.time()),
          query_status = "success"
        ))
      }
    }
    Sys.sleep(attempt)
  }

  empty_result(smiles, "request_failed")
}

normalise_cache <- function(x) {
  required <- c(
    "SMILES", "np_superclass", "np_class", "np_pathway",
    "is_glycoside", "raw", "cached_at", "query_status"
  )
  for (column in setdiff(required, names(x))) x[[column]] <- NA

  x %>%
    transmute(
      SMILES = clean_text(SMILES),
      np_superclass = na_if(clean_text(np_superclass), ""),
      np_class = na_if(clean_text(np_class), ""),
      np_pathway = na_if(clean_text(np_pathway), ""),
      is_glycoside,
      raw = as.character(raw),
      cached_at = as.character(cached_at),
      query_status = as.character(query_status)
    ) %>%
    filter(SMILES != "") %>%
    distinct(SMILES, .keep_all = TRUE)
}

atlas_smiles <- readr::read_csv(
  file.path(project_root, "data", "processed", "paper2_total_list.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    confidence.level = suppressWarnings(as.integer(confidence.level)),
    SMILES = clean_text(smiles)
  ) %>%
  filter(confidence.level <= 3, SMILES != "") %>%
  distinct(SMILES)

source_cache <- readr::read_csv(source_cache_path, show_col_types = FALSE, progress = FALSE) %>%
  normalise_cache()

paper_cache <- if (file.exists(paper_cache_path)) {
  readr::read_csv(paper_cache_path, show_col_types = FALSE, progress = FALSE) %>%
    normalise_cache()
} else {
  source_cache
}
write_csv_stable(paper_cache, paper_cache_path)

cache_for_atlas <- atlas_smiles %>%
  left_join(paper_cache %>% mutate(cache_hit = TRUE), by = "SMILES")
to_query <- cache_for_atlas %>%
  filter(
    is.na(cache_hit) |
      query_status == "request_failed" |
      (is.na(np_superclass) & is.na(np_class) & is.na(np_pathway) & (is.na(raw) | raw == ""))
  ) %>%
  pull(SMILES) %>%
  unique()

message("Level 1-3 unique SMILES: ", nrow(atlas_smiles))
message("SMILES requiring an NPClassifier query: ", length(to_query))

query_missing <- is_true(Sys.getenv("NPC_QUERY_MISSING", unset = "false"))
if (length(to_query) && !query_missing) {
  message(
    "External NPClassifier queries are disabled. Set NPC_QUERY_MISSING=true ",
    "only after approving transmission of the missing SMILES to GNPS."
  )
}

if (length(to_query) && query_missing) {
  batch_size <- 250L
  workers <- min(6L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  batches <- split(to_query, ceiling(seq_along(to_query) / batch_size))
  for (i in seq_along(batches)) {
    message("Querying batch ", i, " of ", length(batches), " (", length(batches[[i]]), " SMILES)")
    batch_results <- future.apply::future_lapply(
      batches[[i]],
      npclassify_one,
      future.seed = TRUE
    ) %>%
      bind_rows()

    paper_cache <- bind_rows(
      paper_cache %>% filter(!SMILES %in% batch_results$SMILES),
      batch_results
    ) %>%
      normalise_cache()

    write_csv_stable(paper_cache, paper_cache_path)
  }
}

mapped <- atlas_smiles %>%
  left_join(paper_cache, by = "SMILES") %>%
  mutate(
    has_superclass = !is.na(np_superclass),
    has_class = !is.na(np_class),
    has_pathway = !is.na(np_pathway),
    has_any_classification = has_superclass | has_class | has_pathway
  )

qc <- tibble(
  metric = c(
    "level_1_3_unique_smiles", "cache_rows_for_level_1_3_smiles",
    "smiles_with_superclass", "smiles_with_class", "smiles_with_pathway",
    "smiles_with_any_classification", "smiles_without_any_classification",
    "request_failures"
  ),
  value = c(
    nrow(mapped), sum(mapped$SMILES %in% paper_cache$SMILES),
    sum(mapped$has_superclass), sum(mapped$has_class), sum(mapped$has_pathway),
    sum(mapped$has_any_classification), sum(!mapped$has_any_classification),
    sum(mapped$query_status == "request_failed", na.rm = TRUE)
  )
)

write_csv_stable(qc, query_log_path)
write_csv_stable(
  mapped %>% filter(!has_any_classification),
  file.path(project_root, "outputs", "qc", "paper2-npclassifier-unclassified-smiles.csv")
)
print(qc)
