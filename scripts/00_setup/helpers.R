suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(yaml)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
paths_file <- file.path(project_root, "config", "paths.yml")
paths_example_file <- file.path(project_root, "config", "paths.example.yml")
if (file.exists(paths_file)) {
  paths <- yaml::read_yaml(paths_file)
} else if (file.exists(paths_example_file)) {
  paths <- yaml::read_yaml(paths_example_file)
  warning(
    "Using config/paths.example.yml. Copy it to config/paths.yml and update ",
    "local paths before running analyses that need external inputs.",
    call. = FALSE
  )
} else {
  stop("Missing config/paths.yml or config/paths.example.yml", call. = FALSE)
}

# University of Melbourne Gen 3 Design System v15.12.0 colour palette.
uom_colours <- c(
  heritage = "#000F46",
  heritage_dark = "#000B34",
  sage_light = "#ABC1A7",
  sage_dark = "#444A40",
  grey_light = "#C8C8C8",
  grey_dark = "#2D2D2D",
  blue_light = "#46C8F0",
  blue_dark = "#003C55",
  pink = "#EB7BBE",
  maroon = "#73234B",
  red_light = "#FF2D3C",
  red_dark = "#78000D",
  yellow = "#FFD629",
  brown = "#A84500",
  green_light = "#9FB825",
  green_dark = "#2C421D"
)

uom_discrete <- c(
  "#000F46", "#404B74", "#8087A2", "#BFC3D1",
  "#000B34", "#003C55", "#406D80", "#809DAA", "#BFCED5",
  "#46C8F0", "#74D6F4", "#A3E4F7", "#D1F1FB",
  "#ABC1A7", "#C0D0BD", "#D5E0D3", "#EAEFE9",
  "#444A40", "#737770", "#A2A49F", "#D0D2CF",
  "#2D2D2D", "#616161", "#969696", "#CACACA",
  "#C8C8C8", "#D6D6D6", "#E4E4E4", "#F1F1F1",
  "#EB7BBE", "#F09CCE", "#F5BDDF", "#FADEEF",
  "#73234B", "#965A78", "#B991A5", "#DCC8D2",
  "#FF2D3C", "#FF616D", "#FF969D", "#FFCACE",
  "#78000D", "#9A4049", "#BB8086", "#DDBFC2",
  "#FFD629", "#FFE05E", "#FFEA94", "#FFF5C9",
  "#A84500", "#BE7440", "#D4A280", "#E9D1BF",
  "#9FB825", "#B7CA5B", "#CFDC92", "#E7EDC8",
  "#2C421D", "#617156", "#95A08E", "#CAD0C6"
)
uom_confidence <- c("1" = uom_colours[["heritage"]], "2" = uom_colours[["blue_light"]], "3" = uom_colours[["yellow"]])

ensure_dirs <- function(...) {
  dirs <- c(...)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

ensure_dirs(
  file.path(project_root, "data", "metadata"),
  file.path(project_root, "data", "intermediate"),
  file.path(project_root, "data", "processed"),
  file.path(project_root, "outputs", "qc"),
  file.path(project_root, "outputs", "tables")
)

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  str_squish(x)
}

is_true <- function(x) {
  clean_text(x) %in% c("1", "TRUE", "True", "true", "Yes", "yes")
}

dataset_annotation_path <- function(dataset_id) {
  flat <- file.path(paths$annotation_root, paste0(dataset_id, ".csv"))
  nested <- file.path(paths$annotation_root, dataset_id, paste0(dataset_id, ".csv"))
  if (file.exists(flat)) return(flat)
  if (file.exists(nested)) return(nested)
  NA_character_
}

dataset_abundance_path <- function(dataset_id) {
  candidates <- c(
    file.path(paths$abundance_root, paste0(dataset_id, "-samples-df.csv")),
    file.path(paths$abundance_root, paste0(dataset_id, "-samples_df.csv")),
    file.path(paths$abundance_root, paste0(dataset_id, "-samples_df.csv"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

write_csv_stable <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
  message("Wrote ", path, " (", nrow(x), " rows)")
}

path_value <- function(name, default = NA_character_) {
  value <- paths[[name]]
  if (is.null(value) || identical(value, "")) default else value
}

read_hgm_sheet <- function(sheet) {
  workbook <- path_value("hgm_workbook")
  if (!is.na(workbook) && file.exists(workbook)) {
    return(readxl::read_excel(workbook, sheet = sheet))
  }

  csv_dir <- path_value("hgm_csv_export_dir")
  if (!is.na(csv_dir) && dir.exists(csv_dir)) {
    csv_path <- file.path(csv_dir, paste0(sheet, ".csv"))
    if (file.exists(csv_path)) {
      return(readr::read_csv(
        csv_path,
        show_col_types = FALSE,
        col_types = readr::cols(.default = "c")
      ))
    }
  }

  stop(
    "Could not find HGM workbook or CSV export for sheet: ", sheet,
    ". Set hgm_workbook or hgm_csv_export_dir in config/paths.yml.",
    call. = FALSE
  )
}
