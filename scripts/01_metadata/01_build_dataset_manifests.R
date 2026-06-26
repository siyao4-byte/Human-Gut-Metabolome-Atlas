source(file.path("scripts", "00_setup", "helpers.R"))

d <- readxl::read_excel(paths$hgm_workbook, sheet = "D - Dataset") %>%
  mutate(
    across(everything(), as.character),
    annotation_path = vapply(HGMD.ID, dataset_annotation_path, character(1)),
    abundance_path = vapply(HGMD.ID, dataset_abundance_path, character(1)),
    annotation_available = !is.na(annotation_path),
    abundance_available = !is.na(abundance_path),
    targeted_flag = is_true(Targeted),
    search_text = str_to_lower(paste(
      clean_text(`column.type`), clean_text(notes),
      clean_text(`processed.data.folder`), clean_text(`book.code`)
    )),
    small_scale_flag = str_detect(search_text, "small[- _]?scale"),
    lipidomic_flag = str_detect(search_text, "lipidomic|lipidomics|lipidome"),
    relevant_platform = str_detect(
      search_text,
      "hilic|phe[- ]?hex|c18|sax"
    ),
    inclusion_reason = case_when(
      HGMD.ID == "HGMD_0294" ~ "excluded: small-scale bile acid dataset",
      small_scale_flag ~ "excluded: small-scale dataset",
      lipidomic_flag ~ "excluded: lipidomic dataset",
      quality != "Good" ~ "excluded: quality is not Good",
      !annotation_available ~ "excluded: annotation CSV unavailable",
      !targeted_flag ~ "included: good untargeted dataset",
      relevant_platform ~ "included: good targeted relevant platform dataset",
      TRUE ~ "excluded: targeted dataset outside specified platforms"
    ),
    include_figure1 = str_starts(inclusion_reason, "included:")
  ) %>%
  select(
    HGMD.ID, acquisition.date, book.code, targeted_flag, QC, quality,
    column.type, ion.mode, raw.data.folder, processed.data.folder, notes,
    annotation_available, annotation_path, abundance_available, abundance_path,
    small_scale_flag, lipidomic_flag, relevant_platform, inclusion_reason,
    include_figure1
  ) %>%
  arrange(HGMD.ID)

figure1 <- d %>% filter(include_figure1)

figure23_ids <- c(sprintf("HGMD_%04d", 314:321), sprintf("HGMD_%04d", 336:351))
figure23 <- d %>%
  filter(HGMD.ID %in% figure23_ids) %>%
  mutate(
    analysis_role = case_when(
      HGMD.ID %in% sprintf("HGMD_%04d", 314:321) ~ "crude",
      str_detect(column.type, regex("HILIC", ignore_case = TRUE)) ~ "polar_fraction",
      str_detect(column.type, regex("Phe[- ]?Hex", ignore_case = TRUE)) ~ "nonpolar_fraction",
      TRUE ~ "review"
    )
  )

write_csv_stable(d, file.path(project_root, "outputs", "qc", "dataset-selection-audit.csv"))
figure1_manifest_path <- file.path(project_root, "data", "metadata", "paper2-dataset-manifest.csv")
figure1_manifest_fallback_path <- file.path(
  project_root, "data", "metadata", "paper2-dataset-manifest.current.csv"
)
tryCatch(
  write_csv_stable(figure1, figure1_manifest_path),
  error = function(e) {
    warning(
      "Canonical Figure 1 manifest is locked; writing the fresh manifest to ",
      figure1_manifest_fallback_path
    )
    write_csv_stable(figure1, figure1_manifest_fallback_path)
  }
)
write_csv_stable(figure23, file.path(project_root, "data", "metadata", "figure23-dataset-manifest.csv"))

summary <- tibble(
  metric = c(
    "Figure 1 included datasets",
    "Figure 1 untargeted datasets",
    "Figure 1 targeted datasets",
    "Small-scale datasets explicitly excluded",
    "Lipidomic datasets explicitly excluded",
    "Figure 2/3 expected datasets",
    "Figure 2/3 annotation files available",
    "Figure 2/3 abundance files available"
  ),
  value = c(
    nrow(figure1),
    sum(!figure1$targeted_flag),
    sum(figure1$targeted_flag),
    sum(d$small_scale_flag),
    sum(d$lipidomic_flag),
    length(figure23_ids),
    sum(figure23$annotation_available),
    sum(figure23$abundance_available)
  )
)
write_csv_stable(summary, file.path(project_root, "outputs", "qc", "dataset-selection-summary.csv"))
print(summary)
