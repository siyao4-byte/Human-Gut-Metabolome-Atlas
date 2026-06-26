source(file.path("scripts", "00_setup", "helpers.R"))

manifest_candidates <- c(
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.csv"),
  file.path(project_root, "data", "metadata", "paper2-dataset-manifest.current.csv")
)
manifest_candidates <- manifest_candidates[file.exists(manifest_candidates)]
if (!length(manifest_candidates)) stop("No Paper 2 dataset manifest found.")
manifest_path <- manifest_candidates[which.max(file.info(manifest_candidates)$mtime)]

manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
paper2_total <- readr::read_csv(
  file.path(project_root, "data", "processed", "paper2_total_list.csv"),
  show_col_types = FALSE
)

output_dir <- file.path(project_root, "outputs", "figure4", "library-inputs")
ensure_dirs(output_dir, file.path(project_root, "outputs", "tables"))

read_rt_evidence <- function(dataset_id, annotation_path, column_type) {
  message("Reading RT evidence from ", dataset_id)
  x <- readr::read_csv(
    annotation_path,
    col_select = any_of(c("compound.name", "rt", "confidence.level", "annotation.type")),
    show_col_types = FALSE,
    progress = FALSE
  )
  missing <- setdiff(c("compound.name", "rt", "confidence.level", "annotation.type"), names(x))
  for (col in missing) x[[col]] <- NA

  x %>%
    transmute(
      compound.name = na_if(clean_text(compound.name), ""),
      rt = suppressWarnings(as.numeric(rt)),
      confidence.level = suppressWarnings(as.numeric(confidence.level)),
      annotation.type = clean_text(annotation.type),
      dataset.ID = dataset_id,
      column.type = clean_text(column_type)
    ) %>%
    filter(
      !is.na(compound.name),
      !is.na(rt),
      !is.na(confidence.level),
      confidence.level <= 3,
      !str_detect(annotation.type, regex("^MSNovelist$", ignore_case = TRUE)),
      !str_detect(compound.name, regex("analogue|candidate|PUBCHEM", ignore_case = TRUE))
    )
}

rt_evidence <- bind_rows(Map(
  read_rt_evidence,
  manifest$HGMD.ID,
  manifest$annotation_path,
  manifest$column.type
))

densest_rt_cluster <- function(x) {
  rts <- x %>%
    distinct(rt, dataset.ID) %>%
    arrange(rt)
  frequency <- n_distinct(rts$dataset.ID)
  if (!nrow(rts) || frequency == 0) {
    return(tibble(
      Frequency = 0L, rt.frequency = 0L,
      cluster_mean_rt = NA_real_, cluster_sd_rt = NA_real_
    ))
  }

  clusters <- lapply(seq_len(nrow(rts)), function(i) {
    in_window <- rts %>% filter(rt >= rts$rt[[i]], rt <= rts$rt[[i]] + 0.2)
    tibble(
      rt.frequency = n_distinct(in_window$dataset.ID),
      cluster_mean_rt = mean(in_window$rt),
      cluster_sd_rt = sd(in_window$rt)
    )
  })
  best <- bind_rows(clusters) %>% slice_max(rt.frequency, n = 1, with_ties = FALSE)
  best %>% mutate(Frequency = frequency, .before = 1)
}

message("Calculating compound/platform RT clusters")
rt_statistics <- rt_evidence %>%
  group_by(compound.name, column.type) %>%
  group_modify(~ densest_rt_cluster(.x)) %>%
  ungroup() %>%
  mutate(rt.confidence.score = rt.frequency / Frequency)

library_rows <- paper2_total %>%
  filter(column.type %in% c("Phe-Hex", "HILIC", "SAX"), ion.mode %in% c("pos", "neg")) %>%
  left_join(rt_statistics, by = c("compound.name", "column.type")) %>%
  mutate(
    rt = suppressWarnings(as.numeric(rt)),
    rt_conf_ge_0_75 = !is.na(rt.confidence.score) & rt.confidence.score >= 0.75,
    rt_for_lib_t075 = if_else(
      rt_conf_ge_0_75 & !is.na(cluster_mean_rt),
      cluster_mean_rt,
      rt
    )
  )

platform_names <- c("Phe-Hex" = "PheHex", "HILIC" = "HILIC", "SAX" = "SAX")
input_groups <- tidyr::expand_grid(
  column.type = names(platform_names),
  ion.mode = c("pos", "neg")
)

input_summary <- bind_rows(lapply(seq_len(nrow(input_groups)), function(i) {
  platform <- input_groups$column.type[[i]]
  mode <- input_groups$ion.mode[[i]]
  output_name <- paste0(platform_names[[platform]], "-", mode, "-paper2-rt075.csv")
  out <- library_rows %>% filter(column.type == platform, ion.mode == mode)
  write_csv_stable(out, file.path(output_dir, output_name))
  tibble(
    column.type = platform,
    ion.mode = mode,
    output_name = output_name,
    rows = nrow(out),
    usable_usi = n_distinct(out$feature.usi[!is.na(out$feature.usi) & out$feature.usi != ""]),
    rt_cluster_mean_used = sum(out$rt_conf_ge_0_75, na.rm = TRUE)
  )
}))

write_csv_stable(
  rt_statistics,
  file.path(project_root, "outputs", "tables", "figure4-rt-cluster-statistics.csv")
)
write_csv_stable(
  input_summary,
  file.path(project_root, "outputs", "tables", "figure4-library-input-summary.csv")
)
print(input_summary)
