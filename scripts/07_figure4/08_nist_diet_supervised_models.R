source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(mixOmics)
  library(randomForest)
})

set.seed(20260618)

nist_dir <- file.path(project_root, "outputs", "mzmine reprocessing", "HGMD_0359")
metadata_path <- file.path(nist_dir, "gnps_metadata.tsv")
annotation_path <- file.path(nist_dir, "data_annotations_NIST.csv")
quant_path <- file.path(nist_dir, "data_iimn_gnps_NIST_quant.csv")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "figures", "figure-4")

theme_supervised <- function(base_size = 15) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#272727"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#272727"),
      axis.text = element_text(colour = "#272727"),
      axis.title = element_text(colour = "#272727"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size - 0.5, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      panel.grid = element_blank()
    )
}

metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    filename = clean_text(filename),
    diet = clean_text(ATTRIBUTE_diet_info),
    sample_type = clean_text(SampleType),
    subject = str_match(filename, "Samp_([0-9]{2})-")[, 2],
    replicate = str_match(filename, "Samp_[0-9]{2}-([0-9]{2})")[, 2],
    quant_col = paste0(str_remove(filename, "\\.mzXML$"), ".mzXML Peak area")
  )

nist_samples <- metadata %>%
  filter(sample_type == "animal", !is.na(subject), diet %in% c("omnivore", "vegetarian")) %>%
  distinct(filename, quant_col, subject, replicate, diet)

annotations <- readr::read_csv(annotation_path, show_col_types = FALSE) %>%
  transmute(
    feature_id = as.character(id),
    compound_name = clean_text(compound_name),
    smiles = clean_text(smiles),
    atlas_score = suppressWarnings(as.numeric(score)),
    precursor_mz = suppressWarnings(as.numeric(precursor_mz)),
    rt = suppressWarnings(as.numeric(rt))
  ) %>%
  filter(compound_name != "") %>%
  arrange(feature_id, desc(atlas_score), compound_name) %>%
  group_by(feature_id) %>%
  summarise(
    top_compound_name = first(compound_name),
    candidate_names = paste(head(unique(compound_name), 3), collapse = "; "),
    n_candidate_names = n_distinct(compound_name),
    best_atlas_score = max(atlas_score, na.rm = TRUE),
    precursor_mz = first(precursor_mz),
    rt = first(rt),
    .groups = "drop"
  ) %>%
  mutate(best_atlas_score = if_else(is.infinite(best_atlas_score), NA_real_, best_atlas_score))

quant <- readr::read_csv(quant_path, show_col_types = FALSE) %>%
  rename(feature_id = `row ID`) %>%
  mutate(feature_id = as.character(feature_id)) %>%
  semi_join(annotations, by = "feature_id")

sample_cols <- intersect(nist_samples$quant_col, names(quant))
if (length(sample_cols) == 0) {
  stop("No NIST sample columns in quant table matched gnps_metadata.tsv filenames.")
}

subject_feature <- quant %>%
  dplyr::select(feature_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "quant_col",
    values_to = "peak_area"
  ) %>%
  left_join(nist_samples, by = "quant_col") %>%
  mutate(
    peak_area = suppressWarnings(as.numeric(peak_area)),
    present = !is.na(peak_area) & peak_area > 0,
    log10_area = log10(coalesce(peak_area, 0) + 1)
  ) %>%
  group_by(quant_col) %>%
  mutate(
    sample_median_log10_area = median(log10_area[present], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    sample_median_log10_area = if_else(
      is.finite(sample_median_log10_area),
      sample_median_log10_area,
      median(sample_median_log10_area[is.finite(sample_median_log10_area)], na.rm = TRUE)
    ),
    global_median_log10_area = median(sample_median_log10_area, na.rm = TRUE),
    normalised_log10_area = log10_area - sample_median_log10_area + global_median_log10_area
  ) %>%
  group_by(subject, diet, feature_id) %>%
  summarise(
    mean_normalised_log10_area = mean(normalised_log10_area, na.rm = TRUE),
    present_replicates = sum(present, na.rm = TRUE),
    subject_present = present_replicates > 0,
    .groups = "drop"
  )

feature_matrix <- subject_feature %>%
  dplyr::select(subject, diet, feature_id, mean_normalised_log10_area) %>%
  pivot_wider(
    names_from = feature_id,
    values_from = mean_normalised_log10_area,
    values_fill = 0
  ) %>%
  arrange(subject)

sample_info <- feature_matrix %>% dplyr::select(subject, diet)
x_raw <- feature_matrix %>% dplyr::select(-subject, -diet)

feature_variance <- vapply(x_raw, stats::var, numeric(1), na.rm = TRUE)
keep_features <- names(feature_variance)[is.finite(feature_variance) & feature_variance > 0]
x_raw <- x_raw[, keep_features, drop = FALSE]

x_scaled <- scale(as.matrix(x_raw))
colnames(x_scaled) <- make.names(colnames(x_raw), unique = TRUE)
feature_name_map <- tibble(
  feature_id = colnames(x_raw),
  model_feature = colnames(x_scaled)
)
y <- factor(sample_info$diet, levels = c("omnivore", "vegetarian"))

plsda_fit <- mixOmics::plsda(x_scaled, y, ncomp = 2)
plsda_scores <- as_tibble(plsda_fit$variates$X[, 1:2], .name_repair = "minimal") %>%
  setNames(c("comp1", "comp2")) %>%
  bind_cols(sample_info) %>%
  mutate(diet = factor(diet, levels = c("omnivore", "vegetarian")))

explained <- plsda_fit$prop_expl_var$X[1:2] * 100

vip_matrix <- tryCatch(mixOmics::vip(plsda_fit), error = function(e) NULL)
if (is.null(vip_matrix)) {
  vip_tbl <- tibble(
    model_feature = rownames(plsda_fit$loadings$X),
    vip_component_1 = abs(plsda_fit$loadings$X[, 1]),
    vip_component_2 = abs(plsda_fit$loadings$X[, 2])
  )
} else {
  vip_tbl <- as_tibble(vip_matrix, rownames = "model_feature")
  names(vip_tbl)[-1] <- paste0("vip_component_", seq_len(ncol(vip_tbl) - 1))
}
vip_tbl <- vip_tbl %>%
  mutate(vip_max = do.call(pmax, c(across(starts_with("vip_component_")), na.rm = TRUE))) %>%
  left_join(feature_name_map, by = "model_feature") %>%
  left_join(annotations, by = "feature_id") %>%
  arrange(desc(vip_max))

plsda_perf <- tryCatch(
  mixOmics::perf(plsda_fit, validation = "Mfold", folds = 3, nrepeat = 50, progressBar = FALSE),
  error = function(e) NULL
)
plsda_error <- if (is.null(plsda_perf)) {
  NA_real_
} else {
  as.numeric(plsda_perf$error.rate$overall[2, "centroids.dist"])
}

rf_data <- as.data.frame(x_scaled) %>%
  mutate(diet = y)

rf_fit <- randomForest::randomForest(
  diet ~ .,
  data = rf_data,
  ntree = 5000,
  mtry = max(1, floor(sqrt(ncol(x_scaled)))),
  importance = TRUE
)

rf_importance_raw <- randomForest::importance(rf_fit, type = 1)
rf_importance <- tibble(
  model_feature = rownames(rf_importance_raw),
  mean_decrease_accuracy = as.numeric(rf_importance_raw[, 1])
) %>%
  left_join(feature_name_map, by = "model_feature") %>%
  left_join(annotations, by = "feature_id") %>%
  arrange(desc(mean_decrease_accuracy))

loo_predictions <- lapply(seq_len(nrow(rf_data)), function(i) {
  train <- rf_data[-i, , drop = FALSE]
  test <- rf_data[i, , drop = FALSE]
  fit <- randomForest::randomForest(
    diet ~ .,
    data = train,
    ntree = 2000,
    mtry = max(1, floor(sqrt(ncol(x_scaled)))),
    importance = FALSE
  )
  prob <- predict(fit, test, type = "prob")
  tibble(
    subject = sample_info$subject[i],
    observed = as.character(y[i]),
    predicted = as.character(predict(fit, test)),
    omnivore_probability = prob[, "omnivore"],
    vegetarian_probability = prob[, "vegetarian"]
  )
}) %>%
  bind_rows() %>%
  mutate(correct = observed == predicted)

model_summary <- tibble(
  metric = c(
    "subjects",
    "omnivore_subjects",
    "vegetarian_subjects",
    "atlas_matched_feature_ids_input",
    "model_features_after_zero_variance_filter",
    "plsda_components",
    "plsda_repeated_3fold_cv_overall_error_centroids_dist",
    "random_forest_ntree",
    "random_forest_oob_error",
    "random_forest_loocv_accuracy"
  ),
  value = c(
    nrow(sample_info),
    sum(y == "omnivore"),
    sum(y == "vegetarian"),
    n_distinct(quant$feature_id),
    ncol(x_scaled),
    2,
    plsda_error,
    rf_fit$ntree,
    tail(rf_fit$err.rate[, "OOB"], 1),
    mean(loo_predictions$correct)
  )
)

diet_cols <- c(omnivore = "#BE7440", vegetarian = "#617156")
subject_label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(size = 4, family = "Arial", show.legend = FALSE, max.overlaps = Inf)
} else {
  geom_text(size = 4, family = "Arial", show.legend = FALSE, vjust = -0.8, check_overlap = TRUE)
}

p_scores <- ggplot(plsda_scores, aes(comp1, comp2, colour = diet, label = subject)) +
  geom_hline(yintercept = 0, colour = "#D8D8D8", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#D8D8D8", linewidth = 0.35) +
  stat_ellipse(aes(group = diet), linewidth = 0.55, alpha = 0.8, show.legend = FALSE) +
  geom_point(size = 4.2, alpha = 0.95) +
  subject_label_layer +
  scale_colour_manual(values = diet_cols, name = "Diet") +
  labs(
    title = "NIST atlas-matched PLS-DA",
    x = paste0("Component 1 (", round(explained[1], 1), "% X variance)"),
    y = paste0("Component 2 (", round(explained[2], 1), "% X variance)")
  ) +
  theme_supervised(15)

top_vip <- vip_tbl %>%
  slice_head(n = 25) %>%
  mutate(label = str_trunc(top_compound_name, 55), label = reorder(label, vip_max))
p_vip <- ggplot(top_vip, aes(vip_max, label, fill = best_atlas_score)) +
  geom_col(width = 0.74, colour = "#303030", linewidth = 0.25) +
  scale_fill_gradient(low = "#D9E2EA", high = "#404B74", na.value = "#BDBDBD", name = "Atlas\nscore") +
  labs(
    title = "Top PLS-DA VIP features",
    x = "VIP score",
    y = NULL
  ) +
  theme_supervised(13)

top_rf <- rf_importance %>%
  filter(is.finite(mean_decrease_accuracy)) %>%
  slice_head(n = 25) %>%
  mutate(label = str_trunc(top_compound_name, 55), label = reorder(label, mean_decrease_accuracy))
p_rf <- ggplot(top_rf, aes(mean_decrease_accuracy, label, fill = best_atlas_score)) +
  geom_col(width = 0.74, colour = "#303030", linewidth = 0.25) +
  scale_fill_gradient(low = "#F0DEC8", high = "#9A4049", na.value = "#BDBDBD", name = "Atlas\nscore") +
  labs(
    title = "Top random-forest diet predictors",
    x = "Mean decrease accuracy",
    y = NULL
  ) +
  theme_supervised(13)

write_csv_stable(plsda_scores, file.path(table_dir, "figure4-nist-plsda-scores.csv"))
write_csv_stable(vip_tbl, file.path(table_dir, "figure4-nist-plsda-vip-scores.csv"))
write_csv_stable(rf_importance, file.path(table_dir, "figure4-nist-random-forest-importance.csv"))
write_csv_stable(loo_predictions, file.path(table_dir, "figure4-nist-random-forest-loocv-predictions.csv"))
write_csv_stable(model_summary, file.path(table_dir, "figure4-nist-supervised-model-summary.csv"))

ggsave(file.path(figure_dir, "figure4-nist-plsda-scores.png"), p_scores, width = 4.32, height = 3.48, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure4-nist-plsda-vip-top25.png"), p_vip, width = 4.92, height = 4.56, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure4-nist-random-forest-top25.png"), p_rf, width = 4.92, height = 4.56, dpi = 300, bg = "white")

message("NIST diet supervised models complete.")
message("PLS-DA CV error (3-fold repeated): ", signif(plsda_error, 3))
message("Random forest OOB error: ", signif(tail(rf_fit$err.rate[, "OOB"], 1), 3))
message("Random forest LOOCV accuracy: ", signif(mean(loo_predictions$correct), 3))
