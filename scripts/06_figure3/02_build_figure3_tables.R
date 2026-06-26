# Figure 3 table build.
# Purpose: convert normalized Figure 3 profiles into PCA scores, crude
# clustering, pool contrasts, prevalence categories, chemical-class summaries,
# and PERMANOVA statistics for plotting.
source(file.path("scripts", "00_setup", "helpers.R"))
suppressPackageStartupMessages({
  library(limma)
  library(vegan)
})

x <- readr::read_csv(file.path(project_root, "data", "intermediate", "figure3-abundance-long.csv"), show_col_types = FALSE)
presence <- readr::read_csv(file.path(project_root, "data", "intermediate", "figure3-pool-presence.csv"), show_col_types = FALSE)
x <- x %>%
  mutate(
    analysis_batch = paste(column.type, ion.mode, processing_level, sep = "::"),
    analysis_feature = paste(analysis_batch, annotation_id, sep = "::")
  ) %>%
  group_by(state, profile_id, combined_pool_id, block_id, processing_level, analysis_batch, analysis_feature) %>%
  summarise(
    value = max(value, na.rm = TRUE),
    annotation_id = first(annotation_id),
    compound.name = first(compound.name),
    smiles = first(smiles),
    confidence.level = min(confidence.level, na.rm = TRUE),
    .groups = "drop"
  )

make_profile_matrix <- function(df) {
  # Rows are profiles and columns are state-specific annotation features. The
  # analysis_feature key keeps column type, ion mode, processing level, and
  # annotation_id distinct before modelling.
  wide <- df %>%
    select(profile_id, analysis_feature, value) %>%
    group_by(profile_id, analysis_feature) %>%
    summarise(value = max(value), .groups = "drop") %>%
    pivot_wider(names_from = analysis_feature, values_from = value, values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$profile_id
  keep <- apply(mat, 2, sd, na.rm = TRUE) > 0
  mat[, keep, drop = FALSE]
}

crude <- x %>% filter(state == "Crude")
fractionated <- x %>% filter(state == "Fractionated")
crude_mat <- make_profile_matrix(crude)
fraction_mat <- make_profile_matrix(fractionated)

crude_meta <- crude %>% distinct(profile_id, combined_pool_id)
fraction_meta <- fractionated %>% distinct(profile_id, combined_pool_id, block_id, processing_level)

crude_pca <- prcomp(crude_mat, center = TRUE, scale. = TRUE)
fraction_pca <- prcomp(fraction_mat, center = TRUE, scale. = TRUE)

pooled_crude_long <- crude %>%
  # For the pooled-crude PCA panel, individual crude profiles are summarized
  # within each pool and crude phase by median abundance before PCA.
  group_by(combined_pool_id, processing_level, annotation_id) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(profile_id = paste(combined_pool_id, processing_level, sep = "::"))
pooled_crude_wide <- pooled_crude_long %>%
  select(profile_id, annotation_id, value) %>%
  pivot_wider(names_from = annotation_id, values_from = value, values_fill = 0)
pooled_crude_mat <- as.matrix(pooled_crude_wide[, -1, drop = FALSE])
rownames(pooled_crude_mat) <- pooled_crude_wide$profile_id
pooled_crude_mat <- pooled_crude_mat[, apply(pooled_crude_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
pooled_crude_meta <- pooled_crude_long %>%
  distinct(profile_id, combined_pool_id, processing_level)
pooled_crude_pca <- prcomp(pooled_crude_mat, center = TRUE, scale. = TRUE)

pooled_fraction_long <- fractionated %>%
  # Candidate standalone PCA: summarize all polar and non-polar fraction
  # profiles within each pool before ordination.
  group_by(combined_pool_id, analysis_feature) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(profile_id = combined_pool_id, processing_level = "all_fractions")
pooled_fraction_wide <- pooled_fraction_long %>%
  select(profile_id, analysis_feature, value) %>%
  pivot_wider(names_from = analysis_feature, values_from = value, values_fill = 0)
pooled_fraction_mat <- as.matrix(pooled_fraction_wide[, -1, drop = FALSE])
rownames(pooled_fraction_mat) <- pooled_fraction_wide$profile_id
pooled_fraction_mat <- pooled_fraction_mat[, apply(pooled_fraction_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
pooled_fraction_meta <- pooled_fraction_long %>%
  distinct(profile_id, combined_pool_id, processing_level)
pooled_fraction_pca <- prcomp(pooled_fraction_mat, center = TRUE, scale. = TRUE)

individual_crude_pca_scores <- as_tibble(crude_pca$x[, 1:2, drop = FALSE], rownames = "profile_id") %>%
  left_join(crude_meta, by = "profile_id") %>%
  mutate(
    state = "Individual crude",
    pc1_variance = crude_pca$sdev[1]^2 / sum(crude_pca$sdev^2),
    pc2_variance = crude_pca$sdev[2]^2 / sum(crude_pca$sdev^2)
  )
fraction_pca_scores <- as_tibble(fraction_pca$x[, 1:2, drop = FALSE], rownames = "profile_id") %>%
  left_join(fraction_meta, by = "profile_id") %>%
  mutate(
    state = "Fractionated",
    pc1_variance = fraction_pca$sdev[1]^2 / sum(fraction_pca$sdev^2),
    pc2_variance = fraction_pca$sdev[2]^2 / sum(fraction_pca$sdev^2)
  )
pca_scores <- bind_rows(
  individual_crude_pca_scores %>% mutate(state = "Crude"),
  fraction_pca_scores
)
pooled_crude_pca_scores <- as_tibble(
  pooled_crude_pca$x[, 1:2, drop = FALSE], rownames = "profile_id"
) %>%
  left_join(pooled_crude_meta, by = "profile_id") %>%
  mutate(
    state = "Pooled crude",
    pc1_variance = pooled_crude_pca$sdev[1]^2 / sum(pooled_crude_pca$sdev^2),
    pc2_variance = pooled_crude_pca$sdev[2]^2 / sum(pooled_crude_pca$sdev^2)
  )
pooled_fraction_pca_scores <- as_tibble(
  pooled_fraction_pca$x[, 1:2, drop = FALSE], rownames = "profile_id"
) %>%
  left_join(pooled_fraction_meta, by = "profile_id") %>%
  mutate(
    state = "Pooled fractions",
    pc1_variance = pooled_fraction_pca$sdev[1]^2 / sum(pooled_fraction_pca$sdev^2),
    pc2_variance = pooled_fraction_pca$sdev[2]^2 / sum(pooled_fraction_pca$sdev^2)
  )

original_bray_path <- file.path(project_root, "data", "intermediate", "figure3-original-hgmd0330-bray-distance.csv")
if (file.exists(original_bray_path)) {
  # Prefer the original donor-level Bray-Curtis distance used to define the
  # four pools; fall back to Euclidean distance only if that cached matrix is
  # unavailable.
  original_bray <- readr::read_csv(original_bray_path, show_col_types = FALSE) %>%
    tibble::column_to_rownames("Sample") %>%
    as.matrix()
  included_donors <- sort(unique(crude_meta$profile_id))
  included_donors <- included_donors[included_donors %in% rownames(original_bray)]
  crude_dist <- original_bray[included_donors, included_donors, drop = FALSE]
  clustering_method <- "Original HGMD_0330 pellet-normalized Bray-Curtis donor distance + Ward.D2"
} else {
  crude_dist <- as.matrix(dist(crude_mat))
  clustering_method <- "Fallback Euclidean donor-profile distance + Ward.D2"
}
crude_hc <- hclust(as.dist(crude_dist), method = "ward.D2")
cluster_order <- tibble(profile_id = crude_hc$labels[crude_hc$order], cluster_order = seq_along(crude_hc$order)) %>%
  left_join(crude_meta, by = "profile_id")
crude_dist_long <- as.data.frame(as.table(crude_dist), stringsAsFactors = FALSE) %>%
  transmute(profile_1 = Var1, profile_2 = Var2, distance = Freq) %>%
  left_join(cluster_order %>% select(profile_id, order_1 = cluster_order), by = c("profile_1" = "profile_id")) %>%
  left_join(cluster_order %>% select(profile_id, order_2 = cluster_order), by = c("profile_2" = "profile_id"))

centroids <- rowsum(crude_mat, group = crude_meta$combined_pool_id[match(rownames(crude_mat), crude_meta$profile_id)]) /
  as.vector(table(crude_meta$combined_pool_id))
pool_dist <- as.matrix(dist(centroids))
balanced <- crude_meta %>% count(combined_pool_id) %>% filter(n >= 3) %>% pull(combined_pool_id)
pool_pairs <- as.data.frame(as.table(pool_dist), stringsAsFactors = FALSE) %>%
  transmute(pool_1 = Var1, pool_2 = Var2, distance = Freq) %>%
  filter(pool_1 < pool_2, pool_1 %in% balanced, pool_2 %in% balanced) %>%
  arrange(desc(distance))
contrast_pools <- pool_pairs %>% slice_head(n = 1)
pool_1 <- contrast_pools$pool_1[[1]]
pool_2 <- contrast_pools$pool_2[[1]]

fit_dataset_contrast <- function(df, state_name, pool_1, pool_2) {
  # Volcano contrasts compare the most separated pools. Fractionated data are
  # blocked by matched fraction sequence; crude data are modelled without a
  # fraction block.
  df %>%
    filter(combined_pool_id %in% c(pool_1, pool_2)) %>%
    group_by(analysis_batch) %>%
    group_modify(~ {
      wide <- .x %>%
        select(profile_id, analysis_feature, value) %>%
        group_by(profile_id, analysis_feature) %>%
        summarise(value = max(value), .groups = "drop") %>%
        pivot_wider(names_from = profile_id, values_from = value, values_fill = 0)
      feature_info <- .x %>% distinct(analysis_feature, annotation_id, compound.name, smiles, confidence.level)
      mat <- as.matrix(wide[, -1, drop = FALSE])
      rownames(mat) <- wide$analysis_feature
      meta <- .x %>% distinct(profile_id, combined_pool_id, block_id) %>% arrange(match(profile_id, colnames(mat)))
      mat <- mat[, meta$profile_id, drop = FALSE]
      keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
      mat <- mat[keep, , drop = FALSE]
      if (!nrow(mat)) return(tibble())

      group <- factor(meta$combined_pool_id, levels = c(pool_1, pool_2))
      if (state_name == "Fractionated") {
        block <- factor(meta$block_id)
        design <- model.matrix(~ block + group)
        coef_name <- tail(colnames(design), 1)
      } else {
        design <- model.matrix(~ group)
        coef_name <- "groupPool_4"
        if (!coef_name %in% colnames(design)) coef_name <- tail(colnames(design), 1)
      }
      fit <- lmFit(mat, design)
      estimable <- is.finite(fit$sigma) & fit$sigma > 0 & is.finite(fit$coefficients[, coef_name])
      if (!any(estimable)) return(tibble())
      fit <- eBayes(fit[estimable, ])
      topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
        as_tibble(rownames = "analysis_feature") %>%
        rename(effect_log10 = logFC, mean_value_log10 = AveExpr, p_value = P.Value, fdr = adj.P.Val) %>%
        left_join(feature_info, by = "analysis_feature") %>%
        mutate(
          effect = effect_log10 * log2(10),
          state = state_name,
          pool_1 = pool_1,
          pool_2 = pool_2
        )
    }) %>%
    ungroup() %>%
    mutate(
      direction = case_when(
        fdr < 0.05 & effect >= 1 ~ paste0("Higher in ", pool_2),
        fdr < 0.05 & effect <= -1 ~ paste0("Higher in ", pool_1),
        TRUE ~ "Not significant"
      )
    )
}

volcano <- bind_rows(
  fit_dataset_contrast(crude, "Crude", pool_1, pool_2),
  fit_dataset_contrast(fractionated, "Fractionated", pool_1, pool_2)
)
top_changes <- volcano %>%
  filter(direction != "Not significant") %>%
  group_by(state) %>%
  arrange(fdr, desc(abs(effect)), .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(label = coalesce(compound.name, annotation_id))

prevalence <- presence %>%
  # These legacy prevalence tables keep crude and fractionated states separate.
  distinct(state, combined_pool_id, annotation_id) %>%
  group_by(state, annotation_id) %>%
  mutate(pool_count = n_distinct(combined_pool_id)) %>%
  ungroup() %>%
  mutate(category = case_when(pool_count == 4 ~ "Common", pool_count == 1 ~ "Pool-unique", TRUE ~ "Shared")) %>%
  count(state, combined_pool_id, category, name = "annotations")

paper_cache_path <- file.path(project_root, "data", "processed", "paper2-npclassifier-cache.csv")
npc <- readr::read_csv(paper_cache_path, show_col_types = FALSE) %>%
  transmute(smiles = clean_text(SMILES), superclass = clean_text(np_superclass)) %>%
  filter(smiles != "", superclass != "") %>%
  distinct(smiles, .keep_all = TRUE)
class_breakdown <- presence %>%
  mutate(smiles = clean_text(smiles)) %>%
  left_join(npc, by = "smiles") %>%
  filter(!is.na(superclass), superclass != "") %>%
  distinct(state, combined_pool_id, annotation_id, superclass) %>%
  group_by(state, annotation_id) %>%
  mutate(pool_count = n_distinct(combined_pool_id)) %>%
  ungroup() %>%
  filter(pool_count %in% c(1, 4)) %>%
  mutate(category = if_else(pool_count == 4, "Common", "Pool-unique")) %>%
  count(state, combined_pool_id, category, superclass, name = "annotations") %>%
  group_by(superclass) %>%
  mutate(total = sum(annotations)) %>%
  ungroup() %>%
  mutate(superclass = if_else(dense_rank(desc(total)) <= 6, superclass, "Other")) %>%
  group_by(state, combined_pool_id, category, superclass) %>%
  summarise(annotations = sum(annotations), .groups = "drop")

combined_presence <- presence %>%
  # Current Figure 3e combines crude and fraction detections first, then counts
  # each annotation once per pool. Category is assigned by how many of the four
  # pools contain that annotation: 4 = common, 2-3 = shared, 1 = pool-unique.
  group_by(annotation_id, combined_pool_id) %>%
  summarise(
    confidence.level = min(confidence.level, na.rm = TRUE),
    compound.name = first(na.omit(compound.name)),
    smiles = first(na.omit(smiles)),
    .groups = "drop"
  ) %>%
  group_by(annotation_id) %>%
  mutate(
    pool_count = n_distinct(combined_pool_id),
    confidence.level = min(confidence.level, na.rm = TRUE),
    category = case_when(
      pool_count == 4 ~ "Common",
      pool_count == 1 ~ "Pool-unique",
      TRUE ~ "Shared"
    )
  ) %>%
  ungroup()

combined_prevalence <- combined_presence %>%
  count(confidence.level, combined_pool_id, category, name = "annotations") %>%
  arrange(confidence.level, combined_pool_id, category)

combined_prevalence_summary <- combined_presence %>%
  distinct(annotation_id, confidence.level, category, pool_count) %>%
  count(confidence.level, category, name = "annotations") %>%
  arrange(confidence.level, category)

combined_unique_by_pool <- combined_presence %>%
  filter(category == "Pool-unique") %>%
  count(confidence.level, combined_pool_id, name = "annotations") %>%
  arrange(confidence.level, combined_pool_id)

combined_common_classes <- combined_presence %>%
  # Figure 3f classifies only annotations common to all four pools. Missing or
  # unmatched NPClassifier calls are kept as "Unclassified" rather than dropped.
  filter(category == "Common") %>%
  group_by(annotation_id) %>%
  summarise(
    confidence.level = min(confidence.level, na.rm = TRUE),
    smiles = first(na.omit(clean_text(smiles))),
    .groups = "drop"
  ) %>%
  left_join(npc, by = "smiles") %>%
  mutate(superclass = if_else(is.na(superclass) | superclass == "", "Unclassified", superclass)) %>%
  arrange(annotation_id, confidence.level, superclass) %>%
  group_by(annotation_id, confidence.level) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(superclass) %>%
  summarise(annotations = n_distinct(annotation_id), .groups = "drop") %>%
  mutate(
    total = sum(annotations),
    fraction = annotations / total
  ) %>%
  arrange(desc(annotations))

combined_pool_unique_classes <- combined_presence %>%
  # Figure 3h classifies annotations detected in only one pool after crude and
  # fraction detections have first been combined. Each annotation contributes
  # once to its unique pool, avoiding double-counting across processing states.
  filter(category == "Pool-unique") %>%
  group_by(annotation_id, combined_pool_id) %>%
  summarise(
    confidence.level = min(confidence.level, na.rm = TRUE),
    smiles = first(na.omit(clean_text(smiles))),
    .groups = "drop"
  ) %>%
  left_join(npc, by = "smiles") %>%
  mutate(superclass = if_else(is.na(superclass) | superclass == "", "Unclassified", superclass)) %>%
  arrange(annotation_id, confidence.level, superclass) %>%
  group_by(annotation_id, combined_pool_id, confidence.level) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(superclass) %>%
  mutate(global_total = n_distinct(annotation_id)) %>%
  ungroup() %>%
  mutate(superclass = if_else(dense_rank(desc(global_total)) <= 8, superclass, "Other")) %>%
  group_by(combined_pool_id, superclass) %>%
  summarise(annotations = n_distinct(annotation_id), .groups = "drop") %>%
  group_by(combined_pool_id) %>%
  mutate(
    pool_total = sum(annotations),
    fraction = annotations / pool_total
  ) %>%
  ungroup() %>%
  arrange(combined_pool_id, desc(annotations))

set.seed(1)
permanova_crude <- vegan::adonis2(dist(crude_mat) ~ combined_pool_id, data = crude_meta, permutations = 999)
set.seed(1)
permanova_fraction <- vegan::adonis2(dist(fraction_mat) ~ combined_pool_id, data = fraction_meta, permutations = 999, strata = fraction_meta$block_id)
set.seed(1)
permanova_pooled_crude <- vegan::adonis2(
  dist(pooled_crude_mat) ~ combined_pool_id + processing_level,
  data = pooled_crude_meta,
  permutations = 999,
  strata = pooled_crude_meta$processing_level,
  by = "margin"
)
permanova <- bind_rows(
  as.data.frame(permanova_crude) %>% tibble::rownames_to_column("term") %>% mutate(state = "Crude"),
  as.data.frame(permanova_fraction) %>% tibble::rownames_to_column("term") %>% mutate(state = "Fractionated"),
  as.data.frame(permanova_pooled_crude) %>% tibble::rownames_to_column("term") %>% mutate(state = "Pooled crude")
)

write_csv_stable(pca_scores, file.path(project_root, "outputs", "tables", "figure3-pca-scores.csv"))
write_csv_stable(pooled_crude_pca_scores, file.path(project_root, "outputs", "tables", "figure3-pooled-crude-pca-scores.csv"))
write_csv_stable(individual_crude_pca_scores, file.path(project_root, "outputs", "tables", "figure3-individual-crude-pca-scores.csv"))
write_csv_stable(pooled_fraction_pca_scores, file.path(project_root, "outputs", "tables", "figure3-pooled-fraction-pca-scores.csv"))
write_csv_stable(cluster_order, file.path(project_root, "outputs", "tables", "figure3-crude-cluster-order.csv"))
write_csv_stable(crude_dist_long, file.path(project_root, "outputs", "tables", "figure3-crude-distance-matrix-long.csv"))
write_csv_stable(tibble(clustering_method = clustering_method), file.path(project_root, "outputs", "qc", "figure3-clustering-method.csv"))
write_csv_stable(pool_pairs, file.path(project_root, "outputs", "tables", "figure3-pool-centroid-distances.csv"))
write_csv_stable(volcano %>% filter(state == "Crude"), file.path(project_root, "outputs", "tables", "figure3-volcano-crude.csv"))
write_csv_stable(volcano %>% filter(state == "Fractionated"), file.path(project_root, "outputs", "tables", "figure3-volcano-fractionated.csv"))
write_csv_stable(top_changes, file.path(project_root, "outputs", "tables", "figure3-top-changing-metabolites.csv"))
write_csv_stable(prevalence, file.path(project_root, "outputs", "tables", "figure3-common-shared-pool-unique.csv"))
write_csv_stable(class_breakdown, file.path(project_root, "outputs", "tables", "figure3-common-unique-class-breakdown.csv"))
write_csv_stable(combined_prevalence, file.path(project_root, "outputs", "tables", "figure3-combined-pool-prevalence-by-confidence.csv"))
write_csv_stable(combined_prevalence_summary, file.path(project_root, "outputs", "tables", "figure3-combined-prevalence-summary.csv"))
write_csv_stable(combined_unique_by_pool, file.path(project_root, "outputs", "tables", "figure3-combined-pool-unique-by-confidence.csv"))
write_csv_stable(combined_common_classes, file.path(project_root, "outputs", "tables", "figure3-combined-common-class-breakdown.csv"))
write_csv_stable(combined_pool_unique_classes, file.path(project_root, "outputs", "tables", "figure3-combined-pool-unique-class-breakdown.csv"))
write_csv_stable(permanova, file.path(project_root, "outputs", "tables", "figure3-permanova.csv"))
