library(patchFPCA)

results_dir <- Sys.getenv("PATCHFPCA_RESULTS_DIR", unset = NA_character_)
if (is.na(results_dir)) {
  stop("Set PATCHFPCA_RESULTS_DIR before running.", call. = FALSE)
}

results_file <- file.path(results_dir, "patch_fpca_results.rds")
if (!file.exists(results_file)) {
  stop("Run 02_run_fpca.R first; results file is missing.", call. = FALSE)
}
results <- readRDS(results_file)

utils::write.csv(
  component_share_table(results),
  file.path(results_dir, "H1_1_component_shares.csv"),
  row.names = FALSE
)

utils::write.csv(
  check_loading_constancy(results, slope_thresh = 0.04, pc = 1L),
  file.path(results_dir, "H1_2_loading_constancy.csv"),
  row.names = FALSE
)

utils::write.csv(
  compare_mapping_invariance(results, pc = 1L),
  file.path(results_dir, "H1_3a_mapping_invariance.csv"),
  row.names = FALSE
)

quality_invariance <- rbind(
  compare_quality_invariance(results, mapping = "broad", pc = 1L),
  compare_quality_invariance(results, mapping = "narrow", pc = 1L)
)
utils::write.csv(
  quality_invariance,
  file.path(results_dir, "H1_3b_quality_invariance.csv"),
  row.names = FALSE
)

utils::write.csv(
  species_centroid_table(results, statistic = "mean", component = 1L),
  file.path(results_dir, "H2_species_centroids.csv"),
  row.names = FALSE
)

utils::write.csv(
  species_rank_correlations(results, statistic = "mean", component = 1L),
  file.path(results_dir, "H2_species_rank_correlations.csv")
)

utils::write.csv(
  rank_species_across_scenarios(
    results,
    statistic = "mean",
    component = 1L,
    digits = 1L
  ),
  file.path(results_dir, "H2_species_ranks.csv"),
  row.names = FALSE
)

