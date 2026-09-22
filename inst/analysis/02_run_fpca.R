library(patchFPCA)

reordered_dir <- Sys.getenv("PATCHFPCA_REORDERED_DIR", unset = NA_character_)
results_dir <- Sys.getenv("PATCHFPCA_RESULTS_DIR", unset = NA_character_)

if (is.na(reordered_dir) || is.na(results_dir)) {
  stop(
    "Set PATCHFPCA_REORDERED_DIR and PATCHFPCA_RESULTS_DIR before running.",
    call. = FALSE
  )
}

reordered <- read_reordered_data(reordered_dir)

# Exact settings used by the supplied manuscript scripts.
results <- run_all_fpca(
  data_rds = reordered$data,
  k_max = 5L,
  seed = 37L,
  center = FALSE,
  scale. = FALSE,
  log_transform = TRUE
)

save_fpca_results(results, results_dir)

# Recommended pre-submission sensitivity comparison for every scenario.
pooled <- pool_distance_profiles(reordered$data, k_max = 5L)
centering_sensitivity <- lapply(names(pooled$Y_by_scheme), function(scenario) {
  comparison <- compare_centering(
    profiles = pooled$Y_by_scheme[[scenario]],
    species = pooled$species_id_by_scheme[[scenario]],
    log_transform = TRUE,
    seed = 37L,
    scale. = FALSE
  )
  data.frame(scenario = scenario, comparison$summary, check.names = FALSE)
})
centering_sensitivity <- do.call(rbind, centering_sensitivity)

utils::write.csv(
  centering_sensitivity,
  file.path(results_dir, "centering_sensitivity.csv"),
  row.names = FALSE
)

