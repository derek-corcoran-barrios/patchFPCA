# patchFPCA

`patchFPCA` contains the reproducible R workflow used for the manuscript
*Species-specific habitat definitions alter inferred habitat connectivity
despite stable drivers of variance in spatial patch geometry*.

It converts raw nearest-neighbour records to ordered distance profiles,
pools and balances profiles across species, fits the six habitat-definition
scenarios, projects every eligible patch into component space, and implements
the manuscript's loading, invariance, and species-rank diagnostics.

## Install

install the package with:

```r
remotes::install_github("derek-corcoran-barrios/"patchFPCA")
```

The package itself uses base R only. `testthat` is needed only to run its
automated tests.

## Manuscript workflow

```r
library(patchFPCA)

# 1. Convert each raw RDS file to a sorted distance-profile RDS and a
#    corresponding neighbour-ID metadata RDS.
reorder_distance_files(
  input_dir = "path/to/Raw_data",
  output_dir = "path/to/Reordered_data",
  min_neighbors = 5
)

# 2. Read the reordered data and run all six manuscript scenarios.
reordered <- read_reordered_data("path/to/Reordered_data")

results <- run_all_fpca(
  data_rds = reordered$data,
  schemes = default_fpca_schemes(),
  k_max = 5,
  seed = 37,
  center = FALSE,
  scale. = FALSE
)

# 3. Save tabular outputs used downstream in the manuscript.
save_fpca_results(results, "path/to/Results")

# 4. Manuscript diagnostics.
component_share_table(results)
check_loading_constancy(results, slope_thresh = 0.04)
compare_mapping_invariance(results)
compare_quality_invariance(results, mapping = "broad")
species_rank_correlations(results)
```

Complete path-safe scripts are in `inst/analysis/`. The five scripts supplied
by Nathalie are retained unchanged in `inst/original-scripts/` for provenance.

## Important centering note

The supplied manuscript code uses `prcomp(..., center = FALSE, scale. = FALSE)`.
That choice is preserved as the package default so results remain
reproducible. It is an **uncentred PCA/SVD**: component shares describe total
squared magnitude about the origin (uncentred inertia), not variance about
column means.

Conventional `center = TRUE` subtracts the mean at each neighbour rank across
patches. It does **not** row-centre each individual profile, and therefore does
not automatically remove all biologically meaningful between-patch spacing.
Run the pre-submission sensitivity check with:

```r
pooled <- pool_distance_profiles(reordered$data, k_max = 5)

sensitivity <- compare_centering(
  profiles = pooled$Y_by_scheme$broad_High,
  species = pooled$species_id_by_scheme$broad_High,
  log_transform = TRUE,
  seed = 37
)

sensitivity$summary
```

For uncentred fits, `pve` is retained as a legacy field because the original
scripts use that name. New code should use `component_share` and inspect
`component_share_type`.

## Mapping from supplied scripts

| Supplied script | Package replacement |
|---|---|
| `Reorder_data.R` | `reorder_distance_data()` and `reorder_distance_files()` |
| `FPCA_function.R` | `fit_patch_fpca()`, `run_fpca_analysis()`, and legacy wrapper `run_FPCA_analysis()` |
| `Run_FPCA_Function.R` | `run_all_fpca()` and `save_fpca_results()` |
| `FPCA_on_distance_rds_No_function_oneliber_script.R` | explicit low-level functions plus `plot()` methods |
| `Hypothesis_testing.R` | diagnostic and invariance functions documented below |

The raw habitat, patch-statistics, genomic, and demographic datasets referenced
by the original scripts are not included. Data-specific downstream RDA/GLM
sections therefore remain in the archived original script rather than being
executed when the package is loaded.
