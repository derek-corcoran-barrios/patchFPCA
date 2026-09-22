# Reproducible manuscript analysis

Run the scripts in numeric order after setting these environment variables:

- `PATCHFPCA_RAW_DIR`: raw per-species distance RDS files.
- `PATCHFPCA_REORDERED_DIR`: destination/source for reordered RDS files.
- `PATCHFPCA_RESULTS_DIR`: destination for model outputs and diagnostics.

The package defaults reproduce the supplied scripts: five ordered neighbours,
`log1p()` transformation, species-balanced training samples, seed 37, and
`center = FALSE`, `scale. = FALSE`.

`02_run_fpca.R` also includes an optional loop for the conventional centred
sensitivity analysis. This should be run before submission because uncentred
component shares are shares of origin-based inertia, whereas centred component
shares are variance shares around column means.

## Spatial-autocorrelation report

`04_spatial_autocorrelation.Rmd` implements the manuscript's spatial robustness
analysis. It:

- rebuilds each final 0 m, 10 m, or 20 m patch from the original shapefile and
  the matching `PolygonToFinal*` lookup;
- joins corrected patch area and patch-level FPCA1 by `final_id`;
- replaces linear `x + y + x:y` terms in variation partitioning with selected
  dbMEM variables;
- fits the same species mixed model without and with an sdmTMB Matern SPDE
  spatial field;
- checks mesh resolution and point-on-surface versus centroid support; and
- compares grouped DHARMa Moran tests before and after spatial modelling.

The report defaults to the 20 m correction. Change `merge_variant` to `zero` or
`10m` only if `02_run_fpca.R` was run from the corresponding corrected distance
files. The report deliberately stops when the patch IDs and selected correction
do not match.

Keep the 0 m, 10 m, and 20 m reordered distance files in separate input
directories. Do not pass a directory containing all three variants to
`read_reordered_data()`, because files from different patch definitions share
the same six scenario tokens and would otherwise be pooled together.

If the package is installed, render from any working directory with the public
wrapper. `project_root` should be the directory containing `Results/` and
`Species_PatchDistances/` (or the parent of whatever paths you supply):

```r
patchFPCA::render_spatial_report(
  project_root = "/path/to/your/analysis-project",
  fpca_results_dir = "Results",
  patch_dir = "Species_PatchDistances",
  spatial_output_dir = "Results/spatial_autocorrelation",
  merge_variant = "20m"
)
```

The wrapper locates the files with `system.file()`, copies the Rmd, helper, and
bibliography to a writable `spatial-report/` directory, and renders the PDF
there. In an installed package, do not use the source-tree path
`inst/analysis/04_spatial_autocorrelation.Rmd` directly.

Package developers working from a source checkout can instead render from the
package root:

```r
rmarkdown::render(
  "inst/analysis/04_spatial_autocorrelation.Rmd",
  params = list(
    project_root = ".",
    fpca_results_dir = "Results",
    patch_dir = "Species_PatchDistances",
    spatial_output_dir = "Results/spatial_autocorrelation",
    merge_variant = "20m"
  ),
  output_file = "patchFPCA_spatial_autocorrelation.pdf"
)
```

The report writes model objects, tables, a run manifest, and the assembled
model data to `spatial_output_dir`. Generated PDF/TeX/cache files and result
folders should remain outside version control.
