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

