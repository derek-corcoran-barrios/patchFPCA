# patchFPCA 0.1.2

- Read the observed and expected Moran's I from DHARMa's `statistic` field;
  reject malformed test results rather than failing after residual simulation.
- Disable knitr's lazy-load cache for the bundled spatial report because
  multi-million-row spatial objects can exceed its serialization limits.

# patchFPCA 0.1.1

- Added a bookdown spatial-robustness report that replaces linear coordinate
  terms with selected dbMEM variables, compares non-spatial and SPDE species
  mixed models, evaluates patch spatial support and mesh resolution, and checks
  grouped DHARMa residual spatial autocorrelation.
- Added `render_spatial_report()` so the bundled report can be rendered from an
  installed package without relying on source-repository `inst/` paths.

# patchFPCA 0.1.0

- Converted the five supplied manuscript scripts into path-safe package
  functions.
- Preserved the manuscript's uncentred and unscaled default FPCA.
- Added explicit uncentred-inertia terminology and a centred sensitivity
  comparison.
- Added tests, manual pages, workflow scripts, and archived source scripts.
