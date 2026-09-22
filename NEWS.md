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
