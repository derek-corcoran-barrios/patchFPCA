#' Render the spatial-autocorrelation robustness report
#'
#' Locates the report bundled with an installed `patchFPCA` package, copies the
#' R Markdown source, helper functions, and bibliography to a writable project
#' directory, and renders the bookdown PDF there. This avoids relying on a
#' source-repository path such as `inst/analysis/...`.
#'
#' @param project_root Existing project directory containing the analysis
#'   inputs, or to which relative input paths refer.
#' @param fpca_results_dir Directory containing `patch_fpca_results.rds`,
#'   relative to `project_root` unless absolute.
#' @param patch_dir Directory containing the patch shapefiles and matching
#'   lookup/area tables, relative to `project_root` unless absolute.
#' @param spatial_output_dir Directory for fitted models and result tables,
#'   relative to `project_root` unless absolute.
#' @param merge_variant Patch correction used to generate the FPCA inputs.
#' @param report_dir Writable directory into which the report sources, cache,
#'   and PDF are copied or generated, relative to `project_root` unless
#'   absolute. When `NULL`, uses `spatial-report/<merge_variant>`.
#' @param output_file PDF filename. Supply a filename, not a directory path.
#'   When `NULL`, includes `merge_variant` in the filename.
#' @param quiet Passed to [rmarkdown::render()].
#' @return Invisibly, the normalized path to the rendered PDF.
#' @export
render_spatial_report <- function(
    project_root = ".",
    fpca_results_dir = "Results",
    patch_dir = "Species_PatchDistances",
    spatial_output_dir = "Results/spatial_autocorrelation",
    merge_variant = c("20m", "10m", "zero"),
    report_dir = NULL,
    output_file = NULL,
    quiet = FALSE) {
  merge_variant <- match.arg(merge_variant)
  if (is.null(report_dir)) {
    report_dir <- file.path("spatial-report", merge_variant)
  }
  if (is.null(output_file)) {
    output_file <- paste0(
      "patchFPCA_spatial_autocorrelation_",
      merge_variant,
      ".pdf"
    )
  }

  required_packages <- c("rmarkdown", "bookdown")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop(
      "Install the report-rendering packages first: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  resolve_project_path <- function(path) {
    path_is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", path)
    normalizePath(
      if (path_is_absolute) path else file.path(project_root, path),
      winslash = "/",
      mustWork = FALSE
    )
  }

  resolved_results_dir <- resolve_project_path(fpca_results_dir)
  expected_results_file <- file.path(
    resolved_results_dir,
    "patch_fpca_results.rds"
  )
  if (!file.exists(expected_results_file)) {
    stop(
      "Cannot render the spatial report because this file is missing: ",
      expected_results_file, ". ",
      "Set fpca_results_dir to the directory containing patch_fpca_results.rds, ",
      "or generate it with patchFPCA::save_fpca_results().",
      call. = FALSE
    )
  }

  resolved_patch_dir <- resolve_project_path(patch_dir)
  if (!dir.exists(resolved_patch_dir)) {
    stop(
      "Cannot render the spatial report because patch_dir does not exist: ",
      resolved_patch_dir,
      call. = FALSE
    )
  }

  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", report_dir)
  report_dir <- if (is_absolute) report_dir else file.path(project_root, report_dir)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  report_dir <- normalizePath(report_dir, winslash = "/", mustWork = TRUE)

  if (basename(output_file) != output_file) {
    stop("output_file must be a filename; use report_dir to choose its directory.", call. = FALSE)
  }

  source_dir <- system.file("analysis", package = "patchFPCA")
  report_files <- c(
    "04_spatial_autocorrelation.Rmd",
    "spatial_analysis_helpers.R",
    "spatial_references.bib"
  )
  source_files <- file.path(source_dir, report_files)
  if (!nzchar(source_dir) || any(!file.exists(source_files))) {
    stop(
      "The installed patchFPCA package does not contain the spatial report. ",
      "Install patchFPCA version 0.1.1 or later.",
      call. = FALSE
    )
  }

  copied <- file.copy(source_files, report_dir, overwrite = TRUE)
  if (!all(copied)) {
    stop(
      "Could not copy the spatial report files to the writable report directory: ",
      report_dir,
      call. = FALSE
    )
  }

  rendered <- rmarkdown::render(
    input = file.path(report_dir, "04_spatial_autocorrelation.Rmd"),
    output_file = output_file,
    output_dir = report_dir,
    params = list(
      project_root = project_root,
      fpca_results_dir = fpca_results_dir,
      patch_dir = patch_dir,
      spatial_output_dir = spatial_output_dir,
      merge_variant = merge_variant
    ),
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )

  invisible(normalizePath(rendered, winslash = "/", mustWork = TRUE))
}
