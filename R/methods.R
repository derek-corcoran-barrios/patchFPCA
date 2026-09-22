#' @export
print.patch_fpca <- function(x, ...) {
  cat("Ordered-distance patch FPCA\n")
  if (!is.null(x$scheme)) cat("  Scenario: ", x$scheme, "\n", sep = "")
  cat("  Profiles: ", x$n_profiles, " total; ", x$n_training,
      " training\n", sep = "")
  cat("  Neighbour ranks: ", x$k, "\n", sep = "")
  cat("  Centred: ", x$center, "; scaled: ", x$scale., "\n", sep = "")
  cat("  Component-share type: ", x$component_share_type, "\n", sep = "")
  cat("  PC1 share: ", format(round(x$component_share[[1L]], 4L)), "\n", sep = "")
  invisible(x)
}

#' @export
summary.patch_fpca <- function(object, ...) {
  data.frame(
    component = seq_along(object$component_share),
    eigenvalue = unname(object$eigenvalues),
    share = unname(object$component_share),
    cumulative_share = unname(object$cumulative_component_share),
    share_type = object$component_share_type,
    stringsAsFactors = FALSE
  )
}

#' @export
plot.patch_fpca <- function(
    x,
    which = c("scores", "loadings", "component_share"),
    components = 1:2,
    sample_fraction = 1,
    seed = x$seed,
    ...) {
  which <- match.arg(which)
  components <- as.integer(components)
  components <- components[components >= 1L & components <= ncol(x$loadings)]
  if (!length(components)) stop("No requested components are available.", call. = FALSE)

  if (which == "scores") {
    if (length(components) < 2L) {
      stop("Score plots require two available components.", call. = FALSE)
    }
    if (!is.numeric(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
      stop("sample_fraction must be in (0, 1].", call. = FALSE)
    }
    set.seed(seed)
    n_keep <- max(1L, floor(nrow(x$scores) * sample_fraction))
    rows <- if (n_keep < nrow(x$scores)) sample(seq_len(nrow(x$scores)), n_keep) else seq_len(nrow(x$scores))
    graphics::plot(
      x$scores[[components[[1L]]]][rows],
      x$scores[[components[[2L]]]][rows],
      xlab = paste0("PC", components[[1L]]),
      ylab = paste0("PC", components[[2L]]),
      main = if (is.null(x$scheme)) "Patch FPCA scores" else x$scheme,
      pch = 16,
      col = grDevices::adjustcolor("grey30", alpha.f = 0.35),
      ...
    )
  } else if (which == "loadings") {
    graphics::matplot(
      seq_len(nrow(x$loadings)),
      x$loadings[, components, drop = FALSE],
      type = "b",
      pch = seq_along(components),
      lty = seq_along(components),
      xlab = "Neighbour rank",
      ylab = "Loading",
      main = "FPCA loading functions",
      ...
    )
    graphics::abline(h = 0, lty = 2, col = "grey")
    graphics::legend(
      "topright",
      legend = paste0("PC", components),
      pch = seq_along(components),
      lty = seq_along(components),
      bty = "n"
    )
  } else {
    graphics::barplot(
      x$component_share,
      names.arg = paste0("PC", seq_along(x$component_share)),
      xlab = "Component",
      ylab = "Component share",
      main = x$component_share_type,
      ...
    )
  }
  invisible(x)
}

#' @export
print.patch_fpca_collection <- function(x, ...) {
  cat("patchFPCA scenario collection\n")
  cat("  Scenarios: ", paste(names(x), collapse = ", "), "\n", sep = "")
  settings <- attr(x, "settings")
  if (!is.null(settings)) {
    cat("  k: ", settings$k_max, "; centred: ", settings$center,
        "; scaled: ", settings$scale., "\n", sep = "")
  }
  invisible(x)
}

