.component_vector <- function(result, pc) {
  pc <- .assert_positive_integer(pc, "pc")
  if (is.null(result$loadings) || ncol(result$loadings) < pc) {
    stop("Requested component is absent from a result.", call. = FALSE)
  }
  result$loadings[, pc]
}

#' Test whether component loadings are approximately constant over rank
#'
#' @param results_list Named list of FPCA results.
#' @param slope_thresh Absolute slope below which a loading vector is labelled
#'   practically constant.
#' @param pc Component number.
#' @return A data frame with slope, p-value, R-squared, coefficient of
#'   variation, range, and the practical-constancy label.
#' @export
check_loading_constancy <- function(results_list, slope_thresh = 0.04, pc = 1L) {
  if (!is.list(results_list) || is.null(names(results_list))) {
    stop("results_list must be a named list.", call. = FALSE)
  }
  rows <- lapply(names(results_list), function(scenario) {
    phi <- .component_vector(results_list[[scenario]], pc)
    rank <- seq_along(phi)
    model <- stats::lm(phi ~ rank)
    model_summary <- summary(model)
    slope <- unname(stats::coef(model)[[2L]])
    coefficient_of_variation <- stats::sd(phi) / abs(mean(phi))
    data.frame(
      scheme = scenario,
      component = pc,
      slope = slope,
      p_value = model_summary$coefficients[2L, 4L],
      r_squared = model_summary$r.squared,
      cv = coefficient_of_variation,
      range = diff(range(phi)),
      practically_constant = abs(slope) < slope_thresh,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Contrast local and regional loading weights
#'
#' @param v Numeric loading vector.
#' @param local Indices assigned to the local group.
#' @param regional Indices assigned to the regional group.
#' @return Named numeric vector containing the two means and their difference.
#' @export
local_regional_test <- function(v, local = 1:2, regional = 3:length(v)) {
  if (!is.numeric(v) || any(!is.finite(v))) {
    stop("v must be a finite numeric vector.", call. = FALSE)
  }
  if (!length(local) || !length(regional) ||
      any(c(local, regional) < 1L) || any(c(local, regional) > length(v))) {
    stop("local and regional must identify valid elements of v.", call. = FALSE)
  }
  local_mean <- mean(v[local])
  regional_mean <- mean(v[regional])
  c(L = local_mean, R = regional_mean, diff = local_mean - regional_mean)
}

#' Local-regional contrast excluding the median rank
#'
#' @param v Numeric loading vector of length at least five.
#' @return Named numeric vector comparing ranks 1-2 with ranks 4-5.
#' @export
local_regional_test_nomedian <- function(v) {
  if (length(v) < 5L) {
    stop("v must contain at least five ranks.", call. = FALSE)
  }
  local_regional_test(v, local = 1:2, regional = 4:5)
}

#' Sign-aligned cosine similarity
#'
#' Principal-component loading signs are indeterminate. By default, the second
#' vector is sign-aligned to the first before similarity is calculated.
#'
#' @param v1,v2 Numeric vectors of equal length.
#' @param align_sign Align the arbitrary component signs before comparison.
#' @return Cosine similarity in `[0, 1]` when sign alignment is requested.
#' @export
cosine_similarity <- function(v1, v2, align_sign = TRUE) {
  if (!is.numeric(v1) || !is.numeric(v2) || length(v1) != length(v2) ||
      any(!is.finite(v1)) || any(!is.finite(v2))) {
    stop("v1 and v2 must be equal-length finite numeric vectors.", call. = FALSE)
  }
  denominator <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))
  if (denominator == 0) stop("Cosine similarity is undefined for zero vectors.", call. = FALSE)
  if (align_sign && sum(v1 * v2) < 0) v2 <- -v2
  similarity <- sum(v1 * v2) / denominator
  max(-1, min(1, similarity))
}

#' Compare broad and narrow mappings within quality levels
#'
#' @param all_results Named list of scenario results.
#' @param pc Component number.
#' @param qualities Quality suffixes present in scenario names.
#' @return A data frame of sign-aligned cosine similarities and angles.
#' @export
compare_mapping_invariance <- function(
    all_results,
    pc = 1L,
    qualities = c("High", "Medium", "Low")) {
  rows <- lapply(qualities, function(quality) {
    broad_name <- paste0("broad_", quality)
    narrow_name <- paste0("narrow_", quality)
    if (!all(c(broad_name, narrow_name) %in% names(all_results))) {
      stop("Missing result pair for quality: ", quality, call. = FALSE)
    }
    similarity <- cosine_similarity(
      .component_vector(all_results[[broad_name]], pc),
      .component_vector(all_results[[narrow_name]], pc)
    )
    data.frame(
      quality = quality,
      cosine_similarity = similarity,
      angle_degrees = acos(similarity) * 180 / pi,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Compare quality levels within one mapping definition
#'
#' @param all_results Named list of scenario results.
#' @param mapping Mapping prefix, normally `"broad"` or `"narrow"`.
#' @param pc Component number.
#' @param qualities Quality suffixes present in scenario names.
#' @return Pairwise sign-aligned cosine similarities and angles.
#' @export
compare_quality_invariance <- function(
    all_results,
    mapping,
    pc = 1L,
    qualities = c("High", "Medium", "Low")) {
  combinations <- utils::combn(qualities, 2L)
  rows <- lapply(seq_len(ncol(combinations)), function(index) {
    q1 <- combinations[1L, index]
    q2 <- combinations[2L, index]
    name1 <- paste0(mapping, "_", q1)
    name2 <- paste0(mapping, "_", q2)
    if (!all(c(name1, name2) %in% names(all_results))) {
      stop("Missing result pair: ", name1, " and ", name2, call. = FALSE)
    }
    similarity <- cosine_similarity(
      .component_vector(all_results[[name1]], pc),
      .component_vector(all_results[[name2]], pc)
    )
    data.frame(
      mapping = mapping,
      q1 = q1,
      q2 = q2,
      cosine = similarity,
      angle = acos(similarity) * 180 / pi,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Mapping-by-quality interaction distance for loading vectors
#'
#' @param all_results Named list of scenario results.
#' @param q1,q2 Two quality suffixes.
#' @param pc Component number.
#' @param align_sign Align all loading vectors to `broad_q1` first.
#' @return Euclidean norm of the difference-in-differences loading vector.
#' @export
interaction_distance <- function(
    all_results,
    q1,
    q2,
    pc = 1L,
    align_sign = TRUE) {
  names_needed <- c(
    paste0("broad_", q1), paste0("broad_", q2),
    paste0("narrow_", q1), paste0("narrow_", q2)
  )
  if (!all(names_needed %in% names(all_results))) {
    stop("Missing one or more required scenarios.", call. = FALSE)
  }
  vectors <- lapply(all_results[names_needed], .component_vector, pc = pc)
  if (align_sign) {
    reference <- vectors[[1L]]
    vectors <- lapply(vectors, function(vector) {
      if (sum(reference * vector) < 0) -vector else vector
    })
  }
  delta <- (vectors[[1L]] - vectors[[2L]]) -
    (vectors[[3L]] - vectors[[4L]])
  sqrt(sum(delta^2))
}

#' Assemble species centroids across scenarios
#'
#' @param results Named list of scenario results.
#' @param statistic Use component means or medians.
#' @param component Component number.
#' @return Wide data frame with one row per species.
#' @export
species_centroid_table <- function(
    results,
    statistic = c("mean", "median"),
    component = 1L) {
  statistic <- match.arg(statistic)
  component <- .assert_positive_integer(component, "component")
  value_column <- paste0("PC", component, "_", statistic)
  scenario_tables <- lapply(names(results), function(scenario) {
    centroids <- results[[scenario]]$species_centroids
    if (is.null(centroids) || !value_column %in% names(centroids)) {
      stop("Missing species-centroid column: ", value_column, call. = FALSE)
    }
    output <- centroids[c("species", value_column)]
    names(output)[[2L]] <- scenario
    output
  })
  Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), scenario_tables)
}

#' Correlate species positions across FPCA scenarios
#'
#' @param results Named list of scenario results.
#' @param statistic Use component means or medians.
#' @param component Component number.
#' @param method Correlation method passed to `stats::cor()`.
#' @return A scenario-by-scenario correlation matrix.
#' @export
species_rank_correlations <- function(
    results,
    statistic = c("mean", "median"),
    component = 1L,
    method = "spearman") {
  table <- species_centroid_table(results, statistic, component)
  stats::cor(table[-1L], method = method, use = "pairwise.complete.obs")
}

#' Rank species across FPCA scenarios
#'
#' @param results Named list of scenario results.
#' @param statistic Use component means or medians.
#' @param component Component number.
#' @param digits Decimal places used before ranking, reproducing the original
#'   manuscript script's one-decimal ranking when `digits = 1`.
#' @param decreasing Rank larger scores first.
#' @return Wide data frame of integer ranks.
#' @export
rank_species_across_scenarios <- function(
    results,
    statistic = c("mean", "median"),
    component = 1L,
    digits = 1L,
    decreasing = TRUE) {
  values <- species_centroid_table(results, statistic, component)
  ranks <- values
  for (column in names(values)[-1L]) {
    x <- round(values[[column]], digits = digits)
    if (decreasing) x <- -x
    ranks[[column]] <- rank(x, ties.method = "min", na.last = "keep")
  }
  ranks
}

#' Compare centred and uncentred FPCA fits
#'
#' @param profiles Numeric ordered-distance profiles.
#' @param species Species identifier per profile.
#' @param log_transform Apply `log1p()`.
#' @param balance Balance species sample sizes.
#' @param seed Sampling seed.
#' @param scale. Scale columns in both fits.
#' @return A list containing both fits and a one-row comparison summary.
#' @export
compare_centering <- function(
    profiles,
    species,
    log_transform = TRUE,
    balance = TRUE,
    seed = 37L,
    scale. = FALSE) {
  uncentred <- fit_patch_fpca(
    profiles, species,
    log_transform = log_transform,
    balance = balance,
    seed = seed,
    center = FALSE,
    scale. = scale.
  )
  centred <- fit_patch_fpca(
    profiles, species,
    log_transform = log_transform,
    balance = balance,
    seed = seed,
    center = TRUE,
    scale. = scale.
  )

  centred_loading <- centred$loadings[, 1L]
  sign_multiplier <- if (sum(uncentred$loadings[, 1L] * centred_loading) < 0) -1 else 1
  centred_loading <- centred_loading * sign_multiplier
  centred_scores <- centred$scores$PC1 * sign_multiplier

  summary_table <- data.frame(
    uncentred_pc1_share = uncentred$component_share[[1L]],
    centred_pc1_share = centred$component_share[[1L]],
    loading_cosine = cosine_similarity(
      uncentred$loadings[, 1L],
      centred_loading,
      align_sign = FALSE
    ),
    score_pearson = stats::cor(uncentred$scores$PC1, centred_scores),
    score_spearman = stats::cor(
      uncentred$scores$PC1,
      centred_scores,
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
  list(uncentred = uncentred, centred = centred, summary = summary_table)
}

