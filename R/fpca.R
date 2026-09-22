.validate_profiles <- function(profiles, min_rows = 2L) {
  if (is.data.frame(profiles)) profiles <- as.matrix(profiles)
  if (!is.matrix(profiles) || !is.numeric(profiles)) {
    stop("profiles must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(profiles) < min_rows || ncol(profiles) < 1L) {
    stop(
      "profiles must contain at least ", min_rows,
      " row(s) and one column.",
      call. = FALSE
    )
  }
  if (any(!is.finite(profiles))) {
    stop("profiles must contain only finite values.", call. = FALSE)
  }
  profiles
}

.balance_profiles <- function(profiles, species, seed) {
  groups <- split(
    seq_len(nrow(profiles)),
    factor(species, levels = unique(species)),
    drop = TRUE
  )
  minimum_n <- min(lengths(groups))
  if (minimum_n < 1L) {
    stop("Every species must have at least one profile.", call. = FALSE)
  }
  set.seed(seed)
  sampled_rows <- unlist(
    lapply(groups, function(index) sample(index, minimum_n, replace = FALSE)),
    use.names = FALSE
  )
  list(
    profiles = profiles[sampled_rows, , drop = FALSE],
    species = species[sampled_rows],
    rows = sampled_rows,
    n_per_species = minimum_n
  )
}

.transform_profiles <- function(profiles, log_transform, profiles_are_logged) {
  if (log_transform && !profiles_are_logged) {
    if (any(profiles < 0)) {
      stop("Distances must be non-negative before log1p transformation.", call. = FALSE)
    }
    return(log1p(profiles))
  }
  profiles
}

#' Fit FPCA to ordered nearest-neighbour distance profiles
#'
#' This implementation treats each row as a discretized function over ordered
#' neighbour rank and uses `stats::prcomp()` on a species-balanced training
#' matrix. All supplied profiles are then projected using the learned loadings.
#'
#' @param profiles Numeric profile matrix; rows are patches and columns are
#'   ordered neighbour ranks.
#' @param species Species identifier for every row. If `NULL`, all rows form
#'   one group.
#' @param log_transform Apply `log1p()` before fitting.
#' @param profiles_are_logged Whether `profiles` have already been transformed.
#' @param balance Downsample every species to the smallest species sample size.
#' @param seed Sampling seed used during balancing.
#' @param center Passed to `stats::prcomp()`. `FALSE` reproduces the supplied
#'   manuscript scripts; `TRUE` performs conventional column-centred PCA.
#' @param scale. Passed to `stats::prcomp()`.
#' @param rank. Optional maximum number of components.
#' @return An object of class `patch_fpca`.
#' @export
fit_patch_fpca <- function(
    profiles,
    species = NULL,
    log_transform = TRUE,
    profiles_are_logged = FALSE,
    balance = TRUE,
    seed = 37L,
    center = FALSE,
    scale. = FALSE,
    rank. = NULL) {
  profiles <- .validate_profiles(profiles)
  if (is.null(species)) species <- rep("all", nrow(profiles))
  if (length(species) != nrow(profiles) || anyNA(species)) {
    stop("species must contain one non-missing value per profile.", call. = FALSE)
  }
  species <- as.character(species)
  if (length(seed) != 1L || is.na(seed)) {
    stop("seed must be one non-missing number.", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(scale.) || length(scale.) != 1L || is.na(scale.)) {
    stop("scale. must be TRUE or FALSE.", call. = FALSE)
  }

  transformed <- .transform_profiles(
    profiles,
    log_transform = log_transform,
    profiles_are_logged = profiles_are_logged
  )

  if (balance) {
    balanced <- .balance_profiles(transformed, species, seed = seed)
  } else {
    balanced <- list(
      profiles = transformed,
      species = species,
      rows = seq_len(nrow(transformed)),
      n_per_species = NA_integer_
    )
  }

  prcomp_arguments <- list(
    x = balanced$profiles,
    center = center,
    scale. = scale.
  )
  if (!is.null(rank.)) {
    prcomp_arguments$rank. <- .assert_positive_integer(rank., "rank.")
  }
  pca_res <- do.call(stats::prcomp, prcomp_arguments)
  scores_matrix <- stats::predict(pca_res, newdata = transformed)
  scores <- as.data.frame(scores_matrix, stringsAsFactors = FALSE)
  names(scores) <- paste0("PC", seq_len(ncol(scores)))

  eigenvalues <- pca_res$sdev^2
  total_eigenvalue <- sum(eigenvalues)
  if (!is.finite(total_eigenvalue) || total_eigenvalue <= 0) {
    stop("The training profiles contain no decomposable variation or magnitude.", call. = FALSE)
  }
  component_share <- eigenvalues / total_eigenvalue
  component_share_type <- if (center) "centred_variance" else "uncentred_inertia"

  result <- list(
    call = match.call(),
    pca_res = pca_res,
    loadings = pca_res$rotation,
    scores = scores,
    species = species,
    training_species = balanced$species,
    training_rows = balanced$rows,
    n_per_species = balanced$n_per_species,
    eigenvalues = eigenvalues,
    component_share = component_share,
    cumulative_component_share = cumsum(component_share),
    component_share_type = component_share_type,
    pve = component_share,
    cum_pve = cumsum(component_share),
    center = center,
    scale. = scale.,
    log_transform = log_transform,
    profiles_are_logged = profiles_are_logged,
    balanced = balance,
    seed = seed,
    n_profiles = nrow(transformed),
    n_training = nrow(balanced$profiles),
    k = ncol(transformed)
  )
  class(result) <- "patch_fpca"
  result
}

#' Project new distance profiles into a fitted FPCA
#'
#' @param object A `patch_fpca` object.
#' @param profiles New numeric profiles in the same neighbour-rank order.
#' @param profiles_are_logged Whether the supplied new profiles have already
#'   received the transformation used by the fit.
#' @return A data frame of component scores named `FPCA1`, `FPCA2`, and so on.
#' @export
project_patch_fpca <- function(object, profiles, profiles_are_logged = FALSE) {
  if (!inherits(object, "patch_fpca")) {
    stop("object must inherit from 'patch_fpca'.", call. = FALSE)
  }
  profiles <- .validate_profiles(profiles, min_rows = 1L)
  if (ncol(profiles) != object$k) {
    stop(
      "New profiles have ", ncol(profiles), " columns; the fit expects ",
      object$k, ".",
      call. = FALSE
    )
  }
  transformed <- .transform_profiles(
    profiles,
    log_transform = object$log_transform,
    profiles_are_logged = profiles_are_logged
  )
  scores <- as.data.frame(
    stats::predict(object$pca_res, newdata = transformed),
    stringsAsFactors = FALSE
  )
  names(scores) <- paste0("FPCA", seq_len(ncol(scores)))
  scores
}

#' Clean species names encoded in source filenames
#'
#' @param x Character vector such as
#'   `Sorted_aglais_urticae_broad_High_distances`.
#' @return Cleaned binomial names.
#' @export
clean_species_name <- function(x) {
  clean_one <- function(value) {
    value <- sub("^Sorted_", "", as.character(value))
    value <- sub("_(broad|narrow)_(High|Medium|Low).*$", "", value)
    value <- gsub("_", " ", value, fixed = TRUE)
    words <- strsplit(trimws(value), "[[:space:]]+")[[1L]]
    if (!length(words) || !nzchar(words[[1L]])) return(NA_character_)
    words[[1L]] <- paste0(
      toupper(substr(words[[1L]], 1L, 1L)),
      tolower(substr(words[[1L]], 2L, nchar(words[[1L]])))
    )
    if (length(words) > 1L) words[-1L] <- tolower(words[-1L])
    paste(words, collapse = " ")
  }
  vapply(x, clean_one, character(1L), USE.NAMES = FALSE)
}

.species_summaries <- function(scores, species) {
  species_levels <- unique(species)
  rows <- lapply(species_levels, function(species_name) {
    index <- which(species == species_name)
    data.frame(
      species = species_name,
      PC1_mean = mean(scores$PC1[index]),
      PC2_mean = if ("PC2" %in% names(scores)) mean(scores$PC2[index]) else NA_real_,
      PC1_median = stats::median(scores$PC1[index]),
      PC2_median = if ("PC2" %in% names(scores)) stats::median(scores$PC2[index]) else NA_real_,
      n_profiles = length(index),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.project_reordered_data <- function(object, data_rds, scheme_name, k_max) {
  token <- paste0("_", scheme_name, "_")
  source_names <- names(data_rds)[grepl(token, names(data_rds), fixed = TRUE)]
  patch_results <- list()

  for (source_name in source_names) {
    df <- data_rds[[source_name]]
    if (!is.data.frame(df) || !nrow(df)) next
    profiles <- extract_distances(df, k_max = k_max)
    if (!nrow(profiles)) next
    retained_rows <- attr(profiles, "row_index")
    score_data <- project_patch_fpca(object, profiles, profiles_are_logged = FALSE)
    patch_id <- if ("final_id_1" %in% names(df)) {
      df$final_id_1[retained_rows]
    } else {
      retained_rows
    }
    patch_results[[source_name]] <- data.frame(
      final_id_1 = patch_id,
      score_data,
      species_full = source_name,
      scheme = scheme_name,
      species = clean_species_name(source_name),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (!length(patch_results)) return(NULL)
  result <- do.call(rbind, patch_results)
  rownames(result) <- NULL
  result
}

#' Run one manuscript FPCA scenario
#'
#' @param scheme_name Name of one habitat-definition scenario.
#' @param profiles_by_scheme Named list of profile matrices.
#' @param species_id_by_scheme Named list of corresponding source/species IDs.
#' @param data_rds Optional reordered data list used to build a patch-level
#'   projection table.
#' @param k_max Number of ordered neighbour ranks.
#' @param thin_frac Retained for compatibility with the supplied function;
#'   plotting subsampling is now handled by `plot()`.
#' @param seed,center,scale. Passed to `fit_patch_fpca()`.
#' @param profiles_are_logged Whether the matrices are already log-transformed.
#' @param log_transform Whether log transformation defines this analysis.
#' @return A `patch_fpca` object augmented with scenario and species summaries.
#' @export
run_fpca_analysis <- function(
    scheme_name,
    profiles_by_scheme,
    species_id_by_scheme,
    data_rds = NULL,
    k_max = 5L,
    thin_frac = 0.05,
    seed = 37L,
    center = FALSE,
    scale. = FALSE,
    profiles_are_logged = TRUE,
    log_transform = TRUE) {
  if (!scheme_name %in% names(profiles_by_scheme)) {
    stop("scheme_name is absent from profiles_by_scheme: ", scheme_name, call. = FALSE)
  }
  if (!scheme_name %in% names(species_id_by_scheme)) {
    stop("scheme_name is absent from species_id_by_scheme: ", scheme_name, call. = FALSE)
  }
  if (!is.numeric(thin_frac) || length(thin_frac) != 1L ||
      is.na(thin_frac) || thin_frac <= 0 || thin_frac > 1) {
    stop("thin_frac must be in (0, 1].", call. = FALSE)
  }

  fit <- fit_patch_fpca(
    profiles = profiles_by_scheme[[scheme_name]],
    species = species_id_by_scheme[[scheme_name]],
    log_transform = log_transform,
    profiles_are_logged = profiles_are_logged,
    balance = TRUE,
    seed = seed,
    center = center,
    scale. = scale.
  )

  cleaned_species <- clean_species_name(species_id_by_scheme[[scheme_name]])
  species_scores <- data.frame(
    species = cleaned_species,
    PC1 = fit$scores$PC1,
    PC2 = if ("PC2" %in% names(fit$scores)) fit$scores$PC2 else NA_real_,
    stringsAsFactors = FALSE
  )
  species_centroids <- .species_summaries(fit$scores, cleaned_species)

  fit$scheme <- scheme_name
  fit$thin_frac <- thin_frac
  fit$phi1 <- fit$loadings[, 1L]
  fit$phi2 <- if (ncol(fit$loadings) >= 2L) fit$loadings[, 2L] else NULL
  fit$species_scores <- species_scores
  fit$species_centroids <- species_centroids
  fit$species_rank_mean <- species_centroids[order(species_centroids$PC1_mean), ]
  fit$species_rank_median <- species_centroids[order(species_centroids$PC1_median), ]
  fit$fpca_patch_df <- if (is.null(data_rds)) NULL else {
    .project_reordered_data(fit, data_rds, scheme_name, k_max)
  }
  class(fit) <- c("patch_fpca_scenario", "patch_fpca")
  fit
}

#' Legacy wrapper for the supplied manuscript function name
#'
#' @param scheme_name,Y_log_list,species_id_by_scheme,data_rds,k_max,thin_frac,seed
#'   Arguments used by the original `run_FPCA_analysis()` function.
#' @param center,scale. PCA settings.
#' @return A `patch_fpca` scenario result.
#' @export
run_FPCA_0_analysis <- function(
    scheme_name,
    Y_log_list,
    species_id_by_scheme,
    data_rds = NULL,
    k_max = 5L,
    thin_frac = 0.05,
    seed = 37L,
    center = FALSE,
    scale. = FALSE) {
  run_fpca_analysis(
    scheme_name = scheme_name,
    profiles_by_scheme = Y_log_list,
    species_id_by_scheme = species_id_by_scheme,
    data_rds = data_rds,
    k_max = k_max,
    thin_frac = thin_frac,
    seed = seed,
    center = center,
    scale. = scale.,
    profiles_are_logged = TRUE,
    log_transform = TRUE
  )
}

#' Run every requested FPCA scenario
#'
#' @param data_rds Named list of reordered distance data frames.
#' @param schemes Scenarios to run.
#' @param k_max,thin_frac,seed,center,scale. Analysis settings.
#' @param log_transform Apply `log1p()` to raw distance profiles.
#' @return A named list of scenario results with class
#'   `patch_fpca_collection`.
#' @export
run_all_fpca <- function(
    data_rds,
    schemes = default_fpca_schemes(),
    k_max = 5L,
    thin_frac = 0.05,
    seed = 37L,
    center = FALSE,
    scale. = FALSE,
    log_transform = TRUE) {
  pooled <- pool_distance_profiles(
    data_rds = data_rds,
    schemes = schemes,
    k_max = k_max
  )
  available_schemes <- intersect(schemes, names(pooled$Y_by_scheme))
  missing_schemes <- setdiff(schemes, available_schemes)
  if (length(missing_schemes)) {
    warning(
      "No complete profiles were found for: ",
      paste(missing_schemes, collapse = ", "),
      call. = FALSE
    )
  }

  results <- lapply(available_schemes, function(scheme_name) {
    run_fpca_analysis(
      scheme_name = scheme_name,
      profiles_by_scheme = pooled$Y_by_scheme,
      species_id_by_scheme = pooled$species_id_by_scheme,
      data_rds = data_rds,
      k_max = k_max,
      thin_frac = thin_frac,
      seed = seed,
      center = center,
      scale. = scale.,
      profiles_are_logged = FALSE,
      log_transform = log_transform
    )
  })
  names(results) <- available_schemes
  class(results) <- c("patch_fpca_collection", "list")
  attr(results, "settings") <- list(
    schemes = available_schemes,
    k_max = k_max,
    seed = seed,
    center = center,
    scale. = scale.,
    log_transform = log_transform
  )
  results
}

#' Tabulate component shares across scenarios
#'
#' @param results Named list returned by `run_all_fpca()`.
#' @return A long data frame of component shares and cumulative shares.
#' @export
component_share_table <- function(results) {
  if (!is.list(results) || is.null(names(results))) {
    stop("results must be a named list.", call. = FALSE)
  }
  rows <- lapply(names(results), function(scenario) {
    result <- results[[scenario]]
    data.frame(
      scenario = scenario,
      component = seq_along(result$component_share),
      share = unname(result$component_share),
      cumulative_share = unname(result$cumulative_component_share),
      share_type = result$component_share_type,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Save manuscript FPCA tables
#'
#' @param results Named list returned by `run_all_fpca()`.
#' @param output_dir Destination directory.
#' @return Invisibly, a character vector of written paths.
#' @export
save_fpca_results <- function(results, output_dir) {
  if (!is.list(results) || is.null(names(results))) {
    stop("results must be a named list.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)

  share_file <- file.path(output_dir, "component_shares.csv")
  utils::write.csv(component_share_table(results), share_file, row.names = FALSE)
  written <- c(written, share_file)

  for (scenario in names(results)) {
    result <- results[[scenario]]
    loadings_file <- file.path(output_dir, paste0(scenario, "_loadings.csv"))
    centroid_file <- file.path(output_dir, paste0(scenario, "_species_centroids.csv"))
    utils::write.csv(
      data.frame(rank = seq_len(nrow(result$loadings)), result$loadings),
      loadings_file,
      row.names = FALSE
    )
    utils::write.csv(result$species_centroids, centroid_file, row.names = FALSE)
    written <- c(written, loadings_file, centroid_file)

    if (!is.null(result$fpca_patch_df)) {
      patch_file <- file.path(output_dir, paste0(scenario, "_fpca_patch_df.csv"))
      utils::write.csv(result$fpca_patch_df, patch_file, row.names = FALSE)
      written <- c(written, patch_file)
    }
  }

  results_file <- file.path(output_dir, "patch_fpca_results.rds")
  saveRDS(results, results_file)
  written <- c(written, results_file)
  invisible(normalizePath(written, winslash = "/", mustWork = TRUE))
}
