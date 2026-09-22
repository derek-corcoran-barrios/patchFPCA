#' Default habitat-definition scenarios
#'
#' @return A character vector containing the six scenarios used in the
#'   manuscript.
#' @export
default_fpca_schemes <- function() {
  c(
    "broad_High", "broad_Medium", "broad_Low",
    "narrow_High", "narrow_Medium", "narrow_Low"
  )
}

.assert_columns <- function(data, columns, object_name = "data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      object_name, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

.assert_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

#' Convert neighbour records to ordered wide profiles
#'
#' Each row of `data` represents the distance between a focal patch and one
#' neighbouring patch. Records are grouped by focal patch, ordered by distance,
#' and widened to `distance_1`, ..., `distance_n` and `id_1`, ..., `id_n`.
#'
#' @param data A data frame of long-form neighbour records.
#' @param id_col Name of the focal-patch identifier column.
#' @param neighbor_col Name of the neighbouring-patch identifier column.
#' @param distance_col Name of the numeric distance column.
#' @param min_neighbors Minimum number of records required for a focal patch.
#' @return A list with `distances` and `metadata` data frames.
#' @export
reorder_distance_data <- function(
    data,
    id_col = "final_id_1",
    neighbor_col = "final_id_2",
    distance_col = "distance_m",
    min_neighbors = 5L) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  .assert_columns(data, c(id_col, neighbor_col, distance_col))
  min_neighbors <- .assert_positive_integer(min_neighbors, "min_neighbors")

  if (!is.numeric(data[[distance_col]])) {
    stop(distance_col, " must be numeric.", call. = FALSE)
  }
  if (anyNA(data[[id_col]])) {
    stop(id_col, " must not contain missing values.", call. = FALSE)
  }

  groups <- split(
    seq_len(nrow(data)),
    as.character(data[[id_col]]),
    drop = TRUE
  )
  groups <- groups[lengths(groups) >= min_neighbors]

  empty_result <- function(prefix, prototype) {
    out <- data.frame(row.names = integer(0))
    out[[id_col]] <- data[[id_col]][FALSE]
    out[[paste0(prefix, "1")]] <- prototype[FALSE]
    out
  }

  if (!length(groups)) {
    return(list(
      distances = empty_result("distance_", data[[distance_col]]),
      metadata = empty_result("id_", data[[neighbor_col]])
    ))
  }

  ordered_groups <- lapply(groups, function(index) {
    index[order(data[[distance_col]][index], na.last = TRUE)]
  })
  first_rows <- vapply(ordered_groups, function(index) index[[1L]], integer(1L))
  max_neighbors <- max(lengths(ordered_groups))

  distances <- data.frame(.row = seq_along(ordered_groups))
  distances[[id_col]] <- data[[id_col]][first_rows]
  distances$.row <- NULL

  metadata <- data.frame(.row = seq_along(ordered_groups))
  metadata[[id_col]] <- data[[id_col]][first_rows]
  metadata$.row <- NULL

  distance_na <- data[[distance_col]][NA_integer_]
  neighbor_na <- data[[neighbor_col]][NA_integer_]

  for (rank in seq_len(max_neighbors)) {
    distance_values <- lapply(ordered_groups, function(index) {
      if (length(index) >= rank) data[[distance_col]][index[[rank]]] else distance_na
    })
    neighbor_values <- lapply(ordered_groups, function(index) {
      if (length(index) >= rank) data[[neighbor_col]][index[[rank]]] else neighbor_na
    })

    distances[[paste0("distance_", rank)]] <- do.call(c, distance_values)
    metadata[[paste0("id_", rank)]] <- do.call(c, neighbor_values)
  }

  rownames(distances) <- NULL
  rownames(metadata) <- NULL
  list(distances = distances, metadata = metadata)
}

#' Reorder one raw distance RDS file
#'
#' @inheritParams reorder_distance_data
#' @param file Path to one raw RDS file.
#' @param output_dir Directory in which to write the two reordered RDS files.
#' @return A one-row manifest describing the written files.
#' @export
reorder_distance_file <- function(
    file,
    output_dir,
    id_col = "final_id_1",
    neighbor_col = "final_id_2",
    distance_col = "distance_m",
    min_neighbors = 5L) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  reordered <- reorder_distance_data(
    readRDS(file),
    id_col = id_col,
    neighbor_col = neighbor_col,
    distance_col = distance_col,
    min_neighbors = min_neighbors
  )

  base_name <- tools::file_path_sans_ext(basename(file))
  distance_file <- file.path(output_dir, paste0("Sorted_", base_name, ".rds"))
  metadata_file <- file.path(
    output_dir,
    paste0("Meta_data_Sorted_", base_name, ".rds")
  )

  saveRDS(reordered$distances, distance_file)
  saveRDS(reordered$metadata, metadata_file)

  data.frame(
    input_file = normalizePath(file, winslash = "/", mustWork = TRUE),
    distance_file = normalizePath(distance_file, winslash = "/", mustWork = TRUE),
    metadata_file = normalizePath(metadata_file, winslash = "/", mustWork = TRUE),
    n_profiles = nrow(reordered$distances),
    stringsAsFactors = FALSE
  )
}

#' Reorder every raw distance RDS file in a directory
#'
#' @inheritParams reorder_distance_file
#' @param input_dir Directory containing raw RDS files.
#' @param pattern Regular expression identifying input files.
#' @return A manifest with one row per input file.
#' @export
reorder_distance_files <- function(
    input_dir,
    output_dir,
    pattern = "\\.rds$",
    id_col = "final_id_1",
    neighbor_col = "final_id_2",
    distance_col = "distance_m",
    min_neighbors = 5L) {
  if (!dir.exists(input_dir)) {
    stop("input_dir does not exist: ", input_dir, call. = FALSE)
  }
  files <- sort(list.files(input_dir, pattern = pattern, full.names = TRUE))
  if (!length(files)) {
    stop("No input RDS files matched in: ", input_dir, call. = FALSE)
  }

  manifests <- lapply(files, function(file) {
    reorder_distance_file(
      file = file,
      output_dir = output_dir,
      id_col = id_col,
      neighbor_col = neighbor_col,
      distance_col = distance_col,
      min_neighbors = min_neighbors
    )
  })
  do.call(rbind, manifests)
}

#' Read reordered distance and metadata RDS files
#'
#' @param data_dir Directory containing reordered RDS files.
#' @return A list containing `data`, `metadata`, and `all` named lists.
#' @export
read_reordered_data <- function(data_dir) {
  if (!dir.exists(data_dir)) {
    stop("data_dir does not exist: ", data_dir, call. = FALSE)
  }
  files <- sort(list.files(data_dir, pattern = "\\.rds$", full.names = TRUE))
  if (!length(files)) {
    stop("No RDS files found in: ", data_dir, call. = FALSE)
  }
  object_names <- tools::file_path_sans_ext(basename(files))
  objects <- stats::setNames(lapply(files, readRDS), object_names)
  is_metadata <- grepl("^Meta_data", names(objects))

  list(
    data = objects[!is_metadata],
    metadata = objects[is_metadata],
    all = objects
  )
}

#' Extract complete ordered distance profiles
#'
#' @param df A reordered distance data frame.
#' @param k_max Number of neighbour ranks to retain.
#' @return A numeric matrix containing complete profiles. Its `row_index`
#'   attribute records the retained rows in `df`.
#' @export
extract_distances <- function(df, k_max = 5L) {
  if (!is.data.frame(df)) {
    stop("df must be a data frame.", call. = FALSE)
  }
  k_max <- .assert_positive_integer(k_max, "k_max")
  distance_columns <- paste0("distance_", seq_len(k_max))
  .assert_columns(df, distance_columns, "df")
  if (!all(vapply(df[distance_columns], is.numeric, logical(1L)))) {
    stop("All selected distance columns must be numeric.", call. = FALSE)
  }

  matrix_data <- as.matrix(df[distance_columns])
  keep <- stats::complete.cases(matrix_data)
  matrix_data <- matrix_data[keep, , drop = FALSE]
  colnames(matrix_data) <- distance_columns
  attr(matrix_data, "row_index") <- which(keep)
  matrix_data
}

#' Pool distance profiles across species within scenarios
#'
#' @param data_rds Named list of reordered distance data frames.
#' @param schemes Scenario names to pool.
#' @param k_max Number of ordered neighbours to retain.
#' @return A list with `Y_by_scheme`, `species_id_by_scheme`, and
#'   `source_row_by_scheme`.
#' @export
pool_distance_profiles <- function(
    data_rds,
    schemes = default_fpca_schemes(),
    k_max = 5L) {
  if (!is.list(data_rds) || is.null(names(data_rds))) {
    stop("data_rds must be a named list of data frames.", call. = FALSE)
  }
  k_max <- .assert_positive_integer(k_max, "k_max")

  Y_by_scheme <- list()
  species_id_by_scheme <- list()
  source_row_by_scheme <- list()

  for (scheme in schemes) {
    token <- paste0("_", scheme, "_")
    source_names <- names(data_rds)[grepl(token, names(data_rds), fixed = TRUE)]
    matrices <- list()
    species_ids <- character(0)
    source_rows <- list()

    for (source_name in source_names) {
      df <- data_rds[[source_name]]
      if (!is.data.frame(df) || !nrow(df)) next
      matrix_data <- extract_distances(df, k_max = k_max)
      if (!nrow(matrix_data)) next

      matrices[[source_name]] <- matrix_data
      species_ids <- c(species_ids, rep(source_name, nrow(matrix_data)))
      source_rows[[source_name]] <- data.frame(
        source = source_name,
        row = attr(matrix_data, "row_index"),
        stringsAsFactors = FALSE
      )
    }

    if (!length(matrices)) next
    Y_by_scheme[[scheme]] <- do.call(rbind, matrices)
    rownames(Y_by_scheme[[scheme]]) <- NULL
    species_id_by_scheme[[scheme]] <- species_ids
    source_row_by_scheme[[scheme]] <- do.call(rbind, source_rows)
    rownames(source_row_by_scheme[[scheme]]) <- NULL
  }

  if (!length(Y_by_scheme)) {
    stop("No complete distance profiles matched the requested schemes.", call. = FALSE)
  }

  list(
    Y_by_scheme = Y_by_scheme,
    species_id_by_scheme = species_id_by_scheme,
    source_row_by_scheme = source_row_by_scheme
  )
}

