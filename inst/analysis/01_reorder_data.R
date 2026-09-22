library(patchFPCA)

raw_dir <- Sys.getenv("PATCHFPCA_RAW_DIR", unset = NA_character_)
reordered_dir <- Sys.getenv("PATCHFPCA_REORDERED_DIR", unset = NA_character_)

if (is.na(raw_dir) || is.na(reordered_dir)) {
  stop(
    "Set PATCHFPCA_RAW_DIR and PATCHFPCA_REORDERED_DIR before running.",
    call. = FALSE
  )
}

manifest <- reorder_distance_files(
  input_dir = raw_dir,
  output_dir = reordered_dir,
  min_neighbors = 5L
)

print(manifest)

