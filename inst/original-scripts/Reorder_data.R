library(dplyr)
library(tidyr)
library(tools)  # for file_path_sans_ext

data_dir <- "C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/Habitats_from_sustainscapes/Habitat_estimates/Raw_data/"
Out_dir  <- "C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/Habitats_from_sustainscapes/Habitat_estimates/Reordered_data/"

rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

for (rds in rds_files) {
  
  df <- readRDS(rds)
  
  # Count number of points per final_id_1
  df_filtered <- df %>%
    group_by(final_id_1) %>%
    filter(n() >= 5) %>%  # keep only rows with 5 or more distances
    ungroup()
  
  # Skip file if nothing left after filtering
  if (nrow(df_filtered) == 0) next
  
  base_name <- file_path_sans_ext(basename(rds))
  
  # Distance matrix (wide)
  distances_wide <- df_filtered %>%
    group_by(final_id_1) %>%
    arrange(distance_m, .by_group = TRUE) %>%
    mutate(col_id = row_number()) %>%
    select(final_id_1, col_id, distance_m) %>%
    pivot_wider(
      names_from = col_id,
      values_from = distance_m,
      names_prefix = "distance_"
    ) %>%
    ungroup()
  
  # Metadata matrix (wide)
  metadata_wide <- df_filtered %>%
    group_by(final_id_1) %>%
    arrange(distance_m, .by_group = TRUE) %>%
    mutate(col_id = row_number()) %>%
    select(final_id_1, col_id, final_id_2) %>%
    pivot_wider(
      names_from = col_id,
      values_from = final_id_2,
      names_prefix = "id_"
    ) %>%
    ungroup()
  
  saveRDS(
    metadata_wide,
    file = file.path(Out_dir, paste0("Meta_data_Sorted_", base_name, ".rds"))
  )
  
  saveRDS(
    distances_wide,
    file = file.path(Out_dir, paste0("Sorted_", base_name, ".rds"))
  )
}
