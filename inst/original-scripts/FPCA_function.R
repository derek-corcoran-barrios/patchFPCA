#Needed packages for running the:  FPCA analysis function 

library(dplyr)   
library(ggplot2)
library(patchwork)
library(tidyr)
library(terra)  # working with spatial data 
library(ggplot2)
library(sf)      # for spatial data 
library(purrr)   # for functional programming (loops)
library(future) 
library(future.apply)  
library(progressr)
library(tidyr)
library(stringr)


##Worker for data prep downstream: 

extract_distances <- function(df, k_max) {
  mat <- as.matrix(df[, paste0("distance_", 1:k_max)])
  mat[complete.cases(mat), ]
}


## ________________________________________________________
## A function of running:
## Functional PCA of nearest-neighbour distance profiles
## Implemented as PCA on discrete functional data with: 
## center = FALSE, scale. = FALSE
## ________________________________________________________


###
### The below function returns a list of dataframes containing many estimates, 
# plot objects and addtional data frames (:D):   


#pca_res = pca residuals
#loadings = PCA loadings extracted  
#scores = all patch scores on all PCs
#phi1 = FPCA1 loadings,
#phi2 = FPCA2 loadings, 
#p_loadings = FPCA1 & FPCA2 loadings 
#pve = Variance explained
#cum_pve = cumulative Variance explained
#p_Varexplained = plot of Variance explained,
#landscape_plot = The subsampled FPCA backround.  
#species_scores = FPCA1 & FPCA2 score for each patch pr. species
#species_centroids = mean of the FPCA1 & FPCA2 score for each species
#combined_plot = combined plot of the overlap in FPCA and species by species plots 
#centroid_plot = plot of FPCA1 & FPCA2 score for each patch pr. species,
#species_rank = species rank o  FPCA1 (the mean),
#fpca_patch_df = a data frame pr species containing patch IDs, FPCA1,FPCA2, FPCA3, FPCA4,FPCA5, species_full (OG data-files name), scheme, species




## FPCA analysis function 

run_FPCA_analysis <- function(scheme_name, Y_log_list, species_id_by_scheme, data_rds = NULL, k_max = 5, thin_frac = 0.05, seed = 37) {
  
  balance_by_species <- function(Y, species_vec, seed = 42) {
    set.seed(seed)
    
    df <- as.data.frame(Y)
    df$species <- species_vec
    
    n_min <- min(table(species_vec))
    
    df_bal <- df %>%
      group_by(species) %>%
      slice_sample(n = n_min) %>%
      ungroup()
    
    # Only numeric columns go into PCA
    Y_bal <- as.matrix(df_bal %>% select(where(is.numeric)))
    species_bal <- df_bal$species
    
    return(list(Y_bal = Y_bal, species_bal = species_bal))
  }
  
  # ---- 1. PCA ----
  Y_log <- Y_log_list[[scheme_name]]
  species_vec <- species_id_by_scheme[[scheme_name]]
  
  # ---- NEW: balance dataset ----
  bal <- balance_by_species(Y_log, species_vec, seed = seed)
  
  Y_bal <- bal$Y_bal
  
  # ---- PCA on balanced data ----
  pca_res <- prcomp(Y_bal, center = FALSE, scale. = FALSE)
  
  scores_all <- predict(pca_res, newdata = Y_log)
  scores <- as.data.frame(scores_all)
  loadings <- pca_res$rotation
  
  # ---- NEW: Patch-level FPCA projection ----
  fpca_patch_df <- NULL
  
  if (!is.null(data_rds)) {
    
    compute_patch_fpca <- function(df, loadings, k_max) {
      
      dist_cols <- paste0("distance_", 1:k_max)
      
      # Ensure columns exist
      dist_cols <- dist_cols[dist_cols %in% colnames(df)]
      
      mat <- as.matrix(df[, dist_cols, drop = FALSE])
      
      keep <- complete.cases(mat)
      mat <- mat[keep, , drop = FALSE]
      
      # Log transform (match PCA input)
      mat_log <- log1p(mat)
      
      # Project
      scores <- mat_log %*% loadings[1:k_max, , drop = FALSE]
      
      scores_df <- as.data.frame(scores)
      colnames(scores_df) <- paste0("FPCA", seq_len(ncol(scores_df)))
      
      out <- df[keep, c("final_id_1"), drop = FALSE]
      out <- cbind(out, scores_df)
      
      return(out)
    }
    
    pattern <- paste0("_", scheme_name, "_")
    scheme_names <- grep(pattern, names(data_rds), value = TRUE)
    
    patch_list <- list()
    
    for (nm in scheme_names) {
      
      df <- data_rds[[nm]]
      
      if (nrow(df) == 0) next
      
      res <- compute_patch_fpca(df, loadings, k_max)
      
      res$species_full <- nm
      res$scheme <- scheme_name
      
      patch_list[[nm]] <- res
    }
    
    fpca_patch_df <- bind_rows(patch_list)
    
    # Optional: clean species names
    fpca_patch_df <- fpca_patch_df %>%
      mutate(
        species = str_extract(species_full, "(?<=Sorted_)[^_]+_[^_]+"),
        species = str_replace(species, "_", " "),
        species = paste0(str_to_title(word(species,1)), " ", word(species,2))
      )
  }
  
  # ---- 2. Eigenfunctions ----
  k <- 1:ncol(Y_log)
  phi1 <- loadings[,1]
  phi2 <- loadings[,2]
  
  # Save eigenfunction plots
  par(mfrow = c(1,2))
  plot(k, phi1, type = "b", pch = 19, ylim = c(-1,1), xlab = "Neighbor rank (k)", 
       ylab = paste0("Loading_", scheme_name), main = "FPCA1")
  abline(h = 0, lty = 2, col = "grey")
  
  plot(k, phi2, type = "b", pch = 19, ylim = c(-1,1), xlab = "Neighbor rank (k)",
       ylab = paste0("Loading_", scheme_name), main = "FPCA2")
  abline(h = 0, lty = 2, col = "grey")
  p_loadings <- recordPlot()
  par(mfrow = c(1,1))
  
  # ---- 3. Variance explained ----
  eigvals <- pca_res$sdev^2
  pve <- eigvals / sum(eigvals)
  cum_pve <- cumsum(pve)
  
  par(mfrow = c(1,2))
  plot(pve, type = "b", pch = 19, xlab = "Component",
       ylab = paste("Proportion of variance explained of:", scheme_name),
       main = paste("FPCA variance explained:", scheme_name))
  
  barplot(cum_pve, xlab = "Component",
          ylab = paste("Cumulative variance explained of:", scheme_name),
          main = paste("Cumulative FPCA variance explained:", scheme_name))
  abline(h = 0.9, lty = 2, col = "grey")
  p_Varexplained <- recordPlot()
  par(mfrow = c(1,1))
  
  # ---- 4. FPCA landscape ----
  set.seed(seed)
  scores_thinned <- scores %>% sample_frac(thin_frac)
  
  landscape_plot <- ggplot(scores_thinned, aes(x = PC1, y = PC2)) +
    geom_point(color = "lightgrey") +
    theme_bw() +
    labs(title = paste("FPCA scores (", scheme_name, ") - thinned", 100*thin_frac, "%", sep=""),
         x = "FPCA 1", y = "FPCA 2") +
    xlim(min(scores_thinned$PC1), max(scores_thinned$PC1)) +
    ylim(min(scores_thinned$PC2), max(scores_thinned$PC2))
  
  # ---- 5. Species projections ----
  
  
  species_scores <- data.frame(
    species = species_id_by_scheme[[scheme_name]],
    PC1 = scores[,1],
    PC2 = scores[,2]
  )
  
  species_scores <- species_scores %>%
    mutate(
      species = str_extract(species, "(?<=Sorted_)[^_]+_[^_]+"),
      species = str_replace(species, "_", " "),
      species = paste0(str_to_title(word(species,1)), " ", word(species,2))
    )
  
  species_centroids <- species_scores %>%
    group_by(species) %>%
    summarise(
      PC1_mean = mean(PC1),
      PC2_mean = mean(PC2),
      PC1_median = median(PC1),
      PC2_median = median(PC2),
      .groups = "drop"
    )
  
  species_scores$species <- factor(species_scores$species)
  species_centroids$species <- factor(species_centroids$species,
                                      levels = levels(species_scores$species))
  
  
  set.seed(seed)
  
  species_scores_thinned <- species_scores %>%
    group_by(species) %>%
    sample_frac(thin_frac) %>%
    ungroup()
  
  species_levels <- levels(species_scores$species)
  cols <- setNames(rainbow(length(species_levels)), species_levels)
  cols_alpha <- paste0(cols, "80")  # 50% opacity for patches
  
  # Top row: All species
  top_plot <- ggplot() +
    geom_point(data = species_scores_thinned, aes(x = PC1, y = PC2, color = species),
               size = 0.5, alpha = 0.5) +
    geom_point(data = species_centroids, aes(x = PC1_mean, y = PC2_mean, fill = species),
               color = "black", pch = 21, size = 3) +
    scale_color_manual(values = c(
      "Coenonympha tullia" = "#dfc27d",
      "Lycaena virgaureae" = "#bf812d",
      "Thymelicus lineola" = "#8c510a",
      "Maniola jurtina" = "darkgrey",                        
      "Pieris napi" = "#74a9cf",
      "Gonepteryx rhamni" = "#0570b0",
      "Aglais urticae" = "#034e7b"
    )) +
    scale_fill_manual(values = c(
      "Coenonympha tullia" = "#dfc27d",
      "Lycaena virgaureae" = "#bf812d",
      "Thymelicus lineola" = "#8c510a",
      "Maniola jurtina" = "darkgrey",                        
      "Pieris napi" = "#74a9cf",
      "Gonepteryx rhamni" = "#0570b0",
      "Aglais urticae" = "#034e7b"
    )) +
    theme_bw() +
    labs(title = "All species", x = "FPCA1", y = "FPCA2")
  

  
  # Desired layout: 3 columns, 3-1-3 species per column
  species_order <- c(
    "Coenonympha tullia", "Lycaena virgaureae", "Thymelicus lineola", # column 1
    NA, # empty to push Maniola jurtina to second column
    "Maniola jurtina", 
    NA, NA, # empty to push last three species to column 3
    "Pieris napi", "Gonepteryx rhamni", "Aglais urticae"
  )
  
  # Add a temporary factor with NA placeholders to force facet_wrap layout
  species_scores_thinned <- species_scores_thinned %>%
    mutate(species_facet = factor(species, levels = species_order))
  
  species_centroids <- species_centroids %>%
    mutate(species_facet = factor(species, levels = species_order))
  
  ##OBS NAT Write as params!! *********!!!!
  
  # Define your species groups
  group1 <- c("Coenonympha tullia", "Lycaena virgaureae", "Thymelicus lineola")
  group2 <- c("Maniola jurtina")
  group3 <- c("Pieris napi", "Gonepteryx rhamni", "Aglais urticae")
  
  # Determine the overall PC1/PC2 limits to keep plots aligned
  xlims <- range(species_scores_thinned$PC1, na.rm = TRUE)
  ylims <- range(species_scores_thinned$PC2, na.rm = TRUE)
  
  # Function to create per-group plot with fixed limits
  plot_group <- function(species_subset, show_legend = FALSE) {
    ggplot(filter(species_scores_thinned, species %in% species_subset), aes(x = PC1, y = PC2)) +
      geom_point(aes(color = species), size = 0.5, alpha = 0.5) +
      geom_point(data = filter(species_centroids, species %in% species_subset),
                 aes(x = PC1_mean, y = PC2_mean, fill = species),
                 color = "black", pch = 21, size = 3) +
      scale_color_manual(values = c(
        "Coenonympha tullia" = "#dfc27d",
        "Lycaena virgaureae" = "#bf812d",
        "Thymelicus lineola" = "#8c510a",
        "Maniola jurtina" = "darkgrey",                        
        "Pieris napi" = "#74a9cf",
        "Gonepteryx rhamni" = "#0570b0",
        "Aglais urticae" = "#034e7b"
      )) +
      scale_fill_manual(values = c(
        "Coenonympha tullia" = "#dfc27d", 
        "Lycaena virgaureae" = "#bf812d",
        "Thymelicus lineola" = "#8c510a",
        "Maniola jurtina" = "darkgrey",                        
        "Pieris napi" = "#74a9cf",
        "Gonepteryx rhamni" = "#0570b0",
        "Aglais urticae" = "#034e7b"
      )) +
      facet_wrap(~species, ncol = 1, scales = "fixed") +
      coord_cartesian(xlim = xlims, ylim = ylims) +  # enforce equal plotting area
      theme_bw() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(angle = 0),
        legend.position = ifelse(show_legend, "right", "none")
      )
  }
  
  # Create the three separate plots
  plot1 <- plot_group(group1, show_legend = FALSE)

  # Middle plot: add dummy rows to make 3 facets
  middle_data <- species_scores_thinned %>%
    filter(species == "Maniola jurtina") %>%
    bind_rows(tibble(
      PC1 = NA, PC2 = NA, species = c("dummy_top", "dummy_bottom")
    ))
  
  middle_centroids <- species_centroids %>%
    filter(species == "Maniola jurtina") %>%
    bind_rows(tibble(
      PC1_mean = NA, PC2_mean = NA, species = c("dummy_top", "dummy_bottom")
    ))
  
  plot2 <- ggplot(middle_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = species), size = 0.5, alpha = 0.5, na.rm = TRUE) +
    geom_point(data = middle_centroids,
               aes(x = PC1_mean, y = PC2_mean, fill = species),
               color = "black", pch = 21, size = 3, na.rm = TRUE) +
    scale_color_manual(values = c(
      "Maniola jurtina" = "darkgrey",
      "dummy_top" = "white",
      "dummy_bottom" = "white"
    )) +
    scale_fill_manual(values = c(
      "Maniola jurtina" = "darkgrey",
      "dummy_top" = "white",
      "dummy_bottom" = "white"
    )) +
    facet_wrap(~species, ncol = 1) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(angle = 0),
      legend.position = "none"
    )
  
  
  plot3 <- plot_group(group3, show_legend = FALSE)
  
  # Combine side by side with patchwork
  bottom_plot <- plot1 + plot2 + plot3 + plot_layout(ncol = 3, guides = "collect")
  
  combined_plot <- top_plot / bottom_plot 
  
  # ---- 6. Species centroid labels plot ----
  centroid_plot <- ggplot(species_centroids, aes(x = PC1_mean, y = PC2_mean, color = species)) +
    geom_point(size = 4) +
    geom_text(aes(label = species), vjust = -1, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(x = "PC1", y = "PC2", color = "Species",
         title = paste("Species mean in FPCA Space:", scheme_name)) +
    scale_color_brewer(palette = "Set2")
  
  # ---- 7. Species ranking by isolation ----
  species_rank_mean <- species_centroids %>%
    arrange(PC1_mean)
  
  species_rank_median <- species_centroids %>%
    arrange(PC1_median)
  
  centroid_plot_median <- ggplot(species_centroids, aes(x = PC1_median, y = PC2_median, color = species)) +
    geom_point(size = 4) +
    geom_text(aes(label = species), vjust = -1, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(x = "PC1", y = "PC2", color = "Species",
         title = paste("Species median in FPCA Space:", scheme_name)) +
    scale_color_brewer(palette = "Set2")
  
  # Return results
  return(list(
    pca_res = pca_res,
    loadings = loadings,
    scores = scores,
    phi1 = phi1,
    phi2 = phi2,
    p_loadings = p_loadings,
    pve = pve,
    cum_pve = cum_pve,
    p_Varexplained = p_Varexplained,
    landscape_plot = landscape_plot,
    species_scores = species_scores,
    species_centroids = species_centroids,
    combined_plot = combined_plot,
    centroid_plot = centroid_plot,
    centroid_plot_median = centroid_plot_median,
    species_rank_mean = species_rank_mean,
    species_rank_median = species_rank_median,
    fpca_patch_df = fpca_patch_df
  ))
} 





## ________________________________________________________
## Functional PCA of nearest-neighbour distance profiles
## Implemented as PCA on discrete functional data
## ________________________________________________________



## ________________________________________________________
## 1. Input: list of species-specific distance files, Read data
## ________________________________________________________

## Obs this will read in both distance data (sorted) and the meta_data (patch IDs)

data_dir ="C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/Habitats_from_sustainscapes/Habitat_estimates/Reordered_data/"



# List all .rds files 
rds_list <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE) 
# Create names without the .rds extension 
file_names <- tools::file_path_sans_ext(basename(rds_list)) 
# Read all .rds files into a named list 
rds_list <- setNames( lapply(rds_list, readRDS), file_names)


head(rds_list)


data_rds <- rds_list[!grepl("^Meta_data", names(rds_list))]

#head(data_rds)

meta_rds <- rds_list[grepl("^Meta_data", names(rds_list))]

#head(meta_rds)




# #### plot the data 
# 
# df_long <- df_counts %>%
#   pivot_longer(-species, names_to = "scheme", values_to = "n")
# 
# ggplot(df_long, aes(x = scheme, y = species, fill = n)) +
#   geom_tile() +
#   scale_fill_viridis_c(option = "C", trans = "log10") +
#   theme_bw() +
#   labs(
#     x = "Scheme",
#     y = "Species",
#     fill = "Row count (log10)"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# 
# 
# 
# 
# 
# #### Join patches pr. scheme across species
# 
# #Currently no downsampling pre building the "background"
# 
# extract_distances <- function(df, k_max) {
#   mat <- as.matrix(df[, paste0("distance_", 1:k_max)])
#   mat[complete.cases(mat), ]
# }




Y_by_scheme <- list()
species_id_by_scheme <- list()





schemes <- c("broad_High", "broad_Medium", "broad_Low",
             "narrow_High", "narrow_Medium", "narrow_Low")

######pool distances profiles cross schemes Only profiles with more than 5 distances 

for (sc in schemes) {
  
  pattern <- paste0("_", sc, "_")
  scheme_names <- grep(pattern, names(data_rds), value = TRUE)
  
  mats <- list()
  species_vec <- c()
  
  for (nm in scheme_names) {
    
    df <- data_rds[[nm]]   # <-- IMPORTANT: this is a tibble, not a list
    
    if (nrow(df) > 0) {
      
      mat <- extract_distances(df, k_max = 5)
      
      mats[[nm]] <- mat
      
      species_vec <- c(species_vec,
                       rep(nm, nrow(mat)))
    }
  }
  
  if (length(mats) > 0) {
    Y_by_scheme[[sc]] <- do.call(rbind, mats)
    species_id_by_scheme[[sc]] <- species_vec
  }
}



str(Y_by_scheme$broad_High)


### A little ekstra QC

#How many patches pr spp pr. scheme

df_counts <- imap_dfr(data_rds, ~{
  tibble(
    name = .y,
    n = nrow(.x)
  )
}) %>%
  mutate(
    # extract species (everything between "Sorted_" and "_broad/_narrow")
    species = str_extract(name, "(?<=Sorted_).*?(?=_broad|_narrow)"),
    
    # extract scheme
    scheme = str_extract(name, "broad_(High|Medium|Low)|narrow_(High|Medium|Low)")
  ) %>%
  select(species, scheme, n) %>%
  pivot_wider(names_from = scheme, values_from = n) %>%
  arrange(species)%>%
  select(species, all_of(schemes))

df_counts



sum(df_counts$broad_High)





#OBS 
## might be Optional(?): log-transform distances 
#(However it is to miniimize noise in the covariance maxtrix)
#c("broad_High", "broad_Medium", "broad_Low",
#  "narrow_High", "narrow_Medium", "narrow_Low")


Y_log_broad_High <- log1p(Y_by_scheme[["broad_High"]])
Y_log_broad_Medium <- log1p(Y_by_scheme[["broad_Medium"]])
Y_log_broad_Low <- log1p(Y_by_scheme[["broad_Low"]])
Y_log_narrow_High <- log1p(Y_by_scheme[["narrow_High"]])
Y_log_narrow_Medium <- log1p(Y_by_scheme[["narrow_Medium"]])
Y_log_narrow_Low <- log1p(Y_by_scheme[["narrow_Low"]])

Y_log_list <- lapply(Y_by_scheme, log1p)

## NEEDS TO BE = 0 

sum(is.na(Y_by_scheme[["broad_High"]]))
sum(is.na(Y_by_scheme[["broad_Medium"]]))
sum(is.na(Y_by_scheme[["broad_Low"]]))

sum(is.na(Y_by_scheme[["narrow_High"]]))
sum(is.na(Y_by_scheme[["narrow_Medium"]]))
sum(is.na(Y_by_scheme[["narrow_Low"]]))




#### Running function on  Mapping schemes (MappingXQuality)

Out_dir <- "../Results"


schemes <- c("broad_High", "broad_Medium", "broad_Low",
             "narrow_High", "narrow_Medium", "narrow_Low")

Y_log_list <- lapply(Y_by_scheme, log1p)

# Loop over schemes
for (sc in schemes) {
  result_name <- paste0(sc, "_results")  # e.g., "broad_High_results"
  
  assign(result_name, run_FPCA_analysis(sc,  Y_log_list, species_id_by_scheme,
                                        data_rds = data_rds))
}



####
# Plot 


broad_High_results$species_rank_mean
broad_Medium_results$species_rank_mean
broad_Low_results$species_rank_mean
narrow_High_results$species_rank_mean
narrow_Medium_results$species_rank_mean
narrow_Low_results$species_rank_mean


broad_High_results$landscape_plot
broad_Medium_results$landscape_plot
broad_Low_results$landscape_plot
narrow_High_results$landscape_plot
narrow_Medium_results$landscape_plot
narrow_Low_results$landscape_plot


broad_High_results$combined_plot
broad_Medium_results$combined_plot
broad_Low_results$combined_plot
narrow_High_results$combined_plot
narrow_Medium_results$combined_plot
narrow_Low_results$combined_plot


#########


broad_High_results$pve
broad_Medium_results$pve
broad_Low_results$pve 
narrow_High_results$pve 
narrow_Medium_results$pve 
narrow_Low_results$pve


broad_High_results$p_loadings 
broad_Medium_results$p_loadings 
broad_Low_results$p_loadings 
narrow_High_results$p_loadings 
narrow_Medium_results$p_loadings 
narrow_Low_results$p_loadings


###
# Check out the loadings !



schemes <- c("broad_High", "broad_Medium", "broad_Low",
             "narrow_High", "narrow_Medium", "narrow_Low")

plot_list <- list()

for (sc in schemes) {
  res <- get(paste0(sc, "_results"))
  k <- 1:length(res$phi1)
  
  # phi1
  p1 <- ggplot(data.frame(k=k, phi=res$phi1), aes(x=k, y=phi)) +
    geom_point() + geom_line() +
    geom_hline(yintercept=0, linetype=2, color="grey") +
    ylim(-1,1) +
    labs(title=paste0(sc, " FPCA1"), x="k", y="phi1") +
    theme_minimal()
  
  # phi2
  p2 <- ggplot(data.frame(k=k, phi=res$phi2), aes(x=k, y=phi)) +
    geom_point() + geom_line() +
    geom_hline(yintercept=0, linetype=2, color="grey") +
    ylim(-1,1) +
    labs(title=paste0(sc, " FPCA2"), x="k", y="phi2") +
    theme_minimal()
  
  plot_list <- c(plot_list, list(p1, p2))
}

# Combine: 2 columns, 6 rows
wrap_plots(plot_list, ncol=2)




#### Test if phi is constant? 

check_phi1_constant <- function(results_list, slope_thresh = 0.1) {
  library(dplyr)
  schemes <- names(results_list)
  
  out <- lapply(schemes, function(sc) {
    res <- results_list[[sc]]
    phi1 <- res$phi1
    k <- 1:length(phi1)
    
    lm_res <- lm(phi1 ~ k)
    slope <- coef(lm_res)[2]
    pval <- summary(lm_res)$coefficients[2,4]
    var_phi <- var(phi1)
    practically_constant <- abs(slope) < slope_thresh
    
    tibble(
      scheme = sc,
      slope = slope,
      p_value = pval,
      variance = var_phi,
      practically_constant = practically_constant
    )
  }) %>% bind_rows()
  
  return(out)
}





# Example usage
check_phi1_constant(list(
  broad_High = broad_High_results,
  broad_Medium = broad_Medium_results,
  broad_Low = broad_Low_results,
  narrow_High = narrow_High_results,
  narrow_Medium = narrow_Medium_results,
  narrow_Low = narrow_Low_results
))




####
# Extract the FPCA pr. species pr. patch pr. mapping scheme 

out_dir <- "../Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (sc in schemes) {
  
  res <- get(paste0(sc, "_results"))
  
  # ---- 1. Save patch-level FPCA ----
#  write.csv(
#    res$fpca_patch_df,
#    file = file.path(out_dir, paste0(sc, "_fpca_patch_df.csv")),
#    row.names = FALSE
#  )
  
  # ---- 2. Save variance explained plot ----
#  png(file.path(out_dir, paste0(sc, "_variance_explained.png")),
#      width = 1200, height = 600, res = 150)
#  replayPlot(res$p_Varexplained)
#  dev.off()
  
  # ---- 3. Save landscape plot ----
  ggsave(
    filename = file.path(out_dir, paste0(sc, "_combined_plot.png")),
    plot = res$combined_plot,
    width = 8, height = 6, dpi = 300
  )
}


