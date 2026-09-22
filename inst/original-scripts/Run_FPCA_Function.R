
## ________________________________________________________
## Functional PCA of nearest-neighbour distance profiles
## Implemented as PCA on discrete functional data
## center = FALSE, scale. = FALSE
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

### Alittle ekstra QC


schemes <- c("broad_High", "broad_Medium", "broad_Low",
             "narrow_High", "narrow_Medium", "narrow_Low")

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


###  Input the mapping schemes 


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






head(Y_by_scheme$broad_High)

str(Y_by_scheme$broad_High)

sum(df_counts$broad_High)




#OBS 
# We log-transform distances 
#(It is to minimize noise in the covariance maxtrix)
#c("broad_High", "broad_Medium", "broad_Low",
#  "narrow_High", "narrow_Medium", "narrow_Low")


Y_log_broad_High <- log1p(Y_by_scheme[["broad_High"]])
Y_log_broad_Medium <- log1p(Y_by_scheme[["broad_Medium"]])
Y_log_broad_Low <- log1p(Y_by_scheme[["broad_Low"]])
Y_log_narrow_High <- log1p(Y_by_scheme[["narrow_High"]])
Y_log_narrow_Medium <- log1p(Y_by_scheme[["narrow_Medium"]])
Y_log_narrow_Low <- log1p(Y_by_scheme[["narrow_Low"]])

Y_log_list <- lapply(Y_by_scheme, log1p)



## Checking for Na
#NEEDS TO BE = 0 

sum(is.na(Y_by_scheme[["broad_High"]]))
sum(is.na(Y_by_scheme[["broad_Medium"]]))
sum(is.na(Y_by_scheme[["broad_Low"]]))

sum(is.na(Y_by_scheme[["narrow_High"]]))
sum(is.na(Y_by_scheme[["narrow_Medium"]]))
sum(is.na(Y_by_scheme[["narrow_Low"]]))





#### Running function on  Mapping schemes (MappingXQuality)

Out_dir <- "../Results"


Y_log_list <- lapply(Y_by_scheme, log1p)

# Loop over schemes
for (sc in schemes) {
  result_name <- paste0(sc, "_results")  # e.g., "broad_High_results"
  
  assign(result_name, run_FPCA_analysis(sc,  Y_log_list, species_id_by_scheme,
                                        data_rds = data_rds))
}


### This funtions rentuns a list of dataframes containing many estimates: 


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


#### Expamples of extraction of objects.


# View objects  
#broad_High_results$ANY OJECT FORM THE ABOVE LIST 


#Eksamples of objects
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

broad_High_results$p_loadings 
broad_Medium_results$p_loadings 
broad_Low_results$p_loadings 
narrow_High_results$p_loadings 
narrow_Medium_results$p_loadings 
narrow_Low_results$p_loadings


###
# Check out the loadings across schemes!



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




####
# Extract the FPCA score for each patch pr. mapping scheme pr. species, 
#Save variance explained as .png
#Save the species projection combined plot.  

out_dir <- "../Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (sc in schemes) {
  
  res <- get(paste0(sc, "_results"))
  
#  ---- 1. Save patch-level FPCA ----
   write.csv(
     res$fpca_patch_df,
     file = file.path(out_dir, paste0(sc, "_fpca_patch_df.csv")),
     row.names = FALSE
   )

#  ---- 2. Save variance explained plot ----
   png(file.path(out_dir, paste0(sc, "_variance_explained.png")),
       width = 1200, height = 600, res = 150)
   replayPlot(res$p_Varexplained)
   dev.off()
  
#  ---- 3. Save landscape plot ----
  ggsave(
    filename = file.path(out_dir, paste0(sc, "_combined_plot.png")),
    plot = res$combined_plot,
    width = 8, height = 6, dpi = 300
  )
}

