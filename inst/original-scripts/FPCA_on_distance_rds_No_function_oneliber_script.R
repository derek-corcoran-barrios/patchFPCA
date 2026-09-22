## ________________________________________________________
## Functional PCA of nearest-neighbour distance profiles
## Implemented as PCA on discrete functional data
## ________________________________________________________


## This scripts is for manually running the FPCA and exploring the effects of the parameteriztion choices. 
# and initial hypothesis testing. 
#The final methods and scripts for hypothesis testing id found in the Hypothesis_testing.R



library(terra)  # working with spatial data 
library(ggplot2)
library(sf)      # for spatial data
library(dplyr)   # for data manipulation
library(purrr)   # for functional programming (loops)
library(future.apply)  # to paralize 
library(progressr)
library(tidyr)
library(stringr)



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



###################Post running the function and saving data.


head(rds_list)


data_rds <- rds_list[!grepl("^Meta_data", names(rds_list))]

head(data_rds)

meta_rds <- rds_list[grepl("^Meta_data", names(rds_list))]

head(meta_rds)





## ________________________________________________________
## 2. QC data,
## ________________________________________________________




Good_counts <- sapply(data_rds, function(df) {
  mat <- as.matrix(df[, paste0("distance_", 1:5)])
  sum(complete.cases(mat))
})

Good_counts

problem_counts <- sapply(data_rds, function(df) {
  mat <- as.matrix(df[, paste0("distance_", 1:5)])
  sum(!complete.cases(mat))
})

problem_counts

proportion_of_problem_distances <- (problem_counts / (Good_counts + problem_counts))*100


df <- data.frame(
  name = names(proportion_of_problem_distances),
  proportion = proportion_of_problem_distances
)

df <- df[order(df$name), ]



# Separate components
df <- df %>%
  separate(name,
           into = c("type","species1", "species2", "niche", "distance", "drop"),
           sep = "_") %>%
  mutate(
    species1 = str_to_title(species1),   # Capitalize first letter
    species = paste(species1, species2, sep = " ")
  ) %>%
  select(species, niche, distance, proportion)

# Alphabetical ordering of species
df$species <- factor(df$species,
                     levels = sort(unique(df$species)))

# Plot
ggplot(df, aes(x = species,
               y = proportion,
               color = species)) +
  geom_point(size = 3) +
  facet_wrap(~ distance) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c(
    "Coenonympha tullia" = "#dfc27d",
    "Lycaena virgaureae" = "#bf812d",
    "Thymelicus lineola" = "#8c510a",
    "Maniola jurtina" = "darkgrey",                        
    "Pieris napi" = "#74a9cf",
    "Gonepteryx rhamni" = "#0570b0",
    "Aglais urticae" = "#034e7b"
  )) +
  guides(color = guide_legend(order=1, title = "Species", nrow = 7)) +
  labs(x = "Species",
       y = "% of problem distances",
       color = "Mapping breadth")


















####3 Below is the oneliner. 
## Use the function instead : FPCA_function.R 













###The full data set have now passed QC, 
#and we can bging to saparate them out for each mapped senario h of H
###   H=(broad/narrow) × (High / H+M / H+M+L).
#                     H=6



## ________________________________________________________
## 3.  build Y of scenerio H, matrix and vector
## ________________________________________________________

#Check out which .rds were are working with ? 

names(data_rds)

head(data_rds$Sorted_aglais_urticae_broad_High_distances_corrected_superpolygons)




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


#### plot the data 

df_long <- df_counts %>%
  pivot_longer(-species, names_to = "scheme", values_to = "n")

ggplot(df_long, aes(x = scheme, y = species, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", trans = "log10") +
  theme_bw() +
  labs(
    x = "Scheme",
    y = "Species",
    fill = "Row count (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )






#### Join patches pr. scheme across species

#Currently no downsampling pre building the "background"

extract_distances <- function(df, k_max) {
  mat <- as.matrix(df[, paste0("distance_", 1:k_max)])
  mat[complete.cases(mat), ]
}

Y_by_scheme <- list()
species_id_by_scheme <- list()





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

sum(is.na(Y_by_scheme[["broad_High"]]))
sum(is.na(Y_by_scheme[["broad_Medium"]]))
sum(is.na(Y_by_scheme[["broad_Low"]]))

sum(is.na(Y_by_scheme[["narrow_High"]]))
sum(is.na(Y_by_scheme[["narrow_Medium"]]))
sum(is.na(Y_by_scheme[["narrow_Low"]]))











































####Fully manual running of the data








## ________________________________________________________
## 4. PCA = FPCA on discrete functional data
## ________________________________________________________


#### NOW PCA ####

## THIS IS IMPORSTANT !!!

### the distance profiles are NOT multivariate data!
#    I
# monotonic functional curve in k-space

#No smoothing (k = 5 is discrete, ordered)

#No centering across k (We are modelling magnitude - Bee! ^_^ )

#No conversion to dense covariance matrices -->
# Now this will:
#Treat each row as one functional observation
#Estimate eigenfunctions in k-space, not patch-space (eg with five covariate) !!! 
#Avoid n × n covariance matrices ^^




pca_res_broad_High <- prcomp(
  Y_log_broad_High,
  center = FALSE,
  scale. = FALSE
)

scores_broad_High <- pca_res_broad_High$x          # FPCA scores
loadings_broad_High <- pca_res_broad_High$rotation # FPCA eigenfunctions



## ________________________________________________________
## 5.a Extract eigen loadings ("functions")
## ________________________________________________________

k <- 1:5

phi1_broad_High <- loadings_broad_High[,1]
phi1_broad_High

phi2_broad_High <- loadings_broad_High[,2]

phi2_broad_High

## Plot eigenfunctions
par(mfrow = c(1,2))


plot(k, phi1_broad_High, type = "b", pch = 19,
     ylim = c(-1,1),
     xlab = "Neighbor rank (k)",
     ylab = "Loading_broad_High",
     main = "FPCA1")
abline(h = 0, lty = 2, col = "grey")

plot(k, phi2_broad_High, type = "b", pch = 19,
     ylim = c(-1,1),
     xlab = "Neighbor rank (k)",
     ylab = "Loading_broad_High",
     main = "FPCA2")
abline(h = 0, lty = 2, col = "grey")

p_loadings_broad_High <- recordPlot()   # ← save plot as object
par(mfrow = c(1,1))





## ________________________________________________________
## 5.b VARIANCE EXPLAINED
## ________________________________________________________



###Variance explained

# Eigenvalues
eigvals_broad_High <- pca_res_broad_High$sdev^2

# Proportion of variance explained
pve_broad_High <- eigvals_broad_High / sum(eigvals_broad_High)

# Cumulative variance explained
cum_pve_broad_High <- cumsum(pve_broad_High)

par(mfrow = c(1,2))

# Inspect

plot(pve_broad_High, type = "b", pch = 19,
     xlab = "Component",
     ylab = "Proportion of variance explained of: broad_High",
     main = "FPCA variance explained: broad_High")


barplot(cum_pve_broad_High,
        xlab = "Component",
        ylab = "Cumulative variance explained of: broad_High",
        main = "Cumulative FPCA variance explained:broad_High")
abline(h = 0.9, lty = 2, col = "grey")

p_Varexplained_broad_High <- recordPlot() 


## ________________________________________________________
## 6. plot the FPCA on the landscape
## ________________________________________________________


# Convert PCA scores to tibble
scores_df <- as.data.frame(scores_broad_High)

# Thin by 90% (keep 10%)
#set.seed(42)  # for reproducibility
scores_thinned <- scores_df %>% sample_frac(0.1)

# Plot
Landscape_broad_High_poit <- ggplot(scores_thinned, aes(x = PC1, y = PC2)) +
  geom_point(color = "lightgrey") +
  theme_bw() +
  labs(
    title = "FPCA scores (broad_High) - thinned 90%",
    x = "FPCA 1",
    y = "FPCA 2"
  )

range(scores_thinned$PC1)
Landscape_broad_High


ggplot(scores_thinned, aes(x = PC1, y = PC2)) +
  geom_density_2d(color = "black", linewidth = 0.6) +
  theme_bw() +
  labs(title = "FPCA scores (broad_High) - thinned 90%", x = "FPCA 1", y = "FPCA 2") +
  xlim(min(scores_thinned$PC1), max(scores_thinned$PC1)) +
  ylim(min(scores_thinned$PC2), max(scores_thinned$PC2))


## ________________________________________________________
## 7. Projecting species paths onto the main axsis of variation. 
## ________________________________________________________


species_scores_broad_High <- data.frame(
  species = species_id_by_scheme[["broad_High"]],
  PC1 = scores_broad_High[,1],
  PC2 = scores_broad_High[,2]
)



## ________________________________________________________
## 7.b Species-level summaries (centroids = the mean FPCA(s))
## ________________________________________________________

species_centroids_broad_High <- species_scores_broad_High %>%
  group_by(species) %>%
  summarise(
    PC1_mean = mean(PC1),
    PC2_mean = mean(PC2),
    .groups = "drop"
  )



## ________________________________________________________
## 8. Alternative Plotting ! ALL species in one window ! 
## ________________________________________________________

# par(mfrow = c(1,1))
# 
# 
# species_levels <- unique(species_centroids_broad_High$species)
# cols <- setNames(rainbow(length(species_levels)), species_levels)
# 
# 
# plot(species_scores_broad_High$PC1, species_scores_broad_High$PC2,
#      col = cols[cols],
#      pch = 16,
#      cex = 0.05,              # SMALL patches
#      xlab = "FPCA1 (Isolation)",
#      ylab = "FPCA2 (arched FPC2)",
#      main = "Patch-level connectivity profiles")
# 
# points(
#   species_centroids_broad_High$PC1_mean,
#   species_centroids_broad_High$PC2_mean,
#   pch = 21,                  # filled circle with border
#   bg = cols[species_centroids_broad_High $species],
#   col = "black",             # black border
#   cex = 2                    # BIG
# )




####Color tv plot! 
# Ensure species is a factor
species_scores_broad_High$species <- factor(species_scores_broad_High$species)
species_centroids_broad_High$species <- factor(species_centroids_broad_High$species,
                                               levels = levels(species_scores_broad_High$species))

# Thin per species
set.seed(42)
species_scores_thinned <- species_scores_broad_High %>%
  group_by(species) %>%
  sample_frac(0.1) %>%
  ungroup()

# Color mapping
species_levels <- levels(species_scores_broad_High$species)
cols <- setNames(rainbow(length(species_levels)), species_levels)
cols_alpha <- paste0(cols, "80")  # ~50% opacity

# Scatter plot for thinned patches (with alpha)
plot(species_scores_thinned$PC1, species_scores_thinned$PC2,
     col = cols_alpha[species_scores_thinned$species],  # alpha applied
     pch = 16,
     cex = 0.05,
     xlab = "FPCA1 (Isolation)",
     ylab = "FPCA2 (Arched FPC2)",
     main = "Patch-level connectivity profiles")

# Overlay centroids (no alpha)
points(
  species_centroids_broad_High$PC1_mean,
  species_centroids_broad_High$PC2_mean,
  pch = 21,
  bg = cols[species_centroids_broad_High$species],  # fully opaque
  col = "black",
  cex = 2
)



###### ONE PLOT PR: SP

library(patchwork)

# Convert to tibble and thin by 90%
set.seed(42)
species_scores_thinned <- species_scores_broad_High_df %>% sample_frac(0.1)

species_scores_thinned$species <- factor(species_scores_thinned$species)
species_centroids_broad_High$species <- factor(species_centroids_broad_High$species)

# Color mapping
species_levels <- levels(species_scores_thinned$species)
cols <- setNames(rainbow(length(species_levels)), species_levels)
cols_alpha <- paste0(cols, "80")  # ~50% opacity for patches

# --- Top row: All species ---
top_plot <- ggplot() +
  geom_point(
    data = species_scores_thinned,
    aes(x = PC1, y = PC2, color = species),
    size = 0.5,
    alpha = 0.5
  ) +
  geom_point(
    data = species_centroids_broad_High,
    aes(x = PC1_mean, y = PC2_mean, fill = species),
    color = "black",
    pch = 21,
    size = 3
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  labs(title = "All species", x = "FPCA1", y = "FPCA2")

# --- Bottom row: One panel per species ---
bottom_plot <- ggplot(species_scores_thinned, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = species), size = 0.5, alpha = 0.5) +
  geom_point(
    data = species_centroids_broad_High,
    aes(x = PC1_mean, y = PC2_mean, fill = species),
    color = "black",
    pch = 21,
    size = 3
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  facet_wrap(~species, ncol = 3) +
  theme_bw() +
  labs(title = "Per-species panels", x = "FPCA1", y = "FPCA2")

# Combine top + bottom
broad_high_final_plot <- top_plot / bottom_plot


## ________________________________________________________
## 9. ARCH???
## ________________________________________________________



species_levels <- unique(species_centroids_broad_High$species)
cols <- setNames(rainbow(length(species_levels)), species_levels)



library(ggplot2)

# Plot is 
ggplot(species_centroids_broad_High, aes(x = PC1_mean, y = PC2_mean, color = species)) +
  geom_point(size = 4) +          # bigger points
  geom_text(aes(label = species),  # optional: add labels
            vjust = -1, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(x = "PC1", y = "PC2", color = "Species",
       title = "Species Centroids in PCA Space") +
  scale_color_brewer(palette = "Set2") # nice qualitative colors





$species_centroids



broad_High_results$species_centroids
broad_Medium_results$species_centroids
broad_Low_results$species_centroids
narrow_High_results$species_centroids
narrow_Medium_results$species_centroids
narrow_Low_results$species_centroids












## ________________________________________________________
## 11. rank species (no hypothesis testing)
## ________________________________________________________

species_rank_isolation_broad_High <- species_centroids_broad_High %>%
  arrange(PC1_mean)


#manually computing the "c" of the mean function. 
phi1_broad_High 

phi2_broad_High 






#manually computing the "c" of the mean function. 
FPCA1_manual <- as.numeric(Y_log %*% phi1_broad_High)
all.equal(FPCA1_manual, scores[,1])
# GIing FPCA1i(s)​=k=1∑5​ϕ1​(k)⋅log(1+dik(s)​)

#manually computing the "c" of the mean function. 
FPCA2_manual <- as.numeric(Y_log %*% phi2)
all.equal(FPCA1_manual, scores[,2])




















