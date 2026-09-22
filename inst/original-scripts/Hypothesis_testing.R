
#### Running the FPCA function on all data
##FPCA_function.R






#Then pull results


all_results <- list(
  broad_High = broad_High_results,
  broad_Medium = broad_Medium_results,
  broad_Low = broad_Low_results,
  narrow_High = narrow_High_results,
  narrow_Medium = narrow_Medium_results,
  narrow_Low = narrow_Low_results)







###############################################################
##### Testing H1.1: 

#Existence of a dominant isolation axis
#evalueted trough variance explained

broad_High_results$pve 
broad_Medium_results$pve 
broad_Low_results$pve 
narrow_High_results$pve 
narrow_Medium_results$pve 
narrow_Low_results$pve

#Ckeck the estimates

broad_High_results$p_Varexplained 
broad_Medium_results$p_Varexplained 
broad_Low_results$p_Varexplained 
narrow_High_results$p_Varexplained 
narrow_Medium_results$p_Varexplained 
narrow_Low_results$p_Varexplained 





###############################################################
##### Testing H1.2: 

#### Test if phi is constant? 



check_loading_constancy <- function(results_list, slope_thresh = 0.04) {
  library(dplyr)
  
  out <- lapply(names(results_list), function(sc) {
    
    phi1 <- results_list[[sc]]$phi1
    k <- seq_along(phi1)
    
    # linear trend
    lm_res <- lm(phi1 ~ k)
    slope <- coef(lm_res)[2]
    pval <- summary(lm_res)$coefficients[2,4]
    r2 <- summary(lm_res)$r.squared
    
    # coefficient of variation (scale-free)
    cv <- sd(phi1) / abs(mean(phi1))
    
    # range
    range_phi <- max(phi1) - min(phi1)
    
    # practical constancy
    practically_constant <- abs(slope) < slope_thresh
    
    data.frame(
      scheme = sc,
      slope = slope,
      p_value = pval,
      r_squared = r2,
      cv = cv,
      range = range_phi,
      practically_constant = practically_constant
    )
    
  }) %>% bind_rows()
  
  return(out)
}



# Example usage
check_loading_constancy(list(
  broad_High = broad_High_results,
  broad_Medium = broad_Medium_results,
  broad_Low = broad_Low_results,
  narrow_High = narrow_High_results,
  narrow_Medium = narrow_Medium_results,
  narrow_Low = narrow_Low_results
))







###############################################################
##### Testing H1.2b — Local vs Regional dominance

#Including the meadian: 

local_regional_test <- function(v){
  
  L <- mean(v[1:2])
  R <- mean(v[3:5])
  
  c(L=L, R=R, diff=L-R)
}

local_regional_test(broad_High_results$loadings[,1])
local_regional_test(broad_Medium_results$loadings[,1])
local_regional_test(broad_Low_results$loadings[,1])

local_regional_test(narrow_High_results$loadings[,1])
local_regional_test(narrow_Medium_results$loadings[,1])
local_regional_test(narrow_Low_results$loadings[,1])




#Excluding the meadian: 

local_regional_test_nomedian <- function(v){
  
  L <- mean(v[1:2])
  R <- mean(v[4:5])
  
  c(L=L, R=R, diff=L-R)
}

local_regional_test_nomedian(broad_High_results$loadings[,1])
local_regional_test_nomedian(broad_Medium_results$loadings[,1])
local_regional_test_nomedian(broad_Low_results$loadings[,1])

local_regional_test_nomedian(narrow_High_results$loadings[,1])
local_regional_test_nomedian(narrow_Medium_results$loadings[,1])
local_regional_test_nomedian(narrow_Low_results$loadings[,1])









###############################################################
##### Testing H1.3: Across mapping scheme and quality there is a stable representation of landscape structure

#Testing H1.3a

### OBS Eigenvectors can flip sign arbitrarily.!!!! 

#Funktion:
cosine_similarity <- function(v1, v2) {
  
  # Sign alignment (PCA indeterminacy)
  if(sum(v1 * v2) < 0) {
    v2 <- -v2
  }
  
  sum(v1 * v2) / 
    (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}



#Compare Broad vs Narrow within quality


compare_mapping_invariance <- function(all_results, pc = 1) {
  
  qualities <- c("High", "Medium", "Low")
  
  out <- data.frame(
    quality = qualities,
    cosine_similarity = NA,
    angle_degrees = NA
  )
  
  for(q in qualities){
    
    v_broad  <- all_results[[paste0("broad_", q)]]$loadings[, pc]
    v_narrow <- all_results[[paste0("narrow_", q)]]$loadings[, pc]
    
    cos_sim <- cosine_similarity(v_broad, v_narrow)
    
    out$cosine_similarity[out$quality == q] <- cos_sim
    out$angle_degrees[out$quality == q] <- acos(cos_sim) * 180 / pi
  }
  
  return(out)
}






H1.2_results <- compare_mapping_invariance(all_results, pc = 1)

H1.2_results


#cosine ≈ 1 → structural invariance (null supported)

#cosine < 0.9 → structural distortion

#cosine ≈ 0 → orthogonal structure (strong alternative)





###############################################################
##### Testing H1.3b 



compare_quality_invariance <- function(all_results, mapping, pc = 1) {
  
  qualities <- c("High","Medium","Low")
  combs <- combn(qualities, 2)
  
  out <- data.frame(
    q1 = combs[1,],
    q2 = combs[2,],
    cosine = NA,
    angle = NA
  )
  
  for(i in 1:ncol(combs)){
    
    v1 <- all_results[[paste0(mapping,"_",combs[1,i])]]$loadings[,pc]
    v2 <- all_results[[paste0(mapping,"_",combs[2,i])]]$loadings[,pc]
    
    cos_sim <- cosine_similarity(v1,v2)
    
    out$cosine[i] <- cos_sim
    out$angle[i]  <- acos(cos_sim)*180/pi
  }
  
  out
}

compare_quality_invariance(all_results,mapping ="broad")

# 
# q1     q2    cosine      angle
# 1   High Medium 0.9999999 0.02613858
# 2   High    Low 0.9998486 0.99701877
# 3 Medium    Low 0.9998409 1.02212291
# 

compare_quality_invariance(all_results,mapping ="narrow")



# 
# q1     q2    cosine     angle
# 1   High Medium 0.9999493 0.5767060
# 2   High    Low 0.9999201 0.7241575
# 3 Medium    Low 0.9999966 0.1496178
# 




###############################################################
##### Testing H1.4 if there is evidence of differences one could have tested for interactiong effects. 

  
# 
# 
# #Δ=​(v1 broad,q1 − v1 broad,q2)−(v1 narrow, q1	−v1narrow,q2)​
# 
# #This is a Euclidean norm in ℝ⁵.
# 
# 
# interaction_distance <- function(all_results, q1, q2, pc = 1){
#   
#   v_b_q1 <- all_results[[paste0("broad_",q1)]]$loadings[,pc]
#   v_b_q2 <- all_results[[paste0("broad_",q2)]]$loadings[,pc]
#   
#   v_n_q1 <- all_results[[paste0("narrow_",q1)]]$loadings[,pc]
#   v_n_q2 <- all_results[[paste0("narrow_",q2)]]$loadings[,pc]
#   
#   delta <- (v_b_q1 - v_b_q2) - (v_n_q1 - v_n_q2)
#   
#   sqrt(sum(delta^2))
# }
# 
# 
# interaction_distance(
#   list(
#     broad_High   = broad_High_results,
#     broad_Medium = broad_Medium_results,
#     broad_Low    = broad_Low_results,
#     narrow_High   = narrow_High_results,
#     narrow_Medium = narrow_Medium_results,
#     narrow_Low    = narrow_Low_results
#   ),
#   q1 = "High",
#   q2 = "Low",
#   pc = 1
# )
# 
# 
# 
# interaction_distance(
#   list(
#     broad_High   = broad_High_results,
#     broad_Medium = broad_Medium_results,
#     broad_Low    = broad_Low_results,
#     narrow_High   = narrow_High_results,
#     narrow_Medium = narrow_Medium_results,
#     narrow_Low    = narrow_Low_results
#   ),
#   q1 = "High",
#   q2 = "Medium",
#   pc = 1
# )
# 
# 
# interaction_distance(
#   list(
#     broad_High   = broad_High_results,
#     broad_Medium = broad_Medium_results,
#     broad_Low    = broad_Low_results,
#     narrow_High   = narrow_High_results,
#     narrow_Medium = narrow_Medium_results,
#     narrow_Low    = narrow_Low_results
#   ),
#   q1 = "Low",
#   q2 = "Medium",
#   pc = 1
# )
# 






#### QUESTION 2   --- Ranks: 
###############################################################
##### Testing H2.1a robustness across mapping


clean_species <- function(x) {
  x |>
    str_remove("Sorted_") |>
    str_remove("_broad_.*|_narrow_.*")
}


df_all <- bind_rows(
  broad_High_results$species_centroids  |> mutate(dataset = "broad_High"),
  broad_Medium_results$species_centroids |> mutate(dataset = "broad_Medium"),
  broad_Low_results$species_centroids   |> mutate(dataset = "broad_Low"),
  narrow_High_results$species_centroids |> mutate(dataset = "narrow_High"),
  narrow_Medium_results$species_centroids |> mutate(dataset = "narrow_Medium"),
  narrow_Low_results$species_centroids  |> mutate(dataset = "narrow_Low")
) |>
  mutate(species = clean_species(species))

df_wide <- df_all |>
  select(species, dataset, PC1_mean) |>
  pivot_wider(names_from = dataset, values_from = PC1_mean)


df_wide_median <- df_all |>
  select(species, dataset, PC1_median) |>
  pivot_wider(names_from = dataset, values_from = PC1_median)


ggplot(df_all, aes(x = dataset, y = PC1_mean, color = species, group = species)) +
  geom_point(size = 2.5) +
  geom_line() +
  scale_color_manual(values = c(
    "Coenonympha tullia" = "#dfc27d",
    "Lycaena virgaureae" = "#bf812d",
    "Thymelicus lineola" = "#8c510a",
    "Maniola jurtina" = "darkgrey",                        
    "Pieris napi" = "#74a9cf",
    "Gonepteryx rhamni" = "#0570b0",
    "Aglais urticae" = "#034e7b"
  )) +
  theme_bw() +
  labs(x = "Dataset", y = "FPCA1 (PC1 mean)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### Corelation: 

df_num <- df_wide |>
  dplyr::select(-species)



cor_mat <- df_wide |>
  select(-species) |>
  cor(method = "spearman")

cor_mat


library(ggplot2)

corr_df <- as.data.frame(as.table(cor_mat))

ggplot(corr_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Correlation")



###2 within Q

###Does niche definition change species positioning?

cor(df_wide$broad_High, df_wide$narrow_High, method = "spearman")
cor(df_wide$broad_Medium, df_wide$narrow_Medium, method = "spearman")
cor(df_wide$broad_Low, df_wide$narrow_Low, method = "spearman")


### With in  board

#“Does habitat quality shift species along FPCA1?”

cor(df_wide$broad_High, df_wide$broad_Medium, method = "spearman")
cor(df_wide$broad_High, df_wide$broad_Low, method = "spearman")
cor(df_wide$broad_Medium, df_wide$broad_Low, method = "spearman")



cor(df_wide$narrow_High, df_wide$narrow_Medium, method = "spearman")
cor(df_wide$narrow_High, df_wide$narrow_Low, method = "spearman")
cor(df_wide$narrow_Medium, df_wide$narrow_Low, method = "spearman")



### First thought

#Species positioning along FPCA1 is somewhat conserved across habitat quality levels, 
#with near-perfect correlation between High and Medium conditions. 
#However, including Low-quality environments induce a restructuring of species rank. 
#Notably, the influence of niche definition (broad vs narrow) is minimal under Medium and Low conditions,
#but becomes pronounced under High-quality conditions, suggesting that ecological differentiation is most
#sensitive to niche definition when environmental constraints are relaxed.

#### Code for ploting ordinal shifts in Rank. 


#Transform into ranks: 

library(dplyr)

df_ranked <- df_wide %>%
  mutate(
    across(
      -species,
      ~ min_rank(desc(round(.x, 1)))
    )
  )

df_ranked


library(dplyr)

df_ranked_medium <- df_wide_median %>%
  mutate(
    across(
      -species,
      ~ min_rank(desc(round(.x, 1)))
    )
  )

df_ranked_medium






library(dplyr)
library(tidyr)
library(ggplot2)

# desired species order
species_order <- c(
  "Coenonympha tullia",
  "Lycaena virgaureae",
  "Thymelicus lineola",
  "Maniola jurtina",
  "Pieris napi",
  "Gonepteryx rhamni",
  "Aglais urticae"
)

species_order <- rev(species_order)

# species colors
color <- c(
  "Coenonympha tullia" = "#dfc27d",
  "Lycaena virgaureae" = "#bf812d",
  "Thymelicus lineola" = "#8c510a",
  "Maniola jurtina" = "darkgrey",
  "Pieris napi" = "#74a9cf",
  "Gonepteryx rhamni" = "#0570b0",
  "Aglais urticae" = "#034e7b"
)

# desired column order
column_order <- c(
  "broad_High",
  "broad_Medium",
  "broad_Low",
  "narrow_High",
  "narrow_Medium",
  "narrow_Low"
)


# long format
df_long <- df_ranked %>%
  pivot_longer(
    -species,
    names_to = "scenario",
    values_to = "rank"
  ) %>%
  mutate(
    species = factor(species, levels = rev(species_order)),
    scenario = factor(scenario, levels = column_order)
  )

# plot
ggplot(df_long, aes(x = scenario, y = species, fill = rank)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = rank), size = 5) +
  
  scale_fill_viridis_c(direction = -1) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "Rank"
  ) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # color species labels
    axis.text.y = element_text(
      color = color[rev(species_order)],
      face = "bold"
    )
  )

sum<-df_wide |> 
  mutate(shift_B_HL = broad_High - broad_Low)

sum<-sum |> 
  mutate(shift_B_HM = broad_High - broad_Medium)

sum<-sum |> 
  mutate(shift_B_ML = broad_Medium - broad_Low)


sum<-sum |> 
  mutate(shift_N_HL = narrow_High - narrow_Low)

sum<-sum |> 
  mutate(shift_N_HM = narrow_High - narrow_Medium)

sum<-sum |> 
  mutate(shift_N_ML = narrow_Medium - narrow_Low)

####


sum<-sum |> 
  mutate(shift_BN_H = broad_High - narrow_High)

sum<-sum |> 
  mutate(shift_BN_M = broad_Medium - narrow_Medium)

sum<-sum |> 
  mutate(shift_BN_L = broad_Low - narrow_Low)



# 
# write.csv(
#   sum,
#   file = file.path(out_dir, "Rank_shift.csv"),
#   row.names = FALSE
# )
# head(sum)




df_long_rank <- sum %>%
  pivot_longer(
    cols = starts_with("shift"),
    names_to = "shift_type",
    values_to = "value"
  )

df_long_rank <- df_long_rank %>%
  mutate(group = case_when(
    grepl("^shift_B_", shift_type) ~ "Broad",
    grepl("^shift_N_", shift_type) ~ "Narrow",
    grepl("^shift_BN_", shift_type) ~ "Within Q, Between mappings"
  ))



ggplot(df_long_rank, aes(x = shift_type, y = species, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~group, scales = "free_x")





# reshape both broad and narrow
df_plot2 <- df_wide |>
  pivot_longer(
    -species,
    names_to = "type_quality",
    values_to = "FPCA1"
  ) |>
  separate(type_quality, into = c("type", "quality"), sep = "_") |>
  mutate(
    type = factor(type, levels = c("broad", "narrow")),
    quality = factor(quality, levels = c("High", "Medium", "Low"))
  )

# simple trajectory plot with both panels
ggplot(df_plot2, aes(x = quality, y = FPCA1, group = species, color = species)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  facet_wrap(~type, nrow = 1) +
  theme_bw() +
  labs(x = "Habitat quality", y = "FPCA1 (PC1)") +
  scale_color_manual(values = c(
    "Coenonympha tullia" = "#dfc27d",
    "Lycaena virgaureae" = "#bf812d",
    "Thymelicus lineola" = "#8c510a",
    "Maniola jurtina" = "darkgrey",                        
    "Pieris napi" = "#74a9cf",
    "Gonepteryx rhamni" = "#0570b0",
    "Aglais urticae" = "#034e7b"
  )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))


#############################################################################################


###########################################           CORRELATION to  summary states


#############################################################################################




###############How well correlated to summary state? 
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(MuMIn)

# -----------------------------
# LOAD + CLEAN
# -----------------------------
stats_clean <- read.csv("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/Habitats_from_sustainscapes/Habitat_estimates/Raw_data/AllSpecies_PatchDistanceStats_CorrectedSuperpolygons.csv") |>
  rename(species = Species) |>
  mutate(
    species = str_replace_all(species, "_", " "),
    species = str_to_sentence(species),
    species = factor(species, levels = c(
      "Coenonympha tullia",
      "Lycaena virgaureae",
      "Thymelicus lineola",
      "Maniola jurtina",
      "Pieris napi",
      "Gonepteryx rhamni",
      "Aglais urticae"
    ))
  )


stats_clean$Total <- stats_clean$n_patches*stats_clean$MeanArea

metrics <- c("n_patches","Total", "MedianArea", "MeanArea", "MedianNN", "MeanNN")

metric_levels <- metrics



#### plot the summary stats: 


#merge data 


library(dplyr)
library(tidyr)

# 1. Prepare df_all: split dataset + rename PC1
df_all2 <- df_all %>%
  separate(dataset, into = c("Definition", "Threshold"), sep = "_") %>%
  rename(FPCA1 = PC1_mean)

# 2. Join with stats_clean
df_joined <- df_all2 %>%
  left_join(stats_clean, by = c("species", "Definition", "Threshold"))

# 3. Optional sanity check (unmatched rows)
unmatched <- anti_join(
  df_all2, stats_clean,
  by = c("species", "Definition", "Threshold")
)

# print if needed
print(unmatched)

# 4. Reshape to long format (INCLUDING Total)
df_plot_all <- df_joined %>%
  pivot_longer(
    cols = c(FPCA1, n_patches, Total, MedianArea, MeanArea, MedianNN, MeanNN),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Definition = factor(Definition, levels = c("broad", "narrow")),
    Threshold  = factor(Threshold, levels = c("High", "Medium", "Low")),
    Metric = factor(Metric, levels = c(
      "FPCA1", "n_patches", "Total", 
      "MedianArea", "MeanArea", "MedianNN", "MeanNN"
    ))
  )
#plot

ggplot(df_plot_all,
       aes(x = Threshold, y = Value, group = species, color = species)) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  
  facet_grid(Metric ~ Definition, scales = "free_y") +
  
  theme_bw() +
  labs(
    x = "Habitat quality",
    y = NULL
  ) +
  
  scale_color_manual(values = c(
    "Coenonympha tullia" = "#dfc27d",
    "Lycaena virgaureae" = "#bf812d",
    "Thymelicus lineola" = "#8c510a",
    "Maniola jurtina" = "darkgrey",
    "Pieris napi" = "#74a9cf",
    "Gonepteryx rhamni" = "#0570b0",
    "Aglais urticae" = "#034e7b"
  )) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )




library(dplyr)
library(tidyr)
library(purrr)
library(MuMIn)
library(scales)
library(ggplot2)



# -----------------------------
# MODEL + METHOD SELECTION
# -----------------------------
choose_method_aicc <- function(x, y) {
  
  df <- data.frame(x = x, y = y) |>
    dplyr::filter(is.finite(x), is.finite(y))
  
  if (nrow(df) < 5) return(NA)
  
  safe_lm <- function(formula) {
    tryCatch(lm(formula, data = df), error = function(e) NULL)
  }
  
  m1 <- safe_lm(y ~ x)
  m2 <- safe_lm(log(y + 1) ~ x)
  m3 <- safe_lm(y ~ poly(x, 2))
  
  get_aicc <- function(model) {
    if (is.null(model)) return(Inf)
    tryCatch(AICc(model), error = function(e) Inf)
  }
  
  aiccs <- c(
    linear = get_aicc(m1),
    exponential = get_aicc(m2),
    quadratic = get_aicc(m3)
  )
  
  best <- names(which.min(aiccs))
  
  p_x <- tryCatch(shapiro.test(df$x)$p.value, error = function(e) NA)
  p_y <- tryCatch(shapiro.test(df$y)$p.value, error = function(e) NA)
  
  approx_normal <- (!is.na(p_x) && p_x > 0.05) &&
    (!is.na(p_y) && p_y > 0.05)
  
  if (best == "quadratic") return(NA)
  if (best == "linear" && approx_normal) return("pearson")
  return("spearman")
}

# -----------------------------
# PAIRWISE CORRELATIONS
# -----------------------------
pairwise_corr <- stats_clean |>
  group_by(Definition, Threshold) |>
  group_modify(~{
    
    df <- .x
    combos <- combn(metrics, 2, simplify = FALSE)
    
    map_dfr(combos, function(vars) {
      
      x <- df[[vars[1]]]
      y <- df[[vars[2]]]
      
      method <- choose_method_aicc(x, y)
      
      cor_val <- if (is.na(method)) NA else
        cor(x, y, method = method, use = "complete.obs")
      
      tibble(
        Metric1 = vars[1],
        Metric2 = vars[2],
        best_method = method,
        correlation = cor_val
      )
    })
  }) |>
  ungroup()

print(pairwise_corr, n =180)

# -----------------------------
# MIRROR MATRIX
# -----------------------------
pairwise_full <- pairwise_corr |>
  bind_rows(
    pairwise_corr |>
      transmute(
        Definition, Threshold,
        Metric1 = Metric2,
        Metric2 = Metric1,
        best_method,
        correlation
      )
  ) |>
  mutate(
    Metric1 = factor(Metric1, levels = metric_levels),
    Metric2 = factor(Metric2, levels = metric_levels)
  )

# -----------------------------
# DIAGONAL
# -----------------------------
diag_data <- expand.grid(
  Definition = unique(stats_clean$Definition),
  Threshold = unique(stats_clean$Threshold),
  Metric1 = metric_levels,
  Metric2 = metric_levels
) |>
  dplyr::filter(Metric1 == Metric2)

# -----------------------------
# RAW DATA FOR SCATTER
# -----------------------------
df_pairs <- stats_clean |>
  group_by(Definition, Threshold) |>
  group_modify(~{
    
    df <- .x
    
    expand.grid(
      Metric1 = metrics,
      Metric2 = metrics,
      species = df$species,
      stringsAsFactors = FALSE
    ) |>
      dplyr::filter(Metric1 != Metric2) |>
      rowwise() |>
      mutate(
        Value1 = df[[Metric1]][match(species, df$species)],
        Value2 = df[[Metric2]][match(species, df$species)]
      ) |>
      ungroup()
    
  }) |>
  ungroup() |>
  group_by(Definition, Threshold, Metric1, Metric2) |>
  mutate(
    x_scaled = scales::rescale(Value1),
    y_scaled = scales::rescale(Value2)
  ) |>
  ungroup() |>
  mutate(
    Metric1 = factor(Metric1, levels = metric_levels),
    Metric2 = factor(Metric2, levels = metric_levels)
  )

# -----------------------------
# ORDER FACETS
# -----------------------------
pairwise_full$Threshold <- factor(pairwise_full$Threshold,
                                  levels = c("High", "Medium", "Low"))
df_pairs$Threshold <- factor(df_pairs$Threshold,
                             levels = c("High", "Medium", "Low"))
diag_data$Threshold <- factor(diag_data$Threshold,
                              levels = c("High", "Medium", "Low"))

pairwise_full$Definition <- factor(pairwise_full$Definition,
                                   levels = c("broad", "narrow"))
df_pairs$Definition <- factor(df_pairs$Definition,
                              levels = c("broad", "narrow"))
diag_data$Definition <- factor(diag_data$Definition,
                               levels = c("broad", "narrow"))

# -----------------------------
# PLOT (UPDATED FOR 7×7)
# -----------------------------
ggplot() +
  
  # DIAGONAL
  geom_tile(
    data = diag_data,
    aes(x = Metric1, y = Metric2),
    fill = "grey85",
    color = "white"
  ) +
  
  # UPPER TRIANGLE
  geom_tile(
    data = pairwise_full |>
      dplyr::filter(match(Metric1, metric_levels) < match(Metric2, metric_levels)),
    aes(x = Metric1, y = Metric2, fill = correlation),
    color = "white"
  ) +
  
  geom_text(
    data = pairwise_full |>
      dplyr::filter(match(Metric1, metric_levels) < match(Metric2, metric_levels)),
    aes(
      x = Metric1,
      y = Metric2,
      label = paste0(round(correlation, 2), "\n(", best_method, ")")
    ),
    size = 2.2   # slightly smaller for 7×7
  ) +
  
  # LOWER TRIANGLE
  geom_point(
    data = df_pairs |>
      dplyr::filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
    aes(
      x = match(Metric1, metric_levels) + (x_scaled - 0.5) * 0.7,
      y = match(Metric2, metric_levels) + (y_scaled - 0.5) * 0.7,
      color = species
    ),
    size = 1,
    alpha = 0.6
  ) +
  
  geom_smooth(
    data = df_pairs |>
      dplyr::filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
    aes(
      x = match(Metric1, metric_levels) + (x_scaled - 0.5) * 0.7,
      y = match(Metric2, metric_levels) + (y_scaled - 0.5) * 0.7,
      group = interaction(Definition, Threshold, Metric1, Metric2)
    ),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_tile(
    data = df_pairs |>
      dplyr::filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
    aes(x = Metric1, y = Metric2),
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  
  facet_grid(Definition ~ Threshold) +
  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey90"
  ) +
  
  scale_color_manual(values = c(
    "Coenonympha tullia" = "#dfc27d",
    "Lycaena virgaureae" = "#bf812d",
    "Thymelicus lineola" = "#8c510a",
    "Maniola jurtina" = "darkgrey",
    "Pieris napi" = "#74a9cf",
    "Gonepteryx rhamni" = "#0570b0",
    "Aglais urticae" = "#034e7b"
  )) +
  
  coord_fixed() +   # 🔥 critical for square tiles
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )






##### CORRELATIONS of FPCA1 to summary stats. 



df_FPCA1 <- df_wide_median |>
  pivot_longer(-species, names_to = "type_threshold", values_to = "FPCA1") |>
  separate(type_threshold, into = c("Definition", "Threshold"), sep = "_") |>
  mutate(
    Definition = factor(Definition, levels = c("broad", "narrow")),
    Threshold = factor(Threshold, levels = c("High", "Medium", "Low"))
  )





# Join FPCA1 with patch metrics


df_joined <- inner_join(
  df_FPCA1,
  stats_clean,
  by = c("species", "Definition", "Threshold")
) 




###PLOT THE DATA: 

df_long <- df_joined |>
  pivot_longer(
    cols = c(n_patches, MeanArea, MedianArea, MeanNN, MedianNN),
    names_to = "Metric",
    values_to = "Value"
  )


ggplot(df_long, aes(x = FPCA1, y = Value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_grid(Metric ~ Definition + Threshold, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "FPCA1",
    y = "Metric value"
  )


ggplot(df_long, aes(x = FPCA1, y = Value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_grid(Metric ~ Definition + Threshold, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "FPCA1",
    y = "Metric value"
  )



library(dplyr)
library(tidyr)

fpca_cor <- df_joined %>%
  group_by(Definition, Threshold) %>%
  summarise(
    across(
      c(n_patches, MedianArea, MeanArea, MedianNN, MeanNN),
      ~ cor(FPCA1, .x, method = "spearman", use = "pairwise.complete.obs")
    ),
    .groups = "drop"
  )



library(ggplot2)

fpca_cor_long <- fpca_cor %>%
  pivot_longer(
    -c(Definition, Threshold),
    names_to = "metric",
    values_to = "correlation"
  )



### order data 

library(dplyr)


metrics <- c("n_patches", "MedianArea", "MeanArea", "MedianNN", "MeanNN")

fpca_cor_long <- fpca_cor_long %>%
  mutate(
    # enforce metric order
    metric = factor(metric, levels = metrics),
    
    # enforce row ordering
    Definition = factor(Definition, levels = c("narrow", "broad")),
    Threshold  = factor(Threshold, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(Definition, Threshold)




ggplot(fpca_cor_long, aes(x = metric, y = 1, fill = correlation)) +
  geom_tile(color = "white") +
  facet_grid(Definition ~ Threshold) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  labs(
    x = "Landscape metric",
    y = "Definition | Threshold",
    fill = "r (FPCA1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




ggplot(fpca_cor_long, aes(x = metric, y = 1, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  facet_grid(Definition ~ Threshold) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  labs(
    x = "Landscape metric",
    y = NULL,
    fill = "r (FPCA1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
















#############################################################################################


###########################################                     Removing C. tullia ! 


#############################################################################################



# 
# 
# 
# # -----------------------------
# # LOAD + CLEAN
# # -----------------------------
# stats_clean <- read.csv("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/Habitats_from_sustainscapes/Habitat_estimates/Raw_data/AllSpecies_PatchDistanceStats_CorrectedSuperpolygons.csv") |>
#   rename(species = Species) |>
#   mutate(
#     species = str_replace_all(species, "_", " "),
#     species = str_to_sentence(species)
#   ) |>
#   filter(species != "Coenonympha tullia") |>
#   mutate(
#     species = factor(species, levels = c(
#       "Lycaena virgaureae",
#       "Thymelicus lineola",
#       "Maniola jurtina",
#       "Pieris napi",
#       "Gonepteryx rhamni",
#       "Aglais urticae"
#     ))
#   )
# 
# metrics <- c("n_patches", "MedianArea", "MeanArea", "MedianNN", "MeanNN")
# metric_levels <- metrics
# 
# # -----------------------------
# # MODEL + METHOD SELECTION
# # -----------------------------
# choose_method_aicc <- function(x, y) {
#   
#   df <- data.frame(x = x, y = y) |>
#     dplyr::filter(is.finite(x), is.finite(y))
#   
#   if (nrow(df) < 5) return(NA)
#   
#   safe_lm <- function(formula) {
#     tryCatch(lm(formula, data = df), error = function(e) NULL)
#   }
#   
#   m1 <- safe_lm(y ~ x)
#   m2 <- safe_lm(log(y + 1) ~ x)
#   m3 <- safe_lm(y ~ poly(x, 2))
#   
#   get_aicc <- function(model) {
#     if (is.null(model)) return(Inf)
#     tryCatch(AICc(model), error = function(e) Inf)
#   }
#   
#   aiccs <- c(
#     linear = get_aicc(m1),
#     exponential = get_aicc(m2),
#     quadratic = get_aicc(m3)
#   )
#   
#   best <- names(which.min(aiccs))
#   
#   p_x <- tryCatch(shapiro.test(df$x)$p.value, error = function(e) NA)
#   p_y <- tryCatch(shapiro.test(df$y)$p.value, error = function(e) NA)
#   
#   approx_normal <- (!is.na(p_x) && p_x > 0.05) &&
#     (!is.na(p_y) && p_y > 0.05)
#   
#   if (best == "quadratic") return(NA)
#   if (best == "linear" && approx_normal) return("pearson")
#   return("spearman")
# }
# 
# # -----------------------------
# # PAIRWISE CORRELATIONS
# # -----------------------------
# pairwise_corr <- stats_clean |>
#   group_by(Definition, Threshold) |>
#   group_modify(~{
#     
#     df <- .x
#     combos <- combn(metrics, 2, simplify = FALSE)
#     
#     map_dfr(combos, function(vars) {
#       
#       x <- df[[vars[1]]]
#       y <- df[[vars[2]]]
#       
#       method <- choose_method_aicc(x, y)
#       
#       cor_val <- if (is.na(method)) NA else
#         cor(x, y, method = method, use = "complete.obs")
#       
#       tibble(
#         Metric1 = vars[1],
#         Metric2 = vars[2],
#         best_method = method,
#         correlation = cor_val
#       )
#     })
#   }) |>
#   ungroup()
# 
# # -----------------------------
# # MIRROR MATRIX (for heatmap)
# # -----------------------------
# pairwise_full <- pairwise_corr |>
#   bind_rows(
#     pairwise_corr |>
#       transmute(
#         Definition, Threshold,
#         Metric1 = Metric2,
#         Metric2 = Metric1,
#         best_method,
#         correlation
#       )
#   ) |>
#   mutate(
#     Metric1 = factor(Metric1, levels = metric_levels),
#     Metric2 = factor(Metric2, levels = metric_levels)
#   )
# 
# # -----------------------------
# # DIAGONAL
# # -----------------------------
# diag_data <- expand.grid(
#   Definition = unique(stats_clean$Definition),
#   Threshold = unique(stats_clean$Threshold),
#   Metric1 = metric_levels,
#   Metric2 = metric_levels
# ) |>
#   filter(Metric1 == Metric2)
# 
# # -----------------------------
# # RAW DATA FOR SCATTER
# # -----------------------------
# df_pairs <- stats_clean |>
#   group_by(Definition, Threshold) |>
#   group_modify(~{
#     
#     df <- .x
#     
#     expand.grid(
#       Metric1 = metrics,
#       Metric2 = metrics,
#       species = df$species,
#       stringsAsFactors = FALSE
#     ) |>
#       filter(Metric1 != Metric2) |>
#       rowwise() |>
#       mutate(
#         Value1 = df[[Metric1]][match(species, df$species)],
#         Value2 = df[[Metric2]][match(species, df$species)]
#       ) |>
#       ungroup()
#     
#   }) |>
#   ungroup() |>
#   group_by(Definition, Threshold, Metric1, Metric2) |>
#   mutate(
#     x_scaled = scales::rescale(Value1),
#     y_scaled = scales::rescale(Value2)
#   ) |>
#   ungroup() |>
#   mutate(
#     Metric1 = factor(Metric1, levels = metric_levels),
#     Metric2 = factor(Metric2, levels = metric_levels)
#   )
# 
# 
# 
# 
# #Order data
# 
# # Order Threshold (columns)
# pairwise_full$Threshold <- factor(
#   pairwise_full$Threshold,
#   levels = c("High", "Medium", "Low")
# )
# 
# df_pairs$Threshold <- factor(
#   df_pairs$Threshold,
#   levels = c("High", "Medium", "Low")
# )
# 
# diag_data$Threshold <- factor(
#   diag_data$Threshold,
#   levels = c("High", "Medium", "Low")
# )
# 
# # Order Definition (rows)
# pairwise_full$Definition <- factor(
#   pairwise_full$Definition,
#   levels = c("broad", "narrow")
# )
# 
# df_pairs$Definition <- factor(
#   df_pairs$Definition,
#   levels = c("broad", "narrow")
# )
# 
# diag_data$Definition <- factor(
#   diag_data$Definition,
#   levels = c("broad", "narrow")
# )
# 
# 
# 
# # -----------------------------
# # PLOT
# # -----------------------------
# ggplot() +
#   
#   # DIAGONAL
#   geom_tile(
#     data = diag_data,
#     aes(x = Metric1, y = Metric2),
#     fill = "grey85",
#     color = "white"
#   ) +
#   
#   # UPPER TRIANGLE
#   geom_tile(
#     data = pairwise_full |>
#       filter(match(Metric1, metric_levels) < match(Metric2, metric_levels)),
#     aes(x = Metric1, y = Metric2, fill = correlation),
#     color = "white"
#   ) +
#   
#   geom_text(
#     data = pairwise_full |>
#       filter(match(Metric1, metric_levels) < match(Metric2, metric_levels)),
#     aes(
#       x = Metric1,
#       y = Metric2,
#       label = paste0(round(correlation, 2), "\n(", best_method, ")")
#     ),
#     size = 2.5
#   ) +
#   
#   # LOWER TRIANGLE (scatter)
#   geom_point(
#     data = df_pairs |>
#       filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
#     aes(
#       x = match(Metric1, metric_levels) + (x_scaled - 0.5) * 0.8,
#       y = match(Metric2, metric_levels) + (y_scaled - 0.5) * 0.8,
#       color = species
#     ),
#     size = 1,
#     alpha = 0.6
#   ) +
#   
#   # SMOOTH
#   geom_smooth(
#     data = df_pairs |>
#       filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
#     aes(
#       x = match(Metric1, metric_levels) + (x_scaled - 0.5) * 0.8,
#       y = match(Metric2, metric_levels) + (y_scaled - 0.5) * 0.8,
#       group = interaction(Definition, Threshold, Metric1, Metric2)
#     ),
#     method = "lm",
#     se = FALSE,
#     color = "black",
#     linewidth = 0.5
#   ) +
#   
#   # BOX AROUND SCATTER
#   geom_tile(
#     data = df_pairs |>
#       filter(match(Metric1, metric_levels) > match(Metric2, metric_levels)),
#     aes(x = Metric1, y = Metric2),
#     fill = NA,
#     color = "black",
#     linewidth = 0.4
#   ) +
#   
#   facet_grid(Definition ~ Threshold) +
#   
#   scale_fill_gradient2(
#     low = "blue", mid = "white", high = "red", midpoint = 0,
#     na.value = "grey90"
#   ) +
#   
#   scale_color_manual(values = c(
#     "Lycaena virgaureae" = "#bf812d",
#     "Thymelicus lineola" = "#8c510a",
#     "Maniola jurtina" = "darkgrey",
#     "Pieris napi" = "#74a9cf",
#     "Gonepteryx rhamni" = "#0570b0",
#     "Aglais urticae" = "#034e7b"
#   )) +
#   
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### CORRELATIONS of FPCA1 to summary stats. 
# 
# 
# 
# df_FPCA1 <- df_wide |>
#   pivot_longer(-species, names_to = "type_threshold", values_to = "FPCA1") |>
#   separate(type_threshold, into = c("Definition", "Threshold"), sep = "_") |>
#   mutate(
#     Definition = factor(Definition, levels = c("broad", "narrow")),
#     Threshold = factor(Threshold, levels = c("High", "Medium", "Low"))
#   )
# 
# 
# 
# 
# 
# # Join FPCA1 with patch metrics
# 
# 
# df_joined <- inner_join(
#   df_FPCA1,
#   stats_clean,
#   by = c("species", "Definition", "Threshold")
# ) 
# 
# 
# 
# 
# ###PLOT THE DATA: 
# 
# df_long <- df_joined |>
#   pivot_longer(
#     cols = c(n_patches, MeanArea, MedianArea, MeanNN, MedianNN),
#     names_to = "Metric",
#     values_to = "Value"
#   )
# 
# 
# ggplot(df_long, aes(x = FPCA1, y = Value)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "loess", se = FALSE, color = "black") +
#   facet_grid(Metric ~ Definition + Threshold, scales = "free_y") +
#   theme_minimal() +
#   labs(
#     x = "FPCA1",
#     y = "Metric value"
#   )
# 
# 
# ggplot(df_long, aes(x = FPCA1, y = Value)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   facet_grid(Metric ~ Definition + Threshold, scales = "free_y") +
#   theme_minimal() +
#   labs(
#     x = "FPCA1",
#     y = "Metric value"
#   )
# 
# 
# library(dplyr)
# library(tidyr)
# 
# fpca_cor <- df_joined %>%
#   group_by(Definition, Threshold) %>%
#   summarise(
#     across(
#       c(n_patches, MedianArea, MeanArea, MedianNN, MeanNN),
#       ~ cor(FPCA1, .x, method = "spearman", use = "pairwise.complete.obs")
#     ),
#     .groups = "drop"
#   )
# 
# 
# 
# library(ggplot2)
# 
# fpca_cor_long <- fpca_cor %>%
#   pivot_longer(
#     -c(Definition, Threshold),
#     names_to = "metric",
#     values_to = "correlation"
#   )
# 
# 
# 
# ### order data 
# 
# library(dplyr)
# 
# 
# metrics <- c("n_patches", "MedianArea", "MeanArea", "MedianNN", "MeanNN")
# 
# fpca_cor_long <- fpca_cor_long %>%
#   mutate(
#     # enforce metric order
#     metric = factor(metric, levels = metrics),
#     
#     # enforce row ordering
#     Definition = factor(Definition, levels = c("narrow", "broad")),
#     Threshold  = factor(Threshold, levels = c("High", "Medium", "Low"))
#   ) %>%
#   arrange(Definition, Threshold)
# 
# 
# 
# 
# 
# 
# 
# ggplot(fpca_cor_long, aes(x = metric, y = 1, fill = correlation)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = round(correlation, 2)), size = 3) +
#   facet_grid(Definition ~ Threshold) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     limits = c(-1, 1)
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Landscape metric",
#     y = NULL,
#     fill = "r^2 (FPCA1)"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )
# 
# 
# 





























##3Prep data for later moddling: 

df_selected_size_metrics <-  df_joined[, c("species", "Definition", "Threshold", "MedianArea", "MeanArea", "Total")]


df_selected_distance_metrics <-  df_joined[, c("species", "Definition", "Threshold", "MedianNN", "MeanNN", "FPCA1", "n_patches")]



all_genomics_indivudial_data_for_Stats <- read_csv("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/Contemporary_Butterflies/Results/Summary_Stats/Individual_based/all_genomics_indivudial_data_for_Stats.csv")



df_selected_genomics <-  all_genomics_indivudial_data_for_Stats[, c("Spp", "ID", "Region", "Het", "Froh_Class_0", "Expressed")]
df_selected_genomics <- df_selected_genomics %>%
  rename(species  = Spp)

#Dispersal 

Dispersal_dimentions <- read_csv("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/DATA/species_data/Dispersal_dimentions.csv")

Dispersal <-  Dispersal_dimentions[, c("Spp", "Nmds_dim1")]

Dispersal <- Dispersal %>%
  rename(species  = Spp)



Dispersal$Nmds_dim1 <- Dispersal$Nmds_dim1*-1



#NE
selected_demograpy_resent <- read_excel("C:/Users/au601132/OneDrive - Aarhus universitet/Skrivebord/PhD/Contemporary_Butterflies/Results/Summary_Stats/Population_based/selected_demograpy_resent.xlsx")

selected_demograpy_resent <-  selected_demograpy_resent[, c("Species", "Ne(t=1)", "Delta_Ne", "Proportion_retained", "1-TD", "Max_Ne", "Current_Ne(t=1)")]


selected_demograpy_resent <- selected_demograpy_resent %>%
  rename(species  = Species)  %>%
  mutate(
    species = str_replace_all(species, "_", " "),
     species = str_to_sentence(species))

#### Build predictors! 

##join landscape data 

df_landscape <- df_selected_distance_metrics %>%
  left_join(df_selected_size_metrics,
            by = c("species", "Definition", "Threshold"))

##Expand to individual level.. eg the species metrics are added: 


### dummy for expressed in T. lineola   OBS remove this from. 

set.seed(42)  # reproducibility

df_selected_genomics <- df_selected_genomics %>%
  mutate(
    Expressed = ifelse(
      species == "Thymelicus lineola" & is.na(Expressed),
      sample(30:50,
             sum(species == "Thymelicus lineola" & is.na(Expressed)),
             replace = TRUE),
      Expressed
    )
  )


df_analysis <- df_selected_genomics %>%
  left_join(Dispersal, by = "species") %>%
  left_join(df_landscape, by = "species") %>%
  left_join(selected_demograpy_resent, by = "species")




df_analysis <- df_analysis %>%
  mutate(
    Definition = factor(Definition),
    Threshold  = factor(Threshold),
    Region     = factor(Region),
    species        = factor(species)
  ) %>%
  drop_na()


df_analysis$FPCA1 <- df_analysis$FPCA1*-1

table(df_analysis$FPCA1)


##Define response matrix

Y <- df_analysis %>%
  select(Het, Froh_Class_0, Expressed)


### building a full RDA 

library(vegan)

rda_full <- rda(
  Y ~ 
    # Distance block
    FPCA1 + MeanNN + MedianNN +
    
    # Size block
    MeanArea + MedianArea + Total +
    
    # Habitat definition (h)
    Definition * Threshold +
    
    # Dispersal
    Nmds_dim1 + Delta_Ne +
    
    # Control
    Condition(Region) +
    Condition(ID),
  
  data = df_analysis
)

table(df_analysis$ID)

anova(rda_full, permutations = 999)
anova(rda_full, by = "terms", permutations = 999)
anova(rda_full, by = "margin", permutations = 999)


RsquareAdj(rda_full)




### now factor species: 

rda_spp <- rda(
  Y ~ FPCA1 + MeanNN + MedianNN +
    MeanArea + MedianArea + Total + 
    species +
    Definition * Threshold +
    Condition(Region),
  data = df_analysis
)

anova(rda_full)
anova(rda_spp)
RsquareAdj(rda_spp)

#STEP 4 — BLOCK TESTING 

## Now partiel RDA

#Distance block

rda_distance <- rda(
  Y ~ FPCA1 + MeanNN + MedianNN +
    Condition(MeanArea + MedianArea + Total +
                Definition + Threshold + Nmds_dim1 + Current_Ne + Region),
  data = df_analysis
)

anova(rda_distance)
RsquareAdj(rda_distance)

#Size block

rda_size <- rda(
  Y ~ MeanArea + MedianArea + Total +
    Condition(FPCA1 + MeanNN + MedianNN +
                Definition + Threshold + Nmds_dim1 + Current_Ne + Region),
  data = df_analysis
)

anova(rda_size)
RsquareAdj(rda_size)


#Dispersal

rda_disp <- rda(
  Y ~ Nmds_dim1 +
    Condition(FPCA1 + MeanNN + MedianNN +
                MeanArea + MedianArea + Total +
                Current_Ne +
                Definition + Threshold + Region),
  data = df_analysis
)

anova(rda_disp)
RsquareAdj(rda_disp)

#Current Ne

rda_Ne <- rda(
  Y ~ Current_Ne +
    Condition(FPCA1 + MeanNN + MedianNN +
                MeanArea + MedianArea + Total +
                Definition + Threshold + Region),
  data = df_analysis
)

anova(rda_disp)
RsquareAdj(rda_disp)


##STEP 5 — Variation partitioning (cleanest ecological result)  
#-- obs can only test between 2-4 "blocks" at a time



X_distance <- df_analysis %>%
  select(FPCA1, MeanNN, MedianNN)

X_size <- df_analysis %>%
  select(MeanArea, MedianArea, Total)

X_disp <- df_analysis %>%
  select(Nmds_dim1)

X_Ne <- df_analysis %>%
  select(Delta_Ne, Max_Ne, `Current_Ne(t=1)`)

X_h <- df_analysis %>%
  select(Definition, Threshold)

varpart_res <- varpart(
  Y,
  X_distance,
  X_size,
  X_disp,
  X_Ne
)


plot(varpart_res)






###Visulisation
# =========================================================
# 0. Libraries
# =========================================================
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(grid)   # needed for arrow()

# =========================================================
# 1. Extract RDA scores (correct scaling)
# =========================================================

scores_sites   <- as.data.frame(scores(rda_full, display = "sites", scaling = 2))
scores_bp      <- as.data.frame(scores(rda_full, display = "bp", scaling = 2))

scores_bp <- scores_bp %>%
  filter(!var %in% c("Definitionnarrow",
                     "ThresholdLow",
                     "ThresholdMedium",
                     "Definitionnarrow:ThresholdLow",
                     "Definitionnarrow:ThresholdMedium"))

scores_species <- as.data.frame(scores(rda_full, display = "species", scaling = 2))

scores_bp$var <- rownames(scores_bp)
scores_species$var <- rownames(scores_species)

# =========================================================
# 2. RDA PLOT (clean, correct geometry)
# =========================================================

p_rda <- ggplot() +
  
  # Individuals
  geom_point(
    data = scores_sites,
    aes(RDA1, RDA2),
    color = "grey40",
    alpha = 0.3,
    size = 1
  ) +
  
  # Environmental arrows (NO rescaling!)
  geom_segment(
    data = scores_bp,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.6,
    color = "#E64B35"
  ) +
  
  # Driver labels (repelled, but NOT rescaled)
  geom_text_repel(
    data = scores_bp,
    aes(RDA1, RDA2, label = var),
    color = "#E64B35",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  # Response variables (genomics)
  geom_text_repel(
    data = scores_species,
    aes(RDA1, RDA2, label = var),
    color = "black",
    fontface = "bold",
    size = 4.5,
    max.overlaps = Inf
  ) +
  
  coord_equal() +
  theme_bw(base_size = 12) +
  
  labs(
    title = "RDA: genomic variation vs landscape drivers",
    x = "RDA1",
    y = "RDA2"
  )

# =========================================================
# 3. VARIANCE PARTITIONING (Venn diagram, fixed)
# =========================================================

par(mar = c(2,2,2,2))

plot(
  varpart_res,
  bg = c("#4DBBD580", "#00A08780", "#E64B3580", "#3C548880"),
  Xnames = c("Distance", "Size", "Dispersal", "Defination (MxQ)"),
  cex = 1.1,
  digits = 2
)

title("Variance partitioning", line = 0.5)

# Capture plot
vp_plot <- recordPlot()

# =========================================================
# 4. Convert varpart to patchwork-compatible object
# =========================================================

p_var <- wrap_elements(
  panel = ~{
    par(mar = c(3,3,3,3))
    plot(
      varpart_res,
      bg = c("#4DBBD580", "#00A08780", "#E64B3580","yellow4", "#3C548880"),
      Xnames = c("Distance", "Size", "Dispersal", "Ne", "Defination (MxQ)"),
      cex = 1.1,
      digits = 2
    )
  }
)

# =========================================================
# 5. Combine panels
# =========================================================

final_plot <- p_var + p_rda + plot_layout(widths = c(1, 2))

final_plot



