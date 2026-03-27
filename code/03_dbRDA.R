################################################################################
# Script: 03_dbRDA.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-03-27
# Description: 
#   Distance-based Redundancy Analysis (dbRDA) to identify environmental drivers
#   of benthic community structure across islands.
#   Uses Bray-Curtis distance (selected via rankindex) and all environmental
#   variables with VIF < 10 (wave, SST, PAR, POC).
#   Model order: wave first (ecologically meaningful), followed by SST, PAR, POC.
# Outputs:
#   - results/figures/Figure_3_dbRDA_ordination.png
#   - results/figures/Figure_3_dbRDA_ordination.tiff
#   - results/tables/Table_4_dbRDA_summary.csv
#   - results/tables/Table_5_dbRDA_axes_summary.csv
#   - results/tables/Table_6_env_variables_significance_sequential.csv
#   - results/tables/Table_6_env_variables_significance_marginal.csv
#   - results/tables/Table_rankindex_results.csv
#   - results/tables/Table_VIF_results.csv
################################################################################

# 1. Load packages ####
required_packages <- c("vegan", "ggplot2", "tidyverse", "dplyr", "ggrepel", "tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
options(scipen = 999)

# 2. Load benthic complete data ####
complete_data <- read.csv("data/raw/benthic_complete_data.csv")

# 3. Load environmental data ####
enviro_raw <- read.csv("data/raw/environment.csv")

# 4. STANDARDIZE ISLAND NAMES ####
complete_data <- complete_data %>%
  mutate(island = recode(island,
                         "trindade" = "TR",
                         "noronha" = "FN",
                         "rocas" = "RA",
                         "stpauls_rocks" = "SP"))

enviro_raw <- enviro_raw %>%
  mutate(island = recode(island,
                         "trindade" = "TR",
                         "noronha" = "FN",
                         "rocas" = "RA",
                         "spspa" = "SP"))

# 5. Recode benthic groups ####
grouped_data <- complete_data %>%
  mutate(
    categoryid = as.character(categoryid),
    group = recode(categoryid, 
                   "BRY" = "INV", "ECH" = "INV", "OTH" = "INV",
                   "CEN" = "MAL", "COR" = "MAL", "COT" = "MAL",
                   "FIL" = "MAL", "FOL" = "MAL", "GLC" = "MAL",
                   "SAR" = "MAL", "STO" = "MAL",
                   "TUR" = "EAM",
                   "BSC" = "SCL", "ODI" = "SCL",
                   "TUN" = "SUS", "POR" = "SUS",
                   .default = categoryid
    )
  ) %>%
  group_by(island, year, group) %>%
  summarise(
    mean_cover = mean(coverpercategory, na.rm = TRUE),
    .groups = "drop"
  )

# 6. Convert benthic data to wide format ####
biotic_wide <- grouped_data %>%
  pivot_wider(names_from = group, values_from = mean_cover, values_fill = 0) %>%
  mutate(site = paste(island, year, sep = "_"))

# Define group order
group_order <- c("EAM", "MAL", "CCA", "ACA", "SCL", "ABI", "INV", "ZOA", "CYA", "SUS")

# Ensure all groups are present
for (grp in group_order) {
  if (!grp %in% colnames(biotic_wide)) {
    biotic_wide[[grp]] <- 0
  }
}

# 7. Prepare environmental data ####
enviro_clean <- enviro_raw %>%
  mutate(site = paste(island, year, sep = "_"))

# 8. Merge datasets ####
df <- biotic_wide %>%
  left_join(enviro_clean, by = c("island", "year", "site"))

# 9. Remove rows with missing data ####
df_clean <- df %>%
  drop_na(SST, PAR, wave, POC)

cat("\n========================================\n")
cat("Data Summary\n")
cat("========================================\n")
cat("Total samples (island × year):", nrow(df_clean), "\n")
df_clean %>%
  group_by(island) %>%
  summarise(n_years = n(), years = paste(sort(year), collapse = ", ")) %>%
  print()

# 10. Create matrices for analysis ####
biotic_mat <- df_clean %>%
  select(all_of(group_order)) %>%
  as.matrix()
rownames(biotic_mat) <- df_clean$site
# Select environmental variables (all with VIF < 10)
env_vars <- c("wave", "SST", "PAR", "POC")
enviro_mat <- df_clean %>%
  select(all_of(env_vars)) %>%
  as.matrix()
rownames(enviro_mat) <- df_clean$site
cat("\nBiotic matrix dimensions:", dim(biotic_mat), "\n")
cat("Environmental matrix dimensions:", dim(enviro_mat), "\n")
cat("\nNames row in biotic_mat:\n")
head(rownames(biotic_mat))

cat("\nNames row in enviro_mat:\n")
head(rownames(enviro_mat))

# 11. Transform data ####
biotic_hell <- decostand(biotic_mat, method = "hellinger")
enviro_std <- decostand(enviro_mat, method = "standardize")
cat("\nNames row in biotic_hell after transformation:\n")
head(rownames(biotic_hell))

# 12. SELECT BEST DISTANCE METRIC USING RANKINDEX ####
cat("\n========================================\n")
cat("Selecting Best Distance Metric\n")
cat("========================================\n")

distance_metrics <- c("euclidean", "manhattan", "gower", "bray", "kulczynski")
rank_scores <- rankindex(enviro_std, biotic_hell, 
                         indices = distance_metrics, 
                         stepacross = FALSE, 
                         method = "spearman")

rank_df <- data.frame(
  Distance = distance_metrics,
  Score = round(rank_scores, 4)
) %>% arrange(desc(Score))

print(rank_df)
write_csv(rank_df, "results/tables/Table_rankindex_results.csv")

best_dist <- distance_metrics[which.max(rank_scores)]
cat("\n✓ Selected distance metric:", best_dist, "\n")

# 13. Check VIF ####
cat("\n========================================\n")
cat("Variance Inflation Factors (VIF)\n")
cat("========================================\n")

cor_matrix <- cor(enviro_std)
vif_values <- diag(solve(cor_matrix))
vif_df <- data.frame(
  Variable = env_vars,
  VIF = round(vif_values, 2)
)
print(vif_df)
write_csv(vif_df, "results/tables/Table_VIF_results.csv")

# 14. Run dbRDA with wave first ####
cat("\n========================================\n")
cat("Running dbRDA - Wave first (wave → SST → PAR → POC)\n")
cat("========================================\n")

dbRDA <- dbrda(biotic_hell ~ wave + SST + PAR + POC, 
               data = as.data.frame(enviro_std), 
               distance = best_dist)

# 15. Model summary ####
r2_adj <- RsquareAdj(dbRDA)$adj.r.squared
cat("\nAdjusted R²:", round(r2_adj, 4), "\n")

# Global test
global_test <- anova(dbRDA, permutations = 999)
cat("\nGlobal dbRDA test:\n")
print(global_test)

# Sequential test (wave first)
terms_test_seq <- anova(dbRDA, by = "terms", permutations = 999)
cat("\nTerms significance (sequential - wave first):\n")
print(terms_test_seq)

# Marginal test (unique effect of each variable)
terms_test_margin <- anova(dbRDA, by = "margin", permutations = 999)
cat("\nTerms significance (marginal - unique effect):\n")
print(terms_test_margin)

# 16. Extract eigenvalues and variance explained ####
all_eig <- eigenvals(dbRDA)
constrained_eig <- all_eig[1:dbRDA$CCA$rank]
total_inertia <- sum(all_eig[all_eig > 0])
constrained_var <- constrained_eig / total_inertia * 100

cat("\nConstrained axes variance explained:\n")
for (i in 1:length(constrained_var)) {
  cat(sprintf("  dbRDA%d: %.1f%% (eigenvalue = %.4f)\n", 
              i, constrained_var[i], constrained_eig[i]))
}

# 17. Save summary tables ####
# Global summary
dbRDA_summary <- data.frame(
  Parameter = c("Distance", "Adjusted R²", "Global F", "Global p"),
  Value = c(best_dist, r2_adj, global_test$F[1], global_test$`Pr(>F)`[1])
)
write_csv(dbRDA_summary, "results/tables/Table_4_dbRDA_summary.csv")

# Axes summary
axes_summary <- data.frame(
  Axis = paste0("dbRDA", 1:length(constrained_eig)),
  Eigenvalue = constrained_eig,
  Proportion = constrained_var / 100,
  Cumulative = cumsum(constrained_var) / 100
)
write_csv(axes_summary, "results/tables/Table_5_dbRDA_axes_summary.csv")

# Sequential test results (wave first)
terms_seq_summary <- data.frame(
  Variable = rownames(terms_test_seq),
  DF = terms_test_seq$Df,
  SumOfSqs = terms_test_seq$SumOfSqs,
  F = terms_test_seq$F,
  P_sequential = terms_test_seq$`Pr(>F)`
)
write_csv(terms_seq_summary, "results/tables/Table_6_env_variables_significance_sequential.csv")

# Marginal test results
terms_margin_summary <- data.frame(
  Variable = rownames(terms_test_margin),
  DF = terms_test_margin$Df,
  SumOfSqs = terms_test_margin$SumOfSqs,
  F = terms_test_margin$F,
  P_marginal = terms_test_margin$`Pr(>F)`
)
write_csv(terms_margin_summary, "results/tables/Table_6_env_variables_significance_marginal.csv")

cat("\n✓ Tables saved successfully!\n")

# 18. Extract scores for plotting ####
# Site scores
sites_scores <- as.data.frame(scores(dbRDA, display = "sites")[, 1:2])
sites_scores$site <- rownames(sites_scores)
colnames(sites_scores)[1:2] <- c("dbRDA1", "dbRDA2")
sites_scores <- sites_scores[, c("site", "dbRDA1", "dbRDA2")]

# Extract island and year from site names
sites_scores <- sites_scores %>%
  separate(site, into = c("island", "year"), sep = "_", remove = FALSE, 
           extra = "merge") %>%
  mutate(island = factor(island, levels = c("SP", "RA", "FN", "TR")))

cat("\nSite scores:\n")
print(head(sites_scores))

# Species scores (benthic groups)
site_scores_axes <- scores(dbRDA, display = "sites")[, 1:2]
species_cor <- cor(biotic_hell, site_scores_axes)

species_scores <- as.data.frame(species_cor)
colnames(species_scores) <- c("dbRDA1", "dbRDA2")
species_scores$group <- rownames(species_scores)
species_scores <- species_scores %>%
  filter(group %in% group_order) %>%
  select(group, dbRDA1, dbRDA2)

species_scores <- species_scores %>%
  mutate(
    dbRDA1 = dbRDA1 * 0.8,
    dbRDA2 = dbRDA2 * 0.8
  )
cat("\nSpecies scores (correlações):\n")
print(species_scores)

# Environmental vector scores
env_scores <- as.data.frame(scores(dbRDA, display = "bp")[, 1:2])
env_scores$variable <- rownames(env_scores)
colnames(env_scores)[1:2] <- c("dbRDA1", "dbRDA2")
env_scores <- env_scores[, c("variable", "dbRDA1", "dbRDA2")]

cat("\nEnvironmental scores:\n")
print(env_scores)


# 19. Calculate scaling factor for vectors ####
max_site <- max(abs(sites_scores[, c("dbRDA1", "dbRDA2")]), na.rm = TRUE)
max_vec <- max(abs(env_scores[, c("dbRDA1", "dbRDA2")]), na.rm = TRUE)
scale_factor <- max_site / max_vec * 0.8

env_scores_scaled <- env_scores %>%
  mutate(
    dbRDA1_scaled = dbRDA1 * scale_factor,
    dbRDA2_scaled = dbRDA2 * scale_factor
  )

cat("\nScaling factor:", round(scale_factor, 3), "\n")
cat("Scaled environmental vectors:\n")
print(env_scores_scaled)

# 20. Create dbRDA plot ####
island_colors <- c("SP" = "deeppink2", "RA" = "chocolate1", "FN" = "blue1", "TR" = "forestgreen")
island_shapes <- c("SP" = 15, "RA" = 16, "FN" = 17, "TR" = 18)

dbRDA_plot <- ggplot() +
  # Sites
  geom_point(data = sites_scores, 
             aes(x = dbRDA1, y = dbRDA2, shape = island, color = island), 
             size = 3, stroke = 1, alpha = 0.9) +
  scale_shape_manual(values = island_shapes) +
  scale_color_manual(values = island_colors) +
  
  # Environmental vectors (scaled)
  geom_segment(data = env_scores_scaled, 
               aes(x = 0, xend = dbRDA1_scaled, y = 0, yend = dbRDA2_scaled),
               color = "grey50", linewidth = 0.8,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_text_repel(data = env_scores_scaled, 
                  aes(x = dbRDA1_scaled * 1.05, y = dbRDA2_scaled * 1.05, 
                      label = variable),
                  color = "grey30", size = 3.5, fontface = "bold") +
  
  # Benthic groups (species scores)
  geom_text_repel(data = species_scores, 
                  aes(x = dbRDA1, y = dbRDA2, label = group),
                  color = "firebrick", size = 3, fontface = "bold") +
  
  # Annotation with model statistics
  annotate("text", x = -1.2, y = 1.3, 
           label = sprintf("dbRDA: %s\nR²Adj. = %.2f\nGlobal p = %.3f", 
                           best_dist, r2_adj, global_test$`Pr(>F)`[1]),
           color = "black", size = 3, fontface = "bold", hjust = 0) +
  
  # Axes
  labs(x = paste0("dbRDA 1 (", round(constrained_var[1], 1), "%)"),
       y = paste0("dbRDA 2 (", round(constrained_var[2], 1), "%)")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  
  # Theme
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.1, 0.15),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(dbRDA_plot)

# 21. Save plots ####
ggsave("results/figures/Figure_3_dbRDA_ordination.png",
       plot = dbRDA_plot, width = 8, height = 6, dpi = 300)

ggsave("results/figures/Figure_3_dbRDA_ordination.tiff",
       plot = dbRDA_plot, width = 8, height = 6, dpi = 300, compression = "lzw")

cat("\n✓ Plot saved: results/figures/Figure_3_dbRDA_ordination.png\n")
cat("✓ Plot saved: results/figures/Figure_3_dbRDA_ordination.tiff\n")

# 22. Print session info for reproducibility ####
sessionInfo()

################################################################################
# End of script
################################################################################

