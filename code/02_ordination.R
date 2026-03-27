################################################################################
# Script: 02_ordination.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-03-27
# Description: 
#   This script performs multivariate analyses of benthic community structure:
#   1. Load benthic complete data (image-level)
#   2. Aggregate to island × year level (mean cover per group)
#   3. Hellinger transformation of benthic cover data
#   4. Calculation of Gower distance matrix
#   5. PERMANOVA to test differences among islands
#   6. Beta dispersion analysis (homogeneity of variances)
#   7. Principal Coordinates Analysis (PCoA) ordination
#   8. envfit to project benthic groups onto ordination space
#   9. Indicator species analysis
# Data structure: each point = island × year combination
# Outputs:
#   - results/figures/Figure_2B_PCoA_ordination.png
#   - results/figures/Figure_2B_PCoA_ordination.tiff
#   - results/tables/Table_1_PERMANOVA_results.csv
#   - results/tables/Table_2_envfit_results.csv
#   - results/tables/Table_3_indicator_species.csv
################################################################################

# 1. Load packages ####
required_packages <- c("vegan", "ggplot2", "tidyverse", "dplyr", "ggrepel", "indicspecies", "tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
options(scipen = 999)

# 2. Load benthic complete data ####
complete_data <- read.csv("data/raw/benthic_complete_data.csv")

# Inspect data
glimpse(complete_data)

# 3. Recode fine-scale categories into benthic groups (same as Script 01) ####
grouped_data <- complete_data %>%
  mutate(
    categoryid = as.character(categoryid),
    group = recode(categoryid, 
                   # Invertebrates
                   "BRY" = "INV", "ECH" = "INV", "OTH" = "INV",
                   # Macroalgae
                   "CEN" = "MAL", "COR" = "MAL", "COT" = "MAL",
                   "FIL" = "MAL", "FOL" = "MAL", "GLC" = "MAL",
                   "SAR" = "MAL", "STO" = "MAL",
                   # Turf/Epilithic Algal Matrix
                   "TUR" = "EAM",
                   # Scleractinian corals
                   "BSC" = "SCL", "ODI" = "SCL",
                   # Suspensivores
                   "TUN" = "SUS", "POR" = "SUS",
                   # Default: keep original categoryid
                   .default = categoryid
    )
  ) %>%
  group_by(island, year, sites, transect, image, group) %>%
  summarise(
    cover_per_group = sum(coverpercategory, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Aggregate to island × year level ####
# Calculate mean cover per group for each island and year
biotic_yearly <- grouped_data %>%
  group_by(island, year, group) %>%
  summarise(
    mean_cover = mean(cover_per_group, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Standardize island names
  mutate(island = recode(island,
                         "trindade" = "TR",
                         "noronha" = "FN",
                         "rocas" = "RA",
                         "stpauls_rocks" = "SP"))

# Check unique combinations
cat("Number of island × year combinations:", 
    n_distinct(paste(biotic_yearly$island, biotic_yearly$year)), "\n")

print("Available combinations:")
biotic_yearly %>%
  distinct(island, year) %>%
  arrange(island, year) %>%
  print()

# 5. Convert to wide format (sites = island_year, columns = benthic groups) ####
biotic_wide <- biotic_yearly %>%
  select(island, year, group, mean_cover) %>%
  pivot_wider(names_from = group, values_from = mean_cover, values_fill = 0) %>%
  mutate(site = paste(island, year, sep = "_")) %>%
  column_to_rownames("site") %>%
  select(-island, -year)

# Create metadata for sites
biotic_meta <- biotic_yearly %>%
  distinct(island, year) %>%
  mutate(site = paste(island, year, sep = "_")) %>%
  arrange(site)

# Define group order for later
group_order <- c("EAM", "MAL", "CCA", "ACA", "SCL", "ABI", "INV", "ZOA", "CYA", "SUS")
# Ensure all groups are present, add missing columns if necessary
for (grp in group_order) {
  if (!grp %in% colnames(biotic_wide)) {
    biotic_wide[[grp]] <- 0
  }
}
# Reorder columns
biotic_wide <- biotic_wide[, group_order]

# 6. Hellinger transformation ####
biotic_hell <- decostand(biotic_wide, method = "hellinger")

# 7. Calculate distance matrix ####
# Gower distance (identified as best in rankindex analysis)
dist_matrix <- vegdist(biotic_hell, method = "gower")

# 8. PERMANOVA ####
permanova <- adonis2(dist_matrix ~ biotic_meta$island, permutations = 999)
print(permanova)

# Save PERMANOVA results
permanova_results <- data.frame(
  Source = rownames(permanova),
  DF = permanova$Df,
  SumOfSqs = permanova$SumOfSqs,
  R2 = permanova$R2,
  F = permanova$F,
  P = permanova$`Pr(>F)`
)
write_csv(permanova_results, "results/tables/Table_1_PERMANOVA_results.csv")

# 9. Beta dispersion (homogeneity of variances) ####
beta_disp <- betadisper(dist_matrix, group = biotic_meta$island)
anova_beta <- anova(beta_disp)
print(anova_beta)

# 10. PCoA ordination ####
pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
pcoa_sites <- data.frame(pcoa$points[, 1:2])
colnames(pcoa_sites) <- c("PCoA1", "PCoA2")
pcoa_sites$island <- biotic_meta$island
pcoa_sites$year <- biotic_meta$year
pcoa_sites$site <- biotic_meta$site

# Define island order
pcoa_sites$island <- factor(pcoa_sites$island, levels = c("SP", "RA", "FN", "TR"))

# Variance explained
var_exp1 <- round(pcoa$eig[1] / sum(pcoa$eig) * 100, 1)
var_exp2 <- round(pcoa$eig[2] / sum(pcoa$eig) * 100, 1)

cat("PCoA variance explained:", var_exp1, "% and", var_exp2, "%\n")

# 11. envfit: project benthic groups onto PCoA ####
env_fit <- envfit(pcoa, biotic_hell, permutations = 999)
envfit_vectors <- as.data.frame(scores(env_fit, display = "vectors"))
envfit_vectors$species <- rownames(envfit_vectors)
envfit_vectors$r <- env_fit$vectors$r
envfit_vectors$p <- env_fit$vectors$pvals

# Select significant vectors (p < 0.05)
sig_vectors <- envfit_vectors %>%
  filter(p < 0.05) %>%
  arrange(desc(r))

print("Significant benthic groups (envfit):")
print(sig_vectors)

# Save envfit results
write_csv(envfit_vectors, "results/tables/Table_2_envfit_results.csv")

# 12. Indicator species analysis ####
# Prepare data for indicator analysis (island × year level)
indicator_data <- biotic_wide %>%
  rownames_to_column("site") %>%
  left_join(biotic_meta, by = "site")

# Calculate indicator values
indicator_values <- multipatt(indicator_data[, group_order], 
                              cluster = indicator_data$island, 
                              func = "r.g",
                              control = how(nperm = 999))

# Extract results
indicator_summary <- data.frame(
  group = rownames(indicator_values$sign),
  stat = indicator_values$sign$stat,
  p.value = indicator_values$sign$p.value,
  island = apply(indicator_values$sign[, 1:4], 1, 
                 function(x) colnames(indicator_values$sign)[which.max(x)])
) %>%
  filter(p.value < 0.05) %>%
  arrange(desc(stat))

print("Indicator species results:")
print(indicator_summary)

write_csv(indicator_summary, "results/tables/Table_3_indicator_species.csv")

# 13. Create PCoA plot ####
# Define colors and shapes
island_colors <- c("SP" = "deeppink2", "RA" = "chocolate1", "FN" = "blue1", "TR" = "forestgreen")
island_shapes <- c("SP" = 15, "RA" = 16, "FN" = 17, "TR" = 18)

pcoa_plot <- ggplot() +
  # Sites (points by island, colored by island)
  geom_point(data = pcoa_sites, 
             aes(x = PCoA1, y = PCoA2, shape = island, color = island), 
             size = 3, stroke = 1, alpha = 0.9) +
  scale_shape_manual(values = island_shapes) +
  scale_color_manual(values = island_colors) +
  
  # Add year labels
  geom_text_repel(data = pcoa_sites, 
                  aes(x = PCoA1, y = PCoA2, label = year),
                  size = 2.5, color = "gray30", max.overlaps = 10) +
  
  # Significant benthic vectors
  geom_segment(data = sig_vectors, 
               aes(x = 0, xend = Dim1 * 0.25, y = 0, yend = Dim2 * 0.25),
               color = "grey50", linewidth = 0.7, linetype = "dashed") +
  geom_text_repel(data = sig_vectors, 
                  aes(x = Dim1 * 0.27, y = Dim2 * 0.27, label = species),
                  color = "grey30", size = 3.5, fontface = "bold",
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0, segment.color = "grey50") +
  
  # PERMANOVA annotation
  annotate("text", x = -0.2, y = 0.25, 
           label = sprintf("PERMANOVA: F = %.1f, R² = %.2f, p = %.3f", 
                           permanova$F[1], permanova$R2[1], permanova$`Pr(>F)`[1]),
           color = "black", size = 3, fontface = "bold") +
  
  # Axes
  labs(x = paste0("PCoA 1 (", var_exp1, "%)"),
       y = paste0("PCoA 2 (", var_exp2, "%)")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  
  # Theme
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.9, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(pcoa_plot)

# 14. Save PCoA plot ####
ggsave("results/figures/Figure_2B_PCoA_ordination.png",
       plot = pcoa_plot, width = 6, height = 6, dpi = 300)

ggsave("results/figures/Figure_2B_PCoA_ordination.tiff",
       plot = pcoa_plot, width = 6, height = 6, dpi = 300, compression = "lzw")

# 15. Print session info for reproducibility ####
sessionInfo()

################################################################################
# End of script
################################################################################