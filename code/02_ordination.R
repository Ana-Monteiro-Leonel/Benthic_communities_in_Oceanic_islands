##############################################################################
# Script: 02_ordination.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-09-03

#   
# Description:
#   Multivariate analyses of benthic community structure using
#   18 benthic categories (17 biological categories + abiotic substrate).
#
#   1. Load image-level benthic cover data
#   2. Recode original category codes into the 18-category dataset
#   3. Aggregate benthic cover to transect level
#   4. Apply Hellinger transformation
#   5. Calculate Bray-Curtis dissimilarities
#   6. PERMANOVA: island, year, and island × year
#   7. Test multivariate dispersion among islands
#   8. Principal Coordinates Analysis (PCoA)
#   9. envfit of benthic categories onto the PCoA
#  10. Indicator analysis of benthic categories among islands
#
# Data structure:
#   Sampling unit = transect
#   Images = subsamples within transects
#   Transects are associated with sites, islands, and survey years
#
# Outputs:
#   - results/figures/Figure_2B_PCoA_ordination.png
#   - results/figures/Figure_2B_PCoA_ordination.tiff
#   - results/tables/Table_1_PERMANOVA_results.csv
#   - results/tables/Table_2_envfit_results.csv
#   - results/tables/Table_3_indicator_analysis.csv
#   - data/processed/ordination_objects.RData
################################################################################

# Set working directory to project root ####
# This script tries to find the project root automatically.
# If it fails, adjust the path below to your local setup.
project_root <- "C:/Users/Ana Monteiro/OneDrive/Documentos/GitHub/Benthic_communities_in_Oceanic_islands"

if (dir.exists(project_root)) {
  setwd(project_root)
  cat("Working directory set to:", getwd(), "\n")
} else {
  # Try to find project root by looking for data/raw directory
  test_dir <- getwd()
  found <- FALSE
  for (i in 1:5) {
    if (file.exists(file.path(test_dir, "data/raw/benthic_complete_data.csv"))) {
      setwd(test_dir)
      found <- TRUE
      cat("Working directory automatically set to:", getwd(), "\n")
      break
    }
    test_dir <- dirname(test_dir)
  }
  if (!found) {
    stop("Could not find project root. Please set 'project_root' manually.\n",
         "Current working directory: ", getwd(), "\n",
         "Expected project path: C:/Users/Ana Monteiro/OneDrive/Documentos/GitHub/Benthic_communities_in_Oceanic_islands")
  }
}

# Verify data file exists
if (!file.exists("data/raw/benthic_complete_data.csv")) {
  stop("benthic_complete_data.csv not found in data/raw/ directory. 
       Please check your working directory.")
}

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

# 3. Recode original codes into the 18 benthic categories ####
grouped_data <- complete_data %>%
  mutate(
    categoryid = as.character(categoryid),
    group = recode(categoryid, 
                   # Invertebrates
                   "OTH" = "INV", "BRY" = "INV", 
                   # Turf/Epilithic Algal Matrix
                   "TUR" = "EAM",
                   # Scleractinian corals
                   "BSC" = "SCL", "ODI" = "SCL",
                   # Sargassum
                   "SAR" = "COR",
                   # Default: keep original categoryid
                   .default = categoryid
    )
  ) %>%
  group_by(island, year, sites, transect, image, group) %>%
  summarise(
    cover_per_group = sum(coverpercategory, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Aggregate to transect level ####
biotic_transect <- grouped_data %>%
  group_by(island, year, sites, transect, group) %>%
  summarise(
    mean_cover = mean(cover_per_group, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    island = recode(island,
                    "trindade" = "TR",
                    "noronha" = "FN",
                    "rocas" = "RA",
                    "stpauls_rocks" = "SP")
  )

cat("Number of transects:",
    n_distinct(paste(biotic_transect$island,
                     biotic_transect$year,
                     biotic_transect$sites,
                     biotic_transect$transect)),
    "\n")

# 5. Convert to wide format (transects × benthic categories) ####
biotic_wide <- biotic_transect %>%
  select(island, year, sites, transect, group, mean_cover) %>%
  pivot_wider(
    names_from = group,
    values_from = mean_cover,
    values_fill = 0
  ) %>%
  mutate(
    sample_id = paste(island,
                      year,
                      sites,
                      transect,
                      sep = "_")
  )

biotic_meta <- biotic_wide %>%
  select(sample_id,
         island,
         year,
         sites,
         transect)

community_matrix <- biotic_wide %>%
  column_to_rownames("sample_id") %>%
  select(-island,
         -year,
         -sites,
         -transect)

# Expected fine-scale benthic categories
category_order <- c(
  "ABI", "ACA", "CCA", "CEN", "COR", "COT",
  "CYA", "EAM", "ECH", "FIL", "FOL", "GLC",
  "INV", "POR", "SCL", "STO", "TUN", "ZOA"
)

# Check for unexpected categories
unexpected_categories <- setdiff(
  colnames(community_matrix),
  category_order
)
if (length(unexpected_categories) > 0) {
  stop(
    "Unexpected benthic categories found: ",
    paste(unexpected_categories, collapse = ", ")
  )
}

# Add categories absent from a particular dataset, if necessary
for (cat in category_order) {
  if (!cat %in% colnames(community_matrix)) {
    community_matrix[[cat]] <- 0
  }
}

# Keep the 18 categories in a fixed order
community_matrix <- community_matrix[, category_order]

cat(
  "Benthic categories used:",
  paste(colnames(community_matrix), collapse = ", "),
  "\nNumber of categories:", ncol(community_matrix), "\n"
)

# 6. Hellinger transformation ####
biotic_hell <- decostand(community_matrix, method = "hellinger")

# 7. Calculate distance matrix ####
# Bray-Curtis distance on Hellinger-transformed data
dist_matrix <- vegdist(biotic_hell, method = "bray")

# Save ordination objects for reproducibility
save(
  community_matrix,
  biotic_meta,
  biotic_hell,
  dist_matrix,
  file = "data/processed/ordination_objects.RData"
)

# 8. PERMANOVA ####
biotic_meta$island <- factor(biotic_meta$island, 
                             levels = c("SP", "RA", "FN", "TR"))
biotic_meta$year <- factor(biotic_meta$year)

# PERMANOVA: spatial effects - Do communities differ among islands?
perm_island <- adonis2(
  dist_matrix ~ island,
  data = biotic_meta,
  permutations = 999
)
print(perm_island)

# PERMANOVA: temporal effects - Are there differences among years?
perm_year <- adonis2(
  dist_matrix ~ year,
  data = biotic_meta,
  permutations = 999
)
print(perm_year)

# PERMANOVA: interaction island/year
perm_interaction <- adonis2(
  dist_matrix ~ island * year,
  data = biotic_meta,
  permutations = 999,
  by = "terms"
)
print(perm_interaction)

# Save PERMANOVA results
permanova_results <- data.frame(
  Source = rownames(perm_interaction),
  DF = perm_interaction$Df,
  SumOfSqs = perm_interaction$SumOfSqs,
  R2 = perm_interaction$R2,
  F = perm_interaction$F,
  P = perm_interaction$`Pr(>F)`
)
write_csv(permanova_results, "results/tables/Table_1_PERMANOVA_results.csv")

# 9. Beta dispersion (homogeneity of variances) ####
beta_disp <- betadisper(dist_matrix, group = biotic_meta$island, add = TRUE)
anova_beta <- anova(beta_disp)
print(anova_beta)

# 10. PCoA ordination ####
pcoa <- wcmdscale(dist_matrix, eig = TRUE, add = TRUE, k = 2)
pcoa_sites <- data.frame(
  pcoa$points[, 1:2]
)
colnames(pcoa_sites) <- c("PCoA1", "PCoA2")

pcoa_sites <- cbind(
  pcoa_sites,
  biotic_meta
)

# Define island order
pcoa_sites$island <- factor(pcoa_sites$island, levels = c("SP", "RA", "FN", "TR"))

# Variance explained (using only positive eigenvalues)
positive_eig <- pcoa$eig[pcoa$eig > 0]
var_exp1 <- round(pcoa$eig[1] / sum(positive_eig) * 100, 1)
var_exp2 <- round(pcoa$eig[2] / sum(positive_eig) * 100, 1)

cat("PCoA variance explained:", var_exp1, "% and", var_exp2, "%\n")

# 11. envfit: project benthic categories  onto PCoA ####
env_fit <- envfit(
  pcoa_sites[, c("PCoA1", "PCoA2")],
  community_matrix,
  permutations = 999
)
envfit_vectors <- as.data.frame(scores(env_fit, display = "vectors"))
envfit_vectors$category <- rownames(envfit_vectors)
envfit_vectors$r2 <- env_fit$vectors$r
envfit_vectors$p <- env_fit$vectors$pvals

# Select significant vectors (p < 0.05)
sig_vectors <- envfit_vectors %>%
  filter(p < 0.05) %>%
  arrange(desc(r2))

# Strongest significant vectors displayed in Figure 2B
plot_vectors <- sig_vectors %>%
  filter(r2 >= 0.20)

print("Significant benthic categories (envfit):")
print(sig_vectors)
print(plot_vectors)

# Save envfit results
write_csv(envfit_vectors, "results/tables/Table_2_envfit_results.csv")

# 12. Indicator analysis ####
# Prepare data for indicator analysis (transect level grouped by island)
indicator_data <- community_matrix %>%
  rownames_to_column("sample_id") %>%
  left_join(biotic_meta, by = "sample_id")

# Multilevel pattern analysis allowing associations
# with individual islands or combinations of islands
set.seed(123)
indicator_values <- multipatt(indicator_data[, category_order], 
                              cluster = indicator_data$island, 
                              func = "r.g",
                              duleg = FALSE,
                              control = how(nperm = 999))

# Extract indicator analysis results
island_cols <- c("s.FN", "s.RA", "s.SP", "s.TR")
indicator_summary <- indicator_values$sign %>%
  as.data.frame() %>%
  mutate(
    category = rownames(.),
    association = apply(
      .[, island_cols],
      1,
      function(x) {
        islands <- sub("s\\.", "", island_cols[x == 1])
        paste(islands, collapse = " + ")
      }
    )
  ) %>%
  select(
    category,
    association,
    IndVal = stat,
    p.value
  ) %>%
  filter(p.value < 0.05) %>%
  arrange(association, desc(IndVal))
print("Indicator analysis results:")
print(indicator_summary)

write_csv(indicator_summary, "results/tables/Table_3_indicator_analysis.csv")

# 13. Create PCoA plot ####
# Define colors and shapes
island_colors <- c("SP" = "deeppink2", "RA" = "chocolate1", "FN" = "blue1", "TR" = "forestgreen")
island_shapes <- c("SP" = 15, "RA" = 16, "FN" = 17, "TR" = 18)

#### Adjust SCL label position for readability ####
plot_vectors_main <- plot_vectors %>%
  filter(category != "SCL")

plot_vector_scl <- plot_vectors %>%
  filter(category == "SCL")

pcoa_plot <- ggplot() +
  # Sites (points by island, colored by island)
  geom_point(data = pcoa_sites, 
             aes(x = PCoA1, y = PCoA2, shape = island, color = island), 
             size = 3, stroke = 1, alpha = 0.9) +
  scale_shape_manual(values = island_shapes) +
  scale_color_manual(values = island_colors) +
  
  # Significant benthic vectors
  geom_segment(data = plot_vectors, 
               aes(x = 0,
                   xend = PCoA1 * 0.30,
                   y = 0,
                   yend = PCoA2 * 0.30),
               color = "grey50",
               linewidth = 0.7,
               linetype = "dashed",
               arrow = arrow(length = grid::unit(0.15, "cm"))
) +
  
  geom_text_repel(data = plot_vectors_main, 
                  aes(x = PCoA1 * 0.33,
                      y = PCoA2 * 0.33,
                      label = category),
                  color = "grey30",
                  size = 3.5,
                  fontface = "bold",
                  box.padding = 0.6,
                  point.padding = 0.4,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  max.overlaps = Inf
                  ) +
    #Adjust SCL label position for readability
    # Manually reposition SCL label to avoid overlap with data points;
    # vector coordinates remain unchanged.
                  geom_text(
                  data = plot_vector_scl,
                  aes(
                  x = PCoA1 * 0.33 + 0.035,
                  y = PCoA2 * 0.33 - 0.015,
                  label = category
                  ),
                  color = "grey30",
                  size = 3.5,
                  fontface = "bold"
                  ) +
  
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
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect( color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(pcoa_plot)

# 14. Save PCoA plot ####
ggsave("results/figures/Figure_2B_PCoA_ordination.png",
       plot = pcoa_plot, width = 7, height = 6, dpi = 300)

ggsave("results/figures/Figure_2B_PCoA_ordination.tiff",
       plot = pcoa_plot, width = 7, height = 6, dpi = 300, compression = "lzw")

# 15. Print session info for reproducibility ####
sessionInfo()

################################################################################
# End of script
################################################################################
