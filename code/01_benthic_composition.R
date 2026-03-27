################################################################################
# Script: 01_benthic_composition.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-03-20
# Description: 
#   This script processes benthic cover data by:
#   1. Recoding fine-scale categories into functional groups
#   2. Aggregating cover per group at image level
#   3. Using transect as the sampling unit to avoid pseudoreplication
#   4. Creating a stacked bar plot of relative abundance of benthic groups
#      with percentage labels on bars (for groups > 5%)
# Outputs:
#   - data/processed/benthic_cover_photos.csv (image-level data)
#   - data/processed/benthic_cover_transects.csv (transect-level data)
#   - data/processed/benthic_cover_summary.csv (island-level summary)
#   - results/figures/Figure_2_benthic_composition.png
#   - results/tables/Table_S1_relative_abundance.csv
################################################################################

# Set working directory (adjust as needed) ####
# Set this to the root of your repository
setwd("C:/Users/Ana Monteiro/OneDrive/Documentos/GitHub/Benthic_communities_in_Oceanic_islands")

# Check current directory
getwd()

# Check and install required packages ####
required_packages <- c("dplyr", "ggplot2", "tidyr", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 1. Load data ####
complete_data <- read_csv("data/raw/benthic_complete_data.csv")

# 2. Inspect data structure ####
glimpse(complete_data)

# 3. Recode fine-scale categories into benthic groups ####
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
  ) %>%
  select(island, year, sites, transect, image, group, cover_per_group)

# 4. Save image-level data (for temporal trends script) ####
write_csv(grouped_data, "data/processed/benthic_cover_photos.csv")


# 5. Calculate mean cover per ISLAND using images as replicates ####
summary_data <- grouped_data %>%
  group_by(island, group) %>%
  summarise(
    mean_cover = mean(cover_per_group, na.rm = TRUE),
    sd_cover = sd(cover_per_group, na.rm = TRUE),
    n_image = n(),
    se_cover = sd_cover / sqrt(n_image),
    .groups = "drop"
  ) %>%
  group_by(island) %>%
  mutate(
    rel_abund = mean_cover / sum(mean_cover) * 100
  ) %>%
  ungroup()

# 6. Standardize island names ####
summary_data <- summary_data %>%
  mutate(
    island = recode(island,
                    "trindade" = "TR",
                    "noronha" = "FN",
                    "rocas" = "RA",
                    "stpauls_rocks" = "SP"
    )
  )

# 7. Define order of groups and convert to factor ####
group_order <- c("EAM", "MAL", "CCA", "ACA", "SCL", "ABI", "INV", "ZOA", "CYA", "SUS")

summary_data <- summary_data %>%
  mutate(
    group = factor(group, levels = group_order)
  ) %>%
  arrange(island, group)

# 8. Save summary data ####
write_csv(summary_data, "data/processed/benthic_cover_summary.csv")

# 9. Define colors for benthic groups ####
benthic_colors <- c(
  EAM = "#F28E2B",   # Orange
  MAL = "#59A14F",   # Green
  CCA = "#E15759",   # Cherry red
  ACA = "#EDC948",   # Yellow
  SCL = "#4E79A7",   # Blue
  ABI = "#595959",   # Gray
  INV = "#9E0142",   # Pink/rose
  ZOA = "#B07AA1",   # Lilac
  CYA = "#17BECF",   # Cyan
  SUS = "#7F9F7F"    # Olive green
)

# 10. Create stacked bar plot ####
stacked_bar <- ggplot(summary_data, aes(x = island, y = rel_abund, fill = group)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  scale_fill_manual(values = benthic_colors, name = "Benthic group", drop = FALSE) +
  scale_x_discrete(limits = c("TR", "FN", "RA", "SP")) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm")
  ) +
  coord_flip() +
  guides(fill = guide_legend(ncol = 1, keyheight = 0.8, keywidth = 0.8))

# 11. Add labels for groups > 5% ####
stacked_bar_with_labels <- stacked_bar +
  geom_text(
    aes(label = ifelse(rel_abund > 5, sprintf("%.1f%%", rel_abund), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white",
    size = 3.5,
    fontface = "bold"
  )

# 12. Display and save plots ####
print(stacked_bar_with_labels)

ggsave("results/figures/Figure_2_benthic_composition.png",
       plot = stacked_bar,
       width = 8,
       height = 6,
       dpi = 300)

ggsave("results/figures/Figure_2_benthic_composition_with_labels.png",
       plot = stacked_bar_with_labels,
       width = 8,
       height = 6,
       dpi = 300)

# 13. Create and save abundance table ####
abundance_table <- summary_data %>%
  select(island, group, rel_abund) %>%
  mutate(rel_abund = round(rel_abund, 1)) %>%
  pivot_wider(names_from = island, values_from = rel_abund, values_fill = 0) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group) %>%
  rename(`Benthic group` = group,
         `Trindade (TR)` = TR,
         `Fernando de Noronha (FN)` = FN,
         `Atol das Rocas (RA)` = RA,
         `São Pedro e São Paulo (SP)` = SP)

write_csv(abundance_table, "results/tables/Table_S1_relative_abundance.csv")

################################################################################
# End of script
################################################################################