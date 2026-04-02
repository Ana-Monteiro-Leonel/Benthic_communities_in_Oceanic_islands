################################################################################
# Script: 07_temporal_trends_TR_transect.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-04-01
# Description: 
#   Temporal trends analysis for Trindade Island (TR)
#   Identified as "trindade" in the dataset.
#   Groups raw categories into functional groups before analysis.
#   Manual selection of groups for trend analysis (ecological prioritization).
# Outputs:
#   - results/figures/Figure_5_temporal_trends_TR.png
#   - results/figures/Figure_5_temporal_trends_TR.tiff
#   - results/tables/Table_TR_summary_statistics.csv
#   - results/tables/Table_TR_trend_summary.csv
################################################################################

# 1. Load packages ####
library(ggplot2)
library(dplyr)
library(zyp)
library(tidyr)
library(patchwork)

# NOTE:
# Set the working directory to the project root before running this script.
# Otherwise, relative paths (e.g., "code/functions_transect.R") will not work.
# Example:
# setwd("path/to/Benthic_communities_in_Oceanic_islands")

# 2. Source functions and global settings ####
source("code/functions_transect.R")

# 3. Load and prepare data with functional grouping ####
df_bio_tr <- read.csv("data/raw/benthic_complete_data.csv") %>%
  filter(island == "trindade") %>%
  mutate(
    categoryid = as.character(categoryid),
    group = case_when(
      categoryid %in% c("CEN", "COR", "COT", "FIL", "FOL", "GLC", "SAR", "STO") ~ "MAL",
      categoryid %in% c("BRY", "ECH") ~ "INV",
      categoryid %in% c("BSC", "ODI") ~ "SCL",
      TRUE ~ categoryid
    )
  ) %>%
  # Sum cover by image
  group_by(island, sites, year, transect, image, group) %>%
  summarise(cover = sum(coverpercategory), .groups = "drop") %>%
  
  # Aggregate by transect
  group_by(island, sites, year, transect, group) %>%
  summarise(
    cover = mean(cover, na.rm = TRUE),
    sd_cover = sd(cover, na.rm = TRUE),   # optional (useful for variability metrics)
    n_images = n(),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    year = as.numeric(year),
    # Coefficient of variation (relative variability)
    cv = ifelse(cover > 0, sd_cover / cover, NA) 
  )

# Check
sampling_effort <- df_bio_tr %>%
  group_by(group, year) %>%
  summarise(n_transects = n(), .groups = "drop")

write.csv(sampling_effort, "results/tables/Table_TR_sampling_effort.csv", 
          row.names = FALSE)

# 4. Summarize median cover by group ####
median_cover <- df_bio_tr %>%
  group_by(group) %>%
  summarise(
    median_cover = median(cover, na.rm = TRUE),
    mean_cover = mean(cover, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_cover))
write.csv(median_cover, "results/tables/Table_TR_median_cover.csv")

cat("\n=== Median cover per group ===\n")
print(median_cover)

# MANUAL SELECTION OF GROUPS FOR TREND ANALYSIS ####
# For Trindade, we prioritize dominant benthic biological groups,
# excluding abiotic cover (ABI) despite high values,
# focusing on ecologically meaningful biotic components.
#
# Based on median cover, top groups are:
#   1. MAL - include (Macroalgae)
#   2. ABI - exclude (Abiotic sand, bare rock, etc)
#   3. EAM - include (Epilithic Algal Matrix) 
#   4. CCA - include (Crustose Coralline Algae) 
#   5. CYA - include (cyanobacteria)
#
# Selected groups: "MAL", "EAM", "CCA", "CYA"

# Manually select groups for trend analysis (Trindade - algal dominated system)
selected_groups <- c("MAL", "EAM", "CCA", "CYA")

# Verify selected groups exist in data
available_groups <- intersect(selected_groups, unique(df_bio_tr$group))
if(length(available_groups) < length(selected_groups)) {
  stop(paste("Missing groups in dataset:", 
             paste(setdiff(selected_groups, available_groups), collapse = ", ")))
}

cat("\n=== Selected groups for trend analysis ===\n")
cat("(Including CYA due to its ecological relevance and observed recent increase)\n")
cat("Selected:", paste(available_groups, collapse = ", "), "\n")

# Show why ABI was excluded
abi_median <- median_cover %>% filter(group == "ABI") %>% pull(median_cover)
cya_median <- median_cover %>% filter(group == "CYA") %>% pull(median_cover)
if(length(abi_median) > 0 && length(cya_median) > 0) {
  cat("\nNOTE: ABI (abiotic substrate) had median cover =", round(abi_median, 1), "%.\n")
  cat("      However, it was excluded as it represents non-biological space (e.g., sand, bare rock),\n")
  cat("      and the analysis focuses on ecologically relevant biotic functional groups.\n")
  cat("      CYA median =", round(cya_median, 1), 
      "% - CYA included due to its ecological relevance in Trindade.\n")
}

# 5. Run trend analysis ####
results <- list()
for(g in available_groups) {
  cat("\nAnalyzing group:", g, "\n")
  result <- analyze_trend(df_bio_tr, g)
  if(!is.null(result)) {
    results[[g]] <- result
    cat("  Tau =", round(result$tau, 3), 
        "p =", format(result$p_value, digits = 4), "\n")
    cat("  Span =", round(result$span_used, 2), "\n")
  } else {
    cat("  Analysis failed\n")
  }
}

# After running trend analysis
for(g in names(results)) {
  print_trend_stats(results[[g]], g)
}

# 6. Create plots (colors are handled automatically by get_color() in functions_transect.R) ####
if(length(results) > 0) {
  plots <- list()
  for(g in names(results)) {
    plots[[g]] <- create_trend_plot(results[[g]], g)
  }
  
  # Arrange layout
  n_plots <- length(plots)
  if(n_plots == 4) {
    combined_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  } else if(n_plots == 3) {
    combined_plot <- (plots[[1]] | plots[[2]]) / plots[[3]]
  } else if(n_plots == 2) {
    combined_plot <- plots[[1]] | plots[[2]]
  } else {
    combined_plot <- plots[[1]]
  }
  
  combined_plot <- combined_plot + 
    plot_annotation(title = NULL,
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  
  # Save outputs
  ggsave("results/figures/Figure_5_temporal_trends_TR.png", combined_plot, 
         width = 10, height = 8, dpi = 300)
  ggsave("results/figures/Figure_5_temporal_trends_TR.tiff", combined_plot, 
         width = 10, height = 8, dpi = 300, compression = "lzw")
  
  print(combined_plot)
}

# 7. Save statistics ####
stats_table <- df_bio_tr %>%
  group_by(group) %>%
  summarise(
    mean_cover = mean(cover, na.rm = TRUE),
    sd_cover = sd(cover, na.rm = TRUE),
    median_cover = median(cover, na.rm = TRUE),
    min_cover = min(cover, na.rm = TRUE),
    max_cover = max(cover, na.rm = TRUE),
    n_obs = n(),
    n_years = n_distinct(year),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_cover))

write.csv(stats_table, "results/tables/Table_TR_summary_statistics.csv", row.names = FALSE)

# 8. Save trend summary ####
trend_summary <- data.frame(
  Island = "trindade",
  Group = names(results),
  Span = sapply(results, function(x) round(x$span_used, 2)),
  Tau = sapply(results, function(x) round(x$tau, 3)),
  P_value = sapply(results, function(x) {
    if(x$p_value < 0.001) "<0.001" else round(x$p_value, 4)
  }),
  Significance = sapply(results, function(x) {
    ifelse(x$p_value < 0.001, "***",
           ifelse(x$p_value < 0.01, "**",
                  ifelse(x$p_value < 0.05, "*", "ns")))
  }),
  Slope = sapply(results, function(x) round(x$trend_slope, 3)),
  CI_lower = sapply(results, function(x) round(x$ci_lower, 3)),
  CI_upper = sapply(results, function(x) round(x$ci_upper, 3)),
  CI = paste0("(", 
              sprintf("%.2f", sapply(results, function(x) x$ci_lower)), ", ",
              sprintf("%.2f", sapply(results, function(x) x$ci_upper)), ")"),
  Trend_percent = sapply(results, function(x) round(x$trend_percent, 1)),
  Range = sapply(results, function(x) round(x$range, 2)),
  
  stringsAsFactors = FALSE
)

write.csv(trend_summary, "results/tables/Table_TR_trend_summary.csv", row.names = FALSE)

# 9. Print summary ####
cat("\n=== TREND ANALYSIS SUMMARY ===\n")
cat("Location: Trindade Island (TR)\n")
cat("Period: 2013-2019\n\n")

for(g in names(results)) {
  res <- results[[g]]
  sig <- ifelse(res$p_value < 0.001, "***",
                ifelse(res$p_value < 0.01, "**",
                       ifelse(res$p_value < 0.05, "*", "ns")))
  direction <- ifelse(res$tau > 0, "Increasing", 
                      ifelse(res$tau < 0, "Decreasing", "Stable"))
  p_val <- ifelse(res$p_value < 0.001, "<0.001", sprintf("%.3f", res$p_value))
  cat(sprintf("%s: %s (τ = %.3f %s, p = %s)\n", 
              g, direction, res$tau, sig, p_val))
  cat(sprintf("  Annual trend: %.1f%% of range\n", res$trend_percent))
}

cat("\nAnalysis complete!\n")

