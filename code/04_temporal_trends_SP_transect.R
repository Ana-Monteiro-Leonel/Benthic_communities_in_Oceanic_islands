################################################################################
# Script: 04_temporal_trends_SP_transect.R
# Author: Monteiro-Leonel, Ana C.
# Date: 2026-03-28
# Description: 
#   Temporal trends analysis for São Pedro and São Paulo Archipelago (SP)
#   Identified as "stpauls_rocks" in the dataset.
#   Automatically selects the 4 groups with highest median cover for trend analysis.
# Outputs:
#   - results/figures/Figure_4_temporal_trends_SP.png
#   - results/figures/Figure_4_temporal_trends_SP.tiff
#   - results/tables/Table_SP_summary_statistics.csv
#   - results/tables/Table_SP_trend_summary.csv
################################################################################

# Set working directory to project root ####
# This script tries to find the project root automatically.
# If it fails, adjust the path below to your local setup.
project_root <- "C:/Users/Ana Leonel/OneDrive/Documentos/GitHub/Benthic_communities_in_Oceanic_islands"

if (dir.exists(project_root)) {
  setwd(project_root)
  cat("Working directory set to:", getwd(), "\n")
} else {
  # Try to find project root by looking for code/functions_transect.R
  test_dir <- getwd()
  found <- FALSE
  for (i in 1:5) {
    if (file.exists(file.path(test_dir, "code/functions_transect.R"))) {
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
         "Expected project path: C:/Users/Ana Leonel/OneDrive/Documentos/GitHub/Benthic_communities_in_Oceanic_islands")
  }
}

# Verify functions file exists
if (!file.exists("code/functions_transect.R")) {
  stop("functions_transect.R not found in code/ directory. 
       Please check your working directory.")
}

# 1. Load packages ####
library(ggplot2)
library(dplyr)
library(zyp)
library(tidyr)
library(patchwork)

# 2. Source functions and global settings ####
source("code/functions_transect.R")

# 3. Load and prepare data with functional grouping ####
df_bio_sp <- read.csv("data/raw/benthic_complete_data.csv") %>%
  filter(island == "stpauls_rocks") %>%
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
    sd_cover = sd(cover, na.rm = TRUE),
    n_images = n(),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    year = as.numeric(year),
    cv = ifelse(cover > 0, sd_cover / cover, NA)
  )

# Check sampling effort
sampling_effort <- df_bio_sp %>%
  group_by(group, year) %>%
  summarise(n_transects = n(), .groups = "drop")

write.csv(sampling_effort, "results/tables/Table_SP_sampling_effort.csv", 
          row.names = FALSE)

# 4. Select top 4 groups with highest median cover ####
median_cover <- df_bio_sp %>%
  group_by(group) %>%
  summarise(
    median_cover = median(cover, na.rm = TRUE),
    mean_cover = mean(cover, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_cover))
write.csv(median_cover, "results/tables/Table_SP_median_cover.csv")

cat("\n=== Median cover per group ===\n")
print(median_cover)

# Select top 4 groups
top_groups <- median_cover %>%
  slice_head(n = 4) %>%
  pull(group)

cat("\n=== Top 4 groups selected for trend analysis ===\n")
print(top_groups)

# 5. Run trend analysis ####
results <- list()
for(g in top_groups) {
  cat("\nAnalyzing group:", g, "\n")
  result <- analyze_trend(df_bio_sp, g)
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
  
  combined_plot <- combined_plot +
    theme(
      plot.margin = margin(5, 5, 5, 5),
      panel.spacing = unit(0.3, "lines")
    )
  
  # Save outputs
  ggsave("results/figures/Figure_4_temporal_trends_SP.png", combined_plot, 
         width = 9, height = 7, dpi = 600)
  ggsave("results/figures/Figure_4_temporal_trends_SP.tiff", combined_plot, 
         width = 9, height = 7, dpi = 600, compression = "lzw")
  
  print(combined_plot)
}

# 7. Save statistics ####
stats_table <- df_bio_sp %>%
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

write.csv(stats_table, "results/tables/Table_SP_summary_statistics.csv", row.names = FALSE)

# 8. Save trend summary ####
trend_summary <- data.frame(
  Island = "stpauls_rocks",
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

write.csv(trend_summary, "results/tables/Table_SP_trend_summary.csv", row.names = FALSE)

# 9. Print summary ####
cat("\n=== TREND ANALYSIS SUMMARY ===\n")
cat("Location: São Pedro and São Paulo Archipelago (stpauls_rocks)\n")
cat("Period: 2013-2019\n\n")

for(g in names(results)) {
  res <- results[[g]]
  sig <- ifelse(res$p_value < 0.001, "***",
                ifelse(res$p_value < 0.01, "**",
                       ifelse(res$p_value < 0.05, "*", "ns")))
  direction <- ifelse(res$tau > 0, "Increasing", 
                      ifelse(res$tau < 0, "Decreasing", "Stable"))
  
  cat(sprintf("%s: %s (τ = %.3f %s, p = %.3f)\n", 
              g, direction, res$tau, sig, res$p_value))
  cat(sprintf("  Annual trend: %.1f%% of range\n", res$trend_percent))
}

cat("\nAnalysis complete!\n")

# 10. Print session info for reproducibility ####
sessionInfo()

################################################################################
# End of script
################################################################################                         
