################################################################################
# File: functions_transect.R
# Author: Monteiro-Leonel, Ana C.
# Description:
#   Custom functions and global settings for benthic time-series analysis
#   using transect-level aggregation (to avoid pseudoreplication, with
#   transects treated as the sampling unit). Includes LOESS smoothing,
#   Mann–Kendall trend testing (Yue & Wang, 2004), Sen’s slope estimation,
#   and visualization utilities.
################################################################################

# 1. Global color palette (standardized across all islands) ####
color_groups <- c(
  SCL = "#4E79A7",    # Blue - Corals/Scleractinians
  EAM = "#F28E2B",    # Orange - Epilithic Algal Matrix
  CCA = "#E15759",    # Red - Crustose Coralline Algae
  CYA = "#76B7B2",    # Teal - Cyanobacteria
  MAL = "#59A14F",    # Green - Macroalgae
  ACA = "#EDC948",    # Yellow - Articulated Coralline Algae
  ZOA = "#B07AA1",    # Purple - Zoanthids
  ABI = "#7F7F7F"     # Gray - Abiotic
)

# Retrieve color for a given benthic group
get_color <- function(group_name) {
  if(group_name %in% names(color_groups)) {
    return(color_groups[group_name])
  } else {
    warning(paste("No color defined for group:", group_name, "using default gray"))
    return("#888888")
  }
}

# 2. Basic graphical theme for all Mann-Kendall plots ####
theme_mk_plots <- function() {
  theme(axis.text.y = element_text(colour = "black", size = 10, angle = 0, 
                                   hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.text.x = element_text(colour = "black", size = 10, angle = 0, 
                                   hjust = 0.5, vjust = 0.0, face = "plain"),
        panel.grid.major = element_line(linewidth = 0.5, colour = "grey95", 
                                        lineend = "butt"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.x.bottom = element_line(linewidth = 0, colour = "black", 
                                          lineend = "butt"),
        axis.line.x.top = element_line(linewidth = 0, colour = "black", 
                                       lineend = "butt"),
        axis.line.y.left = element_line(linewidth = 0, colour = "black", 
                                        lineend = "butt"),
        axis.line.y.right = element_line(linewidth = 0, colour = "black", 
                                         lineend = "butt"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, vjust = -2, colour = "black", 
                                    size = 10, face = "bold"),
        axis.title.y = element_text(hjust = 0.5, vjust = 3, colour = "black", 
                                    size = 10, face = "bold"))
}

find_best_span <- function(data_year, span_values = seq(0.2, 1, by = 0.05)) {
  
  for(span in span_values) {
    
    warning_flag <- FALSE
    # Captura warnings
    tryCatch({
      withCallingHandlers({
        model <- loess(cover ~ year, data = data_year, span = span,
                       control = loess.control(surface = "direct"))
      }, warning = function(w) {
        warning_flag <<- TRUE
        invokeRestart("muffleWarning")
      })
      
    }, error = function(e) {
      warning_flag <<- TRUE
    })
    # Se não deu warning, retorna esse span
    if(!warning_flag) {
      return(span)
    }
  }
  
  # Se TODOS deram warning, retorna NA ou valor padrão
  return(NA)
}

# 3. Analyze Mann-Kendall trend for benthic groups ####
analyze_trend <- function(data, group_name) {
  
  df_raw <- data %>% 
    filter(group == group_name) %>%
    mutate(year = as.numeric(year))
  
  df_year <- df_raw %>%
    group_by(year) %>%
    summarise(cover = mean(cover, na.rm = TRUE), .groups = "drop")
  
  # Check minimum data requirement
  if(nrow(df_year) < 3) {
    warning(paste("Group", group_name, "has less than 3 observations. Analysis skipped."))
    return(NULL)
  }
  
  # Automatic span
  best_span <- find_best_span(df_year)
  
  # Fallback if all spans produce warnings
  if(is.na(best_span)) {
    best_span <- 0.75
  }
  
  # LOESS model
  loess_model <- loess(cover ~ year, data = df_year, span = best_span,
                       control = loess.control(surface = "direct"))
  
  df_year$loess_pred <- predict(loess_model)
  
  # MK on annual mean cover (main inference)
  # Yue-Pilon method accounts for autocorrelation in time series
  mk_result <- zyp.trend.vector(df_year$cover, df_year$year, method = "yuepilon")
  
  # Sen's slope with confidence interval
  trend_slope <- as.numeric(mk_result["trend"])   # Sen's slope (annual trend)
  ci_lower <- as.numeric(mk_result["lbound"])     # Lower CI bound 
  ci_upper <- as.numeric(mk_result["ubound"])     # Upper CI bound
  tau <- as.numeric(mk_result["tau"])             # Kendall's tau statistic
  p_value <- as.numeric(mk_result["sig"])         # p-value
  
  # Residual Standard Error from LOESS
  rse <- sqrt(mean((df_year$cover - df_year$loess_pred)^2, na.rm = TRUE))
  
  # Calculate data range
  data_range <- max(df_year$cover, na.rm = TRUE) - 
                min(df_year$cover, na.rm = TRUE)
  
  # Avoid division by zero
  if(is.na(data_range) || data_range == 0) {
    rse_percent <- 0
    trend_percent <- 0
  } else {
    rse_percent <- (rse / data_range) * 100
    trend_percent <- (trend_slope / data_range) * 100
  }
  # Return results
  return(list(
    data_transect = df_raw,             # transect-level data (for plotting)
    data_year = df_year,                # annual aggregated data (for trend analysis)
    tau = tau,                          # Kendall's tau 
    trend_slope = trend_slope,
    ci_lower = ci_lower,                # Sen's slope (lower)
    ci_upper = ci_upper,                # Sen's slope (upper)
    p_value = p_value,                  # p-value from Mann-Kendall test
    rse = rse,                          # Residual Standard Error
    range = data_range,
    rse_percent = rse_percent,
    trend_percent = trend_percent,
    loess_model = loess_model,
    span_used = best_span,
    mk_result = mk_result               # Keep full MK output for reference
  ))
}
 
# 4. Create trend plot ####
create_trend_plot <- function(result, group_name) {
  
  # Get color automatically from global palette
  group_color <- get_color(group_name)
  
  # Calculate significance label
  sig_label <- ifelse(result$p_value < 0.001, "***",
                      ifelse(result$p_value < 0.01, "**",
                             ifelse(result$p_value < 0.05, "*", "ns")))
  
  # Get axis ranges for annotation placement
  max_y <- max(result$data_transect$cover, na.rm = TRUE)
  min_y <- min(result$data_transect$cover, na.rm = TRUE)
  range_y <- max_y - min_y
  range_y_safe <- ifelse(is.na(range_y) || range_y == 0, 
                         max(max_y * 0.1, 1), 
                         range_y)
  
  min_year <- min(result$data_transect$year, na.rm = TRUE)
  max_year <- max(result$data_transect$year, na.rm = TRUE)
  
  ggplot(data = result$data_transect) +
    geom_boxplot(aes(x = year, y = cover, group = year), 
                 fill = scales::alpha(group_color, 0.5), notch = FALSE) +
    geom_jitter(aes(x = year, y = cover), 
                color = scales::alpha(group_color, 0.8), width = 0.1, height = 0) +
    geom_line(data = result$data_year, 
              aes(x = year, y = loess_pred),
              color = scales::alpha(group_color, 0.9), linewidth = 1) +
    scale_x_continuous(breaks = seq(min(result$data_transect$year, na.rm = TRUE), 
                       max(result$data_transect$year, na.rm = TRUE), by = 1)) +
    labs(x = NULL, y = "Relative benthic cover (%)") +
    theme_mk_plots() +
    annotate("text", x = min_year, y = max_y - 0.02 * range_y_safe, 
             label = paste0("Tau: ", round(result$tau, 2), " ", sig_label, 
                            " | Span: ", round(result$span_used, 2)), 
             color = group_color, fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = min_year, y = max_y - 0.08 * range_y_safe, 
             label = paste0("Slope: ", round(result$trend_slope, 2)),
             color = group_color, fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = max_year, y = max_y - 0.02 * range_y_safe, 
             label = group_name, color = group_color, 
             fontface = "bold", size = 6, hjust = 1)
}

# 5. Print detailed trend statistics ####
print_trend_stats <- function(result, group_name) {
  cat("\n========================================\n")
  cat("Group:", group_name, "\n")
  cat("========================================\n")
  
  cat("LOESS Model:\n")
  cat("  Span:", result$loess_model$pars$span, "\n")
  cat("  Residual Standard Error:", round(result$rse, 4), "\n")
  cat("  Data range:", round(result$range, 4), "\n\n")
  
  cat("Mann-Kendall Test (applied to observed annual mean cover):\n")
  cat("  Tau:", round(result$tau, 6), "\n")
  
  if(result$p_value < 0.001){
    cat("  p-value: < 0.001\n")
  } else {
    cat("  p-value:", round(result$p_value, 4), "\n")
  }
  
  cat("  Sen's slope:", round(result$trend_slope, 6), "cover units per year\n")
  cat("  95% CI:", round(result$ci_lower, 6), "-", round(result$ci_upper, 6), "\n")
  cat("  Annual trend:", round(result$trend_percent, 2), "% of observed range\n")
  cat("  RSE as % of range:", round(result$rse_percent, 1), "%\n")
  trend_direction <- ifelse(result$tau > 0, "increasing", 
                            ifelse(result$tau < 0, "decreasing", "no trend"))
  
  cat("  Trend direction:", trend_direction, "\n")
  # Significance stars
  sig <- ifelse(result$p_value < 0.001, "***",
                ifelse(result$p_value < 0.01, "**",
                       ifelse(result$p_value < 0.05, "*", "ns")))
  cat("\nSignificance:", sig, "\n")
}
