#=============================
# FIRE INSECT CO-OCCURRENCE 3-9 YEARS ANALYSIS
#=============================


# Load required libraries
library(MatchIt)
library(dplyr)
library(marginaleffects)
library(MuMIn)

# Source the file with function definitions
source("/home/goldma34/fire_insect_co-occurence/src/best_model_functions.R")

# Source the file with covariate balance plots
source("/home/goldma34/fire_insect_co-occurence/src/covariate_balance_plots.R")

# Load your data
source("/home/goldma34/fire_insect_co-occurence/src/load_data.R")

# Check if data is loaded
if (!exists("hist_gt90_2_2")) {
  stop("Data not loaded. Make sure hist_gt90_2_2 exists.")
}

# Define parameters
methods <- c("nearest", "optimal")
distances <- c("glm", "mahalanobis")
links <- c("logit", "probit")
m_orders <- c("random")
calipers <- c(0.1, 0.2, 0.3)
replace_list <- c(FALSE)
use_mahvars_list <- c(FALSE, TRUE)

#rename df
intermediate <- hist_gt90_2_2

#============================
# SEVERITY ANALYSIS FOR HIST_GT90_2_2
#============================
cat("\n\n======= SEVERITY ANALYSIS FOR intermediate =======\n")
cat("Running model selection...\n")
best_model_info_sev <- find_best_model(intermediate, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev <- best_model_info_sev$best_model
all_results_sev <- best_model_info_sev$results

# Check if we found a best model
if (!is.null(best_model_sev)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev <- match_data(best_model_sev)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev))

  # Calculate R-squared
  cat("\nR-squared:\n")
  fit.sev_r2 <- r.squaredGLMM(fit.sev)
  print(fit.sev_r2)
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_intermediate_severity.RDS")
  saveRDS(fit.sev, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_intermediate_severity.RDS")
  saveRDS(fit.sev_r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_intermediate_severity_r2.RDS")
  
  # Calculate treatment effects
  cat("\nCalculating treatment effects...\n")
  avg_effects <- avg_comparisons(
    fit.sev, 
    variables = "history",
    newdata = subset(fit.sev$model, history == 1)
  )
  print(avg_effects)
  
  # Save treatment effects
  saveRDS(avg_effects, "/home/goldma34/fire_insect_co-occurence/data/results/treatment_effects_intermediate_severity.RDS")
  
  # Generate balance plots
  cat("\nGenerating balance plots...\n")
  plot_dir <- "/home/goldma34/fire_insect_co-occurence/plots/balance/intermediate/"
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_list <- list("3-9 years - Severity" = best_model_sev)
  balance_plots <- generate_all_balance_plots(model_list, plot_dir)
  
  # Power analysis
  cat("\nPerforming power analysis...\n")
  library(pwr)
  
  # Extract observed effect
  model_summary <- summary(fit.sev)
  treatment_coef <- coef(model_summary)["history", ]
  observed_effect <- treatment_coef["Estimate"]
  observed_se <- treatment_coef["Std. Error"]
  
  # Calculate Cohen's d effect size
  d <- observed_effect / sd(m.data.sev$rbr_w_offset)
  cat("\nObserved effect size (Cohen's d):", round(d, 3), "\n")
  
  # Calculate power for different sample sizes
  sample_sizes <- seq(20, 300, by = 10)
  powers <- sapply(sample_sizes, function(n) {
    pwr.t.test(d = d, n = n/2, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  })
  
  # Create power curve plot
  library(ggplot2)
  power_data <- data.frame(sample_size = sample_sizes, power = powers)
  power_plot <- ggplot(power_data, aes(x = sample_size, y = power)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Total Sample Size", 
      y = "Statistical Power",
      title = paste("Power Analysis for intermediate (d =", round(d, 2), ")")
    ) +
    theme_minimal()
  print(power_plot)
  
  # Save power plot
  power_dir <- "/home/goldma34/fire_insect_co-occurence/plots/power_analysis/intermediate/"
  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(power_dir, "power_curve_severity.png"), power_plot, width = 8, height = 6, dpi = 300)
  
  # Find required sample size for power = 0.8
  min_sample_size <- sample_sizes[min(which(powers >= 0.8))]
  cat("Minimum sample size required for 80% power:", min_sample_size, "\n")
  
  # Calculate minimum detectable effect size with current sample
  current_sample_size <- nrow(m.data.sev)
  min_d <- pwr.t2n.test(
    n1 = sum(m.data.sev$history == 1), 
    n2 = sum(m.data.sev$history == 0), 
    power = 0.8,
    sig.level = 0.05,
    alternative = "two.sided"
  )$d
  
  min_raw_effect <- min_d * sd(m.data.sev$rbr_w_offset)
  cat("\nMinimum detectable effect size (Cohen's d) with current sample size:", round(min_d, 3), "\n")
  cat("Minimum detectable effect size (raw units) with current sample size:", round(min_raw_effect, 3), "\n")
  
} else {
  cat("No suitable severity model found. Try different parameters.\n")
}

#============================
# RECOVERY ANALYSIS FOR intermediate
#============================
cat("\n\n======= RECOVERY ANALYSIS FOR intermediate =======\n")
cat("Running model selection...\n")
best_model_info_rec <- find_best_model(intermediate, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec <- best_model_info_rec$best_model
all_results_rec <- best_model_info_rec$results

# Check if we found a best model
if (!is.null(best_model_rec)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec <- match_data(best_model_rec)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec <- lm(recovery ~ history + host_pct + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec))
  
  # Calculate R-squared
  cat("\nR-squared:\n")
  fit.rec_r2 <- r.squaredGLMM(fit.rec)
  print(fit.rec_r2)
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_intermediate_recovery.RDS")
  saveRDS(fit.rec, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_intermediate_recovery.RDS")
  saveRDS(fit.rec_r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_intermediate_recovery_r2.RDS")
  
  # Calculate treatment effects
  cat("\nCalculating treatment effects...\n")
  avg_effects_rec <- avg_comparisons(
    fit.rec, 
    variables = "history",
    newdata = subset(fit.rec$model, history == 1)
  )
  print(avg_effects_rec)
  
  # Save treatment effects
  saveRDS(avg_effects_rec, "/home/goldma34/fire_insect_co-occurence/data/results/treatment_effects_intermediate_recovery.RDS")
  
  # Generate balance plots
  cat("\nGenerating balance plots...\n")
  plot_dir <- "/home/goldma34/fire_insect_co-occurence/plots/balance/intermediate/"
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_list <- list("3-9 years - Recovery" = best_model_rec)
  balance_plots <- generate_all_balance_plots(model_list, plot_dir)
  
  # Power analysis
  cat("\nPerforming power analysis...\n")
  
  # Extract observed effect
  model_summary <- summary(fit.rec)
  treatment_coef <- coef(model_summary)["history", ]
  observed_effect <- treatment_coef["Estimate"]
  observed_se <- treatment_coef["Std. Error"]
  
  # Calculate Cohen's d effect size
  d <- observed_effect / sd(m.data.rec$recovery)
  cat("\nObserved effect size (Cohen's d):", round(d, 3), "\n")
  
  # Calculate power for different sample sizes
  sample_sizes <- seq(20, 300, by = 10)
  powers <- sapply(sample_sizes, function(n) {
    pwr.t.test(d = d, n = n/2, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  })
  
  # Create power curve plot
  power_data <- data.frame(sample_size = sample_sizes, power = powers)
  power_plot_rec <- ggplot(power_data, aes(x = sample_size, y = power)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Total Sample Size", 
      y = "Statistical Power",
      title = paste("Power Analysis for intermediate Recovery (d =", round(d, 2), ")")
    ) +
    theme_minimal()
  print(power_plot_rec)
  
  # Save power plot
  power_dir <- "/home/goldma34/fire_insect_co-occurence/plots/power_analysis/intermediate/"
  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(power_dir, "power_curve_recovery.png"), power_plot_rec, width = 8, height = 6, dpi = 300)
  
  # Find required sample size for power = 0.8
  min_sample_size <- sample_sizes[min(which(powers >= 0.8))]
  cat("Minimum sample size required for 80% power:", min_sample_size, "\n")
  
  # Calculate minimum detectable effect size with current sample
  current_sample_size <- nrow(m.data.rec)
  min_d <- pwr.t2n.test(
    n1 = sum(m.data.rec$history == 1), 
    n2 = sum(m.data.rec$history == 0), 
    power = 0.8,
    sig.level = 0.05,
    alternative = "two.sided"
  )$d
  
  min_raw_effect <- min_d * sd(m.data.rec$recovery)
  cat("\nMinimum detectable effect size (Cohen's d) with current sample size:", round(min_d, 3), "\n")
  cat("Minimum detectable effect size (raw units) with current sample size:", round(min_raw_effect, 3), "\n")
  
} else {
  cat("No suitable recovery model found. Try different parameters.\n")
}

#============================
# SUMMARY
#============================
cat("\n\n======= SUMMARY OF ANALYSIS FOR intermediate =======\n")

# Create a summary dataframe
summary_df <- data.frame(
  Analysis = c("Severity", "Recovery"),
  Model_Found = c(!is.null(best_model_sev), !is.null(best_model_rec)),
  stringsAsFactors = FALSE
)

# Add sample sizes if models were found
if (!is.null(best_model_sev)) {
  summary_df$Total_Sample_Size[summary_df$Analysis == "Severity"] <- nrow(m.data.sev)
  summary_df$Treated_Sample[summary_df$Analysis == "Severity"] <- sum(m.data.sev$history == 1)
  summary_df$Control_Sample[summary_df$Analysis == "Severity"] <- sum(m.data.sev$history == 0)
}

if (!is.null(best_model_rec)) {
  summary_df$Total_Sample_Size[summary_df$Analysis == "Recovery"] <- nrow(m.data.rec)
  summary_df$Treated_Sample[summary_df$Analysis == "Recovery"] <- sum(m.data.rec$history == 1)
  summary_df$Control_Sample[summary_df$Analysis == "Recovery"] <- sum(m.data.rec$history == 0)
}

# Print summary
print(summary_df)

# Save summary
write.csv(summary_df, "/home/goldma34/fire_insect_co-occurence/data/results/intermediate_analysis_summary.csv", row.names = FALSE)

cat("\nAnalysis complete. Results saved to data/results directory.\n")

#============================
# GENERATE TREATMENT EFFECT PLOTS
#============================
cat("\n\n======= GENERATING TREATMENT EFFECT PLOTS =======\n")

# Source dedicated intermediate treatment effects plotting script
source("/home/goldma34/fire_insect_co-occurence/src/treatment_effect_plots_intermediate.R")
generate_all_intermediate_treatment_effect_plots()

#============================
# GENERATE TREATMENT EFFECT TABLES
#============================
cat("\n\n======= GENERATING TREATMENT EFFECT TABLES =======\n")

# Source dedicated intermediate treatment effects export script
source("/home/goldma34/fire_insect_co-occurence/src/export_treatment_effects_intermediate.R")
export_intermediate_treatment_effects()

cat("\nAnalysis complete. Results saved to data/results directory.\n")
