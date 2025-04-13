# Load required libraries
library(sensemakr)
library(dplyr)

# Define output directory
result_dir <- "/home/goldma34/fire_insect_co-occurence/data/results/"

# Function to safely read RDS files
safe_read_rds <- function(file_path) {
  if (file.exists(file_path)) {
    tryCatch({
      readRDS(file_path)
    }, error = function(e) {
      cat("Error reading", file_path, ":", e$message, "\n")
      NULL
    })
  } else {
    cat("File not found:", file_path, "\n")
    NULL
  }
}

# Load all saved model fits
cat("Loading models...\n")
fit_sev_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_severity.RDS"))
fit_sev_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_severity.RDS"))
fit_sev_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_severity.RDS"))
fit_sev_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_severity.RDS"))

fit_rec_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_recovery.RDS"))
fit_rec_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_recovery.RDS"))
fit_rec_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_recovery.RDS"))
fit_rec_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_recovery.RDS"))

# Define covariates for each model type
covariates_sev <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")
covariates_rec <- c("host_pct", "mean_temperature", "sum_precipitation_mm", "mean_tri")

# Function to analyze one model and return results
analyze_one_model <- function(model, model_name, treatment_var, covariates) {
  if (is.null(model)) {
    cat("Model is NULL for", model_name, "\n")
    return(NULL)
  }
  
  cat("Analyzing", model_name, "...\n")
  
  # Extract model summary
  model_summary <- summary(model)
  
  # Try to extract coefficient info
  tryCatch({
    # Extract coefficient data
    coef_data <- model_summary$coefficients
    
    # Check if treatment variable exists in the model
    if (!(treatment_var %in% rownames(coef_data))) {
      cat("Treatment variable", treatment_var, "not found in model coefficients for", model_name, "\n")
      return(NULL)
    }
    
    # Extract treatment effect stats
    treatment_effect <- coef_data[treatment_var, "Estimate"]
    std_error <- coef_data[treatment_var, "Std. Error"]
    t_value <- coef_data[treatment_var, "t value"]
    p_value <- coef_data[treatment_var, "Pr(>|t|)"]
    
    # Calculate partial R^2 for the treatment variable
    # Formula: t^2 / (t^2 + df)
    df <- model_summary$df[2]  # Residual degrees of freedom
    partial_r2 <- t_value^2 / (t_value^2 + df)
    
    # Calculate robustness values based on formulas from sensemakr paper
    # RV_q1 is the partial R^2 needed with both outcome and treatment to reduce effect to 0
    # A simple approximation is the square root of the partial R^2
    rv_q1 <- sqrt(partial_r2)
    
    # RV_q1_alpha is the partial R^2 needed to reduce to statistical insignificance
    # An approximation is reducing the absolute t-value to the critical value (1.96 for α=0.05)
    critical_t <- qt(0.975, df)  # Two-tailed 5% critical value
    if (abs(t_value) > critical_t) {
      # Calculate how much of a reduction is needed to reach critical value
      reduction_factor <- 1 - (critical_t / abs(t_value))^2
      rv_q1_alpha <- partial_r2 * reduction_factor
    } else {
      # Already insignificant
      rv_q1_alpha <- 0
    }
    
    # Create result row
    result <- data.frame(
      Model = model_name,
      Treatment_Effect = treatment_effect,
      Standard_Error = std_error,
      T_Value = t_value,
      P_Value = p_value,
      Partial_R2 = partial_r2,
      RV_q1 = rv_q1,
      RV_q1_alpha0.05 = rv_q1_alpha,
      Effect_Direction = ifelse(treatment_effect > 0, "Positive", "Negative"),
      Significance = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      stringsAsFactors = FALSE
    )
    
    # Add robustness ratio and assessment
    typical_r2 <- 0.15  # A moderate benchmark partial R2
    result$Robustness_Ratio <- result$RV_q1 / typical_r2
    
    result$Robustness_Assessment <- case_when(
      result$Robustness_Ratio > 2 ~ "Very Strong",
      result$Robustness_Ratio > 1 ~ "Strong",
      result$Robustness_Ratio > 0.5 ~ "Moderate",
      result$Robustness_Ratio > 0.25 ~ "Weak",
      TRUE ~ "Very Weak"
    )
    
    return(result)
  }, error = function(e) {
    cat("Error analyzing", model_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Analyze all models
cat("Starting direct sensitivity analysis...\n")

# Create list of all models, names, etc.
models <- list(fit_sev_1, fit_sev_2, fit_sev_3, fit_sev_4,
              fit_rec_1, fit_rec_2, fit_rec_3, fit_rec_4)

model_names <- c("Subgroup 1 (0-2 years) - Severity",
                "Subgroup 2 (3-5 years) - Severity",
                "Subgroup 3 (6-9 years) - Severity",
                "Subgroup 4 (10+ years) - Severity",
                "Subgroup 1 (0-2 years) - Recovery",
                "Subgroup 2 (3-5 years) - Recovery",
                "Subgroup 3 (6-9 years) - Recovery",
                "Subgroup 4 (10+ years) - Recovery")

treatment_var <- "history"

cov_list <- list(
  covariates_sev, covariates_sev, covariates_sev, covariates_sev,
  covariates_rec, covariates_rec, covariates_rec, covariates_rec
)

# Run analysis for each model and collect results
all_results <- list()

for (i in 1:length(models)) {
  result <- analyze_one_model(models[[i]], model_names[i], treatment_var, cov_list[[i]])
  all_results[[i]] <- result
}

# Combine all non-null results
valid_results <- all_results[!sapply(all_results, is.null)]
if (length(valid_results) > 0) {
  combined_results <- do.call(rbind, valid_results)
  
  # Save full results
  write.csv(combined_results, paste0(result_dir, "sensitivity_analysis_direct.csv"), row.names = FALSE)
  
  # Create interpretable summary
  interpretable_summary <- combined_results %>%
    select(Model, Treatment_Effect, P_Value, Partial_R2, RV_q1, Robustness_Assessment) %>%
    mutate(
      Effect_Size = round(Treatment_Effect, 2),
      P_Value = round(P_Value, 4),
      Partial_R2 = round(Partial_R2, 3),
      Robustness_Value = round(RV_q1, 3)
    ) %>%
    select(Model, Effect_Size, P_Value, Partial_R2, Robustness_Value, Robustness_Assessment)
  
  # Save interpretable summary
  write.csv(interpretable_summary, paste0(result_dir, "sensitivity_analysis_interpretable.csv"), row.names = FALSE)
  
  # Display results
  cat("\nSensitivity Analysis Results:\n")
  print(interpretable_summary)
  
  cat("\nSensitivity analysis complete. Results saved to:",
      paste0(result_dir, "sensitivity_analysis_direct.csv"), "and",
      paste0(result_dir, "sensitivity_analysis_interpretable.csv"), "\n")
} else {
  cat("\nNo valid sensitivity analysis results were generated.\n")
}

# For comparison with sensemakr: Run one model through sensemakr and print raw structure
if (!is.null(fit_sev_1)) {
  cat("\nRunning a test with sensemakr on one model for comparison:\n")
  test_sens <- tryCatch({
    sensemakr(model = fit_sev_1, treatment = "history", benchmark_covariates = covariates_sev)
  }, error = function(e) {
    cat("Error running sensemakr:", e$message, "\n")
    NULL
  })
  
  if (!is.null(test_sens)) {
    cat("sensemakr test successful. Printing structure:\n")
    str(test_sens)
  }
}