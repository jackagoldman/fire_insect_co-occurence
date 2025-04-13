# Load required libraries
library(MatchIt)
library(dplyr)
library(marginaleffects)
library(MuMIn)

# Source the file with function definitions
source("/home/goldma34/fire_insect_co-occurence/src/best_model_functions.R")

#source the file with covariate balance plots
source("/home/goldma34/fire_insect_co-occurence/src/covariate_balance_plots.R")

# Load your data 
source("/home/goldma34/fire_insect_co-occurence/src/load_data.R")  # Replace with actual script that loads hist_gt90_1

# Check if data is loaded
if (!exists("hist_gt90_1")) {
  stop("Data not loaded. Make sure hist_gt90_1 exists.")
}
if (!exists("hist_gt90_2")) {
  stop("Data not loaded. Make sure hist_gt90_2 exists.")
}
if (!exists("hist_gt90_3")) {
  stop("Data not loaded. Make sure hist_gt90_3 exists.")
}
if (!exists("hist_gt90_4")) {
  stop("Data not loaded. Make sure hist_gt90_4 exists.")
}

# Define parameters
methods <- c("nearest", "optimal")
distances <- c("glm", "mahalanobis")
links <- c("logit", "probit")
m_orders <- c("random")
calipers <- c(0.1, 0.2, 0.3)
replace_list <- c(FALSE)
use_mahvars_list <- c(FALSE, TRUE)

#============================
# SUBGROUP 1 - SEVERITY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 1 - SEVERITY ANALYSIS =======\n")
cat("Running model selection for subgroup 1 (severity)...\n")
best_model_info_sev_1 <- find_best_model(hist_gt90_1, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev_1 <- best_model_info_sev_1$best_model
all_results_sev_1 <- best_model_info_sev_1$results

# Check if we found a best model
if (!is.null(best_model_sev_1)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev_1 <- match_data(best_model_sev_1)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev_1 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev_1 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_1, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev_1))

  #rsquared
  cat("\nR-squared:\n")
  fit.sev_1r2 <- r.squaredGLMM(fit.sev_1)
  print(fit.sev_1r2)
  

  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev_1, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup1_severity.RDS")
  saveRDS(fit.sev_1, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup1_severity.RDS")
  saveRDS(fit.sev_1r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup1_severity_r2.RDS")
} else {
  cat("No suitable severity model found for subgroup 1. Try different parameters.\n")
}

#============================
# SUBGROUP 1 - RECOVERY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 1 - RECOVERY ANALYSIS =======\n")
cat("Running model selection for subgroup 1 (recovery)...\n")
best_model_info_rec_1 <- find_best_model(hist_gt90_1, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec_1 <- best_model_info_rec_1$best_model
all_results_rec_1 <- best_model_info_rec_1$results

# Check if we found a best model
if (!is.null(best_model_rec_1)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec_1 <- match_data(best_model_rec_1)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec_1 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec_1 <- lm(rbr_w_offset ~ history + host_pct + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_1, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec_1))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec_1, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup1_recovery.RDS")
  saveRDS(fit.rec_1, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup1_recovery.RDS")
} else {
  cat("No suitable recovery model found for subgroup 1. Try different parameters.\n")
}

#============================
# SUBGROUP 2 - SEVERITY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 2 - SEVERITY ANALYSIS =======\n")
cat("Running model selection for subgroup 2 (severity)...\n")
best_model_info_sev_2 <- find_best_model(hist_gt90_2, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev_2 <- best_model_info_sev_2$best_model
all_results_sev_2 <- best_model_info_sev_2$results

# Check if we found a best model
if (!is.null(best_model_sev_2)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev_2 <- match_data(best_model_sev_2)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev_2 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev_2 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_2, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev_2))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev_2, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup2_severity.RDS")
  saveRDS(fit.sev_2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup2_severity.RDS")
} else {
  cat("No suitable severity model found for subgroup 2. Try different parameters.\n")
}

#============================
# SUBGROUP 2 - RECOVERY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 2 - RECOVERY ANALYSIS =======\n")
cat("Running model selection for subgroup 2 (recovery)...\n")
best_model_info_rec_2 <- find_best_model(hist_gt90_2, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec_2 <- best_model_info_rec_2$best_model
all_results_rec_2 <- best_model_info_rec_2$results

# Check if we found a best model
if (!is.null(best_model_rec_2)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec_2 <- match_data(best_model_rec_2)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec_2 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec_2 <- lm(rbr_w_offset ~ history + host_pct + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_2, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec_2))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec_2, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup2_recovery.RDS")
  saveRDS(fit.rec_2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup2_recovery.RDS")
} else {
  cat("No suitable recovery model found for subgroup 2. Try different parameters.\n")
}

#============================
# SUBGROUP 3 - SEVERITY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 3 - SEVERITY ANALYSIS =======\n")
cat("Running model selection for subgroup 3 (severity)...\n")
best_model_info_sev_3 <- find_best_model(hist_gt90_3, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev_3 <- best_model_info_sev_3$best_model
all_results_sev_3 <- best_model_info_sev_3$results

# Check if we found a best model
if (!is.null(best_model_sev_3)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev_3 <- match_data(best_model_sev_3)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev_3 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev_3 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_3, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev_3))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev_3, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup3_severity.RDS")
  saveRDS(fit.sev_3, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup3_severity.RDS")
} else {
  cat("No suitable severity model found for subgroup 3. Try different parameters.\n")
}

#============================
# SUBGROUP 3 - RECOVERY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 3 - RECOVERY ANALYSIS =======\n")
cat("Running model selection for subgroup 3 (recovery)...\n")
best_model_info_rec_3 <- find_best_model(hist_gt90_3, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec_3 <- best_model_info_rec_3$best_model
all_results_rec_3 <- best_model_info_rec_3$results

# Check if we found a best model
if (!is.null(best_model_rec_3)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec_3 <- match_data(best_model_rec_3)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec_3 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec_3 <- lm(rbr_w_offset ~ history + host_pct + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_3, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec_3))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec_3, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup3_recovery.RDS")
  saveRDS(fit.rec_3, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup3_recovery.RDS")
} else {
  cat("No suitable recovery model found for subgroup 3. Try different parameters.\n")
}

#============================
# SUBGROUP 4 - SEVERITY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 4 - SEVERITY ANALYSIS =======\n")
cat("Running model selection for subgroup 4 (severity)...\n")
best_model_info_sev_4 <- find_best_model(hist_gt90_4, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev_4 <- best_model_info_sev_4$best_model
all_results_sev_4 <- best_model_info_sev_4$results

# Check if we found a best model
if (!is.null(best_model_sev_4)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev_4 <- match_data(best_model_sev_4)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev_4 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev_4 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_4, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev_4))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev_4, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup4_severity.RDS")
  saveRDS(fit.sev_4, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup4_severity.RDS")
} else {
  cat("No suitable severity model found for subgroup 4. Try different parameters.\n")
}

#============================
# SUBGROUP 4 - RECOVERY ANALYSIS
#============================
cat("\n\n======= SUBGROUP 4 - RECOVERY ANALYSIS =======\n")
cat("Running model selection for subgroup 4 (recovery)...\n")
best_model_info_rec_4 <- find_best_model(hist_gt90_4, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec_4 <- best_model_info_rec_4$best_model
all_results_rec_4 <- best_model_info_rec_4$results

# Check if we found a best model
if (!is.null(best_model_rec_4)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec_4 <- match_data(best_model_rec_4)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec_4 %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec_4 <- lm(rbr_w_offset ~ history + host_pct + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_4, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec_4))
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec_4, "/home/goldma34/fire_insect_co-occurence/data/results/best_model_subgroup4_recovery.RDS")
  saveRDS(fit.rec_4, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup4_recovery.RDS")
} else {
  cat("No suitable recovery model found for subgroup 4. Try different parameters.\n")
}

#============================
# COMBINING RESULTS
#============================
cat("\n\n======= COMBINED RESULTS =======\n")

# Create a summary of all subgroups - both severity and recovery
subgroups_summary <- data.frame(
  Subgroup = rep(c("Subgroup 1", "Subgroup 2", "Subgroup 3", "Subgroup 4"), 2),
  Response_Type = c(rep("Severity", 4), rep("Recovery", 4)),
  Time_Period = rep(c("0-2 years", "3-5 years", "6-9 years", "10+ years"), 2),
  Model_Found = c(
    !is.null(best_model_sev_1), !is.null(best_model_sev_2), 
    !is.null(best_model_sev_3), !is.null(best_model_sev_4),
    !is.null(best_model_rec_1), !is.null(best_model_rec_2), 
    !is.null(best_model_rec_3), !is.null(best_model_rec_4)
  ),
  stringsAsFactors = FALSE
)

print(subgroups_summary)

# Save a combined summary of all results
write.csv(subgroups_summary, "/home/goldma34/fire_insect_co-occurence/data/results/subgroups_summary.csv", row.names = FALSE)

#Collect model coefficients for all successful models
model_coefs <- list()
# Severity models
if (!is.null(best_model_sev_1)) model_coefs[["Subgroup 1 - Severity"]] <- coef(summary(fit.sev_1))
if (!is.null(best_model_sev_2)) model_coefs[["Subgroup 2 - Severity"]] <- coef(summary(fit.sev_2))
if (!is.null(best_model_sev_3)) model_coefs[["Subgroup 3 - Severity"]] <- coef(summary(fit.sev_3))
if (!is.null(best_model_sev_4)) model_coefs[["Subgroup 4 - Severity"]] <- coef(summary(fit.sev_4))
# Recovery models
if (!is.null(best_model_rec_1)) model_coefs[["Subgroup 1 - Recovery"]] <- coef(summary(fit.rec_1))
if (!is.null(best_model_rec_2)) model_coefs[["Subgroup 2 - Recovery"]] <- coef(summary(fit.rec_2))
if (!is.null(best_model_rec_3)) model_coefs[["Subgroup 3 - Recovery"]] <- coef(summary(fit.rec_3))
if (!is.null(best_model_rec_4)) model_coefs[["Subgroup 4 - Recovery"]] <- coef(summary(fit.rec_4))

# Save the coefficient results
saveRDS(model_coefs, "/home/goldma34/fire_insect_co-occurence/data/results/all_model_coefficients.RDS")

#============================
# DETAILED MODEL SUMMARY
#============================
cat("\n\n======= DETAILED MODEL SUMMARY =======\n")

# Define paths and names for models
result_dir <- "/home/goldma34/fire_insect_co-occurence/data/results/"
model_paths <- c(
  # Severity models
  paste0(result_dir, "best_model_subgroup1_severity.RDS"),
  paste0(result_dir, "best_model_subgroup2_severity.RDS"),
  paste0(result_dir, "best_model_subgroup3_severity.RDS"),
  paste0(result_dir, "best_model_subgroup4_severity.RDS"),
  # Recovery models
  paste0(result_dir, "best_model_subgroup1_recovery.RDS"),
  paste0(result_dir, "best_model_subgroup2_recovery.RDS"),
  paste0(result_dir, "best_model_subgroup3_recovery.RDS"),
  paste0(result_dir, "best_model_subgroup4_recovery.RDS")
)

fit_paths <- c(
  # Severity models
  paste0(result_dir, "fit_model_subgroup1_severity.RDS"),
  paste0(result_dir, "fit_model_subgroup2_severity.RDS"),
  paste0(result_dir, "fit_model_subgroup3_severity.RDS"),
  paste0(result_dir, "fit_model_subgroup4_severity.RDS"),
  # Recovery models
  paste0(result_dir, "fit_model_subgroup1_recovery.RDS"),
  paste0(result_dir, "fit_model_subgroup2_recovery.RDS"),
  paste0(result_dir, "fit_model_subgroup3_recovery.RDS"),
  paste0(result_dir, "fit_model_subgroup4_recovery.RDS")
)

subgroup_names <- c(
  # Severity models
  "Subgroup 1 (0-2 years) - Severity",
  "Subgroup 2 (3-5 years) - Severity",
  "Subgroup 3 (6-9 years) - Severity",
  "Subgroup 4 (10+ years) - Severity",
  # Recovery models
  "Subgroup 1 (0-2 years) - Recovery",
  "Subgroup 2 (3-5 years) - Recovery",
  "Subgroup 3 (6-9 years) - Recovery",
  "Subgroup 4 (10+ years) - Recovery"
)

# Generate detailed model summary
detailed_summary <- generate_model_summary(
  model_paths,
  fit_paths,
  subgroup_names,
  paste0(result_dir, "detailed_model_summary.csv")
)

# Print summarized version of the detailed summary
cat("\nModel parameters summary:\n")
print(detailed_summary[, c("Subgroup", "Method", "Distance", "Link", "Caliper", "Response", "Sample_Size_Treated", "Sample_Size_Control")])

#============================
# GENERATE MODEL COMPARISON PLOTS
#============================
cat("\n\n======= GENERATING MODEL COMPARISON PLOTS =======\n")

# Try to generate plots if required packages are available
tryCatch({
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(scales)  # For comma formatting
  
  # Create output directory for plots
  plots_dir <- "/home/goldma34/fire_insect_co-occurence/plots/"
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Filter only successful models
  plot_data <- detailed_summary[detailed_summary$Model_Found == TRUE, ]
  
  if(nrow(plot_data) > 0) {
    # Extract response type from Subgroup column
    plot_data$Response <- ifelse(grepl("Recovery", plot_data$Subgroup), "Recovery", "Severity")
    plot_data$Subgroup_Clean <- gsub(" - (Severity|Recovery)", "", plot_data$Subgroup)
    plot_data$Subgroup_Short <- gsub("Subgroup ([0-9]) \\(([0-9]+-[0-9]+) years\\).*", "\\2 yrs", plot_data$Subgroup_Clean)
    
    # Function to add value labels to bars
    add_bar_labels <- function(p) {
      p + geom_text(aes(label = comma(after_stat(y))), 
                   stat = "identity", 
                   vjust = -0.5, 
                   size = 3.5)
    }
    
    #-----------------------------------------
    # 1. TOTAL SAMPLE SIZE BY RESPONSE TYPE
    #-----------------------------------------
    p1 <- ggplot(plot_data, aes(x = Subgroup_Short, y = Total_Sample_Size, fill = Response)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = comma(Total_Sample_Size)), 
                position = position_dodge(width = 0.9), 
                vjust = -0.5, 
                size = 3.5) +
      theme_minimal() +
      labs(title = "Total Sample Size by Subgroup and Response Type", 
           x = "Time Since Defoliation", 
           y = "Number of Fires",
           fill = "Model Type") +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        legend.position = "top"
      ) +
      scale_fill_manual(values = c("Recovery" = "#2ca25f", "Severity" = "#8856a7"))
    
    #-----------------------------------------
    # 2. TREATMENT GROUPS - SEVERITY ONLY
    #-----------------------------------------
    # Filter data for severity models
    severity_data_long <- melt(
      plot_data[plot_data$Response == "Severity", ], 
      id.vars = c("Subgroup_Clean", "Subgroup_Short"),
      measure.vars = c("Sample_Size_Treated", "Sample_Size_Control"),
      variable.name = "Group",
      value.name = "Count"
    )
    
    # Only create plot if we have severity data
    if(nrow(severity_data_long) > 0) {
      p_sev_groups <- ggplot(severity_data_long, aes(x = Subgroup_Short, y = Count, fill = Group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        geom_text(aes(label = comma(Count)), 
                  position = position_dodge(width = 0.9), 
                  vjust = -0.5, 
                  size = 3.5) +
        theme_minimal() +
        labs(title = "Severity Models: Sample Sizes by Treatment Group", 
             subtitle = "Comparison between defoliated and non-defoliated areas",
             x = "Time Since Defoliation", 
             y = "Number of Fires") +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid.major.x = element_blank(),
          legend.position = "top"
        ) +
        scale_fill_manual(
          values = c("Sample_Size_Treated" = "#e31a1c", "Sample_Size_Control" = "#1f78b4"),
          labels = c("Sample_Size_Treated" = "Defoliated", "Sample_Size_Control" = "Non-Defoliated")
        )
    } else {
      p_sev_groups <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "No severity models found") + 
        theme_void()
    }
    
    #-----------------------------------------
    # 3. TREATMENT GROUPS - RECOVERY ONLY
    #-----------------------------------------
    # Filter data for recovery models
    recovery_data_long <- melt(
      plot_data[plot_data$Response == "Recovery", ], 
      id.vars = c("Subgroup_Clean", "Subgroup_Short"),
      measure.vars = c("Sample_Size_Treated", "Sample_Size_Control"),
      variable.name = "Group",
      value.name = "Count"
    )
    
    # Only create plot if we have recovery data
    if(nrow(recovery_data_long) > 0) {
      p_rec_groups <- ggplot(recovery_data_long, aes(x = Subgroup_Short, y = Count, fill = Group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        geom_text(aes(label = comma(Count)), 
                  position = position_dodge(width = 0.9), 
                  vjust = -0.5, 
                  size = 3.5) +
        theme_minimal() +
        labs(title = "Recovery Models: Sample Sizes by Treatment Group", 
             subtitle = "Comparison between defoliated and non-defoliated areas",
             x = "Time Since Defoliation", 
             y = "Number of Fires") +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid.major.x = element_blank(),
          legend.position = "top"
        ) +
        scale_fill_manual(
          values = c("Sample_Size_Treated" = "#33a02c", "Sample_Size_Control" = "#a6cee3"),
          labels = c("Sample_Size_Treated" = "Defoliated", "Sample_Size_Control" = "Non-Defoliated")
        )
    } else {
      p_rec_groups <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "No recovery models found") + 
        theme_void()
    }
    
    #-----------------------------------------
    # 4. SEVERITY ONLY - TOTAL SAMPLE SIZE
    #-----------------------------------------
    severity_data <- plot_data[plot_data$Response == "Severity", ]
    if(nrow(severity_data) > 0) {
      p_sev_total <- ggplot(severity_data, aes(x = Subgroup_Short, y = Total_Sample_Size)) +
        geom_bar(stat = "identity", fill = "#8856a7") +
        geom_text(aes(label = comma(Total_Sample_Size)), 
                  vjust = -0.5, 
                  size = 3.5) +
        theme_minimal() +
        labs(title = "Severity Models: Total Sample Size by Subgroup", 
             x = "Time Since Defoliation", 
             y = "Number of Fires") +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid.major.x = element_blank()
        )
    }
    
    #-----------------------------------------
    # 5. RECOVERY ONLY - TOTAL SAMPLE SIZE
    #-----------------------------------------
    recovery_data <- plot_data[plot_data$Response == "Recovery", ]
    if(nrow(recovery_data) > 0) {
      p_rec_total <- ggplot(recovery_data, aes(x = Subgroup_Short, y = Total_Sample_Size)) +
        geom_bar(stat = "identity", fill = "#2ca25f") +
        geom_text(aes(label = comma(Total_Sample_Size)), 
                  vjust = -0.5, 
                  size = 3.5) +
        theme_minimal() +
        labs(title = "Recovery Models: Total Sample Size by Subgroup", 
             x = "Time Since Defoliation", 
             y = "Number of Fires") +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid.major.x = element_blank()
        )
    }
    
    #-----------------------------------------
    # SAVE ALL PLOTS
    #-----------------------------------------
    # Save combined plots
    pdf(paste0(plots_dir, "sample_size_summary.pdf"), width = 11, height = 12)
    if(nrow(plot_data) > 0) {
      grid.arrange(
        p1,
        arrangeGrob(p_sev_groups, p_rec_groups, ncol = 1),
        nrow = 2,
        heights = c(1, 2)
      )
    }
    dev.off()
    
    # Save individual plots
    ggsave(paste0(plots_dir, "total_sample_size_by_response.pdf"), p1, width = 10, height = 6)
    
    if(nrow(severity_data) > 0) {
      ggsave(paste0(plots_dir, "severity_total_sample_sizes.pdf"), p_sev_total, width = 8, height = 6)
      ggsave(paste0(plots_dir, "severity_treatment_groups.pdf"), p_sev_groups, width = 8, height = 6)
    }
    
    if(nrow(recovery_data) > 0) {
      ggsave(paste0(plots_dir, "recovery_total_sample_sizes.pdf"), p_rec_total, width = 8, height = 6)
      ggsave(paste0(plots_dir, "recovery_treatment_groups.pdf"), p_rec_groups, width = 8, height = 6)
    }
    
    cat("Plots saved to", plots_dir, "\n")
  } else {
    cat("No successful models to plot\n")
  }
}, error = function(e) {
  cat("Could not generate plots:", e$message, "\n")
  cat("Make sure ggplot2, gridExtra, reshape2, and scales packages are installed\n")
})

cat("\nAnalysis complete. Results saved to data/results directory.\n")


#============================
# GENERATE BALANCE ASSESSMENT PLOTS
#============================
cat("\n\n======= GENERATING BALANCE ASSESSMENT PLOTS =======\n")

plot_dir <- "/home/goldma34/fire_insect_co-occurence/plots/balance/"
result_dir <- "/home/goldma34/fire_insect_co-occurence/data/results/"

# Create a list of available models
model_list <- list()

# Add severity models if available
if (!is.null(best_model_sev_1)) model_list[["Subgroup 1 (0-2 years) - Severity"]] <- best_model_sev_1
if (!is.null(best_model_sev_2)) model_list[["Subgroup 2 (3-5 years) - Severity"]] <- best_model_sev_2
if (!is.null(best_model_sev_3)) model_list[["Subgroup 3 (6-9 years) - Severity"]] <- best_model_sev_3
if (!is.null(best_model_sev_4)) model_list[["Subgroup 4 (10+ years) - Severity"]] <- best_model_sev_4

# Add recovery models if available
if (!is.null(best_model_rec_1)) model_list[["Subgroup 1 (0-2 years) - Recovery"]] <- best_model_rec_1
if (!is.null(best_model_rec_2)) model_list[["Subgroup 2 (3-5 years) - Recovery"]] <- best_model_rec_2
if (!is.null(best_model_rec_3)) model_list[["Subgroup 3 (6-9 years) - Recovery"]] <- best_model_rec_3
if (!is.null(best_model_rec_4)) model_list[["Subgroup 4 (10+ years) - Recovery"]] <- best_model_rec_4

# Generate balance plots for all available models
if (length(model_list) > 0) {
  cat("Generating balance plots for", length(model_list), "models...\n")
  balance_plots <- generate_all_balance_plots(model_list, plot_dir)
  cat("Balance plot generation complete. Plots saved to", plot_dir, "\n")
} else {
  cat("No models available to generate balance plots.\n")
}

#============================
# GENERATE AND SAVE R2SQUAREDGLMM
#============================

# Add R-squared GLMM for Subgroup 2 severity
if (!is.null(best_model_sev_1)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.sev_1r2 <- r.squaredGLMM(fit.sev_1)
  print(fit.sev_1r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.sev_1r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup1_severity_r2.RDS")
}


# Add R-squared GLMM for Subgroup 2 severity
if (!is.null(best_model_sev_2)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.sev_2r2 <- r.squaredGLMM(fit.sev_2)
  print(fit.sev_2r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.sev_2r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup2_severity_r2.RDS")
}

# Add R-squared GLMM for Subgroup 3 severity
if (!is.null(best_model_sev_3)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.sev_3r2 <- r.squaredGLMM(fit.sev_3)
  print(fit.sev_3r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.sev_3r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup3_severity_r2.RDS")
}

# Add R-squared GLMM for Subgroup 4 severity
if (!is.null(best_model_sev_4)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.sev_4r2 <- r.squaredGLMM(fit.sev_4)
  print(fit.sev_4r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.sev_4r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup4_severity_r2.RDS")
}

# Add R-squared GLMM for Subgroup 1 recovery
if (!is.null(best_model_rec_1)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.rec_1r2 <- r.squaredGLMM(fit.rec_1)
  print(fit.rec_1r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.rec_1r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup1_recovery_r2.RDS")
}

# Add R-squared GLMM for Subgroup 2 recovery
if (!is.null(best_model_rec_2)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.rec_2r2 <- r.squaredGLMM(fit.rec_2)
  print(fit.rec_2r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.rec_2r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup2_recovery_r2.RDS")
}

# Add R-squared GLMM for Subgroup 3 recovery
if (!is.null(best_model_rec_3)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.rec_3r2 <- r.squaredGLMM(fit.rec_3)
  print(fit.rec_3r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.rec_3r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup3_recovery_r2.RDS")
}

# Add R-squared GLMM for Subgroup 4 recovery
if (!is.null(best_model_rec_4)) {
  # After printing model summary
  cat("\nR-squared:\n")
  fit.rec_4r2 <- r.squaredGLMM(fit.rec_4)
  print(fit.rec_4r2)
  
  # Add to the saveRDS statements
  saveRDS(fit.rec_4r2, "/home/goldma34/fire_insect_co-occurence/data/results/fit_model_subgroup4_recovery_r2.RDS")
}