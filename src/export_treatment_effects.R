#' Export Treatment Effect Results to CSV
#'
#' This function extracts treatment effect estimates, confidence intervals, and p-values
#' from model objects and exports them to a CSV file.
#'
#' @param model_list List of fitted model objects
#' @param model_names Character vector of model names (must match model_list)
#' @param response_type Character either "Severity" or "Recovery"
#' @param output_dir Directory to save CSV file
#' @return Data frame with compiled treatment effect results
#' @examples
#' # For severity models
#' sev_models <- list(fit.sev_1, fit.sev_2, fit.sev_3, fit.sev_4)
#' names <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
#' export_treatment_effects(sev_models, names, "Severity")
export_treatment_effects <- function(model_list, model_names, response_type,
                                   output_dir = "/home/goldma34/fire_insect_co-occurence/data/results/") {
  
  # Load required libraries
  require(dplyr)
  require(marginaleffects)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Validate inputs
  if (length(model_list) != length(model_names)) {
    stop("Length of model_list and model_names must be the same")
  }
  
  if (!response_type %in% c("Severity", "Recovery")) {
    stop("response_type must be either 'Severity' or 'Recovery'")
  }
  
  # Create empty data frame to store all results
  all_effects <- data.frame()
  
  # Process each model
  for (i in 1:length(model_list)) {
    # Skip if model is NULL
    if (is.null(model_list[[i]])) {
      cat("Skipping", model_names[i], "- model is NULL\n")
      next
    }
    
    tryCatch({
      # Calculate average treatment effects
      effects <- marginaleffects::avg_comparisons(
        model_list[[i]], 
        variables = "history",
        newdata = subset(model_list[[i]]$model, history == 1)
      ) 
      
      # Convert to data frame and add metadata
      effects_df <- as.data.frame(effects) %>%
        mutate(
          Model = model_names[i],
          Response = response_type,
          Subgroup = model_names[i]
        ) %>%
        # Rename columns for clarity
        rename(
          Estimate = estimate,
          Std_Error = std.error,
          Conf_Low = conf.low,
          Conf_High = conf.high,
          P_Value = p.value
        ) %>%
        # Select and reorder columns
        select(Response, Subgroup, Estimate, Std_Error, Conf_Low, Conf_High, P_Value, 
               term, contrast, everything())
      
      # Append to main results
      all_effects <- bind_rows(all_effects, effects_df)
      
      cat("Processed treatment effects for", model_names[i], "\n")
      
    }, error = function(e) {
      cat("Error processing", model_names[i], ":", e$message, "\n")
    })
  }
  
  # Save results to CSV
  if (nrow(all_effects) > 0) {
    csv_filename <- paste0(tolower(response_type), "_treatment_effects.csv")
    csv_path <- file.path(output_dir, csv_filename)
    
    write.csv(all_effects, csv_path, row.names = FALSE)
    cat("Saved treatment effects to", csv_path, "\n")
  } else {
    cat("No effects data to save for", response_type, "\n")
  }
  
  return(all_effects)
}

#' Export Treatment Effects for All Subgroup Models
#'
#' This function loads all the subgroup models and exports their treatment effects to CSV files.
export_all_treatment_effects <- function() {
  # Result directory
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
  
  # Load all model fits
  cat("Loading models...\n")
  fit_sev_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_severity.RDS"))
  fit_sev_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_severity.RDS"))
  fit_sev_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_severity.RDS"))
  fit_sev_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_severity.RDS"))
  
  fit_rec_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_recovery.RDS"))
  fit_rec_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_recovery.RDS"))
  fit_rec_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_recovery.RDS"))
  fit_rec_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_recovery.RDS"))
  
  # Define model groups
  sev_models <- list(fit_sev_1, fit_sev_2, fit_sev_3, fit_sev_4)
  rec_models <- list(fit_rec_1, fit_rec_2, fit_rec_3, fit_rec_4)
  
  model_names <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
  
  # Export severity effects
  cat("Exporting severity treatment effects...\n")
  sev_effects <- export_treatment_effects(
    sev_models, 
    model_names, 
    "Severity", 
    result_dir
  )
  
  # Export recovery effects
  cat("Exporting recovery treatment effects...\n")
  rec_effects <- export_treatment_effects(
    rec_models, 
    model_names, 
    "Recovery", 
    result_dir
  )
  
  # Combine all effects into one comprehensive file
  all_effects <- bind_rows(sev_effects, rec_effects)
  
  if (nrow(all_effects) > 0) {
    write.csv(all_effects, file.path(result_dir, "all_treatment_effects.csv"), row.names = FALSE)
    cat("Saved combined treatment effects to", file.path(result_dir, "all_treatment_effects.csv"), "\n")
  }
  
  cat("All treatment effects exported!\n")
}

# Run the function to export all treatment effects when this script is sourced
if (interactive()) {
  cat("To export all treatment effects, run:\n")
  cat("export_all_treatment_effects()\n")
} else {
  # Auto-run when sourced non-interactively
  export_all_treatment_effects()
}