# Load required libraries
library(dplyr)

# Define results directory
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

# Function to extract R-squared values from r.squaredGLMM results
extract_r2 <- function(r2_obj, model_name) {
  if (is.null(r2_obj)) {
    return(NULL)
  }
  
  tryCatch({
    # Extract R-squared values (format may vary depending on model type)
    if (is.data.frame(r2_obj)) {
      # Some r.squaredGLMM results are data frames
      marginal <- r2_obj$R2m[1]  # Marginal R-squared (fixed effects only)
      conditional <- r2_obj$R2c[1]  # Conditional R-squared (fixed + random effects)
    } else if (is.numeric(r2_obj) && length(r2_obj) >= 2) {
      # Some results are numeric vectors
      marginal <- r2_obj[1]
      conditional <- r2_obj[2]
    } else {
      # Handle other cases
      marginal <- NA
      conditional <- NA
      cat("Unknown r.squaredGLMM result format for", model_name, "\n")
    }
    
    # Create a data frame row
    return(data.frame(
      Model = model_name,
      R2_Marginal = marginal,
      R2_Conditional = conditional,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat("Error extracting R² for", model_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Define model names and their corresponding R² file paths
models <- c(
  "Subgroup 1 (0-2 years) - Severity",
  "Subgroup 2 (3-5 years) - Severity",
  "Subgroup 3 (6-9 years) - Severity",
  "Subgroup 4 (10+ years) - Severity",
  "Subgroup 1 (0-2 years) - Recovery",
  "Subgroup 2 (3-5 years) - Recovery", 
  "Subgroup 3 (6-9 years) - Recovery",
  "Subgroup 4 (10+ years) - Recovery"
)

file_paths <- c(
  paste0(result_dir, "fit_model_subgroup1_severity_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup2_severity_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup3_severity_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup4_severity_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup1_recovery_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup2_recovery_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup3_recovery_r2.RDS"),
  paste0(result_dir, "fit_model_subgroup4_recovery_r2.RDS")
)

# Load all R² results and combine into a table
cat("Loading R² results...\n")
r2_results <- list()
for (i in 1:length(models)) {
  r2_obj <- safe_read_rds(file_paths[i])
  r2_row <- extract_r2(r2_obj, models[i])
  r2_results[[i]] <- r2_row
}

# Combine non-null results
r2_results <- r2_results[!sapply(r2_results, is.null)]
if (length(r2_results) > 0) {
  combined_r2 <- bind_rows(r2_results)
  
  # Add model type and time period columns for easier analysis
  combined_r2 <- combined_r2 %>%
    mutate(
      Model_Type = ifelse(grepl("Severity", Model), "Severity", "Recovery"),
      Time_Period = gsub(".*\\(([0-9]+-[0-9]+).*\\).*", "\\1", Model),
      Time_Period = factor(Time_Period, levels = c("0-2", "3-5", "6-9", "10+"))
    )
  
  # Format marginal and conditional R² values (round to 4 decimal places)
  combined_r2$R2_Marginal <- round(combined_r2$R2_Marginal, 4)
  combined_r2$R2_Conditional <- round(combined_r2$R2_Conditional, 4)
  
  # Print the combined table
  cat("\nR-squared GLMM Results Table:\n")
  print(combined_r2)
  
  # Save as CSV
  output_csv <- paste0(result_dir, "r_squared_summary.csv")
  write.csv(combined_r2, output_csv, row.names = FALSE)
  cat("\nR-squared summary saved to:", output_csv, "\n")
  
  # Create a formatted HTML table for reports
  if (requireNamespace("knitr", quietly = TRUE) && 
      requireNamespace("kableExtra", quietly = TRUE)) {
    
    cat("Creating formatted HTML table...\n")
    
    library(knitr)
    library(kableExtra)
    
    # Create a nicely formatted table
    formatted_table <- combined_r2 %>%
      select(Model, R2_Marginal, R2_Conditional) %>%
      arrange(factor(Model, levels = models))
      
    html_table <- kable(formatted_table, format = "html", 
                          caption = "R-squared Results for all Models", 
                          col.names = c("Model", 
                                        "Marginal R² (Fixed Effects)", 
                                        "Conditional R² (Fixed + Random Effects)")) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                    full_width = FALSE) %>%
      add_header_above(c(" " = 1, "R-squared GLMM Values" = 2))
    
    # Save the HTML table
    html_output <- paste0(result_dir, "r_squared_table.html")
    cat(html_table, file = html_output)
    cat("Formatted HTML table saved to:", html_output, "\n")
  }
  
  # Create summary statistics by model type
  cat("\nSummary Statistics by Model Type:\n")
  model_type_summary <- combined_r2 %>%
    group_by(Model_Type) %>%
    summarise(
      Mean_Marginal_R2 = mean(R2_Marginal, na.rm = TRUE),
      Median_Marginal_R2 = median(R2_Marginal, na.rm = TRUE),
      Min_Marginal_R2 = min(R2_Marginal, na.rm = TRUE),
      Max_Marginal_R2 = max(R2_Marginal, na.rm = TRUE),
      Mean_Conditional_R2 = mean(R2_Conditional, na.rm = TRUE),
      Median_Conditional_R2 = median(R2_Conditional, na.rm = TRUE),
      Min_Conditional_R2 = min(R2_Conditional, na.rm = TRUE),
      Max_Conditional_R2 = max(R2_Conditional, na.rm = TRUE)
    )
  
  print(model_type_summary)
  write.csv(model_type_summary, paste0(result_dir, "r_squared_summary_by_model_type.csv"), row.names = FALSE)
  
} else {
  cat("No R-squared results were found.\n")
}