#' Generate Balance Assessment Plots
#'
#' This function creates a comprehensive set of balance assessment plots for a MatchIt model,
#' including a love plot and individual balance plots for each covariate.
#'
#' @param match_obj A MatchIt object to assess
#' @param response_type Character string: "severity" or "recovery" to determine which variables to plot
#' @param plot_dir Directory to save plots (will be created if it doesn't exist)
#' @param subgroup_name Name of the subgroup (e.g., "Subgroup 1 (0-2 years)")
#' @return A list containing all generated plots
#' @examples
#' # For a severity model:
#' plots_sev_1 <- generate_balance_plots(best_model_sev_1, "severity", 
#'                  "/home/goldma34/fire_insect_co-occurence/plots/balance/", 
#'                  "Subgroup 1 (0-2 years)")
#'                  
#' # For a recovery model:
#' plots_rec_1 <- generate_balance_plots(best_model_rec_1, "recovery",
#'                  "/home/goldma34/fire_insect_co-occurence/plots/balance/",
#'                  "Subgroup 1 (0-2 years)")
generate_balance_plots <- function(match_obj, response_type = "severity", plot_dir, subgroup_name) {
  # Required packages
  require(cobalt)
  require(ggplot2)
  require(cowplot)
  require(gridExtra)
  
  # Create directory if it doesn't exist
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define proper variable names based on response type
  if (response_type == "severity") {
    var_df <- data.frame(
      old = c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri"),
      new = c("Host Species Percentage", "Initial Spread Index", "Drought Code", 
              "Duff Moisture Code", "Fine Fuel Moisture Code", "Build Up Index", 
              "Fire Weather Index", "Mean Terrain Ruggedness Index")
    )
    potential_vars <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")
  } else if (response_type == "recovery") {
    var_df <- data.frame(
      old = c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri"),
      new = c("Host Species Percentage", "Burn Severity (RBR)", "Mean Temperature", 
              "Total Precipitation (mm)", "Mean Terrain Ruggedness Index")
    )
    potential_vars <- c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri")
  } else {
    stop("response_type must be either 'severity' or 'recovery'")
  }
  
  # Check which variables are actually available in the model
  avail_vars <- c()
  if (!is.null(match_obj$X)) {
    avail_vars <- intersect(potential_vars, names(match_obj$X))
    cat("Available variables for plotting:", paste(avail_vars, collapse=", "), "\n")
  } else {
    avail_vars <- intersect(potential_vars, names(match_obj$data))
    cat("Using data from model$data, available variables:", paste(avail_vars, collapse=", "), "\n")
  }
  
  # Filter variable mapping to only include available variables
  var_df <- var_df[var_df$old %in% avail_vars, ]
  plot_vars <- avail_vars
  
  if (length(plot_vars) == 0) {
    warning("No matching variables found in the model data. Cannot create balance plots.")
    return(list(love_plot = NULL, balance_plots = list(), combined_plot = NULL))
  }
  
  # Clean subgroup name for filenames
  clean_name <- gsub(" |\\(|\\)", "_", subgroup_name)
  clean_name <- gsub("__", "_", clean_name)
  clean_name <- gsub("_+$", "", clean_name)
  
  # Create the love plot with error handling
  love_plt <- tryCatch({
    cobalt::love.plot(match_obj, 
                     vars = avail_vars,
                     stats = c("mean.diffs"), 
                     threshold = c(m = .25), 
                     binary = "std",
                     abs = TRUE,
                     var.order = "unadjusted",
                     var.names = setNames(var_df$new, var_df$old),
                     limits = c(0, 1),
                     grid = FALSE,
                     wrap = 20,
                     sample.names = c("Unmatched", "Matched"),
                     position = "top",
                     title = paste0("Balance for ", subgroup_name, " (", response_type, ")"),
                     shapes = c("circle", "triangle"),
                     colors = c("#FF8C00A0", "#8B0000A0"))
  }, error = function(e) {
    cat("Error creating love plot:", e$message, "\n")
    NULL
  })
  
  # Save the love plot if successful
  if (!is.null(love_plt)) {
    love_plot_path <- file.path(plot_dir, paste0(clean_name, "_", response_type, "_love_plot.pdf"))
    tryCatch({
      ggsave(love_plot_path, love_plt, width = 10, height = 8)
      cat("Saved love plot to", love_plot_path, "\n")
    }, error = function(e) {
      cat("Error saving love plot:", e$message, "\n")
    })
  }
  
  # Create individual balance plots for each variable
  bal_plots <- list()
  
  for (var in plot_vars) {
    cat("Creating balance plot for variable:", var, "\n")
    
    # Get display name for variable
    var_display_name <- var_df$new[var_df$old == var]
    if (length(var_display_name) == 0) var_display_name <- var
    
    # Create balance plot with error handling
    bal_plt <- tryCatch({
      p <- cobalt::bal.plot(match_obj, 
                          var.name = var, 
                          which = "both",
                          colors = c("#FF8C00A0", "#8B0000A0"),
                          sample.names = c("Unmatched", "Matched"))
      
      # Add customizations if the plot was created successfully
      if (!is.null(p)) {
        p <- p + scale_fill_manual(
                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
          xlab(var_display_name) + 
          ggtitle(NULL) +
          theme_bw() + 
          theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
      }
      p
    }, error = function(e) {
      cat("Error in bal.plot for", var, ":", e$message, "\n")
      NULL
    })
    
    # Save the plot if successful
    if (!is.null(bal_plt)) {
      bal_plots[[var]] <- bal_plt
      
      # Save individual plot
      indiv_plot_path <- file.path(plot_dir, paste0(clean_name, "_", response_type, "_", var, "_balance.pdf"))
      tryCatch({
        ggsave(indiv_plot_path, bal_plt, width = 8, height = 6)
        cat("Saved balance plot for", var, "to", indiv_plot_path, "\n")
      }, error = function(e) {
        cat("Error saving balance plot for", var, ":", e$message, "\n")
      })
    }
  }
  
  # Create a combined plot if any balance plots were created
  combo_plot <- NULL
  if (length(bal_plots) > 0) {
    tryCatch({
      combo_plot <- cowplot::plot_grid(plotlist = bal_plots, ncol = 2, align = 'v')
      combo_plot_path <- file.path(plot_dir, paste0(clean_name, "_", response_type, "_all_vars_balance.pdf"))
      ggsave(combo_plot_path, combo_plot, width = 14, height = 5 * ceiling(length(bal_plots)/2))
      cat("Saved combined balance plot to", combo_plot_path, "\n")
    }, error = function(e) {
      cat("Error creating combined plot:", e$message, "\n")
    })
  } else {
    cat("No balance plots were created successfully, skipping combined plot.\n")
  }
  
  # Return all plots
  return(list(
    love_plot = love_plt,
    balance_plots = bal_plots,
    combined_plot = combo_plot
  ))
}

#' Generate Balance Plots for Multiple Models
#'
#' This function applies the balance plot generation to multiple models for different subgroups
#' and response types.
#'
#' @param model_list List of MatchIt models, named by subgroup and response type
#' @param plot_dir Directory to save plots
#' @return A nested list of all generated plots
#' @examples
#' models <- list(
#'   "Subgroup 1 - Severity" = best_model_sev_1,
#'   "Subgroup 1 - Recovery" = best_model_rec_1,
#'   "Subgroup 2 - Severity" = best_model_sev_2
#' )
#' all_plots <- generate_all_balance_plots(models, "/home/goldma34/fire_insect_co-occurence/plots/balance/")
generate_all_balance_plots <- function(model_list, plot_dir) {
  results <- list()
  
  for (name in names(model_list)) {
    # Extract subgroup name and response type from the model name
    parts <- strsplit(name, " - ")[[1]]
    if (length(parts) != 2) {
      warning(paste("Skipping", name, "- expected format 'Subgroup X - ResponseType'"))
      next
    }
    
    subgroup_name <- parts[1]
    response_type <- tolower(parts[2])
    
    # Only process if model exists
    if (!is.null(model_list[[name]])) {
      cat("Generating balance plots for", name, "...\n")
      results[[name]] <- generate_balance_plots(
        model_list[[name]], 
        response_type, 
        plot_dir,
        subgroup_name
      )
    } else {
      cat("Skipping", name, "- model is NULL\n")
    }
  }
  
  return(results)
}

#' Process All Subgroups for Balance Assessment
#'
#' This function loads all saved MatchIt models and generates balance assessment plots
#' for all available models across subgroups and response types.
#'
#' @param result_dir Directory containing saved model RDS files
#' @param plot_dir Directory to save generated plots
#' @return A nested list of all generated plots
#' @examples
#' process_all_balance_plots(
#'   "/home/goldma34/fire_insect_co-occurence/data/results/", 
#'   "/home/goldma34/fire_insect_co-occurence/plots/balance/"
#' )
process_all_balance_plots <- function(result_dir, plot_dir) {
  # Define file patterns for both response types
  severity_files <- list.files(result_dir, pattern = "best_model_subgroup[0-9]+_severity\\.RDS", full.names = TRUE)
  recovery_files <- list.files(result_dir, pattern = "best_model_subgroup[0-9]+_recovery\\.RDS", full.names = TRUE)
  
  # Combine all files and extract subgroup info
  all_files <- c(severity_files, recovery_files)
  
  if (length(all_files) == 0) {
    cat("No model files found in the specified directory.\n")
    return(NULL)
  }
  
  # Load models and build the model list
  models <- list()
  subgroup_pattern <- "best_model_subgroup([0-9]+)_(severity|recovery)\\.RDS"
  time_periods <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
  
  for (file in all_files) {
    # Extract subgroup number and response type
    matches <- regexec(subgroup_pattern, basename(file))
    if (length(matches[[1]]) > 0) {
      subgroup_num <- as.numeric(regmatches(basename(file), matches)[[1]][2])
      response_type <- regmatches(basename(file), matches)[[1]][3]
      
      # Create proper name
      if (subgroup_num <= length(time_periods)) {
        time_period <- time_periods[subgroup_num]
        model_name <- paste0("Subgroup ", subgroup_num, " (", time_period, ") - ", 
                            toupper(substr(response_type, 1, 1)), 
                            substr(response_type, 2, nchar(response_type)))
      } else {
        model_name <- paste0("Subgroup ", subgroup_num, " - ", 
                            toupper(substr(response_type, 1, 1)), 
                            substr(response_type, 2, nchar(response_type)))
      }
      
      # Load the model
      tryCatch({
        models[[model_name]] <- readRDS(file)
        cat("Loaded model:", model_name, "\n")
      }, error = function(e) {
        warning(paste("Failed to load", file, ":", e$message))
      })
    }
  }
  
  # Generate plots for all loaded models
  if (length(models) > 0) {
    return(generate_all_balance_plots(models, plot_dir))
  } else {
    cat("No models could be loaded.\n")
    return(NULL)
  }
}

