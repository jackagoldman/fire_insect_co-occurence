#' Generate Treatment Effect Plots for Subgroup Models
#'
#' This function creates standardized treatment effect plots showing
#' the difference between defoliated and non-defoliated areas for
#' all your subgroup models.
#'
#' @param model_list List of fitted model objects
#' @param model_names Character vector of model names (must match model_list)
#' @param response_type Character either "Severity" or "Recovery"
#' @param output_dir Directory to save plots
#' @param combined_plot Logical; whether to also create a combined plot of all models
#' @return List of ggplot objects
#' @examples
#' # For severity models
#' sev_models <- list(fit.sev_1, fit.sev_2, fit.sev_3, fit.sev_4)
#' names <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
#' plot_treatment_effects(sev_models, names, "Severity", "plots/treat_effects/")
generate_treatment_effect_plots <- function(model_list, model_names, response_type,
                                           output_dir = "/home/goldma34/fire_insect_co-occurence/plots/treat_effects/",
                                           combined_plot = TRUE) {
  
  # Load required libraries
  require(ggplot2)
  require(dplyr)
  require(marginaleffects)
  require(patchwork)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Validate inputs
  if (length(model_list) != length(model_names)) {
    stop("Length of model_list and model_names must be the same")
  }
  
  if (!response_type %in% c("Severity", "Recovery")) {
    stop("response_type must be either 'Severity' or 'Recovery'")
  }
  
  # Set y-axis label based on response type
  y_label <- if(response_type == "Severity") {
    "Burn Severity (RBR)"
  } else {
    "Recovery Magnitude (%)"
  }
  
  # Set plot colors based on response type
  colors <- if(response_type == "Severity") {
    c("Defoliated" = "#D55E00", "Non-Defoliated" = "#0072B2")  # Orange and Blue
  } else {
    c("Defoliated" = "#009E73", "Non-Defoliated" = "#CC79A7")  # Green and Pink
  }
  
  # List to store plot objects
  plot_list <- list()
  all_predictions <- data.frame()
  
  # Generate plots for each model
  for (i in 1:length(model_list)) {
    # Skip if model is NULL
    if (is.null(model_list[[i]])) {
      cat("Skipping", model_names[i], "- model is NULL\n")
      next
    }
    
    tryCatch({
      # Create predictions
      preds <- plot_predictions(model_list[[i]], condition = c("history"), draw = FALSE)
      
      # Add model info
      preds <- preds %>% 
        mutate(Model = model_names[i]) %>%
        mutate(history = case_when(history == 1 ~ "Defoliated",
                                   history == 0 ~ "Non-Defoliated"))
      
      # Store for combined plot
      all_predictions <- bind_rows(all_predictions, preds)
      
      # Create individual plot
      p <- ggplot(preds, aes(x = history, y = estimate, color = history)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
        labs(
          title = paste(response_type, "-", model_names[i]),
          y = y_label, 
          x = "Defoliation History", 
          color = "History"
        ) +
        scale_color_manual(values = colors) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 10)
        )
      
      # Store plot object
      plot_list[[model_names[i]]] <- p
      
      # Save individual plot
      plot_filename <- paste0("fig_", 
                             tolower(response_type), 
                             "_treat_effect_", 
                             gsub(" ", "_", tolower(model_names[i])), 
                             ".png")
      
      ggsave(
        plot = p,
        filename = file.path(output_dir, plot_filename),
        width = 6, 
        height = 4, 
        dpi = 300
      )
      
      cat("Saved plot for", model_names[i], "to", file.path(output_dir, plot_filename), "\n")
      
    }, error = function(e) {
      cat("Error creating plot for", model_names[i], ":", e$message, "\n")
    })
  }
  
  # Create combined plot if requested and if we have data
  if (combined_plot && nrow(all_predictions) > 0) {
    # Create faceted plot
    p_combined <- ggplot(all_predictions, aes(x = history, y = estimate, color = history)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
      facet_wrap(~ Model, scales = "free_y", ncol = 2) +
      labs(
        title = paste(response_type, "by Time Since Defoliation"),
        y = y_label, 
        x = "Defoliation History", 
        color = "History"
      ) +
      scale_color_manual(values = colors) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    # Save combined plot
    combined_filename <- paste0("fig_", tolower(response_type), "_treat_effects_combined.png")
    ggsave(
      plot = p_combined,
      filename = file.path(output_dir, combined_filename),
      width = 10, 
      height = 8, 
      dpi = 300
    )
    
    cat("Saved combined plot to", file.path(output_dir, combined_filename), "\n")
    
    # Add to return list
    plot_list[["combined"]] <- p_combined
  }
  
  return(plot_list)
}

#' Run Treatment Effect Plots for All Subgroup Models
#'
#' This function loads all the subgroup models and generates treatment effect plots.
generate_all_treatment_effect_plots <- function() {
  # Result directory
  result_dir <- "/home/goldma34/fire_insect_co-occurence/data/results/"
  output_dir <- "/home/goldma34/fire_insect_co-occurence/plots/treat_effects/"
  
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
  
  # Generate severity plots
  cat("Generating severity treatment effect plots...\n")
  sev_plots <- generate_treatment_effect_plots(
    sev_models, 
    model_names, 
    "Severity", 
    output_dir
  )
  
  # Generate recovery plots
  cat("Generating recovery treatment effect plots...\n")
  rec_plots <- generate_treatment_effect_plots(
    rec_models, 
    model_names, 
    "Recovery", 
    output_dir
  )
  
  # Generate comprehensive plot with both responses if all models loaded
  if (all(c(!is.null(fit_sev_1), !is.null(fit_rec_1)))) {
    # Compare severity and recovery for each time period
    for (i in 1:length(model_names)) {
      if (!is.null(sev_models[[i]]) && !is.null(rec_models[[i]])) {
        comparison_plot <- sev_plots[[model_names[i]]] + rec_plots[[model_names[i]]] +
          plot_layout(ncol = 2, guides = "collect") & 
          theme(legend.position = "bottom")
        
        comparison_filename <- paste0("fig_comparison_", 
                                    gsub(" ", "_", tolower(model_names[i])), 
                                    ".png")
        
        ggsave(
          plot = comparison_plot,
          filename = file.path(output_dir, comparison_filename),
          width = 12, 
          height = 5, 
          dpi = 300
        )
      }
    }
  }
  cat("All treatment effect plots generated!\n")
}

# Run the function to generate all plots when this script is sourced
if (interactive()) {
  cat("To generate all treatment effect plots, run:\n")
  cat("generate_all_treatment_effect_plots()\n")
} else {
  # Auto-run when sourced non-interactively
  generate_all_treatment_effect_plots()
}