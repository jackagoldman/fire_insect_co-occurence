#' Adjust MatchIt Parameters
#'
#' This function adjusts the parameters of the MatchIt model based on the provided inputs.
#'
#' @param data A data frame containing the dataset.
#' @param response A character string specifying the response type ("severity" or "recovery").
#' @param method A character string specifying the matching method. Default is "nearest".
#' @param distance A character string specifying the distance measure. Default is "glm".
#' @param link A character string specifying the link function. Default is "logit".
#' @param m_order A character string specifying the matching order. Default is "random".
#' @param caliper A numeric value specifying the maximum distance for matches. Default is 0.2.
#' @param replace A logical value indicating whether matches should be with replacement. Default is TRUE.
#' @param use_mahvars A logical value indicating whether to use mahvars. Default is FALSE.
#' @return A MatchIt object containing the matching results.
adjust_matchit_params <- function(data, response = "severity", method = "nearest", distance = "glm", link = "logit", m_order = "random", caliper = 0.2, replace = TRUE, use_mahvars = FALSE) {
  # Add error handling for undefined variables
  tryCatch({
    # Define covariates based on response type
    if (response == "severity") {
      covariates <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")
    } else if (response == "recovery") {
      covariates <- c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri")
    } else {
      stop("Invalid response argument. Use 'severity' or 'recovery'.")
    }
    
    # Create model formula dynamically
    model_formula <- as.formula(paste("history ~", paste(covariates, collapse = " + ")))
    
    # Define mahvars based on use_mahvars and distance and response type
    mahvars <- if (use_mahvars && distance != "mahalanobis") {
      as.formula(paste("~", paste(covariates, collapse = " + ")))
    } else {
      NULL
    }
    
    # When using mahalanobis distance and a caliper, we need named calipers
    if (distance == "mahalanobis" && !is.null(caliper)) {
      # Handle the case when mahvars is NULL but distance is mahalanobis
      if (is.null(mahvars)) {
        # Use the covariates directly since mahvars is NULL
        caliper <- setNames(rep(caliper, length(covariates)), covariates)
      } else {
        # Use mahvars
        caliper <- setNames(rep(caliper, length(all.vars(mahvars))), all.vars(mahvars))
      }
    }
    
    # Fit the model with the dynamic formula
    model <- matchit(model_formula,
                    data = data,
                    method = method,
                    distance = distance,
                    link = link,
                    m.order = m_order,
                    caliper = caliper,
                    replace = replace,
                    mahvars = mahvars)
    return(model)
  }, error = function(e) {
    message("Error in adjust_matchit_params: ", e$message)
    return(NULL)
  })
}

#' Find Best MatchIt Model
#'
#' This function iterates over all combinations of the given parameters, stores the results in a list, and selects the best model based on the standard mean difference (SMD) being below 0.25.
#'
#' @param data A data frame containing the dataset.
#' @param response A character string specifying the response type ("severity" or "recovery").
#' @param methods A character vector specifying the matching methods to try.
#' @param distances A character vector specifying the distance measures to try.
#' @param links A character vector specifying the link functions to try.
#' @param m_orders A character vector specifying the matching orders to try.
#' @param calipers A numeric vector specifying the maximum distances for matches to try.
#' @param replace_list A list of logical values indicating whether matches should be with replacement.
#' @param use_mahvars_list A list of logical values indicating whether to use mahvars.
#' @return A list containing the best MatchIt model and all results.
find_best_model <- function(data, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list) {
  # Add validation for inputs
  if (missing(data) || is.null(data) || nrow(data) == 0) {
    stop("Input data is missing, NULL, or empty")
  }
  
  # Validate response parameter
  if (!response %in% c("severity", "recovery")) {
    stop("Response must be either 'severity' or 'recovery'")
  }
  
  # Initialize results
  results <- list()
  best_model <- NULL
  best_smd <- Inf
  
  # Add progress tracking
  total_combos <- length(methods) * length(distances) * length(links) * 
                 length(m_orders) * length(calipers) * length(replace_list) * 
                 length(use_mahvars_list)
  cat("Testing", total_combos, "parameter combinations for", response, "response\n")
  count <- 0
  
  # Nested loops for all parameter combinations
  for (method in methods) {
    for (distance in distances) {
      for (link in links) {
        for (m_order in m_orders) {
          for (caliper in calipers) {
            for (replace in replace_list) {
              for (use_mahvars in use_mahvars_list) {
                count <- count + 1
                cat("Testing combination", count, "of", total_combos, ":\n")
                cat("  method =", method, ", distance =", distance, ", link =", link, "\n")
                cat("  m_order =", m_order, ", caliper =", caliper, ", replace =", replace, ", use_mahvars =", use_mahvars, "\n")
                
                # Try to fit model with error handling
                model <- adjust_matchit_params(data, response, method, distance, link, m_order, caliper, replace, use_mahvars)
                
                # Skip to next iteration if model is NULL
                if (is.null(model)) {
                  cat("  Failed to fit model with this combination\n")
                  next
                }
                
                # Try to get model summary with error handling
                tryCatch({
                  model_summary <- summary(model)
                  
                  if (is.list(model_summary) && "sum.matched" %in% names(model_summary)) {
                    smd <- as.data.frame(model_summary$sum.matched)$`Std. Mean Diff.`
                    
                    cat("  Mean SMD:", mean(smd), "\n")
                    
                    if (all(smd < 0.25)) {
                      cat("  All SMD < 0.25, adding to results\n")
                      model_key <- paste(method, distance, link, m_order, caliper, replace, use_mahvars, sep = "_")
                      results[[model_key]] <- model
                      
                      # Check if this is the best model so far
                      if (mean(smd) < best_smd) {
                        best_smd <- mean(smd)
                        best_model <- model
                        # Add model info for easier reference later
                        best_model$info <- list(
                          response = response,
                          method = method,
                          distance = distance,
                          link = link,
                          m_order = m_order,
                          caliper = caliper,
                          replace = replace,
                          mahalanobis = as.character(use_mahvars)
                        )
                        cat("  New best model! Mean SMD =", best_smd, "\n")
                      }
                    } else {
                      cat("  Some SMD values exceed 0.25, skipping\n")
                    }
                  } else {
                    cat("  Model summary does not contain 'sum.matched', skipping\n")
                  }
                }, error = function(e) {
                  cat("  Error in model summary:", e$message, "\n")
                })
              }
            }
          }
        }
      }
    }
  }
  
  # Report results
  if (is.null(best_model)) {
    cat("No suitable model found! Try different parameters.\n")
  } else {
    cat("\nBest model found with mean SMD:", best_smd, "\n")
    cat("Best model parameters:\n")
    print(best_model$info)
  }
  
  return(list(best_model = best_model, results = results))
}

#' Extract Model Information
#'
#' This function extracts key information from a MatchIt model and its corresponding
#' linear model fit, including parameters used and sample sizes.
#'
#' @param model A MatchIt model object, or path to a saved MatchIt model RDS file.
#' @param model_fit An lm model object fit on matched data, or path to a saved model fit RDS file.
#' @param subgroup_name A character string identifying the subgroup (e.g., "Subgroup 1 (0-2 years)").
#' @return A data frame containing model parameters and sample size information.
#' @examples
#' # With model objects
#' info <- extract_model_info(best_model_1, fit.w1, "Subgroup 1 (0-2 years)")
#' 
#' # With file paths
#' info <- extract_model_info(
#'   "/path/to/best_model_1.RDS", 
#'   "/path/to/fit_model_1.RDS", 
#'   "Subgroup 1 (0-2 years)"
#' )
extract_model_info <- function(model, model_fit, subgroup_name) {
  require(dplyr)
  require(MatchIt)
  
  # Check if model is a file path and load it
  if (is.character(model)) {
    if (!file.exists(model)) {
      return(data.frame(
        Subgroup = subgroup_name,
        Model_Found = FALSE,
        Method = NA,
        Distance = NA,
        Link = NA,
        M_Order = NA,
        Caliper = NA,
        Replace = NA,
        Use_Mahvars = NA,
        Response = NA,
        Sample_Size_Treated = NA,
        Sample_Size_Control = NA,
        Total_Sample_Size = NA,
        stringsAsFactors = FALSE
      ))
    }
    model <- readRDS(model)
  }
  
  # Check if model_fit is a file path and load it
  if (is.character(model_fit)) {
    if (file.exists(model_fit)) {
      model_fit <- readRDS(model_fit)
    } else {
      model_fit <- NULL
    }
  }
  
  # Get matched data
  m.data <- match_data(model)
  
  # Extract sample sizes
  sample_sizes <- m.data %>%
    group_by(history) %>%
    summarise(Count = n()) %>%
    ungroup()
  
  # Get sample sizes for treated (history=1) and control (history=0)
  treated_size <- sample_sizes %>% filter(history == 1) %>% pull(Count)
  control_size <- sample_sizes %>% filter(history == 0) %>% pull(Count)
  
  if (length(treated_size) == 0) treated_size <- 0
  if (length(control_size) == 0) control_size <- 0
  
  total_size <- nrow(m.data)
  
  # Extract model parameters from info if available
  if (!is.null(model$info)) {
    # Use the stored info
    response <- if (is.null(model$info$response)) NA else model$info$response
    method <- if (is.null(model$info$method)) NA else model$info$method
    distance <- if (is.null(model$info$distance)) NA else model$info$distance
    link <- if (is.null(model$info$link)) NA else model$info$link
    m_order <- if (is.null(model$info$m_order)) NA else model$info$m_order
    caliper <- if (is.null(model$info$caliper)) NA else model$info$caliper
    replace <- if (is.null(model$info$replace)) NA else model$info$replace
    use_mahvars <- if (is.null(model$info$mahalanobis)) NA else model$info$mahalanobis
  } else {
    # Try to extract from model object directly
    response <- NA  # Cannot determine reliably without info
    method <- if (is.null(model$method)) NA else model$method
    distance <- if (is.null(model$distance)) NA else model$distance
    link <- if (is.null(model$link)) NA else model$link
    m_order <- if (is.null(model$m.order)) NA else model$m.order
    caliper <- if (is.null(model$caliper)) NA else {
      if (length(model$caliper) > 1) model$caliper[1] else model$caliper
    }
    replace <- if (is.null(model$replace)) NA else model$replace
    use_mahvars <- if (is.null(model$mahvars)) "FALSE" else "TRUE"
  }
  
  # Create result data frame
  result <- data.frame(
    Subgroup = subgroup_name,
    Model_Found = TRUE,
    Method = method,
    Distance = distance,
    Link = link,
    M_Order = m_order,
    Caliper = caliper,
    Replace = replace,
    Use_Mahvars = use_mahvars,
    Response = response,
    Sample_Size_Treated = treated_size,
    Sample_Size_Control = control_size,
    Total_Sample_Size = total_size,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Generate Model Summary Report
#'
#' This function creates a detailed summary report of all models from different subgroups.
#'
#' @param model_paths A vector of paths to saved MatchIt model RDS files.
#' @param fit_paths A vector of paths to saved lm model fit RDS files.
#' @param subgroup_names A vector of names identifying each subgroup.
#' @param output_file Path where the summary CSV file should be saved.
#' @return A data frame containing the combined model information for all subgroups.
#' @examples
#' model_paths <- c(
#'   "data/results/best_model_subgroup1.RDS",
#'   "data/results/best_model_subgroup2.RDS"
#' )
#' fit_paths <- c(
#'   "data/results/fit_model_subgroup1.RDS",
#'   "data/results/fit_model_subgroup2.RDS"
#' )
#' subgroup_names <- c("Subgroup 1 (0-2 years)", "Subgroup 2 (3-5 years)")
#' summary <- generate_model_summary(model_paths, fit_paths, subgroup_names, "data/results/model_summary.csv")
generate_model_summary <- function(model_paths, fit_paths, subgroup_names, output_file) {
  # Validate inputs
  if (length(model_paths) != length(fit_paths) || length(model_paths) != length(subgroup_names)) {
    stop("model_paths, fit_paths, and subgroup_names must have the same length")
  }
  
  # Extract model info for all subgroups
  model_info <- list()
  for (i in 1:length(model_paths)) {
    model_info[[i]] <- extract_model_info(model_paths[i], fit_paths[i], subgroup_names[i])
  }
  
  # Combine all model info into a single data frame
  all_model_info <- do.call(rbind, model_info)
  
  # Save to CSV if output_file is provided
  if (!missing(output_file)) {
    write.csv(all_model_info, output_file, row.names = FALSE)
    cat("\nModel summary saved to:", output_file, "\n")
  }
  
  return(all_model_info)
}