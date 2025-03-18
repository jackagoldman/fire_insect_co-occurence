#' Adjust MatchIt Parameters
#'
#' This function adjusts the parameters of the MatchIt model based on the provided inputs.
#'
#' @param data A data frame containing the dataset.
#' @param method A character string specifying the matching method. Default is "nearest".
#' @param distance A character string specifying the distance measure. Default is "glm".
#' @param link A character string specifying the link function. Default is "logit".
#' @param m_order A character string specifying the matching order. Default is "random".
#' @param caliper A numeric value specifying the maximum distance for matches. Default is 0.2.
#' @param replace A logical value indicating whether matches should be with replacement. Default is TRUE.
#' @param mahvars A formula specifying the covariates for Mahalanobis distance matching. Default is NULL.
#' @return A MatchIt object containing the matching results.
#' @examples
#' adjusted_model <- adjust_matchit_params(hist_gt90_1, method = "nearest", caliper = 0.25)


adjust_matchit_params <- function(data, method = "nearest", distance = "glm", link = "logit", m_order = "random", caliper = 0.2, replace = TRUE, use_mahvars = FALSE) {
  mahvars <- if (use_mahvars && distance != "mahalanobis") {
    ~host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90
  } else {
    NULL
  }
  
  if (distance == "mahalanobis" && !is.null(caliper)) {
    caliper <- setNames(rep(caliper, length(all.vars(mahvars))), all.vars(mahvars))
  }
  
  model <- matchit(history ~ host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90,
                   data = data,
                   method = method,
                   distance = distance,
                   link = link,
                   m.order = m_order,
                   caliper = caliper,
                   replace = replace,
                   mahvars = mahvars)
  return(model)
}



#' Find Best MatchIt Model
#'
#' This function iterates over all combinations of the given parameters, stores the results in a list, and selects the best model based on the standard mean difference (SMD) being below 0.25.
#'
#' @param data A data frame containing the dataset.
#' @param methods A character vector specifying the matching methods to try.
#' @param distances A character vector specifying the distance measures to try.
#' @param links A character vector specifying the link functions to try.
#' @param m_orders A character vector specifying the matching orders to try.
#' @param calipers A numeric vector specifying the maximum distances for matches to try.
#' @param replace A logical value indicating whether matches should be with replacement.
#' @param mahvars_list A list of formulas specifying the covariates for Mahalanobis distance matching to try.
#' @return A list containing the best MatchIt model and all results.
#' @examples
#' methods <- c("nearest", "optimal")
#' distances <- c("glm", "mahalanobis")
#' links <- c("logit", "probit")
#' m_orders <- c("random", "largest")
#' calipers <- c(0.1, 0.2, 0.3)
#' replace <- TRUE
#' mahvars_list <- list(NULL, ~host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90)
#' best_model_info <- find_best_model(hist_gt90_1, methods, distances, links, m_orders, calipers, replace, mahvars_list)
#' best_model <- best_model_info$best_model
#' all_results <- best_model_info$results
#' 
#' 
#' 
find_best_model <- function(data, methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list) {
  results <- list()
  best_model <- NULL
  best_smd <- Inf
  
  for (method in methods) {
    for (distance in distances) {
      for (link in links) {
        for (m_order in m_orders) {
          for (caliper in calipers) {
            for (replace in replace_list) {
              for (use_mahvars in use_mahvars_list) {
                model <- adjust_matchit_params(data, method, distance, link, m_order, caliper, replace, use_mahvars)
                model_summary <- summary(model)
                
                if (is.list(model_summary) && "sum.matched" %in% names(model_summary)) {
                  smd <- as.data.frame(model_summary$sum.matched)$`Std. Mean Diff.`
                  
                  if (all(smd < 0.25)) {
                    results[[paste(method, distance, link, m_order, caliper, replace, use_mahvars, sep = "_")]] <- model
                    
                    if (mean(smd) < best_smd) {
                      best_smd <- mean(smd)
                      best_model <- model
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(list(best_model = best_model, results = results))
}


#' Extract Parameters from Best Models
#'
#' This function extracts parameters from up to three best models and combines them into a single data frame. 
#' It also includes a subgrouping string and a time series descriptor (TSD) extracted from the model names.
#'
#' @param best_model_1 A list representing the first best model.
#' @param best_model_2 A list representing the second best model.
#' @param best_model_3 A list representing the third best model.
#' @param subgrouping A string representing the subgrouping information.
#' @return A data frame containing the extracted parameters from the best models, including the subgrouping and TSD.
#' @examples
#' best_model_1 <- list(info = list(method = "nearest", mahalanobis = "TRUE", distance = "glm", link = "logit"), caliper = 0.2)
#' best_model_2 <- list(info = list(method = "optimal", mahalanobis = "FALSE", distance = "mahalanobis", link = "probit"), caliper = 0.3)
#' best_model_3 <- list(info = list(method = "full", mahalanobis = "TRUE", distance = "glm", link = "logit"), caliper = 0.1)
#' subgrouping <- "Group A"
#' extract_best_model_params(best_model_1, best_model_2, best_model_3)
extract_best_model_params <- function(best_model_1, best_model_2, best_model_3) {
  
  # Helper function to extract the suffix and remove it from the name
  extract_suffix <- function(model_name) {
    parts <- strsplit(model_name, "_")[[1]]
    suffix <- tail(parts, 1)
    base_name <- paste(head(parts, -1), collapse = "_")
    return(list(base_name = base_name, suffix = suffix))
  }
  
  # Extract suffixes and update best_model arguments
  best_model_1_info <- extract_suffix(deparse(substitute(best_model_1)))
  best_model_2_info <- extract_suffix(deparse(substitute(best_model_2)))
  best_model_3_info <- extract_suffix(deparse(substitute(best_model_3)))
  
  best_model_1$subgroup <- best_model_1_info$suffix
  best_model_2$subgroup <- best_model_2_info$suffix
  best_model_3$subgroup <- best_model_3_info$suffix
  
  # Combine the three best models into a list
  best_models <- list(best_model_1, best_model_2, best_model_3)

  # Function to extract parameters from a single model
  extract_params <- function(model) {
    if (!is.list(model) || !is.list(model$info)) {
      stop("Each best_model must be a list with an 'info' element that is also a list.")
    }
    
    params <- list(
      Subgroup = if(model$subgroup == 1) "0-2 years after defoliation" else if(model$subgroup == 2) "3-9 years after defoliation" else "10+ years after defoliation",
      Method = if (model$info$method == "nearest") "1:1 nearest neighbor matching without replacement" else model$info$method,
      Matching = if (model$info$mahalanobis == "TRUE" & model$info$distance == "mahalanobis") "Mahalanobis" else if (model$info$mahalanobis == "TRUE" & model$info$distance == "glm") "Mahalanobis distance for matching and a propensity score caliper to ensure close matches" else "Propensity Score",
      Distance = if (model$info$distance == "glm") "estimated with logistic regression (glm)" else if (model$info$distance == "mahalanobis") "estimated with mahalanobis distances" else model$info$distance,
      Link = model$info$link,
      Caliper = model$caliper
    )
    
    # Replace NULL values with NA
    params <- lapply(params, function(x) {
      if (is.null(x)) {
        x <- NA
      }
      return(x)
    })
    
    # Ensure all elements have the same length
    max_length <- max(sapply(params, length))
    params <- lapply(params, function(x) {
      if (length(x) < max_length) {
        x <- rep(x, length.out = max_length)
      }
      return(x)
    })
    
    params_df <- as.data.frame(params, stringsAsFactors = FALSE)
    return(params_df)
  }
  
  # Extract parameters for each model and combine into a single data frame
  params_list <- lapply(best_models, extract_params)
  combined_params_df <- do.call(rbind, params_list)
  
  return(combined_params_df)
}

#' Calculate R-squared for GLS Model
#'
#' This function calculates the R-squared value for a Generalized Least Squares (GLS) model.
#'
#' @param model A GLS model object of class `gls`.
#' @param formula A formula specifying the model.
#' @return The R-squared value for the GLS model.
#' @examples
#' library(nlme)
#' gls_model <- gls(y ~ x, data = your_data, correlation = corExp(form = ~ x))
#' rsquared.gls(gls_model, y ~ x)
#' @export
rsquared.gls <- function(model, formula) {
  # Extract the data used to fit the model
  data <- getData(model)
  
  # Extract the response variable
  response <- model$terms[[2]]
  y <- data[[as.character(response)]]
  
  # Create the model matrix
  X <- model.matrix(formula, data = data)
  
  # Calculate the variance of the fitted values
  sigmaF <- var(as.vector(model$coefficients %*% t(X)))
  
  # Calculate the variance of the residuals
  sigmaE <- var(resid(model))
  
  # Calculate R-squared
  R_squared <- sigmaF / (sigmaF + sigmaE)
  
  print(R_squared)
}


#' Extract P-Values from GLS Model Summary and Convert to Kable
#'
#' This function extracts the p-value table from a GLS model summary and converts it to a `kable` object for easy visualization.
#'
#' @param model_gls A GLS model object of class `gls`.
#' @return A `kable` object containing the p-values from the GLS model summary.
#' @examples
#' \dontrun{
#' library(nlme)
#' library(knitr)
#' library(dplyr)
#'
#' # Fit a sample GLS model
#' model_gls <- gls(distance ~ age + Sex, data = Orthodont, correlation = corGaus(form = ~ age))
#'
#' # Get the p-value table as a kable object
#' gls_pvalues_to_kable(model_gls)
#' }
#' @export
gls_pvalues_to_kable <- function(model_gls) {
  require(kableExtra)
  # Get the summary of the GLS model
  summary_gls <- summary(model_gls)
  
  # Extract the p-value table
  p_values <-summary_gls$tTable %>%
    as.data.frame() %>% 
    select(c(Value, Std.Error ,`t-value` ,`p-value`)) %>% 
    rename(p = `p-value`) %>% 
    rownames_to_column(., var = "Covariates") %>% 
    rename("se ± " = "Std.Error") %>%
    rename(estimate = Value) %>% 
    mutate(Covariates = recode(Covariates, "host_pct" = "Host %" , 
                               "isi_90"= "Initial Spread Index",
                               "dc_90"="Drought Code" ,
                               "dmc_90" = "Duff Moisture Code",
                               "ffmc_90" = "Fine Fuel Moisture Code",
                               "bui_90" = "Build Up Index",
                               "fwi_90" = "Fire Weather Index",
                               "x" = "Longitude",
                               "y" = "Latitude",
                               "Cumulative_Years_Defol:window_opp1" = "Years Defoliated × 0-2 Years Since Defoliation",
                               "Cumulative_Years_Defol:window_opp2" = "Years Defoliated × 3-10 Years Since Defoliation",
                               "Cumulative_Years_Defol:window_opp3" = "Years Defoliated × 10-15 Years Since Defoliation")) %>% 
    mutate(`s` = sapply(p, add_significance_stars)) %>% 
    mutate(across(where(is.numeric), round, 2)) %>% 
    mutate(p = if_else(`p` < 0.0001, "< 0.0001", as.character(`p`))) %>% 
    rename(" " = s)
  
  # Convert the p-value table to a kable object
  kable(p_values)
}



