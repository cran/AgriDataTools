#' @importFrom utils write.table
NULL

#' Comprehensive Descriptive and Summary Statistics Engine for Phenotypic Traits
#'
#' @description
#' The \code{compute_summary_stats} function performs an exhaustive descriptive statistical sweep
#' across multiple numeric traits in an agricultural dataset. It calculates central tendency,
#' dispersion, and distribution shape metrics (skewness and kurtosis) for line screening.
#'
#' @details
#' Before executing hypothesis testing models like ANOVA, establishing dataset distribution profiles
#' is critical. This engine parses target numerical vectors to extract metrics: Standard Error of the Mean
#' is calculated as \eqn{SE = \frac{SD}{\sqrt{n}}}, Skewness measures distribution asymmetry, and
#' Kurtosis indicates tail weight relative to a normal curve. It dynamically filters out environmental factors
#' like 'Genotype' or 'Replication' and targets purely phenotypic observations.
#'
#' @param data A verified \code{data.frame} containing the experimental trial records.
#' @param traits A character vector specifying the exact column names to analyze. If \code{NULL},
#'   the system automatically discovers and evaluates all numeric columns. Defaults to \code{NULL}.
#' @param reporting_level An integer vector flag defining console trace settings: \code{0} for silent,
#'   \code{1} for descriptive summary grids. Defaults to \code{1}.
#'
#' @return A detailed structured \code{data.frame} where rows represent traits and columns contain calculated metrics.
#'
#' @seealso \code{\link{validate_agri_data}}
#'
#' @importFrom stats na.omit var sd
#' @export
#'
#' @examples
#' # Generate standard summary profiles across all phenotypic traits
#' descriptive_grid <- compute_summary_stats(data = gv_data)
#'
compute_summary_stats <- function(data, traits = NULL, reporting_level = 1) {
  
  # =========================================================================
  # BLOCK 1: STARTUP HEADER
  # =========================================================================
  if (reporting_level >= 1) {
    cat("\nDESCRIPTIVE SUMMARY STATISTICS\n\n")
  }
  
  # =========================================================================
  # BLOCK 2: SYSTEM DEFENSE AND AUTOMATIC DISCOVERY PIPELINE
  # =========================================================================
  if (missing(data) || !is.data.frame(data)) {
    stop("CRITICAL DATA FAULT: Input must be a valid data frame structure.", call. = FALSE)
  }
  
  if (is.null(traits)) {
    ignore_fields <- c("Genotype", "Genotypes", "Replication", "Replications", "Rep", "Reps", 
                       "Block", "Blocks", "Line", "Lines", "Cultivar", "Cultivars", 
                       "Variety", "Varieties", "Env", "Environment", "Loc", "Location", "Year")
    all_cols      <- colnames(data)
    numeric_cols  <- all_cols[sapply(data, is.numeric)]
    traits        <- setdiff(numeric_cols, ignore_fields)
    
    if (length(traits) == 0) {
      stop("DISCOVERY FAULT: Failed to capture numeric columns for phenotypic trait profiling.", call. = FALSE)
    }
  } else {
    missing_fields <- setdiff(traits, colnames(data))
    if (length(missing_fields) > 0) {
      stop(paste0("REGISTRATION FAULT: Columns not found in dataset: [",
                  paste(missing_fields, collapse = ", "), "]."), call. = FALSE)
    }
  }
  
  # =========================================================================
  # BLOCK 3: MATHEMATICAL MOMENTS COMPUTATION ENGINE
  # =========================================================================
  summary_records <- list()
  
  for (trait in traits) {
    vector_clean <- stats::na.omit(as.numeric(data[[trait]]))
    n_obs        <- length(vector_clean)
    
    if (n_obs < 3) {
      warning(paste0("DATA SHORTAGE WARNING: Trait '", trait, "' has insufficient observations. Skipping distribution shape moments."), call. = FALSE)
      next
    }
    
    mean_val <- mean(vector_clean)
    var_val  <- stats::var(vector_clean)
    sd_val   <- stats::sd(vector_clean)
    se_val   <- sd_val / sqrt(n_obs)
    min_val  <- min(vector_clean)
    max_val  <- max(vector_clean)
    
    # Adjusted Sample Moments (Unbiased Skewness and Kurtosis)
    deviations <- vector_clean - mean_val
    m2         <- sum(deviations^2)
    m3         <- sum(deviations^3)
    m4         <- sum(deviations^4)
    
    if (n_obs > 2 && var_val > 0) {
      skewness_val <- (n_obs * m3) / ((n_obs - 1) * (n_obs - 2) * (sd_val^3))
    } else {
      skewness_val <- 0
    }
    
    if (n_obs > 3 && var_val > 0) {
      kurtosis_val <- ((n_obs * (n_obs + 1) * m4) / ((n_obs - 1) * (n_obs - 2) * (n_obs - 3) * (var_val^2))) - 
        ((3 * ((n_obs - 1)^2)) / ((n_obs - 2) * (n_obs - 3)))
    } else {
      kurtosis_val <- 0
    }
    
    summary_records[[trait]] <- data.frame(
      "Trait"     = trait,
      "N"         = n_obs,
      "Mean"      = mean_val,
      "Variance"  = var_val,
      "Std_Dev"   = sd_val,
      "Std_Error" = se_val,
      "Min"       = min_val,
      "Max"       = max_val,
      "Skewness"  = skewness_val,
      "Kurtosis"  = kurtosis_val,
      stringsAsFactors = FALSE
    )
  }
  
  compiled_summary <- do.call(rbind, summary_records)
  rownames(compiled_summary) <- NULL
  
  # =========================================================================
  # BLOCK 4: FORMATTED CONSOLE TABLE WITH MINIMAL LINES
  # =========================================================================
  if (reporting_level >= 1) {
    cat("SUMMARY STATISTICS FOR ALL TRAITS:\n")
    cat(rep("-", 90), "\n", sep = "")
    print(compiled_summary, row.names = FALSE)
    cat(rep("-", 90), "\n\n", sep = "")
  }
  
  return(compiled_summary)
}