#' Analysis of Covariance (ANCOVA) Engine for Plant Breeding Trials
#'
#' @description
#' The `compute_ancova` function performs an exhaustive Analysis of Covariance
#' across agricultural experimental records. It adjusts the treatment (genotypic)
#' means for differences in a concomitant variable (covariate) to improve experimental
#' precision and error control.
#'
#' @details
#' In agricultural experiments, variation in a dependent phenotypic trait can sometimes
#' be partially controlled by measuring a secondary auxiliary property (covariate, e.g., initial
#' stand count or flowering duration). This engine fits a linear model incorporating both
#' the categorical genotypic design structure and the continuous covariate vector, extracting
#' the adjusted sums of squares, F-statistics, and adjusted marginal means.
#'
#' @param data A verified \code{data.frame} containing the experimental trial records.
#' @param response_trait A character string specifying the dependent phenotypic trait column.
#' @param covariate_trait A character string specifying the auxiliary covariate column.
#' @param reporting_level An integer vector flag defining console trace settings:
#'       \code{0} for silent execution, \code{1} for printing the variance partitioning table, and \code{2} for exhaustive diagnostic tracking logs. Defaults to \code{2}.
#'
#' @return Invisibly returns a structured named \code{list} of class \code{"list"} containing 5 computational components:
#' \item{response_trait}{A character string indicating the target phenotypic response variable evaluated.}
#' \item{covariate_trait}{A character string indicating the concomitant covariate variable utilized for model adjustment.}
#' \item{model}{The underlying fitted linear model object of class \code{\link[stats]{aov}}.}
#' \item{ancova_table}{A summary list structure of class \code{"summary.aov"} containing the Analysis of Covariance table with degrees of freedom, sums of squares, mean squares, F-values, and p-values.}
#' \item{adjusted_means}{A \code{data.frame} containing the covariate-adjusted genotypic/cultivar means (least-squares means) calculated at the mean value of the covariate.}
#' 
#' If \code{reporting_level >= 1}, the compiled ANCOVA table is printed directly to the console before returning the output object.
#'
#' @seealso \code{\link{validate_agri_data}}, \code{\link{compute_lsd}}
#'
#' @export
#'
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' 
#' # Perform ANCOVA using Plant Height (PH) as response and Spike Length (SL) as covariate
#' ancova_res <- compute_ancova(
#'   data = gv_data,
#'   response_trait = "PH",
#'   covariate_trait = "SL",
#'   reporting_level = 2
#' )
#' print(ancova_res$ancova_table)
#'
compute_ancova <- function(data, response_trait, covariate_trait, reporting_level = 2) {
  
  # =========================================================================
  # BLOCK 1: STARTUP PARAMETERS AND SYSTEM DEFENSE
  # =========================================================================
  timestamp_start <- Sys.time()
  
  if (reporting_level >= 1) {
    cat(rep("=", 85), "\n", sep = "")
    cat("AGRIDATATOOLS PACKAGED ENGINE v0.1.0 - ANALYSIS OF COVARIANCE (ANCOVA)\n")
    cat("Computation Inception: ", as.character(timestamp_start), "\n")
    cat(rep("-", 85), "\n", sep = "")
  }
  
  if (missing(data) || !is.data.frame(data)) {
    stop("CRITICAL DATA FAULT: Input must be a valid data frame structure.", call. = FALSE)
  }
  
  if (missing(response_trait) || missing(covariate_trait)) {
    stop("CRITICAL INPUT FAULT: Both response_trait and covariate_trait must be provided.", call. = FALSE)
  }
  
  required_cols <- c("Genotype", "Replication", response_trait, covariate_trait)
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste0("REGISTRATION FAULT: Columns not found in dataset: [",
                paste(missing_cols, collapse = ", "), "]."), call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 2: MODEL FITTING (ANCOVA LINEAR MODEL)
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - ANCOVA]: Fitting linear model with covariate adjustment...\n")
  }
  
  # Ensure Genotype and Replication are factors
  data$Genotype <- as.factor(data$Genotype)
  data$Replication <- as.factor(data$Replication)
  
  # Construct formula dynamically
  formula_str <- paste0(response_trait, " ~ Replication + Genotype + ", covariate_trait)
  ancova_model <- stats::aov(as.formula(formula_str), data = data)
  
  # Extract Type I or standard ANOVA table for ANCOVA
  ancova_summary <- summary(ancova_model)
  
  # =========================================================================
  # BLOCK 3: EXTRACTING ADJUSTED MEANS (LSMEANS)
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - ANCOVA]: Computing covariate-adjusted genotypic means...\n")
  }
  
  # Calculate adjusted means using standard linear model predictions at covariate mean
  covar_mean <- mean(data[[covariate_trait]], na.rm = TRUE)
  genotypes <- levels(data$Genotype)
  reps_first <- levels(data$Replication)[1]
  
  pred_data <- data.frame(
    Genotype = factor(genotypes, levels = genotypes),
    Replication = factor(rep(reps_first, length(genotypes)), levels = levels(data$Replication))
  )
  pred_data[[covariate_trait]] <- covar_mean
  
  pred_data$Adjusted_Mean <- stats::predict(ancova_model, newdata = pred_data)
  
  adjusted_means_df <- pred_data[, c("Genotype", "Adjusted_Mean")]
  rownames(adjusted_means_df) <- NULL
  
  # =========================================================================
  # BLOCK 4: CONSOLE PRESENTATION PIPELINE
  # =========================================================================
  if (reporting_level >= 1) {
    cat("\n", rep("-", 75), "\n", sep = "")
    cat(" COMPILED ANCOVA VARIANCE PARTITIONING TABLE\n")
    cat(rep("-", 75), "\n", sep = "")
    print(ancova_summary)
    cat(rep("=", 85), "\n\n", sep = "")
  }
  
  # =========================================================================
  # BLOCK 5: CLOSURE AND RETURN STRUCTURE
  # =========================================================================
  timestamp_end <- Sys.time()
  if (reporting_level >= 2) {
    cat("[LOG - FINALIZE]: ANCOVA execution completed in ", round(as.numeric(difftime(timestamp_end, timestamp_start, units = "secs")), 5), " seconds.\n")
  }
  
  result_payload <- list(
    response_trait = response_trait,
    covariate_trait = covariate_trait,
    model = ancova_model,
    ancova_table = ancova_summary,
    adjusted_means = adjusted_means_df
  )
  
  return(invisible(result_payload))
}

utils::globalVariables(c("Genotype", "Replication"))