#' Genetic Variability and Quantitative Inheritance Parameters Estimation Engine
#'
#' @description
#' The `estimate_variability` function calculates comprehensive biometric genetic profiles
#' from replication-based agricultural trial datasets. It partitions phenotypic variance into
#' Genotypic Variance (Vg), Phenotypic Variance (Vp), and Environmental Variance (Ve), and computes
#' critical breeding metrics including Genotypic Coefficient of Variation (GCV), Phenotypic
#' Coefficient of Variation (PCV), Broad-Sense Heritability (H2), Genetic Advance (GA), and Genetic Advance as % of Mean (GAM).
#'
#' @details
#' In quantitative genetics, phenotypic variance must be dissected into its components to
#' determine the role of genetic factors versus environmental noise. This function extracts
#' the Error Mean Square (EMS) and Genotypic Mean Square (GMS) directly from completed ANOVA matrices:
#' \deqn{V_g = \frac{GMS - EMS}{r}}
#' \deqn{V_p = V_g + EMS}
#' \deqn{H^2 = \frac{V_g}{V_p}}
#' Where \eqn{r} represents the total absolute replication or block count.
#'
#' If high environmental variations cause the computed Genotypic Variance to become negative,
#' the system automatically applies a mathematical lower boundary floor at \code{0.0001} to preserve
#' downstream pipeline integrity and issues a detailed structural warning message.
#'
#' @param anova_results A structured \code{list} returned by either the \code{anova_rcbd}
#'   or \code{anova_crd} analysis pipelines within this package.
#' @param total_replications An integer specifying the total number of replications/blocks used
#'   in the experimental trial layout.
#' @param reporting_level An integer flag defining console output details: \code{0} for silent,
#'   \code{1} for summary parameter tables, and \code{2} for comprehensive descriptive logs. Defaults to \code{2}.
#'
#' @return A structured named \code{list} containing calculated genetic variability components:
#' \item{genotypic_variance}{Estimated genotypic variance component (\eqn{V_g}).}
#' \item{phenotypic_variance}{Total phenotypic variance component (\eqn{V_p}).}
#' \item{environmental_variance}{Environmental variance/Error Mean Square (\eqn{V_e}).}
#' \item{gcv}{Genotypic Coefficient of Variation percentage (GCV\%).}
#' \item{pcv}{Phenotypic Coefficient of Variation percentage (PCV\%).}
#' \item{heritability_percentage}{Broad-Sense Heritability percentage (\eqn{H^2 \%}).}
#' \item{genetic_advance}{Expected Genetic Advance (GA) at 5\% selection intensity (\eqn{k = 2.06}).}
#' \item{gam_percentage}{Genetic Advance as a percentage of the Grand Mean (GAM\%).}
#'
#' @seealso \code{\link{anova_rcbd}}, \code{\link{anova_crd}}
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'    # Execute complete genetic variability partitioning
#'    rcbd_out <- anova_rcbd(data = gv_data, trait = "PH", reporting_level = 0)
#'    var_metrics <- estimate_variability(anova_results = rcbd_out, total_replications = 3)
#'    print(var_metrics$heritability_percentage)
#' }
estimate_variability <- function(anova_results, total_replications, reporting_level = 2) {
  
  # =========================================================================
  # BLOCK 1: ENVIRONMENTAL REGISTRATION AND DIAGNOSTIC TRACE SETUP
  # =========================================================================
  timestamp_start <- Sys.time()
  
  if (reporting_level >= 1) {
    cat(rep("=", 85), "\n", sep = "")
    cat("AGRIDATATOOLS PACKAGED ENGINE v0.1.0 - VARIABILITY ANALYSIS ENGINE\n")
    cat("Parameter Sweep Inception: ", as.character(timestamp_start), "\n")
    cat(rep("-", 85), "\n", sep = "")
  }
  
  # =========================================================================
  # BLOCK 2: SYSTEM DEFENSE AND STRUCTURAL GATEKEEPER
  # =========================================================================
  if (missing(anova_results)) {
    stop("CRITICAL METRIC FAULT: Argument 'anova_results' is missing. Cannot estimate variability without an ANOVA model.", call. = FALSE)
  }
  
  if (missing(total_replications)) {
    stop("CRITICAL CONFIGURATION FAULT: Argument 'total_replications' is missing. Replication count is required.", call. = FALSE)
  }
  
  # Verify structural integrity of input list
  required_elements <- c("anova_table", "mean_square_error", "grand_mean")
  missing_elements  <- setdiff(required_elements, names(anova_results))
  
  if (length(missing_elements) > 0) {
    stop(paste0("STRUCTURAL DEVIATION FAULT: Input 'anova_results' is missing required nodes: [",
                paste(missing_elements, collapse = ", "), "]. Process aborted."), call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 3: VARIANCE EXTRACTION AND PARTITIONING
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - VARIABILITY]: Mapping data variables out of the ANOVA payload structures...\n")
  }
  
  anova_matrix   <- anova_results$anova_table
  mean_sq_error  <- as.numeric(anova_results$mean_square_error)
  grand_mean_val <- as.numeric(anova_results$grand_mean)
  
  # Dynamically capture Genotype row index variations safely (Genotype / Genotypes)
  genotype_row_idx <- grep("Genotype", anova_matrix$Source)
  if (length(genotype_row_idx) == 0) {
    stop("METADATA SYNTAX CRASH: Could not discover 'Genotype' or 'Genotypes' row in the provided variance analysis matrix.", call. = FALSE)
  }
  
  mean_sq_genotype <- as.numeric(anova_matrix$MS[genotype_row_idx[1]])
  
  # Core mathematical calculations
  environmental_variance_val <- mean_sq_error
  genotypic_variance_val     <- (mean_sq_genotype - mean_sq_error) / total_replications
  
  # Defensive floor handling for negative variance anomalies
  if (genotypic_variance_val < 0) {
    warning(paste0("QUANTITATIVE BIOMETRIC WARNING: Calculated Genotypic Variance is negative (",
                   round(genotypic_variance_val, 5), "). Bounding to system floor boundary of 0.0001."), call. = FALSE)
    genotypic_variance_val <- 0.0001
  }
  
  phenotypic_variance_val <- genotypic_variance_val + environmental_variance_val
  
  # =========================================================================
  # BLOCK 4: BREEDING COEFFICIENTS AND GENETIC ADVANCE ESTIMATION
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - VARIABILITY]: Computing GCV, PCV, GA, and GAM metrics...\n")
  }
  
  # Basic Coefficients
  gcv_val <- (sqrt(genotypic_variance_val) / grand_mean_val) * 100
  pcv_val <- (sqrt(phenotypic_variance_val) / grand_mean_val) * 100
  h2_perc <- (genotypic_variance_val / phenotypic_variance_val) * 100
  
  # Genetic Advance (GA) and GAM Calculation (k = 2.06 for 5% selection intensity)
  k <- 2.06
  ga <- k * sqrt(phenotypic_variance_val) * (genotypic_variance_val / phenotypic_variance_val)
  gam_val <- (ga / grand_mean_val) * 100
  
  # =========================================================================
  # BLOCK 5: CONSOLE PRESENTATION PIPELINE
  # =========================================================================
  if (reporting_level >= 1) {
    cat("\n", rep("-", 75), "\n", sep = "")
    cat(" QUANTITATIVE BIOMETRIC GENETIC PARAMETERS AND VARIABILITY SUMMARY\n")
    cat(rep("-", 75), "\n", sep = "")
    cat("  Genotypic Variance (Vg)             : ", round(genotypic_variance_val, 5), "\n")
    cat("  Phenotypic Variance (Vp)            : ", round(phenotypic_variance_val, 5), "\n")
    cat("  Environmental Variance (Ve)         : ", round(environmental_variance_val, 5), "\n")
    cat("  Genotypic Coeff of Variation (GCV%) : ", round(gcv_val, 2), "%\n")
    cat("  Phenotypic Coeff of Variation (PCV%): ", round(pcv_val, 2), "%\n")
    cat("  Broad-Sense Heritability (H2 %)     : ", round(h2_perc, 2), "%\n")
    cat("  Genetic Advance (GA)                : ", round(ga, 4), "\n")
    cat("  Genetic Advance as % of Mean (GAM)  : ", round(gam_val, 2), "%\n")
    cat(rep("=", 85), "\n\n", sep = "")
  }
  
  # =========================================================================
  # BLOCK 6: SYSTEM PAYMENT AND EXIT PAYLOAD
  # =========================================================================
  timestamp_end      <- Sys.time()
  processing_seconds <- as.numeric(difftime(timestamp_end, timestamp_start, units = "secs"))
  
  if (reporting_level >= 2) {
    cat("[LOG - FINALIZE]: Variability evaluation completed in ", round(processing_seconds, 5), " seconds.\n")
  }
  
  return(list(
    genotypic_variance      = genotypic_variance_val,
    phenotypic_variance     = phenotypic_variance_val,
    environmental_variance  = environmental_variance_val,
    gcv                     = gcv_val,
    pcv                     = pcv_val,
    heritability_percentage = h2_perc,
    genetic_advance         = ga,
    gam_percentage          = gam_val
  ))
}