#' Comprehensive Analysis of Variance (ANOVA) Engine for Randomized Complete Block Design (RCBD)
#'
#' @description
#' The `anova_rcbd` function executes a complete, high-precision linear model analysis for
#' agricultural trials laid out under an RCBD framework. It computes partition sums of squares,
#' hypothesis testing statistics, significance flags, and the Coefficient of Variation (CV%).
#'
#' @details
#' In plant breeding and agronomy trials, isolating block variance from the true experimental
#' error is vital to properly evaluate lines, cultivars, or treatments. This function uses
#' standard least-squares projection to build the classic orthogonal ANOVA matrix:
#' \deqn{Y_{ij} = \mu + G_i + R_j + e_{ij}}
#' Where \eqn{G_i} represents the genotype effect, \eqn{R_j} is the replication block effect,
#' and \eqn{e_{ij}} is the residual experimental error.
#'
#' @param data A verified \code{data.frame} containing the columns \code{Genotype}, \code{Replication},
#'    and the target phenotypic trait response.
#' @param trait A single character string specifying the exact column name of the numeric trait to analyze.
#' @param reporting_level An integer vector flag defining console trace settings: \code{0} for silent,
#'    \code{1} for basic summary tables, and \code{2} for comprehensive descriptive metrics. Defaults to \code{2}.
#'
#' @return A structured named \code{list} containing 4 computational components:
#' \item{anova_table}{A \code{data.frame} acting as the standard ANOVA source matrix table containing Df, SS, MS, F_value, and p_value.}
#' \item{cv_percentage}{The computed Coefficient of Variation percentage scalar (\eqn{CV\%}).}
#' \item{mean_square_error}{The isolated Residual Error Mean Square (EMS), ready for genetic parameter engines.}
#' \item{grand_mean}{The general mean arithmetic value of the evaluated trait.}
#'
#' @seealso \code{\link{compute_lsd}}
#'
#' @importFrom stats lm anova as.formula
#' @export
#'
#' @examples
#'    # Execute complete RCBD partition on gv_data asset for Plant Height (PH)
#'    rcbd_results <- anova_rcbd(data = gv_data, trait = "PH")
#' 
anova_rcbd <- function(data, trait, reporting_level = 2) {
  
  # -------------------------------------------------------------------------
  # LAYER 1: DISCOVERY AND INPUT GUARDRAILS
  # -------------------------------------------------------------------------
  if (missing(data) || missing(trait)) {
    stop("ANALYSIS ERROR: Both 'data' and 'trait' arguments must be provided.", call. = FALSE)
  }
  
  if (!trait %in% colnames(data)) {
    stop(paste0("TRAIT REGISTRATION FAULT: Target column '", trait, "' was not found."), call. = FALSE)
  }
  
  data$Genotype    <- as.factor(data$Genotype)
  data$Replication <- as.factor(data$Replication)
  response_vector  <- as.numeric(data[[trait]])
  
  if (any(is.na(response_vector))) {
    stop("DATA TRUNCATION FAULT: Missing values (NA) detected. Please patch records.", call. = FALSE)
  }
  
  # -------------------------------------------------------------------------
  # LAYER 2: LINEAR MODEL MATRIX ALGEBRA
  # -------------------------------------------------------------------------
  formula_string <- paste0("`", trait, "` ~ Replication + Genotype")
  fitted_model   <- stats::lm(stats::as.formula(formula_string), data = data)
  raw_anova      <- stats::anova(fitted_model)
  
  df_rep <- raw_anova["Replication", "Df"]
  df_gen <- raw_anova["Genotype", "Df"]
  df_err <- raw_anova["Residuals", "Df"]
  df_tot <- df_rep + df_gen + df_err
  
  ss_rep <- raw_anova["Replication", "Sum Sq"]
  ss_gen <- raw_anova["Genotype", "Sum Sq"]
  ss_err <- raw_anova["Residuals", "Sum Sq"]
  ss_tot <- ss_rep + ss_gen + ss_err
  
  ms_rep <- raw_anova["Replication", "Mean Sq"]
  ms_gen <- raw_anova["Genotype", "Mean Sq"]
  ms_err <- raw_anova["Residuals", "Mean Sq"]
  
  f_rep <- raw_anova["Replication", "F value"]
  f_gen <- raw_anova["Genotype", "F value"]
  p_rep <- raw_anova["Replication", "Pr(>F)"]
  p_gen <- raw_anova["Genotype", "Pr(>F)"]
  
  grand_mean_val <- mean(response_vector)
  cv_val         <- (sqrt(ms_err) / grand_mean_val) * 100
  
  compiled_anova <- data.frame(
    "Source"  = c("Replications (Blocks)", "Genotypes (Lines)", "Error (Residual)", "Total"),
    "Df"      = c(df_rep, df_gen, df_err, df_tot),
    "SS"      = c(ss_rep, ss_gen, ss_err, ss_tot),
    "MS"      = c(ms_rep, ms_gen, ms_err, NA),
    "F_value" = c(f_rep, f_gen, NA, NA),
    "p_value" = c(p_rep, p_gen, NA, NA),
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------------------
  # LAYER 3: CONSOLE PRESENTATION PIPELINE
  # -------------------------------------------------------------------------
  if (reporting_level >= 1) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat(" ANALYSIS OF VARIANCE (ANOVA) FOR RCBD - TRAIT: ", trait, "\n")
    cat(rep("-", 70), "\n", sep = "")
    
    print_table <- compiled_anova
    print_table$SS <- round(print_table$SS, 4)
    print_table$MS <- round(print_table$MS, 4)
    print_table$F_value <- round(print_table$F_value, 3)
    print_table$p_value <- format.pval(print_table$p_value, digits = 4, eps = 0.001)
    
    print(print_table, row.names = FALSE)
    
    cat(rep("-", 70), "\n", sep = "")
    cat(" Trial Grand Mean         : ", round(grand_mean_val, 4), "\n")
    cat(" Coefficient of Var (CV%): ", round(cv_val, 2), "%\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }
  
  return(invisible(list(
    anova_table = compiled_anova,
    cv_percentage = cv_val,
    mean_square_error = ms_err,
    grand_mean = grand_mean_val
  )))
}