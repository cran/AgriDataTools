#' Comprehensive Analysis of Variance (ANOVA) Engine for Randomized Complete Block Design (RCBD)
#'
#' @description
#' The \code{anova_rcbd} function executes a complete, high-precision linear model analysis for
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
#' @param data A verified \code{data.frame} containing genotype, replication block, and phenotypic response vectors.
#' @param trait A single character string specifying the exact column name of the numeric trait to analyze.
#' @param genotype_col An optional character string specifying the column name for genotypes. If \code{NULL}, automatic column detection is performed. Defaults to \code{NULL}.
#' @param rep_col An optional character string specifying the column name for replications/blocks. If \code{NULL}, automatic column detection is performed. Defaults to \code{NULL}.
#' @param reporting_level An integer flag defining console trace settings: \code{0} for silent execution, and \code{1} for printing formatted summary tables. Defaults to \code{1}.
#'
#' @return Invisibly returns a structured named \code{list} of class \code{"list"} containing 4 computational components:
#' \itemize{
#'   \item \strong{anova_table}: A \code{data.frame} acting as the standard ANOVA source matrix table containing Df, SS, MS, F_value, p_value, and significance flags.
#'   \item \strong{cv_percentage}: The computed Coefficient of Variation percentage scalar (CV%).
#'   \item \strong{mean_square_error}: The isolated Residual Error Mean Square (EMS), ready for genetic parameter engines.
#'   \item \strong{grand_mean}: The general mean arithmetic value of the evaluated trait.
#' }
#'
#' @seealso \code{\link{anova_crd}}
#'
#' @importFrom stats lm anova as.formula complete.cases
#' @export
#'
#' @examples
#' # Load your own dataset
#' data(gv_data, package = "AgriDataTools")
#' # Execute complete RCBD partition on target trait Plant Height(PH)
#' rcbd_results <- anova_rcbd(data = gv_data, trait = "PH")
#'
anova_rcbd <- function(data, trait, genotype_col = NULL, rep_col = NULL, reporting_level = 1) {
  
  if (missing(data) || missing(trait)) {
    stop("CRITICAL EXECUTION FAULT: Both 'data' and 'trait' arguments must be provided.", call. = FALSE)
  }
  
  if (!trait %in% colnames(data)) {
    stop(paste0("TRAIT REGISTRATION FAULT: Target column '", trait, "' was not found in dataset."), call. = FALSE)
  }
  
  # Dynamic Genotype Column Auto-detection
  if (is.null(genotype_col)) {
    possible_geno_cols <- c("Genotype", "genotype", "Genotypes", "genotypes", "Gen", "gen", "Cultivar", "cultivar", "Line", "line", "Treatment", "treatment")
    matched_col <- intersect(possible_geno_cols, colnames(data))
    if (length(matched_col) > 0) {
      genotype_col <- matched_col[1]
    } else {
      stop("GENOTYPE COLUMN FAULT: Could not automatically detect Genotype/Treatment column. Please specify 'genotype_col'.", call. = FALSE)
    }
  }
  
  # Dynamic Replication Column Auto-detection
  if (is.null(rep_col)) {
    possible_rep_cols <- c("Replication", "replication", "Replications", "replications", "Rep", "rep", "Block", "block", "Reps", "reps")
    matched_rep <- intersect(possible_rep_cols, colnames(data))
    if (length(matched_rep) > 0) {
      rep_col <- matched_rep[1]
    } else {
      stop("REPLICATION COLUMN FAULT: Could not automatically detect Replication/Block column. Please specify 'rep_col'.", call. = FALSE)
    }
  }
  
  if (!genotype_col %in% colnames(data)) {
    stop(paste0("GENOTYPE COLUMN FAULT: Specified genotype column '", genotype_col, "' not found in dataset."), call. = FALSE)
  }
  
  if (!rep_col %in% colnames(data)) {
    stop(paste0("REPLICATION COLUMN FAULT: Specified replication column '", rep_col, "' not found in dataset."), call. = FALSE)
  }
  
  # Handle NA Values Gracefully
  clean_data <- data[!is.na(data[[trait]]) & !is.na(data[[genotype_col]]) & !is.na(data[[rep_col]]), , drop = FALSE]
  if (nrow(clean_data) < nrow(data)) {
    warning(paste0("MISSING DATA NOTICE: Removed ", nrow(data) - nrow(clean_data), " rows containing NA values in target trait or factor columns."), call. = FALSE)
  }
  
  clean_data[[genotype_col]] <- as.factor(clean_data[[genotype_col]])
  clean_data[[rep_col]]      <- as.factor(clean_data[[rep_col]])
  response_vector            <- as.numeric(clean_data[[trait]])
  
  # Fit RCBD Linear Model
  formula_string <- paste0("`", trait, "` ~ `", rep_col, "` + `", genotype_col, "`")
  fitted_model   <- stats::lm(stats::as.formula(formula_string), data = clean_data)
  raw_anova      <- stats::anova(fitted_model)
  
  df_rep <- raw_anova[1, "Df"]
  df_gen <- raw_anova[2, "Df"]
  df_err <- raw_anova["Residuals", "Df"]
  df_tot <- df_rep + df_gen + df_err
  
  ss_rep <- raw_anova[1, "Sum Sq"]
  ss_gen <- raw_anova[2, "Sum Sq"]
  ss_err <- raw_anova["Residuals", "Sum Sq"]
  ss_tot <- ss_rep + ss_gen + ss_err
  
  ms_rep <- raw_anova[1, "Mean Sq"]
  ms_gen <- raw_anova[2, "Mean Sq"]
  ms_err <- raw_anova["Residuals", "Mean Sq"]
  
  f_rep <- raw_anova[1, "F value"]
  f_gen <- raw_anova[2, "F value"]
  
  p_rep <- raw_anova[1, "Pr(>F)"]
  p_gen <- raw_anova[2, "Pr(>F)"]
  
  # Significance flags
  assign_stars <- function(p_val) {
    if (is.na(p_val)) return("")
    if (p_val < 0.001) return("***")
    if (p_val < 0.01)  return("**")
    if (p_val < 0.05)  return("*")
    return("ns")
  }
  
  signif_rep <- assign_stars(p_rep)
  signif_gen <- assign_stars(p_gen)
  
  grand_mean_val <- mean(response_vector)
  cv_val         <- (sqrt(ms_err) / grand_mean_val) * 100
  
  compiled_anova <- data.frame(
    "Source"    = c("Replications", "Genotypes", "Error", "Total"),
    "Df"        = c(df_rep, df_gen, df_err, df_tot),
    "SS"        = c(ss_rep, ss_gen, ss_err, ss_tot),
    "MS"        = c(ms_rep, ms_gen, ms_err, NA),
    "F_value"   = c(f_rep, f_gen, NA, NA),
    "p_value"   = c(p_rep, p_gen, NA, NA),
    "Signif"    = c(signif_rep, signif_gen, "", ""),
    stringsAsFactors = FALSE
  )
  
  # Console Presentation
  if (reporting_level >= 1) {
    cat("\n", rep("=", 85), "\n", sep = "")
    cat(" ANALYSIS OF VARIANCE (ANOVA) FOR RCBD - TRAIT: ", trait, "\n")
    cat(rep("=", 85), "\n", sep = "")
    
    print_table <- compiled_anova
    print_table$SS <- ifelse(is.na(print_table$SS), "", sprintf("%.4f", print_table$SS))
    print_table$MS <- ifelse(is.na(print_table$MS), "", sprintf("%.4f", print_table$MS))
    print_table$F_value <- ifelse(is.na(print_table$F_value), "", sprintf("%.3f", print_table$F_value))
    print_table$p_value <- ifelse(is.na(print_table$p_value), "", format.pval(print_table$p_value, digits = 4, eps = 0.001))
    
    print(print_table, row.names = FALSE)
    
    cat(rep("-", 85), "\n", sep = "")
    cat(" Grand Mean                        : ", sprintf("%.4f", grand_mean_val), "\n")
    cat(" Coefficient of Variation (CV%)    : ", sprintf("%.2f", cv_val), "%\n")
    cat(rep("=", 85), "\n\n", sep = "")
  }
  
  return(invisible(list(
    anova_table         = compiled_anova,
    cv_percentage       = cv_val,
    mean_square_error   = ms_err,
    grand_mean          = grand_mean_val
  )))
}