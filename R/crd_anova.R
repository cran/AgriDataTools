#' Comprehensive Analysis of Variance (ANOVA) Engine for Completely Randomized Design (CRD)
#'
#' @description
#' The \code{anova_crd} function executes a complete, high-precision linear model analysis for
#' agricultural, laboratory, or greenhouse trials laid out under a Completely Randomized Design (CRD).
#' It computes partition sums of squares, hypothesis testing statistics, treatment variances,
#' significance flags, and the Coefficient of Variation (CV%).
#'
#' @details
#' In laboratory experiments, growth chamber studies, or field trials with completely homogeneous
#' environments, blocking is unnecessary. This function utilizes standard least-squares projection
#' to build the classic orthogonal CRD ANOVA matrix, modeling the response vector as a function
#' of treatment effects without blocking constraints:
#' \deqn{Y_{ij} = \mu + T_i + \varepsilon_{ij}}
#' Where \eqn{T_i} represents the treatment/genotype effect, and \eqn{\varepsilon_{ij}} is the
#' residual experimental error. The function handles both balanced and unbalanced data structures
#' perfectly, ensuring proper adjustments to degrees of freedom if replication numbers vary across lines.
#'
#' @param data A verified \code{data.frame} containing genotype/treatment identifiers and phenotypic responses.
#' @param trait A single character string specifying the exact column name of the numeric trait to analyze.
#' @param genotype_col An optional character string specifying the column name for genotypes. If \code{NULL}, automatic column detection is performed. Defaults to \code{NULL}.
#' @param reporting_level An integer flag defining console trace settings: \code{0} for silent execution and \code{1} for printing formatted ANOVA summary tables. Defaults to \code{1}.
#'
#' @return Invisibly returns a structured named \code{list} containing 4 computational components:
#' \itemize{
#'   \item \strong{anova_table}: A \code{data.frame} acting as the standard ANOVA source matrix table for CRD.
#'   \item \strong{cv_percentage}: A numeric scalar representing the computed Coefficient of Variation percentage.
#'   \item \strong{mean_square_error}: A numeric scalar representing the isolated Residual Error Mean Square.
#'   \item \strong{grand_mean}: A numeric scalar representing the overall arithmetic mean value.
#' }
#'
#' @seealso \code{\link{anova_rcbd}}
#' @importFrom stats anova as.formula lm complete.cases
#' @export
#'
#' @examples
#' # Load your own dataset
#' data(gv_data, package = "AgriDataTools")
#' # Execute complete CRD partition on target trait Plant Height(PH)
#' crd_results <- anova_crd(data = gv_data, trait = "PH", genotype_col = "Genotype")
#' 
anova_crd <- function(data, trait, genotype_col = NULL, reporting_level = 1) {
  
  if (missing(data) || missing(trait)) {
    stop("CRITICAL EXECUTION FAULT: Both 'data' and 'trait' arguments must be supplied.", call. = FALSE)
  }
  
  if (!trait %in% colnames(data)) {
    stop(paste0("TRAIT REGISTRATION FAULT: Target column '", trait, "' was not discovered in dataset."), call. = FALSE)
  }
  
  # Dynamic Genotype Column Detection
  if (is.null(genotype_col)) {
    possible_geno_cols <- c("Genotype", "genotype", "Genotypes", "genotypes", "Gen", "gen", "Variety", "variety", "Line", "line", "Treatment", "treatment")
    matched_col <- intersect(possible_geno_cols, colnames(data))
    if (length(matched_col) > 0) {
      genotype_col <- matched_col[1]
    } else {
      stop("GENOTYPE COLUMN FAULT: Could not automatically detect Genotype/Treatment column. Please specify 'genotype_col'.", call. = FALSE)
    }
  }
  
  if (!genotype_col %in% colnames(data)) {
    stop(paste0("GENOTYPE COLUMN FAULT: Specified genotype column '", genotype_col, "' not found in dataset."), call. = FALSE)
  }
  
  # Clean NA values safely
  clean_data <- data[!is.na(data[[trait]]) & !is.na(data[[genotype_col]]), , drop = FALSE]
  if (nrow(clean_data) < nrow(data)) {
    warning(paste0("MISSING DATA NOTICE: Removed ", nrow(data) - nrow(clean_data), " rows containing NA values in trait or genotype column."), call. = FALSE)
  }
  
  clean_data[[genotype_col]] <- as.factor(clean_data[[genotype_col]])
  response_vector <- as.numeric(clean_data[[trait]])
  
  # Fit Linear Model
  formula_string <- paste0("`", trait, "` ~ `", genotype_col, "`")
  fitted_model   <- stats::lm(stats::as.formula(formula_string), data = clean_data)
  raw_anova      <- stats::anova(fitted_model)
  
  df_gen <- raw_anova[1, "Df"]
  df_err <- raw_anova["Residuals", "Df"]
  df_tot <- df_gen + df_err
  
  ss_gen <- raw_anova[1, "Sum Sq"]
  ss_err <- raw_anova["Residuals", "Sum Sq"]
  ss_tot <- ss_gen + ss_err
  
  ms_gen <- raw_anova[1, "Mean Sq"]
  ms_err <- raw_anova["Residuals", "Mean Sq"]
  
  f_gen  <- raw_anova[1, "F value"]
  p_gen  <- raw_anova[1, "Pr(>F)"]
  
  # Significance stars flag assignment
  signif_flag <- ifelse(is.na(p_gen), "",
                        ifelse(p_gen < 0.001, "***",
                               ifelse(p_gen < 0.01,  "**",
                                      ifelse(p_gen < 0.05,  "*", "ns"))))
  
  grand_mean_val <- mean(response_vector)
  cv_val         <- (sqrt(ms_err) / grand_mean_val) * 100
  
  compiled_anova <- data.frame(
    "Source"    = c("Genotypes", "Error", "Total"),
    "Df"        = c(df_gen, df_err, df_tot),
    "SS"        = c(ss_gen, ss_err, ss_tot),
    "MS"        = c(ms_gen, ms_err, NA),
    "F_value"   = c(f_gen, NA, NA),
    "p_value"   = c(p_gen, NA, NA),
    "Signif"    = c(signif_flag, "", ""),
    stringsAsFactors = FALSE
  )
  
  # Console Output Pipeline
  if (reporting_level >= 1) {
    cat("\n", rep("=", 85), "\n", sep = "")
    cat(" ANALYSIS OF VARIANCE (ANOVA) FOR CRD - TRAIT: ", trait, "\n")
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