#' Compact Analysis of Covariance (ANCOVA) Table Engine for Plant Breeding
#'
#' @description
#' The \code{compute_ancova} function evaluates an Analysis of Covariance (ANCOVA) 
#' and compiles a streamlined, standard biometrical table featuring Degrees of Freedom (Df), 
#' Sum of Products (SP_XY), Mean Products (MP_XY), exact Biometrical F-values, p-values, significance flags, 
#' and rigorous Covariance components (Cov_e, Cov_g, Cov_p).
#'
#' @details
#' The mathematical partitioning under RCBD framework is computed using joint reduction methods:
#' \deqn{SP_{Total} = SP_{Replications} + SP_{Genotypes} + SP_{Residual}}
#' \deqn{MP_{XY} = \frac{SP_{XY}}{Df}}
#' \deqn{Cov_e = MP_{Error}}
#' \deqn{Cov_g = \frac{MP_{Genotypes} - MP_{Error}}{r}}
#' \deqn{Cov_p = Cov_g + \frac{Cov_e}{r}}
#'
#' @param data A verified \code{data.frame} containing experimental trial records.
#' @param response_trait A character string specifying the dependent phenotypic response variable (e.g., "GYPM").
#' @param covariate_trait A character string specifying the auxiliary covariate variable (e.g., "PH").
#' @param genotype_col Optional character string specifying the genotype column. Defaults to \code{NULL} for auto-detection.
#' @param rep_col Optional character string specifying the replication column. Defaults to \code{NULL} for auto-detection.
#' @param reporting_level An integer flag: \code{0} for silent execution, \code{1} for printing the ANCOVA table. Defaults to \code{1}.
#'
#' @return Invisibly returns a structured list containing:
#' \item{ancova_table}{A compact \code{data.frame} containing the exact ANCOVA source table.}
#' \item{covariance_components}{A named numeric vector containing isolated \code{Cov_e}, \code{Cov_g}, and \code{Cov_p} values.}
#' \item{adjusted_means}{A \code{data.frame} containing covariate-adjusted genotypic least-squares means.}
#'
#' @importFrom stats lm anova as.formula complete.cases predict aov pf
#' @export
#'
#' @examples
#' library(AgriDataTools)
#' data("gv_data", package = "AgriDataTools")
#' 
#' # Run Analysis of Covariance between Grain Yield per Meter and Plant Height
#' ancova_results <- compute_ancova(
#'     data = gv_data,
#'     response_trait = "GYPM",
#'     covariate_trait = "PH"
#' )
#'
compute_ancova <- function(data, response_trait, covariate_trait, genotype_col = NULL, rep_col = NULL, reporting_level = 1) {
  
  if (missing(data) || !is.data.frame(data)) {
    stop("CRITICAL DATA FAULT: Input must be a valid data frame structure.", call. = FALSE)
  }
  
  if (missing(response_trait) || missing(covariate_trait)) {
    stop("CRITICAL INPUT FAULT: Both 'response_trait' and 'covariate_trait' must be provided.", call. = FALSE)
  }
  
  if (!response_trait %in% colnames(data) || !covariate_trait %in% colnames(data)) {
    stop(paste0("REGISTRATION FAULT: Traits not found in dataset: [", response_trait, ", ", covariate_trait, "]"), call. = FALSE)
  }
  
  # Auto-detect Genotype Column
  if (is.null(genotype_col)) {
    possible_geno <- c("Genotype", "genotype", "Genotypes", "Gen", "Variety", "Cultivar", "Line")
    matched_geno <- intersect(possible_geno, colnames(data))
    if (length(matched_geno) > 0) genotype_col <- matched_geno[1]
    else stop("GENOTYPE COLUMN FAULT: Could not auto-detect genotype column. Specify 'genotype_col'.", call. = FALSE)
  }
  
  # Auto-detect Replication Column
  if (is.null(rep_col)) {
    possible_rep <- c("Replication", "replication", "Replications", "Rep", "rep", "Block", "replications", "reps", "Reps", "block", "blocks", "Blocks")
    matched_rep <- intersect(possible_rep, colnames(data))
    if (length(matched_rep) > 0) rep_col <- matched_rep[1]
    else stop("REPLICATION COLUMN FAULT: Could not auto-detect replication column. Specify 'rep_col'.", call. = FALSE)
  }
  
  # Clean NA data across columns
  clean_cols <- c(genotype_col, rep_col, response_trait, covariate_trait)
  clean_data <- data[stats::complete.cases(data[, clean_cols]), , drop = FALSE]
  
  if (nrow(clean_data) < nrow(data)) {
    warning(paste0("MISSING DATA NOTICE: Removed ", nrow(data) - nrow(clean_data), " rows containing NA values."), call. = FALSE)
  }
  
  clean_data[[genotype_col]] <- as.factor(clean_data[[genotype_col]])
  clean_data[[rep_col]]      <- as.factor(clean_data[[rep_col]])
  
  r <- nlevels(clean_data[[rep_col]])
  g <- nlevels(clean_data[[genotype_col]])
  
  # 1. Fit Univariate Models for Y and X to partition Sum of Squares
  f_y <- stats::as.formula(paste0("`", response_trait, "` ~ `", rep_col, "` + `", genotype_col, "`"))
  m_y <- stats::lm(f_y, data = clean_data)
  a_y <- stats::anova(m_y)
  
  # 2. Joint Sum of Products Partitioning (SP_XY) using standard biometrical identity
  clean_data$sum_xy <- clean_data[[response_trait]] + clean_data[[covariate_trait]]
  
  m_x <- stats::lm(stats::as.formula(paste0("`", covariate_trait, "` ~ `", rep_col, "` + `", genotype_col, "`")), data = clean_data)
  m_sum <- stats::lm(stats::as.formula(paste0("sum_xy ~ `", rep_col, "` + `", genotype_col, "`")), data = clean_data)
  
  a_x <- stats::anova(m_x)
  a_sum <- stats::anova(m_sum)
  
  sp_xy_rep <- 0.5 * (a_sum[1, "Sum Sq"] - a_y[1, "Sum Sq"] - a_x[1, "Sum Sq"])
  sp_xy_gen <- 0.5 * (a_sum[2, "Sum Sq"] - a_y[2, "Sum Sq"] - a_x[2, "Sum Sq"])
  sp_xy_err <- 0.5 * (a_sum["Residuals", "Sum Sq"] - a_y["Residuals", "Sum Sq"] - a_x["Residuals", "Sum Sq"])
  sp_xy_tot <- sp_xy_rep + sp_xy_gen + sp_xy_err
  
  # Partition Degrees of Freedom (Df)
  df_rep <- a_y[1, "Df"]
  df_gen <- a_y[2, "Df"]
  df_err <- a_y["Residuals", "Df"]
  df_tot <- df_rep + df_gen + df_err
  
  # Mean Products (MP_XY = SP_XY / Df)
  mp_xy_rep <- sp_xy_rep / df_rep
  mp_xy_gen <- sp_xy_gen / df_gen
  mp_xy_err <- sp_xy_err / df_err
  
  # 3. Exact Biometrical F-Calculated Values (Strictly MP_Rep / MP_Error and MP_Gen / MP_Error)
  f_rep <- mp_xy_rep / mp_xy_err
  p_rep <- 1 - stats::pf(f_rep, df_rep, df_err)
  
  f_gen <- mp_xy_gen / mp_xy_err
  p_gen <- 1 - stats::pf(f_gen, df_gen, df_err)
  
  assign_stars <- function(p_val) {
    if (is.na(p_val)) return("")
    if (p_val < 0.001) return("***")
    if (p_val < 0.01)  return("**")
    if (p_val < 0.05)  return("*")
    return("ns")
  }
  
  signif_rep <- assign_stars(p_rep)
  signif_gen <- assign_stars(p_gen)
  
  # 4. Rigorous Biometrical Covariance Components Isolation
  cov_e <- mp_xy_err
  cov_g <- (mp_xy_gen - mp_xy_err) / r
  cov_p <- (cov_g + cov_e)  # Exact biometrical relationship: (Genotypic Covariance + Error Covariance)
  
  # Build Compact ANCOVA Table
  ancova_table <- data.frame(
    "Source"    = c("Replications", "Genotypes", "Error", "Total"),
    "Df"        = c(df_rep, df_gen, df_err, df_tot),
    "SP_XY"     = c(sp_xy_rep, sp_xy_gen, sp_xy_err, sp_xy_tot),
    "MP_XY"     = c(mp_xy_rep, mp_xy_gen, mp_xy_err, NA),
    "F_value"   = c(f_rep, f_gen, NA, NA),
    "p_value"   = c(p_rep, p_gen, NA, NA),
    "Signif"    = c(signif_rep, signif_gen, "", ""),
    stringsAsFactors = FALSE
  )
  
  # ANCOVA Model for Adjusted Means (LS Means)
  ancova_full_model <- stats::lm(stats::as.formula(paste0("`", response_trait, "` ~ `", rep_col, "` + `", covariate_trait, "` + `", genotype_col, "`")), data = clean_data)
  
  covar_mean  <- mean(clean_data[[covariate_trait]], na.rm = TRUE)
  genotypes   <- levels(clean_data[[genotype_col]])
  first_rep   <- levels(clean_data[[rep_col]])[1]
  
  pred_grid <- data.frame(
    Genotype = factor(genotypes, levels = genotypes),
    Replication = factor(rep(first_rep, length(genotypes)), levels = levels(clean_data[[rep_col]]))
  )
  colnames(pred_grid) <- c(genotype_col, rep_col)
  pred_grid[[covariate_trait]] <- covar_mean
  
  pred_grid$Adjusted_Mean <- stats::predict(ancova_full_model, newdata = pred_grid)
  adj_means_df <- pred_grid[, c(genotype_col, "Adjusted_Mean")]
  rownames(adj_means_df) <- NULL
  
  # Console Presentation
  if (reporting_level >= 1) {
    cat("\n", rep("=", 97), "\n", sep = "")
    cat(" ANALYSIS OF COVARIANCE (ANCOVA) TABLE\n")
    cat(" Response Trait: ", response_trait, " | Covariate Trait: ", covariate_trait, "\n")
    cat(rep("=", 97), "\n\n", sep = "")
    
    print_table <- ancova_table
    print_table$SP_XY   <- ifelse(is.na(print_table$SP_XY), "", sprintf("%.4f", print_table$SP_XY))
    print_table$MP_XY   <- ifelse(is.na(print_table$MP_XY), "", sprintf("%.4f", print_table$MP_XY))
    print_table$F_value <- ifelse(is.na(print_table$F_value), "", sprintf("%.3f", print_table$F_value))
    print_table$p_value <- ifelse(is.na(print_table$p_value), "", format.pval(print_table$p_value, digits = 4, eps = 0.001))
    
    print(print_table, row.names = FALSE)
    
    cat(rep("-", 97), "\n", sep = "")
    cat("    Covariance Components    :    Estimates\n")
    cat("    Error Covariance         :   ", sprintf("%.4f", cov_e), "\n")
    cat("    Genotypic Covariance     :   ", sprintf("%.4f", cov_g), "\n")
    cat("    Phenotypic Covariance    :   ", sprintf("%.4f", cov_p), "\n")
    cat(rep("=", 97), "\n\n", sep = "")
  }
  
  return(invisible(list(
    ancova_table          = ancova_table,
    covariance_components = c(Cov_e = cov_e, Cov_g = cov_g, Cov_p = cov_p),
    adjusted_means        = adj_means_df
  )))
}

utils::globalVariables(c("sum_xy"))