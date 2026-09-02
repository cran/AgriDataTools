#' Multi-Level Genetic, Phenotypic, and Environmental Correlation Engine
#'
#' @description
#' The \code{compute_correlation} function calculates genotypic (\eqn{r_g}), phenotypic (\eqn{r_p}),
#' and environmental (\eqn{r_e}) correlation coefficient matrices across quantitative traits using
#' analysis of variance (ANOVA) and covariance (ANCOVA) partitions, with integrated significance flags.
#'
#' @details
#' The engine partitions variance and covariance components using mean squares (MS) and mean cross-products (MCP):
#' \deqn{r_g = \frac{Cov_g}{\sqrt{\sigma^2_{g1} \cdot \sigma^2_{g2}}}}
#' \deqn{r_p = \frac{Cov_p}{\sqrt{\sigma^2_{p1} \cdot \sigma^2_{p2}}}}
#' \deqn{r_e = \frac{Cov_e}{\sqrt{\sigma^2_{e1} \cdot \sigma^2_{e2}}}}
#'
#' @param data A \code{data.frame} containing experimental phenotypic records with \code{Genotype} and \code{Replication} (or \code{Rep}) factors.
#' @param traits A character vector specifying numeric trait columns to evaluate. Defaults to \code{NULL} for automatic detection.
#' @param genotype_col Optional character string specifying the genotype column. Defaults to \code{NULL} for auto-detection.
#' @param rep_col Optional character string specifying the replication column. Defaults to \code{NULL} for auto-detection.
#' @param reporting_level An integer flag defining console trace settings: \code{0} for silent execution and \code{1} for displaying summary matrices with significance codes. Defaults to \code{1}.
#'
#' @return Invisibly returns a structured named \code{list} containing 9 correlation and significance matrices:
#' \item{genotypic_correlation}{A \code{data.frame} matrix of genotypic correlation coefficients (\eqn{r_g}).}
#' \item{phenotypic_correlation}{A \code{data.frame} matrix of phenotypic correlation coefficients (\eqn{r_p}).}
#' \item{environmental_correlation}{A \code{data.frame} matrix of environmental correlation coefficients (\eqn{r_e}).}
#' \item{genotypic_significance}{A \code{data.frame} matrix of genotypic correlations formatted with significance stars.}
#' \item{phenotypic_significance}{A \code{data.frame} matrix of phenotypic correlations formatted with significance stars.}
#' \item{environmental_significance}{A \code{data.frame} matrix of environmental correlations formatted with significance stars.}
#' \item{genotypic_p_values}{A \code{data.frame} matrix containing raw calculated p-values for genotypic correlations.}
#' \item{phenotypic_p_values}{A \code{data.frame} matrix containing raw calculated p-values for phenotypic correlations.}
#' \item{environmental_p_values}{A \code{data.frame} matrix containing raw calculated p-values for environmental correlations.}
#'
#' @importFrom stats aov pt complete.cases
#' @export
#'
#' @examples
#' # Load your own dataset
#' data(gv_data, package = "AgriDataTools")
#' 
#' # Specify trait columns matching your dataset structure
#' traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
#' 
#' # Run correlation engine
#' corr_results <- compute_correlation(
#'   data = gv_data,
#'   traits = traits,
#'   reporting_level = 1
#' )
compute_correlation <- function(data, traits = NULL, genotype_col = NULL, rep_col = NULL, reporting_level = 1) {
  
  if (missing(data) || !is.data.frame(data)) {
    stop("CRITICAL ERROR: Input 'data' must be a valid structured data frame.", call. = FALSE)
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
    possible_rep <- c("Replication", "replication", "Replications", "Rep", "rep", "Block")
    matched_rep <- intersect(possible_rep, colnames(data))
    if (length(matched_rep) > 0) rep_col <- matched_rep[1]
    else stop("REPLICATION COLUMN FAULT: Could not auto-detect replication column. Specify 'rep_col'.", call. = FALSE)
  }
  
  # Auto-detect Traits
  if (is.null(traits)) {
    ignore_fields <- c(genotype_col, rep_col, "Block", "Line", "Cultivar")
    all_cols      <- colnames(data)
    numeric_cols  <- all_cols[sapply(data, is.numeric)]
    traits        <- setdiff(numeric_cols, ignore_fields)
  }
  
  missing_traits <- traits[!traits %in% colnames(data)]
  if (length(missing_traits) > 0) {
    stop(paste("DATASET ERROR: Traits missing from provided frame:", paste(missing_traits, collapse = ", ")), call. = FALSE)
  }
  
  num_traits <- length(traits)
  if (num_traits < 2) {
    stop("DISCOVERY ERROR: At least two numeric traits are required for correlation analysis.", call. = FALSE)
  }
  
  # Filter complete cases across selected traits and design factors
  clean_cols <- c(genotype_col, rep_col, traits)
  clean_data <- data[stats::complete.cases(data[, clean_cols]), , drop = FALSE]
  
  clean_data[[genotype_col]] <- as.factor(clean_data[[genotype_col]])
  clean_data[[rep_col]]      <- as.factor(clean_data[[rep_col]])
  
  g <- nlevels(clean_data[[genotype_col]])
  r <- nlevels(clean_data[[rep_col]])
  
  r_g <- matrix(1, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  r_p <- matrix(1, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  r_e <- matrix(1, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  
  p_g <- matrix(NA, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  p_p <- matrix(NA, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  p_e <- matrix(NA, nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
  
  df_g <- g - 1
  df_e <- (g - 1) * (r - 1)
  
  for (i in 1:num_traits) {
    for (j in i:num_traits) {
      if (i == j) next
      
      t1 <- traits[i]
      t2 <- traits[j]
      
      fit1 <- stats::aov(clean_data[[t1]] ~ clean_data[[rep_col]] + clean_data[[genotype_col]])
      fit2 <- stats::aov(clean_data[[t2]] ~ clean_data[[rep_col]] + clean_data[[genotype_col]])
      
      sum_fit1 <- summary(fit1)[[1]]
      sum_fit2 <- summary(fit2)[[1]]
      
      ms_g1 <- sum_fit1[2, "Mean Sq"]
      ms_g2 <- sum_fit2[2, "Mean Sq"]
      ms_e1 <- sum_fit1["Residuals", "Mean Sq"]
      ms_e2 <- sum_fit2["Residuals", "Mean Sq"]
      
      # Cross-product ANOVA for Covariance
      cp_fit <- stats::aov((clean_data[[t1]] + clean_data[[t2]]) ~ clean_data[[rep_col]] + clean_data[[genotype_col]])
      sum_cp <- summary(cp_fit)[[1]]
      
      mcp_g_total <- sum_cp[2, "Mean Sq"]
      mcp_e_total <- sum_cp["Residuals", "Mean Sq"]
      
      mcp_g <- (mcp_g_total - ms_g1 - ms_g2) / 2
      mcp_e <- (mcp_e_total - ms_e1 - ms_e2) / 2
      
      sigma2_e1 <- ms_e1
      sigma2_e2 <- ms_e2
      cov_e     <- mcp_e
      
      sigma2_g1 <- max(0.00001, (ms_g1 - ms_e1) / r)
      sigma2_g2 <- max(0.00001, (ms_g2 - ms_e2) / r)
      cov_g     <- (mcp_g - mcp_e) / r
      
      sigma2_p1 <- sigma2_g1 + (sigma2_e1 / r)
      sigma2_p2 <- sigma2_g2 + (sigma2_e2 / r)
      cov_p     <- cov_g + (cov_e / r)
      
      val_rg <- cov_g / sqrt(sigma2_g1 * sigma2_g2)
      val_rp <- cov_p / sqrt(sigma2_p1 * sigma2_p2)
      val_re <- cov_e / sqrt(sigma2_e1 * sigma2_e2)
      
      r_g[i, j] <- r_g[j, i] <- max(-1, min(1, val_rg))
      r_p[i, j] <- r_p[j, i] <- max(-1, min(1, val_rp))
      r_e[i, j] <- r_e[j, i] <- max(-1, min(1, val_re))
      
      # Standard Error & t-tests
      t_g <- (r_g[i, j] * sqrt(df_g)) / sqrt(max(0.00001, 1 - r_g[i, j]^2))
      t_p <- (r_p[i, j] * sqrt(df_g)) / sqrt(max(0.00001, 1 - r_p[i, j]^2))
      t_e <- (r_e[i, j] * sqrt(df_e)) / sqrt(max(0.00001, 1 - r_e[i, j]^2))
      
      p_g[i, j] <- p_g[j, i] <- 2 * stats::pt(abs(t_g), df = df_g, lower.tail = FALSE)
      p_p[i, j] <- p_p[j, i] <- 2 * stats::pt(abs(t_p), df = df_g, lower.tail = FALSE)
      p_e[i, j] <- p_e[j, i] <- 2 * stats::pt(abs(t_e), df = df_e, lower.tail = FALSE)
    }
  }
  
  format_with_stars <- function(r_mat, p_mat) {
    formatted_mat <- matrix("", nrow = num_traits, ncol = num_traits, dimnames = list(traits, traits))
    for (i in 1:num_traits) {
      for (j in 1:num_traits) {
        if (i == j) {
          formatted_mat[i, j] <- "1.0000"
        } else {
          p <- p_mat[i, j]
          stars <- if (is.na(p)) "ns" else if (p <= 0.001) "***" else if (p <= 0.01) "**" else if (p <= 0.05) "*" else "ns"
          formatted_mat[i, j] <- paste0(sprintf("%.4f", r_mat[i, j]), " ", stars)
        }
      }
    }
    return(as.data.frame(formatted_mat, stringsAsFactors = FALSE))
  }
  
  sig_g_df <- format_with_stars(r_g, p_g)
  sig_p_df <- format_with_stars(r_p, p_p)
  sig_e_df <- format_with_stars(r_e, p_e)
  
  if (reporting_level >= 1) {
    cat("\n", rep("=", 85), "\n", sep = "")
    cat(" MULTI-LEVEL CORRELATION ANALYSIS (Genotypic, Phenotypic & Environmental)\n")
    cat(" Significance Codes: *** p<=0.001, ** p<=0.01, * p<=0.05, ns = non-significant\n")
    cat(rep("=", 85), "\n\n", sep = "")
    
    cat("\nTABLE 1: GENOTYPIC CORRELATION MATRIX (rg)\n")
    print(sig_g_df, quote = FALSE)
    
    cat("\nTABLE 2: PHENOTYPIC CORRELATION MATRIX (rp)\n")
    print(sig_p_df, quote = FALSE)
    
    cat("\nTABLE 3: ENVIRONMENTAL CORRELATION MATRIX (re)\n")
    print(sig_e_df, quote = FALSE)
    cat("\n", rep("=", 85), "\n\n", sep = "")
  }
  
  return(invisible(list(
    genotypic_correlation      = as.data.frame(r_g),
    phenotypic_correlation     = as.data.frame(r_p),
    environmental_correlation  = as.data.frame(r_e),
    genotypic_significance     = sig_g_df,
    phenotypic_significance    = sig_p_df,
    environmental_significance = sig_e_df,
    genotypic_p_values         = as.data.frame(p_g),
    phenotypic_p_values        = as.data.frame(p_p),
    environmental_p_values     = as.data.frame(p_e)
  )))
}