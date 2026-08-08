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
#' @param reporting_level An integer flag defining console trace settings: \code{0} for silent execution and \code{1} for displaying summary matrices with significance codes. Defaults to \code{1}.
#'
#' @return Invisibly returns a structured named \code{list} of class \code{"list"} containing 9 correlation and significance matrices:
#' \item{genotypic_correlation}{A \code{data.frame} matrix of genotypic correlation coefficients (\eqn{r_g}) between evaluated traits.}
#' \item{phenotypic_correlation}{A \code{data.frame} matrix of phenotypic correlation coefficients (\eqn{r_p}) between evaluated traits.}
#' \item{environmental_correlation}{A \code{data.frame} matrix of environmental correlation coefficients (\eqn{r_e}) between evaluated traits.}
#' \item{genotypic_significance}{A \code{data.frame} matrix of genotypic correlations formatted with significance stars (\code{***}, \code{**}, \code{*}, or \code{ns}).}
#' \item{phenotypic_significance}{A \code{data.frame} matrix of phenotypic correlations formatted with significance stars.}
#' \item{environmental_significance}{A \code{data.frame} matrix of environmental correlations formatted with significance stars.}
#' \item{genotypic_p_values}{A \code{data.frame} matrix containing raw calculated two-tailed p-values for genotypic correlations.}
#' \item{phenotypic_p_values}{A \code{data.frame} matrix containing raw calculated two-tailed p-values for phenotypic correlations.}
#' \item{environmental_p_values}{A \code{data.frame} matrix containing raw calculated two-tailed p-values for environmental correlations.}
#' 
#' If \code{reporting_level >= 1}, formatted correlation matrices with significance flags are printed directly to the console before returning the output object.
#'
#' @importFrom stats aov pt
#' @export
#'
#' @examples
#' # Load benchmark breeding dataset
#' data(gv_data, package = "AgriDataTools")
#' 
#' # Define trait columns
#' my_traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
#' 
#' # Run multi-level correlation engine
#' corr_results <- compute_correlation(
#'   data = gv_data,
#'   traits = my_traits,
#'   reporting_level = 1
#' )
compute_correlation <- function(data, traits = NULL, reporting_level = 1) {
  
  if (missing(data) || !is.data.frame(data)) {
    stop("CRITICAL ERROR: Input 'data' must be a valid structured data frame.", call. = FALSE)
  }
  
  if (!all(c("Genotype", "Replication") %in% colnames(data)) && !all(c("Genotype", "Rep") %in% colnames(data))) {
    stop("CRITICAL ERROR: Data must contain 'Genotype' and 'Replication' (or 'Rep') columns for genetic partitioning.", call. = FALSE)
  }
  
  rep_col <- if ("Replication" %in% colnames(data)) "Replication" else "Rep"
  
  if (is.null(traits)) {
    ignore_fields <- c("Genotype", rep_col, "Block", "Line", "Cultivar")
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
  
  data$Genotype <- as.factor(data$Genotype)
  data[[rep_col]] <- as.factor(data[[rep_col]])
  
  g <- length(levels(data$Genotype))
  r <- length(levels(data[[rep_col]]))
  
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
      
      fit1 <- stats::aov(data[[t1]] ~ data[[rep_col]] + data$Genotype)
      fit2 <- stats::aov(data[[t2]] ~ data[[rep_col]] + data$Genotype)
      
      ms_g1 <- summary(fit1)[[1]]["data$Genotype", "Mean Sq"]
      ms_g2 <- summary(fit2)[[1]]["data$Genotype", "Mean Sq"]
      ms_e1 <- summary(fit1)[[1]]["Residuals", "Mean Sq"]
      ms_e2 <- summary(fit2)[[1]]["Residuals", "Mean Sq"]
      
      cp_fit <- stats::aov((data[[t1]] + data[[t2]]) ~ data[[rep_col]] + data$Genotype)
      mcp_g_total <- summary(cp_fit)[[1]]["data$Genotype", "Mean Sq"]
      mcp_e_total <- summary(cp_fit)[[1]]["Residuals", "Mean Sq"]
      
      mcp_g <- (mcp_g_total - ms_g1 - ms_g2) / 2
      mcp_e <- (mcp_e_total - ms_e1 - ms_e2) / 2
      
      sigma2_e1 <- ms_e1
      sigma2_e2 <- ms_e2
      cov_e    <- mcp_e
      
      sigma2_g1 <- (ms_g1 - ms_e1) / r
      sigma2_g2 <- (ms_g2 - ms_e2) / r
      cov_g    <- (mcp_g - mcp_e) / r
      
      sigma2_p1 <- sigma2_g1 + (sigma2_e1 / r)
      sigma2_p2 <- sigma2_g2 + (sigma2_e2 / r)
      cov_p    <- cov_g + (cov_e / r)
      
      val_rg <- cov_g / sqrt(max(0.00001, sigma2_g1 * sigma2_g2))
      val_rp <- cov_p / sqrt(max(0.00001, sigma2_p1 * sigma2_p2))
      val_re <- cov_e / sqrt(max(0.00001, sigma2_e1 * sigma2_e2))
      
      r_g[i, j] <- r_g[j, i] <- max(-1, min(1, val_rg))
      r_p[i, j] <- r_p[j, i] <- max(-1, min(1, val_rp))
      r_e[i, j] <- r_e[j, i] <- max(-1, min(1, val_re))
      
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
          stars <- if (p <= 0.001) "***" else if (p <= 0.01) "**" else if (p <= 0.05) "*" else "ns"
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
    cat("\n", rep("=", 80), "\n", sep = "")
    cat(" MULTI-LEVEL CORRELATION ANALYSIS (Genotypic, Phenotypic & Environmental)\n")
    cat(" Significance Codes: *** p<=0.001, ** p<=0.01, * p<=0.05, ns = non-significant\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    cat("--- GENOTYPIC CORRELATION MATRIX (rg) WITH SIGNIFICANCE ---\n")
    print(sig_g_df, quote = FALSE)
    
    cat("\n--- PHENOTYPIC CORRELATION MATRIX (rp) WITH SIGNIFICANCE ---\n")
    print(sig_p_df, quote = FALSE)
    
    cat("\n--- ENVIRONMENTAL CORRELATION MATRIX (re) WITH SIGNIFICANCE ---\n")
    print(sig_e_df, quote = FALSE)
    cat("\n", rep("=", 80), "\n\n", sep = "")
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