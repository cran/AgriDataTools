#' Principal Component Analysis for Agronomic Traits
#'
#' Performs Principal Component Analysis (PCA) on targeted quantitative agronomic 
#' parameters. Supports dynamic trait mapping to convert trait abbreviations into 
#' full descriptive names, dynamic genotype column recognition, automatic replication 
#' aggregation to genotypic means, and computes modern multivariate metrics including 
#' Kaiser-Guttman retention rules, eigenvector loadings, and percentage trait contributions.
#'
#' @param data A data frame containing genotype information and trait columns.
#' @param traits A character vector specifying the exact trait column names to include.
#' @param genotype_col An optional character string specifying the column name for genotypes. If \code{NULL}, automatic column detection is performed. Defaults to \code{NULL}.
#' @param scale Logical. If TRUE (default), variables are standardized to unit variance.
#' @param reporting_level Integer. Control output verbosity: 0 (silent), 1 (summary), or 2 (exhaustive full reporting). Defaults to \code{2}.
#' @param trait_lookup An optional named character vector for mapping trait abbreviations 
#'   to full descriptive labels (e.g., c("PH" = "Plant Height")).
#'
#' @return A structured named \code{list} containing 6 multivariate components:
#' \item{pca_object}{The raw \code{prcomp} output object based on genotypic means.}
#' \item{eigenvalues}{Data frame of eigenvalues, variance percentages, cumulative variance, and Kaiser retention decision.}
#' \item{loadings}{Data frame of eigenvector loadings matrix.}
#' \item{contributions}{Data frame of percentage contributions of each trait across components.}
#' \item{cos2}{Data frame representing quality of representation (\eqn{Cos^2}) for each trait.}
#' \item{scores}{Data frame of principal component scores assigned to unique genotypes.}
#'
#' @importFrom stats complete.cases prcomp aggregate na.omit
#' @export
#'
#' @examples
#' library(AgriDataTools)
#' data("gv_data", package = "AgriDataTools")
#' 
#' # Define custom trait mapping before running (Edit names as needed)
#' custom_traits_map <- c(
#'     "PH"   = "Plant Height",
#'     "SL"   = "Spike Length",
#'     "PL"   = "Peduncle Length",
#'     "NOT"  = "Number of Tillers",
#'     "NOSS" = "Number of Spikelets per Spike",
#'     "TGW"  = "Thousand Grain Weight",
#'     "GYPM" = "Grain Yield per Meter"
#' )
#' 
#' # Run Modern PCA with full trait names
#' pca_results <- analyze_pca(
#'     data = gv_data,
#'     traits = names(custom_traits_map),
#'     genotype_col = "Genotype",
#'     trait_lookup = custom_traits_map
#' )
#'
analyze_pca <- function(data, traits, genotype_col = NULL, scale = TRUE, reporting_level = 2, trait_lookup = NULL) {
  
  if (missing(data) || missing(traits)) {
    stop("CRITICAL PCA FAULT: Missing core arguments.", call. = FALSE)
  }
  
  missing_traits <- traits[!traits %in% colnames(data)]
  if (length(missing_traits) > 0) {
    stop(paste("DATASET ERROR: Traits missing from provided data frame:", paste(missing_traits, collapse = ", ")), call. = FALSE)
  }
  
  # Dynamic Genotype Column Detection (Zero Hard-coding)
  if (is.null(genotype_col)) {
    possible_geno_cols <- c("Genotype", "genotype", "Genotypes", "genotypes", "Gen", "gen", "Variety", "variety", "Line", "line", "Cultivar", "cultivar")
    matched_col <- intersect(possible_geno_cols, colnames(data))
    if (length(matched_col) > 0) {
      genotype_col <- matched_col[1]
    } else {
      stop("CRITICAL DATA FAULT: Could not automatically detect a genotype column. Please specify 'genotype_col'.", call. = FALSE)
    }
  }
  
  if (!genotype_col %in% colnames(data)) {
    stop(paste("DATASET ERROR: Specified genotype column '", genotype_col, "' not found in data frame."), call. = FALSE)
  }
  
  # Extract relevant columns (Genotype + Traits)
  sub_data <- data[, c(genotype_col, traits), drop = FALSE]
  
  # Handle missing values safely
  complete_rows <- stats::complete.cases(sub_data)
  if (sum(complete_rows) < nrow(sub_data)) {
    sub_data <- sub_data[complete_rows, , drop = FALSE]
  }
  
  # Ensure all trait columns are numeric
  for (col in traits) {
    sub_data[[col]] <- as.numeric(sub_data[[col]])
  }
  
  # Automatic Replication Aggregation: Compute mean values per unique genotype
  # This guarantees that if a user has 40 genotypes with 3 replications each, 
  # it correctly aggregates into exactly 40 unique genotypic rows.
  agg_formula <- stats::as.formula(paste("~", genotype_col))
  genotypic_means <- stats::aggregate(sub_data[, traits, drop = FALSE], 
                                      by = list(Genotype = sub_data[[genotype_col]]), 
                                      FUN = mean, na.rm = TRUE)
  
  unique_genotypes <- as.character(genotypic_means$Genotype)
  clean_matrix <- genotypic_means[, traits, drop = FALSE]
  rownames(clean_matrix) <- unique_genotypes
  
  # Dynamic Trait Lookup Mapping (Zero background hardcoding)
  if (!is.null(trait_lookup) && is.vector(trait_lookup) && !is.null(names(trait_lookup))) {
    current_cols <- colnames(clean_matrix)
    mapped_cols <- sapply(current_cols, function(code) {
      code_clean <- trimws(toupper(code))
      lookup_names <- trimws(toupper(names(trait_lookup)))
      if (code_clean %in% lookup_names) {
        return(as.character(trait_lookup[which(lookup_names == code_clean)[1]]))
      } else {
        return(code)
      }
    })
    colnames(clean_matrix) <- mapped_cols
  }
  
  # Core PCA Computation on Genotypic Means Matrix
  pca_execution <- stats::prcomp(clean_matrix, center = TRUE, scale. = scale)
  
  calculated_eigenvalues <- pca_execution$sdev^2
  total_variance_sum <- sum(calculated_eigenvalues)
  variance_proportions <- (calculated_eigenvalues / total_variance_sum) * 100
  cumulative_variance <- cumsum(variance_proportions)
  
  # Kaiser-Guttman Rule Retention Check (Eigenvalue >= 1.0)
  retained_kaiser <- calculated_eigenvalues >= 1.0
  
  eigen_summary_table <- data.frame(
    "Component"          = paste0("PC", 1:length(calculated_eigenvalues)),
    "Eigenvalue"         = calculated_eigenvalues,
    "Variance_Percent"   = variance_proportions,
    "Cumulative_Percent" = cumulative_variance,
    "Retain_Kaiser"      = ifelse(retained_kaiser, "Yes (EV >= 1)", "No"),
    stringsAsFactors     = FALSE
  )
  
  # Loadings (Eigenvectors)
  factor_loadings_matrix <- as.data.frame(pca_execution$rotation)
  
  # Trait Contribution Percentages
  trait_contributions <- as.data.frame(apply(factor_loadings_matrix^2, 2, function(x) (x / sum(x)) * 100))
  
  # Cos2 Quality of Representation
  cos2_matrix <- as.data.frame(factor_loadings_matrix^2)
  
  # Scores Matrix and Genotype Label Assignment for Unique Genotypes
  individual_scores <- as.data.frame(pca_execution$x)
  rownames(individual_scores) <- make.unique(unique_genotypes)
  individual_scores <- cbind(
    "Genotype" = as.factor(unique_genotypes),
    individual_scores,
    stringsAsFactors = FALSE
  )
  
  # Console Reporting
  if (reporting_level >= 1) {
    cat("\n", rep("=", 85), "\n", sep = "")
    cat(" PRINCIPAL COMPONENT ANALYSIS (PCA) SUMMARY (Genotypic Means)\n")
    cat(rep("=", 85), "\n\n", sep = "")
    
    cat(rep("-", 75), "\n", sep = "")
    cat(" 1. PRINCIPAL COMPONENT EIGENVALUE & VARIANCE SUMMARY (Kaiser Rule)\n")
    cat(rep("-", 75), "\n", sep = "")
    
    eigen_output <- eigen_summary_table
    numeric_cols <- c("Eigenvalue", "Variance_Percent", "Cumulative_Percent")
    eigen_output[, numeric_cols] <- round(eigen_output[, numeric_cols], 4)
    print(eigen_output, row.names = FALSE)
    
    cat("\n", rep("-", 75), "\n", sep = "")
    cat(" 2. TRAIT EIGENVECTOR LOADINGS MATRIX\n")
    cat(rep("-", 75), "\n", sep = "")
    print(round(factor_loadings_matrix, 4))
    
    cat("\n", rep("-", 75), "\n", sep = "")
    cat(" 3. TRAIT CONTRIBUTIONS TO COMPONENTS (% Contribution)\n")
    cat(rep("-", 75), "\n", sep = "")
    print(round(trait_contributions, 2))
    
    cat(rep("=", 85), "\n\n", sep = "")
  }
  
  return(invisible(list(
    pca_object    = pca_execution,
    eigenvalues   = eigen_summary_table,
    loadings      = factor_loadings_matrix,
    contributions = trait_contributions,
    cos2          = cos2_matrix,
    scores        = individual_scores
  )))
}