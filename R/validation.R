#' Rigid Multi-Layered Structural Data Validation and Matrix Integrity Engine
#'
#' @description
#' The \code{validate_agri_data} function serves as the primary data defense matrix of the
#' \code{AgriDataTools} package. It performs automated multi-dimensional quality control,
#' type alignment checks, and semantic structural audits on agricultural research datasets
#' before passing them to downstream quantitative genetic workflows.
#'
#' @details
#' In quantitative genetics, downstream models like ANOVA, Heritability estimations, and Path analysis
#' are highly vulnerable to layout irregularities. Silent formatting issues can bias variance
#' components or trigger execution failures.
#'
#' \code{validate_agri_data} runs a dynamic defensive pipeline:
#' \enumerate{
#'   \item \strong{Object Matrix Verification:} Ensures input is a valid data.frame with observations.
#'   \item \strong{Dynamic Header Discovery:} Auto-detects Genotype, Replication, and Trait columns without rigid hardcoding.
#'   \item \strong{Factor Balance Audits:} Verifies minimum levels for lines and blocks.
#'   \item \strong{Numeric Type & Range Compliance:} Validates numeric integrity and flags biological anomalies (e.g., negative values).
#'   \item \strong{Orthogonality Sweep:} Evaluates design layout cell frequencies for balance.
#' }
#'
#' @param data A non-null \code{data.frame} containing experimental field trial records.
#' @param reporting_level An integer mapping scale for console trace output: \code{0} for silent,
#'   \code{1} for clean summary status report. Defaults to \code{1}.
#'
#' @return A logical scalar \code{TRUE} if the dataset completely satisfies the operational constraints
#'   of the biometrical pipeline. Throws informative errors if critical structural faults are detected.
#'
#' @references
#' \itemize{
#'   \item Cochran, W.G. and Cox, G.M. (1957). Experimental Designs. 2nd Edition, John Wiley & Sons, New York.
#'   \item Falconer, D.S. and Mackay, T.F.C. (1996). Introduction to Quantitative Genetics. 4th Edition, Longman, Essex.
#' }
#'
#' @export
#'
#' @examples
#' # Trigger dynamic validation sweep
#' validation_status <- validate_agri_data(data = gv_data)
#'
validate_agri_data <- function(data, reporting_level = 1) {
  
  # =========================================================================
  # BLOCK 1: OBJECT TYPE AND NULL VALUE DEFENSE LAYER
  # =========================================================================
  if (missing(data)) {
    stop("CRITICAL VALIDATION FAULT: Argument 'data' is missing. You must supply an active agricultural data frame.", call. = FALSE)
  }
  
  if (is.null(data)) {
    stop("CRITICAL DATA POINTER FAULT: Input object 'data' resolves to NULL reference.", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop(paste0("CRITICAL CLASS TYPE CONFLICT: Input object must be a 'data.frame'. Received class: [",
                paste(class(data), collapse = ", "), "]. Cast your object using as.data.frame()."), call. = FALSE)
  }
  
  total_rows <- nrow(data)
  total_cols <- ncol(data)
  
  if (total_rows == 0) {
    stop("CRITICAL MATRIX EMPTY FAULT: The data frame contains zero rows.", call. = FALSE)
  }
  
  if (total_cols == 0) {
    stop("CRITICAL ARCHITECTURE FAULT: The data frame contains zero columns.", call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 2: DYNAMIC COLUMN HEADER IDENTIFICATION (NO HARDCODING)
  # =========================================================================
  active_headers <- colnames(data)
  
  # Dynamic matching for Genotype factor column
  g_patterns <- "(?i)^genotype$|^genotypes$|^line$|^lines$|^cultivar$|^cultivars$|^variety$|^varieties$"
  g_col_match <- grep(g_patterns, active_headers, value = TRUE)
  
  if (length(g_col_match) == 0) {
    stop(paste0("CRITICAL STRUCTURAL ATTRIBUTE MISMATCH: Could not discover 'Genotype' identification column.\n",
                "Available Columns: [ ", paste(active_headers, collapse = ", "), " ]"), call. = FALSE)
  }
  genotype_col <- g_col_match[1]
  
  # Dynamic matching for Replication factor column
  r_patterns <- "(?i)^replication$|^replications$|^rep$|^reps$|^block$|^blocks$"
  r_col_match <- grep(r_patterns, active_headers, value = TRUE)
  
  if (length(r_col_match) == 0) {
    stop(paste0("CRITICAL STRUCTURAL ATTRIBUTE MISMATCH: Could not discover 'Replication' column.\n",
                "Available Columns: [ ", paste(active_headers, collapse = ", "), " ]"), call. = FALSE)
  }
  replication_col <- r_col_match[1]
  
  # Dynamic discovery of phenotypic trait columns (all numeric columns excluding G and R)
  non_trait_cols <- c(genotype_col, replication_col, "Env", "Environment", "Loc", "Location", "Year")
  numeric_cols   <- active_headers[sapply(data, is.numeric)]
  trait_cols     <- setdiff(numeric_cols, non_trait_cols)
  
  if (length(trait_cols) == 0) {
    stop("CRITICAL TRAIT INVENTORY FAULT: Failed to discover any numeric phenotypic traits in the dataset.", call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 3: CATEGORICAL FACTOR BALANCING AUDIT
  # =========================================================================
  raw_genotypes    <- data[[genotype_col]]
  raw_replications <- data[[replication_col]]
  
  unique_genotypes    <- unique(stats::na.omit(raw_genotypes))
  unique_replications <- unique(stats::na.omit(raw_replications))
  
  total_unique_g <- length(unique_genotypes)
  total_unique_r <- length(unique_replications)
  
  if (total_unique_g < 2) {
    stop(paste0("INTEGRITY COMPLIANCE FAULT: Variance partitioning requires a minimum of 2 distinct genotypes. Discovered: ", total_unique_g), call. = FALSE)
  }
  
  if (total_unique_r < 2) {
    stop(paste0("DESIGN INFRASTRUCTURE FAULT: Replicated designs require a minimum of 2 replication blocks. Discovered: ", total_unique_r), call. = FALSE)
  }
  
  # Check for hidden leading/trailing whitespace in genotype names
  if (is.character(raw_genotypes) || is.factor(raw_genotypes)) {
    char_g <- as.character(raw_genotypes)
    if (any(char_g != trimws(char_g))) {
      warning("DATA QUALITY ANOMALY: Hidden leading or trailing spaces detected inside 'Genotype' factor entries.", call. = FALSE)
    }
  }
  
  # =========================================================================
  # BLOCK 4: QUANTITATIVE VECTOR INTEGRITY & BIOLOGICAL RANGE INSPECTIONS
  # =========================================================================
  for (trait in trait_cols) {
    target_vector <- data[[trait]]
    
    # Verify numeric class
    if (!is.numeric(target_vector)) {
      stop(paste0("CRITICAL QUANTITATIVE TYPE VIOLATION: Trait column [ ", trait, " ] is not numeric. Class: ", paste(class(target_vector), collapse = ", ")), call. = FALSE)
    }
    
    # Missing data ratio profiling
    missing_count <- sum(is.na(target_vector))
    missing_ratio <- missing_count / total_rows
    
    if (missing_count > 0) {
      if (missing_ratio > 0.50) {
        stop(paste0("CRITICAL DATA LOSS: Trait '", trait, "' has exceeded the 50% missing data limit (NA ratio: ", round(missing_ratio * 100, 2), "%)."), call. = FALSE)
      }
      if (reporting_level >= 1) {
        warning(paste0("DATA LOSS ALERT: Trait '", trait, "' has ", missing_count, " missing observations (", round(missing_ratio * 100, 2), "%)."), call. = FALSE)
      }
    }
    
    # Negative biological bound inspection
    clean_elements <- target_vector[!is.na(target_vector)]
    if (length(clean_elements) > 0) {
      if (any(clean_elements < 0)) {
        negative_offsets <- which(target_vector < 0)
        stop(paste0("BIOLOGICAL BOUND EXCEPTION: Negative numeric values detected in trait [ ", trait,
                    " ] at row indices: [ ", paste(negative_offsets, collapse = ", "), " ]."), call. = FALSE)
      }
    }
  }
  
  # =========================================================================
  # BLOCK 5: CROSS-TABULATION ORTHOGONALITY ANALYSIS
  # =========================================================================
  frequency_table <- table(as.factor(data[[genotype_col]]), as.factor(data[[replication_col]]))
  
  if (any(frequency_table == 0)) {
    warning("DESIGN METRIC WARNING: Unbalanced structure detected. Some Genotype x Replication combinations contain 0 observations.", call. = FALSE)
  }
  
  if (any(frequency_table > 1)) {
    warning("DESIGN DEFICIENCY WARNING: Multiple entries detected for identical Genotype-Replication cells (subsampling records).", call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 6: CONSOLE REPORTING (CLEAN & MINIMAL)
  # =========================================================================
  if (reporting_level >= 1) {
    cat("\nDATA VALIDATION REPORT\n")
    cat(rep("-", 70), "\n", sep = "")
    cat(" SUCCESS: Dataset passed all structural and biometrical integrity audits.\n")
    cat(sprintf(" Verified Rows: %d | Genotypes: %d | Replications: %d | Phenotypic Traits: %d\n",
                total_rows, total_unique_g, total_unique_r, length(trait_cols)))
    cat(rep("-", 70), "\n\n", sep = "")
  }
  
  return(TRUE)
}