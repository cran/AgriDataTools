#' Rigid Multi-Layered Structural Data Validation and Matrix Integrity Engine
#'
#' @description
#' The `validate_agri_data` function serves as the primary data defense matrix of the
#' \code{AgriDataTools} package. It is engineered to perform exhaustive, multi-dimensional
#' quality control, type alignment checks, and semantic structural audits on user-provided
#' datasets before they are passed to sensitive breeding and quantitative genetic workflows.
#'
#' @details
#' In biometric computing, agricultural research, and plant breeding quantitative data analysis,
#' downstream models like ANOVA, Heritability estimations, and Path analysis are highly vulnerable
#' to layout irregularities. Silent formatting issues inside spreadsheets can bias variance
#' components or trigger generic engine failures.
#'
#' To solve this, `validate_agri_data` runs an automated, multi-tiered defensive pipeline:
#' \enumerate{
#'   \item \strong{Object and Class Matrix Verification:} Ensures the input object is a true dataframe.
#'   \item \strong{Dimension Capability Evaluation:} Checks for absolute row/column structural thresholds.
#'   \item \strong{Strict Column Header Matching:} Verifies the existence of structural keys and core traits.
#'   \item \strong{Categorical Factor Integrity Audits:} Verifies grouping layouts for Genotypes and Replications.
#'   \item \strong{Quantitative Type Compliance Testing:} Validates phenotypic trait vectors for numeric compliance.
#'   \item \strong{Biological Range and Bound Inspections:} Screens for mathematical anomalies like negative values.
#'   \item \strong{Missing Value (NA) Variance Profiling:} Analyzes missing data distribution across blocks.
#' }
#'
#' @param data A non-null \code{data.frame} containing the experimental field trial records.
#' @param strict_mode A logical scalar. If \code{TRUE}, the validation engine demands an absolute structural
#'   match with the core dataset archetype: exactly 120 rows, 40 distinct genotypes, 3 replications,
#'   and all 7 mandatory phenotypic traits (\code{PH}, \code{SL}, \code{PL}, \code{NOT}, \code{NOSS},
#'   \code{TGW}, \code{GYPM}). Defaults to \code{TRUE}.
#' @param reporting_level An integer mapping scale indicating console log verbosity: \code{0} for dead silent,
#'   \code{1} for critical system steps, and \code{2} for exhaustive vector diagnostics. Defaults to \code{2}.
#'
#' @return A logical scalar \code{TRUE} if the dataset completely satisfies the operational constraints
#'   of the biometrical pipeline. If any structural non-compliance is detected, it throws a highly structured,
#'   informative exception detailing the exact coordinates of the failure.
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
#' # Assuming package environments are active and datasets are loaded via data(gv_data)
#' if (interactive()) {
#'     # Trigger standard validation sweep
#'     validation_status <- validate_agri_data(data = gv_data, strict_mode = TRUE)
#'     message("Is dataset operational? ", validation_status)
#' }
validate_agri_data <- function(data,
                               strict_mode = TRUE,
                               reporting_level = 2) {
  
  # =========================================================================
  # BLOCK 1: ENVIRONMENT DIAGNOSTICS & SYSTEM INITIALIZATION
  # =========================================================================
  timestamp_init <- Sys.time()
  
  if (reporting_level >= 1) {
    cat(rep("=", 85), "\n", sep = "")
    cat("AGRIDATATOOLS PACKAGED ENGINE v0.1.0 - MULTI-LAYER DATA AUDIT MATRIX\n")
    cat("System Validation Sweep Timestamp: ", as.character(timestamp_init), "\n")
    cat("Operation Mode Setting: ", ifelse(strict_mode, "STRICT ARCHETYPE MANDATE", "FLEXIBLE LAYOUT"), "\n")
    cat(rep("-", 85), "\n", sep = "")
  }
  
  # =========================================================================
  # BLOCK 2: OBJECT TYPE AND NULL VALUE DEFENSE LAYER
  # =========================================================================
  if (missing(data)) {
    stop("CRITICAL VALIDATION FAULT: The target dataset argument 'data' is missing from the function call. You must supply an active agricultural data frame to proceed.", call. = FALSE)
  }
  
  if (is.null(data)) {
    stop("CRITICAL DATA POINTER FAULT: The input object passed to 'data' resolves to a NULL reference. The analysis pipeline cannot extract matrices from empty memory vectors.", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop(paste0("CRITICAL CLASS TYPE CONFLICT: Input object must be a standard base 'data.frame'. The received object possesses class attributes: [",
                paste(class(data), collapse = ", "), "]. Please cast your object using as.data.frame() before validation."), call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 3: DIMENSION CAPABILITY AND EMPTY ROW EVALUATION
  # =========================================================================
  total_rows <- nrow(data)
  total_cols <- ncol(data)
  
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - DIMENSION]: Analyzing initial structure matrix indices...\n")
    cat("[DIAGNOSTIC - DIMENSION]: Total Rows Detected:     ", total_rows, "\n")
    cat("[DIAGNOSTIC - DIMENSION]: Total Columns Detected: ", total_cols, "\n")
  }
  
  if (total_rows == 0) {
    stop("CRITICAL MATRIX EMPTY FAULT: The data frame contains exactly zero rows. There are no agronomic field records available to compile models.", call. = FALSE)
  }
  
  if (total_cols == 0) {
    stop("CRITICAL ARCHITECTURE FAULT: The data frame contains exactly zero columns. No variable vectors were isolated.", call. = FALSE)
  }
  
  if (strict_mode) {
    if (total_rows != 120) {
      stop(paste0("STRICT DEVIATION EXCEPTION: Strict mode requires the exact 'gv_data' profile composition of exactly 120 rows. Your dataset contains [",
                  total_rows, "] rows. Adjust 'strict_mode = FALSE' if processing custom dimensions."), call. = FALSE)
    }
  }
  
  # =========================================================================
  # BLOCK 4: ABSOLUTE COLUMN HEADER IDENTIFICATION AND SANITY SWEEP
  # =========================================================================
  active_headers <- colnames(data)
  
  # Define mandatory experimental structural design anchors
  mandatory_anchors <- c("Genotype", "Replication")
  missing_anchors <- mandatory_anchors[!(mandatory_anchors %in% active_headers)]
  
  if (length(missing_anchors) > 0) {
    stop(paste0("\n", rep("*", 85), "\n",
                "CRITICAL STRUCTURAL ATTRIBUTE MISMATCH\n",
                "The data validation matrix failed to discover the mandatory experimental variables.\n",
                "Missing Identification Columns: [ ", paste(missing_anchors, collapse = " | "), " ]\n",
                "Available Header Roster in Dataset: [ ", paste(active_headers, collapse = ", "), " ]\n",
                "Correction Required: Re-align variable assignments to match capitalization conventions exactly.\n",
                rep("*", 85), "\n"), call. = FALSE)
  }
  
  # Define target plant breeding phenotypic traits requested by user specifications
  mandatory_traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
  missing_traits <- mandatory_traits[!(mandatory_traits %in% active_headers)]
  
  if (strict_mode) {
    if (length(missing_traits) > 0) {
      stop(paste0("STRICT TRAIT INVENTORY FAULT: Strict validation mode is active. The dataset lacks core plant traits:\n",
                  "Missing Target Traits: [ ", paste(missing_traits, collapse = ", "), " ]\n",
                  "Ensure your dataset provides all 7 indices for downstream multivariate models."), call. = FALSE)
    }
  } else {
    # If not in strict mode, ensure at least one trait exists
    available_traits_subset <- mandatory_traits[mandatory_traits %in% active_headers]
    if (length(available_traits_subset) == 0) {
      stop("FLEXIBLE LAYOUT FAULT: The dataset does not contain any recognizable phenotypic traits from the core specification roster.", call. = FALSE)
    }
  }
  
  # =========================================================================
  # BLOCK 5: CATEGORICAL FACTOR BALANCING AND COMPLIANCE AUDIT
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - CATEGORICAL]: Auditing grouping factors and experimental alignment...\n")
  }
  
  raw_genotypes <- data$Genotype
  raw_replications <- data$Replication
  
  unique_genotypes <- unique(raw_genotypes[!is.na(raw_genotypes)])
  unique_replications <- unique(raw_replications[!is.na(raw_replications)])
  
  total_unique_g <- length(unique_genotypes)
  total_unique_r <- length(unique_replications)
  
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - CATEGORICAL]: Unique Genotypes Isolated:     ", total_unique_g, "\n")
    cat("[DIAGNOSTIC - CATEGORICAL]: Unique Replications Isolated: ", total_unique_r, "\n")
  }
  
  if (total_unique_g < 2) {
    stop(paste0("INTEGRITY COMPLIANCE FAULT: Genetic variance partitioning requires a minimum of 2 distinct genotypes. Found: ", total_unique_g), call. = FALSE)
  }
  
  if (total_unique_r < 2) {
    stop(paste0("DESIGN INFRASTRUCTURE FAULT: Error component decomposition requires a minimum of 2 replication blocks. Found: ", total_unique_r), call. = FALSE)
  }
  
  if (strict_mode) {
    if (total_unique_g != 40) {
      stop(paste0("STRICT INVENTORY DEVIATION: Expected exactly 40 genotypes (G1 to G40) under strict mode guidelines. Discovered: ", total_unique_g), call. = FALSE)
    }
    if (total_unique_r != 3) {
      stop(paste0("STRICT DESIGN DEVIATION: Expected exactly 3 replication structural blocks under strict mode guidelines. Discovered: ", total_unique_r), call. = FALSE)
    }
  }
  
  # Verify if genotype character expressions contain structural spaces that distort names
  if (is.character(raw_genotypes) || is.factor(raw_genotypes)) {
    string_converted_g <- as.character(raw_genotypes)
    trimmed_converted_g <- trimws(string_converted_g)
    if (any(string_converted_g != trimmed_converted_g)) {
      warning("DATA QUALITY ANOMALY: Hidden leading or trailing spaces detected inside the 'Genotype' factor entries. This can duplicate line names during string sorting.", call. = FALSE)
    }
  }
  
  # =========================================================================
  # BLOCK 6: QUANTITATIVE VECTOR INTEGRITY AND NUMERIC COMPLIANCE
  # =========================================================================
  traits_to_validate <- if (strict_mode) mandatory_traits else mandatory_traits[mandatory_traits %in% active_headers]
  
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - NUMERIC]: Initiating quantitative scanning on phenotypic data columns...\n")
  }
  
  for (trait in traits_to_validate) {
    target_vector <- data[[trait]]
    
    # Absolute numeric verification test
    if (!is.numeric(target_vector)) {
      stop(paste0("\n", rep("!", 85), "\n",
                  "CRITICAL QUANTITATIVE TYPE VIOLATION\n",
                  "The phenotypic trait column vector [ ", trait, " ] is not recognized as a numeric data type.\n",
                  "Current R Object Storage Class: ", paste(class(target_vector), collapse = ", "), "\n",
                  "Downstream matrix models will fail. Remove text annotations or currency signs from values.\n",
                  rep("!", 85), "\n"), call. = FALSE)
    }
    
    # Isolate missing metrics safely
    missing_count <- sum(is.na(target_vector))
    missing_ratio <- missing_count / total_rows
    
    if (missing_count > 0) {
      msg <- paste0("DATA LOSS ALERT: Trait vector '", trait, "' has ", missing_count,
                    " uninitialized observations (NA). Missing proportion: ", round(missing_ratio * 100, 2), "%.")
      if (reporting_level >= 1) warning(msg, call. = FALSE)
      
      # Block if data loss exceeds processing limits (50%)
      if (missing_ratio > 0.50) {
        stop(paste0("CRITICAL DATA LOSS: Trait '", trait, "' has exceeded the 50% data loss threshold. Linear model components cannot evaluate missing matrices."), call. = FALSE)
      }
    }
    
    # =========================================================================
    # BLOCK 7: BIOLOGICAL RANGE AND BOUND INSPECTIONS
    # =========================================================================
    # Isolate non-missing elements for range testing
    clean_numeric_elements <- target_vector[!is.na(target_vector)]
    
    if (length(clean_numeric_elements) > 0) {
      # Screening for negative values
      if (any(clean_numeric_elements < 0)) {
        negative_offsets <- which(target_vector < 0)
        stop(paste0("BIOLOGICAL BOUND EXCEPTION: Negative numeric limits detected in column [ ", trait,
                    " ] at row index offsets: [ ", paste(negative_offsets + 1, collapse = ", "), " ]. ",
                    "Agronomic metrics like plant height or grain mass cannot fall below absolute zero."), call. = FALSE)
      }
      
      # Extreme outlier alert limits (e.g., zero value validations)
      if (any(clean_numeric_elements == 0)) {
        zero_offsets <- which(target_vector == 0)
        warning(paste0("BIOLOGICAL ZERO WARNING: Observation exactly equal to zero isolated in column [ ", trait,
                       " ] at row offsets: [ ", paste(zero_offsets + 1, collapse = ", "), " ]. ",
                       "Verify if this represents a missing record or actual biological expression."), call. = FALSE)
      }
    }
  }
  
  # =========================================================================
  # BLOCK 8: CROSS-TABULATION ORTHOGONALITY ANALYSIS
  # =========================================================================
  if (reporting_level >= 2) {
    cat("[DIAGNOSTIC - ORTHOGONALITY]: Testing design layout cell frequencies...\n")
  }
  
  factor_g <- as.factor(data$Genotype)
  factor_r <- as.factor(data$Replication)
  
  frequency_table <- table(factor_g, factor_r)
  
  # Detect empty cells or missing layout nodes
  if (any(frequency_table == 0)) {
    warning("DESIGN METRIC WARNING: Serious structural imbalance detected. Certain Genotype x Replication blocks have zero records. Downstream models will default to unbalanced configurations.", call. = FALSE)
  }
  
  # Detect replicate duplication cells
  if (any(frequency_table > 1)) {
    warning("DESIGN DEFICIENCY WARNING: Multiple entries discovered for identical Genotype-Replication coordinates. Ensure your dataset does not contain duplicate subsampling records.", call. = FALSE)
  }
  
  # =========================================================================
  # BLOCK 9: REPORT COMPILATION AND EXIT SIGNATURE
  # =========================================================================
  timestamp_complete <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(timestamp_complete, timestamp_init, units = "secs"))
  
  if (reporting_level >= 1) {
    cat("[LOG - COMPLIANCE]: Integrity check completed successfully.\n")
    cat("[LOG - COMPLIANCE]: Total data records verified as stable: ", total_rows, "\n")
    cat("[LOG - COMPLIANCE]: Resource compute time required: ", round(elapsed_seconds, 5), " seconds.\n")
    cat(rep("=", 85), "\n", sep = "")
  }
  
  return(TRUE)
}