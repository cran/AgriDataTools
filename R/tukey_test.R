#' Tukey's Honestly Significant Difference (HSD) Post-Hoc Mean Comparison Engine
#'
#' @description
#' The \code{compute_tukey} function executes a high-precision pairwise post-hoc mean separation
#' analysis using Tukey's Honestly Significant Difference framework. It calculates family-wise error
#' control thresholds, pairwise differences, and outputs comprehensive ranking tables with group letters.
#'
#' @param data A verified \code{data.frame} containing treatment lines/cultivars and the target phenotypic trait.
#' @param trait A single character string specifying the column name of the target trait.
#' @param anova_results A structured \code{list} derived from upstream ANOVA layouts containing an \code{anova_table}.
#' @param total_replications An integer specifying the absolute number of replication blocks (\eqn{r}).
#' @param geno_col A character string specifying the treatment/genotype column name. Defaults to auto-detecting \code{"Genotype"}.
#' @param alpha A numeric value defining the adjusted family-wise error rate threshold. Defaults to \code{0.05}.
#' @param reporting_level An integer vector flag defining console trace settings: \code{0} for silent execution, \code{1} for summary, and \code{2} for exhaustive tracking. Defaults to \code{2}.
#'
#' @return A structured named \code{list} of class \code{"list"} containing 4 computational components:
#' \item{tukey_value}{The absolute calculated scalar value of Tukey's Honestly Significant Difference threshold.}
#' \item{se_mean}{The isolated Standard Error of the treatment mean scalar (\eqn{SE_{\bar{y}}}).}
#' \item{comparison_matrix}{A detailed data frame containing pairwise differences, q-statistics, p-values, and significance markers.}
#' \item{ranked_means}{A structured data frame containing sorted treatment means and significance group letters (\code{Tukey_Letters}).}
#'
#' @importFrom stats qtukey ptukey
#' @export
#'
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' reps <- length(unique(gv_data$Replication))
#' rcbd_results <- anova_rcbd(data = gv_data, trait = "PH", reporting_level = 0)
#' 
#' tukey_output <- compute_tukey(
#'   data = gv_data,
#'   trait = "PH",
#'   anova_results = rcbd_results,
#'   total_replications = reps
#' )
compute_tukey <- function(data, trait, anova_results, total_replications, geno_col = NULL, alpha = 0.05, reporting_level = 2) {
  
  if (missing(data) || missing(trait) || missing(anova_results) || missing(total_replications)) {
    stop("CRITICAL PAIRWISE FAULT: Missing core arguments ('data', 'trait', 'anova_results', or 'total_replications').", call. = FALSE)
  }
  
  # Dynamic Genotype Column Resolution
  if (is.null(geno_col)) {
    possible_cols <- c("Genotype", "genotype", "Genotypes", "genotypes", "Line", "line", "Cultivar", "cultivar", "Variety", "variety", "Entry", "entry")
    found_col <- intersect(possible_cols, colnames(data))
    if (length(found_col) > 0) {
      geno_col <- found_col[1]
    } else {
      stop("COLUMN FAULT: Genotype/Line column not detected automatically. Please specify 'geno_col'.", call. = FALSE)
    }
  }
  
  if (!geno_col %in% colnames(data)) {
    stop(paste0("COLUMN FAULT: Specified genotype column '", geno_col, "' not found in input data."), call. = FALSE)
  }
  
  if (!trait %in% colnames(data)) {
    stop(paste0("TRAIT FAULT: Specified trait '", trait, "' not found in input data."), call. = FALSE)
  }
  
  # Clean Data
  clean_data <- data[!is.na(data[[trait]]) & !is.na(data[[geno_col]]), , drop = FALSE]
  clean_data[[geno_col]] <- as.factor(clean_data[[geno_col]])
  response_vector  <- as.numeric(clean_data[[trait]])
  
  # Robust Extraction of Error MS and Df from ANOVA Table
  anova_table <- anova_results$anova_table
  error_row <- anova_table[grep("Error|Residual", anova_table$Source, ignore.case = TRUE), ]
  
  if (nrow(error_row) == 0) {
    stop("ANOVA FAULT: Could not locate 'Error' or 'Residual' row in the provided anova_table.", call. = FALSE)
  }
  
  df_error <- as.numeric(error_row$Df[1])
  ems_val  <- as.numeric(error_row$MS[1])
  
  if (is.na(ems_val) || is.na(df_error) || df_error <= 0 || ems_val <= 0) {
    stop("COMPUTATIONAL FAULT: Invalid Mean Square Error or Degrees of Freedom retrieved from ANOVA.", call. = FALSE)
  }
  
  genotype_levels <- levels(clean_data[[geno_col]])
  num_genotypes   <- length(genotype_levels)
  
  # Tukey's HSD Calculations
  se_mean_val   <- sqrt(ems_val / total_replications)
  q_critical    <- stats::qtukey(1 - alpha, nmeans = num_genotypes, df = df_error)
  tukey_hsd_val <- q_critical * se_mean_val
  
  # Mean Computation
  mean_records <- tapply(response_vector, clean_data[[geno_col]], mean, na.rm = TRUE)
  ranked_df    <- data.frame(
    Genotype = names(mean_records),
    Mean     = as.numeric(mean_records),
    stringsAsFactors = FALSE
  )
  colnames(ranked_df)[1] <- geno_col
  
  # Sort Means Descending
  ranked_df            <- ranked_df[order(-ranked_df$Mean), ]
  rownames(ranked_df) <- NULL
  
  # Pairwise Difference Matrix Generation
  n <- nrow(ranked_df)
  comp_list <- list()
  comp_counter <- 1
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      g1 <- ranked_df[[geno_col]][i]
      g2 <- ranked_df[[geno_col]][j]
      diff_val <- ranked_df$Mean[i] - ranked_df$Mean[j]
      
      q_calc   <- abs(diff_val) / se_mean_val
      p_val    <- 1 - stats::ptukey(q_calc, nmeans = num_genotypes, df = df_error)
      sig_flag <- if (abs(diff_val) >= tukey_hsd_val) "*" else "ns"
      
      comp_list[[comp_counter]] <- data.frame(
        Comparison   = paste(g1, "vs", g2),
        Difference   = diff_val,
        q_Calculated = q_calc,
        p_value      = p_val,
        Significance = sig_flag,
        stringsAsFactors = FALSE
      )
      comp_counter <- comp_counter + 1
    }
  }
  comparison_matrix <- do.call(rbind, comp_list)
  
  # Group Lettering Algorithm
  M <- matrix(FALSE, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (abs(ranked_df$Mean[i] - ranked_df$Mean[j]) <= tukey_hsd_val) {
        M[i, j] <- TRUE
      }
    }
  }
  
  groups <- list()
  for (i in 1:n) {
    current_group <- i
    for (j in 1:n) {
      if (i != j && all(M[j, current_group])) {
        current_group <- c(current_group, j)
      }
    }
    groups[[i]] <- sort(unique(current_group))
  }
  
  unique_groups <- list()
  for (g in groups) {
    if (!any(sapply(unique_groups, function(ug) all(g %in% ug)))) {
      unique_groups <- unique_groups[!sapply(unique_groups, function(ug) all(ug %in% g))]
      unique_groups[[length(unique_groups) + 1]] <- g
    }
  }
  
  unique_groups  <- unique_groups[order(sapply(unique_groups, function(x) min(x)))]
  letters_vector <- rep("", n)
  
  for (g_idx in seq_along(unique_groups)) {
    current_letter <- letters[g_idx]
    if (g_idx > 26) current_letter <- paste0("z", g_idx - 26)
    for (member in unique_groups[[g_idx]]) {
      letters_vector[member] <- paste0(letters_vector[member], current_letter)
    }
  }
  
  ranked_df$Tukey_Letters <- letters_vector
  
  # Console Reporting
  if (reporting_level >= 1) {
    cat("\n", rep("=", 90), "\n", sep = "")
    cat("            TUKEY'S HONESTLY SIGNIFICANT DIFFERENCE (HSD) TEST\n")
    cat("            Target Trait: ", trait, " | Significance Level (Alpha): ", alpha, "\n", sep = "")
    cat(rep("=", 90), "\n", sep = "")
    
    cat("\n   TUKEY'S HSD CRITICAL METRICS         :  Estimates \n")
    cat("   * Error Mean Square (EMS)            : ", sprintf("%.5f", ems_val), "\n")
    cat("   * Error Degrees of Freedom (df)      : ", df_error, "\n")
    cat("   * Number of Treatments/Genotypes     : ", num_genotypes, "\n")
    cat("   * Standard Error of Mean (SE_mean)   : ", sprintf("%.5f", se_mean_val), "\n")
    cat("   * Critical q-value (k=", num_genotypes, ", df=", df_error, ")     :  ", sprintf("%.5f", q_critical), "\n", sep = "")
    cat("   * Tukey HSD Value (Alpha = ", alpha, ")     :  ", sprintf("%.5f", tukey_hsd_val), "\n", sep = "")
    
    cat("\n", rep("-", 90), "\n", sep = "")
    cat(" RANKED MEANS AND STATISTICAL SIGNIFICANCE GROUPS:\n")
    cat(rep("-", 90), "\n", sep = "")
    print(ranked_df, row.names = FALSE)
    cat(rep("-", 90), "\n", sep = "")
    cat(" Note: Means sharing the same letter are not significantly different at p <= ", alpha, ".\n", sep = "")
    cat(rep("=", 90), "\n\n", sep = "")
  }
  
  return(invisible(list(
    tukey_value       = tukey_hsd_val,
    se_mean           = se_mean_val,
    comparison_matrix = comparison_matrix,
    ranked_means      = ranked_df
  )))
}