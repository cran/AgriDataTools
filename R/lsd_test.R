#' Fisher's Least Significant Difference (LSD) Post-Hoc Mean Comparison Engine
#'
#' @description
#' The `compute_lsd` function executes a rigorous, high-precision pairwise post-hoc mean separation
#' analysis using Fisher's Least Significant Difference protocol. It isolates critical differences,
#' evaluates pairwise significance metrics, and outputs comprehensive ranking tables with group letters.
#'
#' @param data A verified \code{data.frame} containing the columns \code{Genotype} and the
#'    phenotypic trait under investigation.
#' @param trait A single character string specifying the column name of the target trait.
#' @param anova_results A structured \code{list} derived from upstream ANOVA layouts containing an \code{anova_table}.
#' @param total_replications An integer specifying the absolute number of replication blocks (\eqn{r}).
#' @param alpha A numeric value defining the Type-I error rate probability threshold. Defaults to \code{0.05}.
#' @param reporting_level An integer vector flag defining console trace settings: \code{0} for silent execution, \code{1} for summary, and \code{2} for exhaustive tracking. Defaults to \code{2}.
#'
#' @return A structured named \code{list} of class \code{"list"} containing 4 computational components:
#' \item{lsd_value}{A numeric scalar representing the absolute calculated value of Fisher's Least Significant Difference at the specified alpha level.}
#' \item{sed}{A numeric scalar representing the isolated Standard Error of Difference (SED) between two treatment means.}
#' \item{comparison_matrix}{A \code{data.frame} or \code{NULL} layout reserved for detailed pairwise lines differences.}
#' \item{ranked_means}{A \code{data.frame} containing lines/cultivars sorted by mean performance alongside their assigned statistical significance group letters (\code{LSD_Letters}).}
#'
#' @importFrom stats qt
#' @export
#'
#' @examples
#' # Load integrated data matrices
#' data(gv_data, package = "AgriDataTools")
#' 
#' # Total replications matching the experimental design blocks
#' reps <- length(unique(gv_data$Replication))
#' 
#' # Generating upstream ANOVA structure from gv_data for the target trait (PH)
#' model_fit <- aov(PH ~ Genotype + Replication, data = gv_data)
#' anova_summary <- summary(model_fit)[[1]]
#' 
#' mock_anova <- list(
#'   anova_table = data.frame(
#'     Source = c("Genotype", "Replication", "Error"),
#'     Df = anova_summary$Df,
#'     MS = anova_summary$`Mean Sq`,
#'     stringsAsFactors = FALSE
#'   )
#' )
#' 
#' # Running post-hoc separation sweeps with clean letters group layout
#' lsd_output <- compute_lsd(
#'   data = gv_data, 
#'   trait = "PH", 
#'   anova_results = mock_anova, 
#'   total_replications = reps, 
#'   alpha = 0.05
#' )
#' print(lsd_output$ranked_means)
compute_lsd <- function(data, trait, anova_results, total_replications, alpha = 0.05, reporting_level = 2) {
  timestamp_start <- Sys.time()
  
  if (missing(data) || missing(trait) || missing(anova_results) || missing(total_replications)) {
    stop("CRITICAL PAIRWISE FAULT: Missing core arguments.", call. = FALSE)
  }
  
  data$Genotype   <- as.factor(data$Genotype)
  response_vector <- as.numeric(data[[trait]])
  
  anova_matrix <- anova_results$anova_table
  error_row_idx <- grep("Error", anova_matrix$Source)
  df_error <- as.numeric(anova_matrix$Df[error_row_idx])
  ems_val  <- as.numeric(anova_matrix$MS[error_row_idx])
  
  sed_val        <- sqrt((2 * ems_val) / total_replications)
  t_critical     <- stats::qt(1 - (alpha / 2), df = df_error)
  lsd_scalar_val <- t_critical * sed_val
  
  genotype_levels <- levels(data$Genotype)
  mean_records <- sapply(genotype_levels, function(g) mean(response_vector[data$Genotype == g], na.rm = TRUE))
  
  ranked_df <- data.frame(Genotype = names(mean_records), Mean = as.numeric(mean_records), stringsAsFactors = FALSE)
  ranked_df <- ranked_df[order(-ranked_df$Mean), ]
  rownames(ranked_df) <- NULL
  
  # --- EXPERT GROUP LETTERING CALCULATOR ENGINE ---
  n <- nrow(ranked_df)
  M <- matrix(FALSE, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(abs(ranked_df$Mean[i] - ranked_df$Mean[j]) <= lsd_scalar_val) M[i, j] <- TRUE
    }
  }
  groups <- list()
  for(i in 1:n) {
    current_group <- i
    for(j in 1:n) {
      if(i != j && all(M[j, current_group])) current_group <- c(current_group, j)
    }
    groups[[i]] <- sort(unique(current_group))
  }
  unique_groups <- list()
  for(g in groups) {
    if(!any(sapply(unique_groups, function(ug) all(g %in% ug)))) {
      unique_groups <- unique_groups[!sapply(unique_groups, function(ug) all(ug %in% g))]
      unique_groups[[length(unique_groups) + 1]] <- g
    }
  }
  unique_groups <- unique_groups[order(sapply(unique_groups, function(x) min(x)))]
  letters_vector <- rep("", n)
  for(g_idx in seq_along(unique_groups)) {
    current_letter <- letters[g_idx]
    if(g_idx > 26) current_letter <- paste0("z", g_idx - 26)
    for(member in unique_groups[[g_idx]]) letters_vector[member] <- paste0(letters_vector[member], current_letter)
  }
  ranked_df$LSD_Letters <- letters_vector
  
  return(list(lsd_value = lsd_scalar_val, sed = sed_val, comparison_matrix = NULL, ranked_means = ranked_df))
}