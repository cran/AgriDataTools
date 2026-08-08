#' Tukey's Honestly Significant Difference (HSD) Post-Hoc Mean Comparison Engine
#'
#' @description
#' The `compute_tukey` function executes a high-precision pairwise post-hoc mean separation
#' analysis using Tukey's Honestly Significant Difference framework.
#'
#' @param data A verified \code{data.frame} containing the columns \code{Genotype} and the
#'    phenotypic trait under evaluation.
#' @param trait A single character string specifying the column name of the target trait.
#' @param anova_results A structured \code{list} derived from upstream ANOVA layouts.
#' @param total_replications An integer specifying the absolute number of replication blocks (\eqn{r}).
#' @param alpha A numeric value defining the adjusted family-wise error rate threshold. Defaults to \code{0.05}.
#' @param reporting_level An integer vector flag defining console trace settings. Defaults to \code{2}.
#'
#' @return A structured named \code{list} of class \code{"list"} containing 4 computational components:
#' \item{tukey_value}{The absolute calculated scalar value of Tukey's Honestly Significant Difference threshold.}
#' \item{se_mean}{The isolated Standard Error of the treatment mean scalar (\eqn{SE_{\bar{y}}}).}
#' \item{comparison_matrix}{A detailed data frame or \code{NULL} layout reserved for pairwise lines differences.}
#' \item{ranked_means}{A structured data frame containing sorted treatment means and significance group letters (\code{Tukey_Letters}).}
#'
#' @author Faheem Khan (\email{2022ag94@@uaf.edu.pk})
#'
#' @importFrom stats qtukey
#' @export
#'
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' reps <- length(unique(gv_data$Replication))
#' 
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
#' tukey_output <- compute_tukey(
#'   data = gv_data,
#'   trait = "PH",
#'   anova_results = mock_anova,
#'   total_replications = reps
#' )
#' print(tukey_output$ranked_means)
compute_tukey <- function(data, trait, anova_results, total_replications, alpha = 0.05, reporting_level = 2) {
  if (missing(data) || missing(trait) || missing(anova_results) || missing(total_replications)) {
    stop("CRITICAL PAIRWISE FAULT: Missing core arguments.", call. = FALSE)
  }
  
  data$Genotype   <- as.factor(data$Genotype)
  response_vector <- as.numeric(data[[trait]])
  
  anova_matrix <- anova_results$anova_table
  error_row_idx <- grep("Error", anova_matrix$Source)
  df_error <- as.numeric(anova_matrix$Df[error_row_idx])
  ems_val  <- as.numeric(anova_matrix$MS[error_row_idx])
  
  genotype_levels <- levels(data$Genotype)
  num_genotypes   <- length(genotype_levels)
  
  se_mean_val    <- sqrt(ems_val / total_replications)
  q_critical     <- stats::qtukey(1 - alpha, nmeans = num_genotypes, df = df_error)
  tukey_hsd_val  <- q_critical * se_mean_val
  
  mean_records <- sapply(genotype_levels, function(g) mean(response_vector[data$Genotype == g], na.rm = TRUE))
  ranked_df <- data.frame(Genotype = names(mean_records), Mean = as.numeric(mean_records), stringsAsFactors = FALSE)
  ranked_df <- ranked_df[order(-ranked_df$Mean), ]
  rownames(ranked_df) <- NULL
  
  # --- EXPERT GROUP LETTERING CALCULATOR ENGINE ---
  n <- nrow(ranked_df)
  M <- matrix(FALSE, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(abs(ranked_df$Mean[i] - ranked_df$Mean[j]) <= tukey_hsd_val) M[i, j] <- TRUE
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
  ranked_df$Tukey_Letters <- letters_vector
  
  return(list(tukey_value = tukey_hsd_val, se_mean = se_mean_val, comparison_matrix = NULL, ranked_means = ranked_df))
}