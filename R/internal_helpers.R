#' Internal Package Helpers and Matrix Formatting Utilities
#'
#' @description
#' The functions documented here are internal utility structures meant to handle mathematical transformation,
#' table cleanups, and defensive checks across multiple high-level analysis engines within AgriDataTools.
#'
#' @details
#' These helper tools are not exported to the user environment to maintain a clean package namespace.
#' They provide centralized, low-level string processing, matrix layout adjustments, and basic scalar transformations.
#'
#' @param input_table A raw data frame layout containing numerical extraction columns.
#' @param decimal_places Integer value specifying the rounding precision. Defaults to \code{4}.
#'
#' @return Transformed vectors or data frames adjusted to system presentation requirements.
#'
#' @author Faheem Khan (\email{2022ag94@@uaf.edu.pk})
#'
#' @keywords internal
#' @noRd
.clean_output_matrix <- function(input_table, decimal_places = 4) {
  
  # Safeguard matrix input structural properties
  if (!is.data.frame(input_table)) {
    return(input_table)
  }
  
  # Loop through columns and safely round numeric representations
  for (col_name in colnames(input_table)) {
    if (is.numeric(input_table[[col_name]])) {
      input_table[[col_name]] <- round(input_table[[col_name]], decimal_places)
    }
  }
  
  return(input_table)
}

#' Generate Significance Stars for Probability Vectors
#'
#' @param p_value_vector A numeric vector of calculated probability values.
#'
#' @return A character vector containing statistical significance indicators.
#'
#' @keywords internal
#' @noRd
.assign_significance_stars <- function(p_value_vector) {
  
  if (is.null(p_value_vector) || !is.numeric(p_value_vector)) {
    return(character(0))
  }
  
  stars_vector <- vapply(p_value_vector, function(p) {
    if (is.na(p)) return("")
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    if (p <= 0.1)   return(".")
    return("ns")
  }, FUN.VALUE = character(1))
  
  return(stars_vector)
}

#' Validate General Matrix Dimensions and Uniform Replications
#'
#' @param data Input data frame layout.
#' @param expected_cols Character vector of mandatory column names.
#'
#' @return Logical \code{TRUE} if layout requirements are fulfilled.
#'
#' @keywords internal
#' @noRd
.verify_orthogonal_dimensions <- function(data, expected_cols) {
  
  # Assert structural column presence
  missing_fields <- setdiff(expected_cols, colnames(data))
  if (length(missing_fields) > 0) {
    stop(paste0("CRITICAL ENGINE ERROR: Internal structure verification failed. Missing columns: ",
                paste(missing_fields, collapse = ", ")), call. = FALSE)
  }
  
  return(TRUE)
}

# =========================================================================
# END OF FILE: R/internal_helpers.R
# =========================================================================