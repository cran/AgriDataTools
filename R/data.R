#' Genotypic Variability and Agricultural Research Dataset
#'
#' Evaluated performance parameters across multiple line and cultivar iterations 
#' under a randomized complete block layout.
#'
#' @docType data
#' @name gv_data
#' @keywords datasets
#'
#' @format A \code{data.frame} containing phenotypic records with wheat morphological traits:
#' \describe{
#'   \item{Genotype}{Factor or character vector identifying evaluated breeding germplasm lines or cultivars.}
#'   \item{Replication}{Factor or integer vector indicating experimental replications or blocks within the layout.}
#'   \item{PH}{Plant Height measured in centimeters (cm).}
#'   \item{PL}{Peduncle Length measured in centimeters (cm).}
#'   \item{SL}{Spike Length measured in centimeters (cm).}
#'   \item{NOT}{Number of tillers per plant.}
#'   \item{NOSS}{Number of spikelets per spike.}
#'   \item{TGW}{Thousand grain weight measured in grams (g).}
#'   \item{GYPM}{Grain yield per meter measured in grams (g).}
#' }
#' 
#' @source Experimental plant breeding field trial records.
#' 
#' @usage data(gv_data)
#' 
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' head(gv_data)
"gv_data"