#' AgriDataTools: Automated Statistical Analysis and Tools for Agricultural Research
#'
#' @description
#' A comprehensive, high-precision biometrical computing toolkit engineered specifically for plant breeding,
#' agronomic trial evaluations, and quantitative genetic research. Provides end-to-end processing pipelines
#' for completely randomized designs (CRD), randomized complete block designs (RCBD), and analysis of covariance (ANCOVA),
#' including ANOVA, variance component partitioning (Vg, Vp, Ve), broad-sense heritability (H2), genetic advance (GA),
#' and post-hoc pairwise mean separation tests (Tukey's HSD, LSD).
#'
#' The package provides advanced quantitative tools including:
#' \itemize{
#'   \item \strong{Descriptive Statistical Summaries:} Automated calculation of overall grand mean, range (minimum to maximum limits), standard deviation, standard error, and phenotypic Coefficient of Variation (CV%).
#'   \item \strong{Analysis of Covariance (ANCOVA):} Advanced experimental error control and adjusted mean estimations across treatment layouts.
#'   \item \strong{Genetic Variability Analysis:} Comprehensive estimation of Phenotypic Coefficient of Variation (PCV), Genotypic Coefficient of Variation (GCV), Heritability (H2), and Genetic Advance (GA).
#'   \item \strong{Multi-Level Correlation Analysis:} Partitioning of association matrices into Genotypic, Phenotypic, and Environmental correlation components.
#'   \item \strong{Path Coefficient Analysis:} Dissection of simple correlation vectors into direct and indirect contribution effects across both Genotypic and Phenotypic levels.
#'   \item \strong{Principal Component Analysis (PCA):} Multivariate dimensionality reduction, eigenvalue extraction, and variance proportion profiling across phenotypic traits.
#'   \item \strong{Cluster Analysis:} Hierarchical agglomerative and k-means clustering strategies for genetic diversity pattern classification.
#' }
#'
#' @author Faheem Khan (\email{2022ag94@@uaf.edu.pk})
#'
#' @docType package
#' @name AgriDataTools-package
#' @aliases AgriDataTools
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats lm aov anova as.formula cor pt qt qf pf var sd qtukey ptukey prcomp dist hclust cutree kmeans aggregate qqnorm qqline
#' @importFrom graphics abline arrows axis image par text plot
#' @importFrom grDevices heat.colors
#' @importFrom utils head packageVersion
## usethis namespace: end
NULL