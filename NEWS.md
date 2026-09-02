# AgriDataTools 0.2.1

* **Comprehensive Graphical Engine Refactoring (`plot_agri_graphics`)**:
  * Fully restored and optimized all plotting modules (`residual`, `correlation`, `mean`, `pca`, `cluster`, and `path`) with zero data loss.
  * Resolved hardcoding issues, label bugs (such as G-1, G-2 naming conflicts), and color palette limitations across all visualization types.
  * Implemented dynamic text sizing (`cex`) and automated genotype aggregation in PCA biplots to eliminate text overlapping and ensure clean rendering for large datasets.
  * Added unconstrained dynamic color palettes in circular dendrograms to support unlimited cluster groups smoothly.
  * Configured robust graphics state safety guards (`on.exit` and `par` restoration) to prevent graphical device leaks during residual and circular cluster plotting.

* **Package Infrastructure & CRAN Compliance**:
  * Updated version number to `0.2.0` to reflect major functional enhancements and critical bug fixes.
  * Integrated official `BugReports` and GitHub repository metadata links into the `DESCRIPTION` file.
  * Enhanced internal error-handling and input validation routines for all analytical payloads across statistical modules.

# AgriDataTools 0.1.2

* Initial structural deployment of the core infrastructure.
* Integrated complete analytical engines for CRD and RCBD Analysis of Variance (ANOVA).
* Added multiple group comparison post-hoc matrices (LSD, Tukey, and Scheffe testing blocks).
* Integrated core statistical functions including Summary Statistics, Correlation Analysis, Path Coefficients, and Principal Component Analysis (PCA).
* Finalized multivariate clustering modules (Hierarchical Agglomerative Clustering) and multi-trait plotting systems.