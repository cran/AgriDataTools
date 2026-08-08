#' Advanced High-Precision Publication-Ready Graphics Suite
#'
#' @description
#' The `plot_agri_graphics` function serves as the unified visualization hub for
#' the AgriDataTools package. It handles basic statistical diagnostics
#' (residuals, correlations) alongside modern publication-grade graphical
#' representations for mean performance, PCA space, hierarchical dendrograms, 
#' and path analysis direct effects.
#'
#' @import ggplot2
#' @import dplyr
#' @import factoextra
#' @importFrom dendextend color_branches circlize_dendrogram
#' @importFrom circlize circos.clear circos.initialize
#' @importFrom stats as.dendrogram dist hclust reorder aov lm residuals fitted qqnorm qqline aggregate
#' @importFrom graphics par abline plot title
#' @importFrom utils globalVariables
#' @export
#' 
#' @details
#' Visualizing high-dimensional screening metrics across diverse lines or
#' cultivars requires balancing diagnostic model validation checks with
#' advanced multivariate aesthetics. This engine supports base diagnostic
#' rendering as well as optimized ggplot2 geometries featuring dynamic color
#' palettes, non-overlapping labels, and geometric layout vector mapping fields.
#'
#' @param type A single character string specifying the target chart module:
#'         \code{"residual"}, \code{"correlation"}, \code{"mean"}, \code{"pca"},
#'         \code{"cluster"}, or \code{"path"}.
#' @param payload A structured analysis \code{list} derived from computational
#'         engines (e.g., \code{compute_lsd}, \code{analyze_pca},
#'         \code{analyze_clustering}, \code{compute_path_analysis}).
#' @param trait_name A character string defining the target phenotypic trait
#'         title label. Used primarily in \code{"mean"} and \code{"path"} layouts.
#' @param reporting_level An integer vector flag defining console trace settings:
#'         \code{0} for silent, \code{1} for structural updates, and \code{2} for
#'         exhaustive analytical tracing. Defaults to \code{2}.
#' @param num_clusters An integer specifying the number of cluster groups to color 
#'         in the circular dendrogram module. Defaults to \code{4}.
#'
#' @return Invisibly returns a logical scalar \code{TRUE} upon successful execution.
#'    This function is primarily invoked for its side effect of rendering publication-grade
#'    graphical plots (e.g., residual diagnostic plots, correlation heatmaps, mean performance
#'    barcharts, PCA biplots, circular dendrograms, or path analysis plots) to the active 
#'    graphics device.
#'
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
#' 
#' # Define custom mapping or number of clusters beforehand
#' k_groups <- 4
#'
#' # 1. Mean performance: Genotypic performance with LSD
#' reps <- length(unique(gv_data$Replication))
#' fit <- aov(PH ~ Genotype + Replication, data = gv_data)
#' m_anova <- list(anova_table = data.frame(
#'    Source = c("Genotype", "Replication", "Error"),
#'    Df = summary(fit)[[1]]$Df,
#'    MS = summary(fit)[[1]][[3]]
#' ))
#' lsd_res <- compute_lsd(gv_data, "PH", m_anova, reps)
#' plot_agri_graphics(type = "mean", payload = lsd_res,
#'                    trait_name = "Plant Height")
#'
#' # 2. PCA: Multivariate variation
#' custom_traits_map <- c(
#'      "PH"    = "Plant Height",
#'      "SL"    = "Spike Length",
#'      "PL"    = "Peduncle Length",
#'      "NOT"   = "Number of Tillers",
#'      "NOSS"  = "Number of Spikelets per Spike",
#'      "TGW"   = "Thousand Grain Weight",
#'      "GYPM"  = "Grain Yield per Meter"
#' )
#' pca_res <- analyze_pca(data = gv_data, traits = traits, trait_lookup = custom_traits_map)
#' plot_agri_graphics(type = "pca", payload = pca_res,
#'                    trait_name = "PCA Plot")
#'
#' # 3. Clustering: Dendrogram with flexible cluster parameter option
#' cl_res <- analyze_clustering(data = gv_data, traits = traits, k = k_groups)
#' plot_agri_graphics(type = "cluster", payload = cl_res,
#'                    trait_name = "Clustering", num_clusters = k_groups)
#'
#' # 4. Residuals: Diagnostic plots
#' fit <- lm(PH ~ Genotype, data = gv_data)
#' res_pl <- list(residuals = residuals(fit),
#'                fitted_values = fitted(fit))
#' plot_agri_graphics(type = "residual", payload = res_pl,
#'                    trait_name = "Residuals")
#'
#' # 5. Correlations: Phenotypic matrix
#' cor_m <- cor(gv_data[, traits], use = "pairwise.complete.obs")
#' plot_agri_graphics(type = "correlation",
#'                    payload = list(correlation_matrix = cor_m),
#'                    trait_name = "Correlation")
#'
#' # 6. Path Analysis: Direct Effects Plot for Grain Yield per Meter (GYPM)
#' corr_res <- compute_correlation(data = gv_data, traits = traits, reporting_level = 0)
#' all_predictors <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW")
#' path_res <- compute_path_analysis(
#'    correlation_payload = corr_res, 
#'    response_trait = "GYPM", 
#'    predictor_traits = all_predictors,
#'    reporting_level = 0
#' )
#' plot_agri_graphics(type = "path", payload = path_res, trait_name = "Grain Yield per Meter (GYPM)")
#' 
plot_agri_graphics <- function(type, payload, trait_name = "Target Character Matrix", reporting_level = 2, num_clusters = 4) {
  
  # =========================================================================
  # BLOCK 1: STARTUP PARAMETERS AND SYSTEM DEFENSE ROUTINES
  # =========================================================================
  timestamp_start <- Sys.time()
  
  Var1 <- Var2 <- value <- Effect <- Impact <- Trait <- Mean <- Genotype <- NULL
  
  if (missing(type) || missing(payload)) {
    stop("CRITICAL INPUT FAULT: Missing arguments. Provide both visualization type and payload.", call. = FALSE)
  }
  
  chart_type <- match.arg(tolower(trimws(type)), c("residual", "correlation", "mean", "pca", "cluster", "path"))
  
  if (reporting_level >= 1) {
    cat(rep("=", 85), "\n", sep = "")
    cat("AGRIDATATOOLS PACKAGED ENGINE v0.1.0 - INTEGRATED GRAPHICS VISUALIZATION\n")
    cat("Plot Generation Inception: ", as.character(timestamp_start), "\n")
    cat(rep("-", 85), "\n", sep = "")
  }
  
  # --- Theme Configuration ---
  agri_custom_theme <- function() {
    ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15, color = "#1a252f", hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.major = ggplot2::element_line(color = "gray90", linetype = "dashed")
      )
  }
  
  # =========================================================================
  # MODULES EXECUTION CHAIN (Pure Plotting Rendering)
  # =========================================================================
  
  if (chart_type == "correlation") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("Package 'reshape2' required.")
    corr_mat <- payload$correlation_matrix
    corr_mat[upper.tri(corr_mat)] <- NA
    melted_cormat <- reshape2::melt(corr_mat, na.rm = TRUE)
    g <- ggplot2::ggplot(melted_cormat, ggplot2::aes(Var2, Var1, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limit = c(-1, 1)) +
      agri_custom_theme() + ggplot2::labs(title = "Pearson Correlation Matrix Heatmap", x = "", y = "", fill = "r Value")
    print(g)
    
  } else if (chart_type == "residual") {
    resids <- payload$residuals
    fitted_vals <- payload$fitted_values
    
    # CRAN GUARD: Save and restore graphics settings immediately
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    
    graphics::par(mfrow = c(1, 2))
    graphics::plot(fitted_vals, resids, main = "Residuals vs Fitted", pch = 19, col = "darkgreen")
    graphics::abline(h = 0, lty = 2, col = "red")
    stats::qqnorm(resids, main = "Normal Q-Q Distribution", pch = 19, col = "darkgreen")
    stats::qqline(resids, col = "red")
    
  } else if (chart_type == "mean") {
    target_df <- payload$ranked_means
    g <- ggplot2::ggplot(target_df, ggplot2::aes(x = reorder(Genotype, -Mean), y = Mean, fill = Mean)) +
      ggplot2::geom_bar(stat = "identity", color = "black", width = 0.7) + 
      ggplot2::scale_fill_viridis_c(option = "mako") +
      agri_custom_theme() +
      ggplot2::labs(title = paste("Mean Performance:", trait_name), x = "Genotype / Cultivar", y = "Mean Value")
    print(g)
    
  } else if (chart_type == "pca") {
    if (!requireNamespace("factoextra", quietly = TRUE)) stop("Package 'factoextra' required.")
    
    pca_obj <- if (is.list(payload) && "pca_object" %in% names(payload)) payload$pca_object else payload
    res.pca <- pca_obj
    
    # Dynamic Aggregation for prcomp / PCA scores to ensure unique genotypes
    if (!is.null(res.pca$x)) {
      scores_df <- as.data.frame(res.pca$x)
      orig_rows <- nrow(scores_df)
      
      if (orig_rows > 1) {
        orig_rownames <- rownames(scores_df)
        clean_names <- sub("(_[Rr]ep?[0-9]+|\\.[0-9]+|_R[0-9]+)$", "", orig_rownames)
        scores_df$Genotype_Group <- clean_names
        
        if (length(unique(scores_df$Genotype_Group)) == orig_rows && orig_rows >= 3) {
          inferred_genos <- ceiling(orig_rows / 3)
          scores_df$Genotype_Group <- rep(1:inferred_genos, each = 3, length.out = orig_rows)
        }
        
        if (orig_rows > length(unique(scores_df$Genotype_Group))) {
          pc_cols <- grep("^PC|^Dim", colnames(scores_df), value = TRUE)
          if (length(pc_cols) == 0) pc_cols <- colnames(scores_df)
          
          agg_scores <- stats::aggregate(scores_df[, pc_cols, drop = FALSE], 
                                         by = list(Genotype = scores_df$Genotype_Group), 
                                         FUN = mean, na.rm = TRUE)
          
          rownames(agg_scores) <- paste0("G_", agg_scores$Genotype)
          agg_scores$Genotype <- NULL
          res.pca$x <- as.matrix(agg_scores)
        }
      }
    }
    
    var_exp_1 <- round((res.pca$sdev[1]^2 / sum(res.pca$sdev^2)) * 100, 1)
    var_exp_2 <- round((res.pca$sdev[2]^2 / sum(res.pca$sdev^2)) * 100, 1)
    
    p <- factoextra::fviz_pca_biplot(res.pca,
                                     repel = TRUE,
                                     col.ind = "steelblue", col.var = "darkred",
                                     arrowsize = 1.2, geom.ind = c("point", "text"),
                                     labelsize = 4, pointsize = 3, ggtheme = ggplot2::theme_minimal()) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray85", linetype = "dashed"),
                     axis.line = ggplot2::element_line(color = "black", linewidth = 1),
                     plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(face = "bold", size = 12)) +
      ggplot2::labs(title = "Principal Component Analysis (GxT Biplot)",
                    subtitle = "Multivariate Trait Interrelationships (Genotype Means)",
                    x = paste("PC1 (", var_exp_1, "%)"),
                    y = paste("PC2 (", var_exp_2, "%)"))
    print(p)
    
  } else if (chart_type == "cluster") {
    hc_object <- if (is.list(payload) && "hc_object" %in% names(payload)) payload$hc_object else payload
    if (is.null(hc_object) && inherits(payload, "hclust")) hc_object <- payload
    if (is.null(hc_object)) stop("Cluster payload must contain a valid 'hclust' object.", call. = FALSE)
    
    dend <- stats::as.dendrogram(hc_object)
    
    cluster_colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")
    k_use <- min(num_clusters, length(cluster_colors))
    
    dend <- dendextend::color_branches(dend, k = k_use, col = cluster_colors[1:k_use])
    
    # CRAN GUARD: Reset par and circos graphics state safely
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
      graphics::par(old_par)
      circlize::circos.clear()
    }, add = TRUE)
    
    graphics::par(mar = c(2, 2, 4, 2))
    
    circlize::circos.initialize(factors = rep("a", length(dend)), xlim = c(0, 1))
    dendextend::circlize_dendrogram(dend, labels_cex = 0.75, dend_track_height = 0.55)
    
    graphics::title(
      main = paste("Circular Dendrogram:", trait_name),
      sub = paste0("Divided into ", k_use, " distinct clusters"),
      col.main = "#1a252f", font.main = 2, cex.main = 1.3
    )
    circlize::circos.clear()
    
  } else if (chart_type == "path") {
    direct_eff <- NULL
    
    if (is.list(payload)) {
      if ("Phenotypic" %in% names(payload) && is.list(payload$Phenotypic) && "direct_effects" %in% names(payload$Phenotypic)) {
        direct_eff <- payload$Phenotypic$direct_effects
      } else if ("Genotypic" %in% names(payload) && is.list(payload$Genotypic) && "direct_effects" %in% names(payload$Genotypic)) {
        direct_eff <- payload$Genotypic$direct_effects
      } else if ("direct_effects" %in% names(payload)) {
        direct_eff <- payload$direct_effects
      } else {
        first_elem <- payload[[1]]
        if (is.list(first_elem) && "direct_effects" %in% names(first_elem)) {
          direct_eff <- first_elem$direct_effects
        }
      }
    }
    
    if (is.null(direct_eff)) {
      direct_eff <- payload
    }
    
    if (is.null(direct_eff) || (is.null(names(direct_eff)) && !is.data.frame(direct_eff))) {
      stop("Path payload must be a named vector or valid structured list containing direct effects.", call. = FALSE)
    }
    
    path_df <- data.frame(
      Trait = names(direct_eff), 
      Effect = as.numeric(direct_eff),
      stringsAsFactors = FALSE
    )
    path_df$Impact <- ifelse(path_df$Effect >= 0, "Positive Direct Effect", "Negative Direct Effect")
    
    g <- ggplot2::ggplot(path_df, ggplot2::aes(x = reorder(Trait, Effect), y = Effect, fill = Impact)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6, color = "black") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("Positive Direct Effect" = "#2ca02c", "Negative Direct Effect" = "#d62728")) +
      agri_custom_theme() + 
      ggplot2::labs(title = paste("Path Analysis Direct Effects on", trait_name), x = "Predictor Traits", y = "Path Coefficient Value")
    print(g)
  }
  
  # =========================================================================
  # BLOCK 4: CLOSURE
  # =========================================================================
  timestamp_end <- Sys.time()
  if (reporting_level >= 2) {
    cat("[LOG - FINALIZE]: Stream closed in ", round(as.numeric(difftime(timestamp_end, timestamp_start, units = "secs")), 5), " seconds.\n")
  }
  return(invisible(TRUE))
}

utils::globalVariables(c("PH", "GYPM", "NOSS", "NOT", "SL", "TGW", "PL", "gv_data", 
                         "across", "where", "select", "summarise", "group_by", 
                         "fviz_pca_biplot", "theme_minimal", "theme", 
                         "element_line", "element_text", "labs", "Effect", "Impact", "Trait", "Mean", "Genotype", "Genotype_Group", "Var1", "Var2", "value"))