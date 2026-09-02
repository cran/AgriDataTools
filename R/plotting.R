#' Advanced High-Precision Publication-Ready Graphics Suite
#'
#' @description
#' The \code{plot_agri_graphics} function serves as the unified visualization hub for
#' the \code{AgriDataTools} package. It handles basic statistical diagnostics
#' (residuals, correlations) alongside modern publication-grade graphical
#' representations for mean performance, PCA space, hierarchical dendrograms, 
#' and path analysis direct effects.
#'
#' @import ggplot2
#' @import dplyr
#' @import factoextra
#' @import ggrepel
#' @importFrom dendextend color_branches circlize_dendrogram color_labels
#' @importFrom circlize circos.clear circos.initialize circos.trackPlotRegion get.cell.meta.data circos.text circos.rect draw.sector
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dendrogram dist hclust reorder aov lm residuals fitted qqnorm qqline aggregate prcomp as.formula cutree
#' @importFrom graphics par abline plot title legend
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
#'          \code{"residual"}, \code{"correlation"}, \code{"mean"}, \code{"pca"},
#'          \code{"cluster"}, or \code{"path"}.
#' @param payload A structured analysis \code{list} derived from computational
#'          engines (e.g., \code{compute_lsd}, \code{analyze_pca},
#'          \code{analyze_clustering}, \code{compute_path_analysis}).
#' @param trait_name A character string defining the target phenotypic trait
#'          title label. Used primarily in \code{"mean"} and \code{"path"} layouts.
#' @param reporting_level An integer vector flag defining console trace settings:
#'          \code{0} for silent, \code{1} for structural updates, and \code{2} for
#'          exhaustive analytical tracing. Defaults to \code{0} (Clean Console).
#' @param num_clusters An integer specifying the number of cluster groups to color 
#'          in the circular dendrogram module. Defaults to \code{4}.
#'
#' @return Invisibly returns a logical scalar \code{TRUE} upon successful execution.
#'      This function is primarily invoked for its side effect of rendering publication-grade
#'      graphical plots (e.g., residual diagnostic plots, correlation heatmaps, mean performance
#'      barcharts, PCA biplots, circular dendrograms, or path analysis plots) to the active 
#'      graphics device.
#'
#' @examples
#' data(gv_data, package = "AgriDataTools")
#' traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
#' 
#' # Define custom mapping or number of clusters beforehand
#' k_groups <- 4
#'
#' # 1. Mean performance: Genotypic performance with LSD
#' if (interactive()) {
#'   reps <- length(unique(gv_data$Replication))
#'   fit <- aov(PH ~ Genotype + Replication, data = gv_data)
#'   m_anova <- list(anova_table = data.frame(
#'     Source = c("Genotype", "Replication", "Error"),
#'     Df = summary(fit)[[1]]$Df,
#'     MS = summary(fit)[[1]][[3]]
#'   ))
#'   lsd_output <- compute_lsd(gv_data, "PH", m_anova, reps, reporting_level = 0)
#'   plot_agri_graphics(type = "mean", payload = lsd_output,
#'                      trait_name = "Plant Height")
#' }
#'
#' # 2. PCA: Multivariate variation with unique actual genotype names
#' if (interactive()) {
#'    custom_traits_map <- c(
#'        "PH"    = "Plant Height",
#'        "SL"    = "Spike Length",
#'        "PL"    = "Peduncle Length",
#'        "NOT"   = "Number of Tillers",
#'        "NOSS"  = "Number of Spikelets per Spike",
#'        "TGW"   = "Thousand Grain Weight",
#'        "GYPM"  = "Grain Yield per Meter"
#'    )
#'    pca_results <- analyze_pca(data = gv_data, traits = traits, 
#'                               trait_lookup = custom_traits_map, 
#'                               reporting_level = 0)
#'    plot_agri_graphics(type = "pca", payload = pca_results,
#'                       trait_name = "PCA Plot")
#' }
#'
#' # 3. Clustering: Dendrogram with flexible cluster parameter option
#' if (interactive()) {
#'    cluster_results <- analyze_clustering(data = gv_data, traits = traits,
#'     k = k_groups, reporting_level = 0)
#'    plot_agri_graphics(type = "cluster", payload = cluster_results,
#'                       trait_name = "Clustering", num_clusters = k_groups)
#' }
#'
#' # 4. Residuals: Diagnostic plots
#' if (interactive()) {
#'   fit <- lm(PH ~ Genotype, data = gv_data)
#'   res_pl <- list(residuals = residuals(fit),
#'                  fitted_values = fitted(fit))
#'   plot_agri_graphics(type = "residual", payload = res_pl,
#'                      trait_name = "Residuals")
#' }
#'
#' # 5. Correlations: Phenotypic matrix
#' if (interactive()) {
#'   cor_m <- cor(gv_data[, traits], use = "pairwise.complete.obs")
#'   plot_agri_graphics(type = "correlation",
#'                      payload = list(correlation_matrix = cor_m),
#'                      trait_name = "Correlation")
#' }
#'
#' # 6. Path Analysis: Direct Effects Plot for Grain Yield per Meter (GYPM)
#' if (interactive()) {
#'   corr_results <- compute_correlation(data = gv_data, traits = traits, reporting_level = 0)
#'   all_predictors <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW")
#'   path_results <- compute_path_analysis(
#'     correlation_payload = corr_results, 
#'     response_trait = "GYPM", 
#'     predictor_traits = all_predictors,
#'     reporting_level = 0
#'   )
#'   plot_agri_graphics(type = "path", payload = path_results, 
#'                      trait_name = "Grain Yield per Meter (GYPM)")
#' }
#'
plot_agri_graphics <- function(type, payload, trait_name = "Target Character Matrix", reporting_level = 0, num_clusters = 4) {
  
  # =========================================================================
  # BLOCK 1: STARTUP PARAMETERS AND SYSTEM DEFENSE ROUTINES
  # =========================================================================
  
  Var1 <- Var2 <- value <- Effect <- Impact <- Trait <- Mean <- Genotype <- PC1 <- PC2 <- PC1_scaled <- PC2_scaled <- Name <- Type <- x <- y <- NULL
  
  if (missing(type) || missing(payload)) {
    stop("CRITICAL INPUT FAULT: Missing arguments. Provide both visualization type and payload.", call. = FALSE)
  }
  
  chart_type <- match.arg(tolower(trimws(type)), c("residual", "correlation", "mean", "pca", "cluster", "path"))
  
  if (reporting_level >= 1) {
    message("AGRI-GRAPHICS ENGINE: Initializing rendering pipeline for module -> [", toupper(chart_type), "]")
  }
  
  # --- Theme Configuration ---
  agri_custom_theme <- function() {
    ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15, color = "#1a252f", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11, color = "#34495e", hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.major = ggplot2::element_line(color = "gray90", linetype = "dashed"),
        legend.position = "right"
      )
  }
  
  # =========================================================================
  # BLOCK 2: MODULES EXECUTION CHAIN & RENDERING LOGIC
  # =========================================================================
  
  if (chart_type == "correlation") {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("CRITICAL DEPENDENCY FAULT: Package 'reshape2' is required for correlation heatmaps.", call. = FALSE)
    }
    
    corr_mat <- if (is.list(payload) && "correlation_matrix" %in% names(payload)) {
      payload$correlation_matrix
    } else {
      payload
    }
    
    corr_mat <- as.matrix(corr_mat)
    corr_mat[upper.tri(corr_mat)] <- NA
    melted_cormat <- reshape2::melt(corr_mat, na.rm = TRUE)
    
    g <- ggplot2::ggplot(melted_cormat, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient2(
        low = "#d73027", mid = "#ffffbf", high = "#4575b4",
        midpoint = 0, limit = c(-1, 1), space = "Lab",
        name = "Pearson\nCorrelation (r)"
      ) +
      agri_custom_theme() +
      ggplot2::labs(
        title = paste("Pearson Correlation Matrix Heatmap"),
        subtitle = "Multi-trait Linear Associations across Examined Traits",
        x = "", y = ""
      )
    
    print(g)
    
  } else if (chart_type == "residual") {
    resids <- payload$residuals
    fitted_vals <- payload$fitted_values
    
    if (is.null(resids) || is.null(fitted_vals)) {
      stop("RESIDUAL FAULT: Payload must contain both 'residuals' and 'fitted_values' vectors.", call. = FALSE)
    }
    
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5))
    
    graphics::plot(
      fitted_vals, resids,
      main = "Residuals vs Fitted Values",
      xlab = "Fitted Values", ylab = "Residuals",
      pch = 19, col = "#1b4f72", cex = 1.1
    )
    graphics::abline(h = 0, lty = 2, col = "#c0392b", lwd = 2)
    
    stats::qqnorm(
      resids,
      main = "Normal Q-Q Distribution",
      xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
      pch = 19, col = "#1b4f72", cex = 1.1
    )
    stats::qqline(resids, col = "#c0392b", lwd = 2)
    
  } else if (chart_type == "mean") {
    target_df <- if (is.list(payload) && "ranked_means" %in% names(payload)) {
      payload$ranked_means
    } else if (is.data.frame(payload)) {
      payload
    } else {
      as.data.frame(payload)
    }
    
    geno_col_candidates <- c("Genotype", "Genotypes", "Cultivar", "Cultivars", "Line", "Lines", "Gen", "Accession", "Accessions", "Row.names")
    matched_geno_col <- geno_col_candidates[geno_col_candidates %in% colnames(target_df)]
    
    if (length(matched_geno_col) > 0) {
      target_df$Genotype <- target_df[[matched_geno_col[1]]]
    } else {
      target_df$Genotype <- rownames(target_df)
    }
    
    if (!"Mean" %in% colnames(target_df)) {
      num_cols <- sapply(target_df, is.numeric)
      if (any(num_cols)) {
        target_df$Mean <- target_df[[which(num_cols)[1]]]
      } else {
        stop("MEAN PLOT FAULT: Could not locate numeric mean performance column in payload.", call. = FALSE)
      }
    }
    
    g <- ggplot2::ggplot(target_df, ggplot2::aes(x = reorder(Genotype, -Mean), y = Mean, fill = Mean)) +
      ggplot2::geom_bar(stat = "identity", color = "black", width = 0.75, linewidth = 0.4) +
      ggplot2::scale_fill_viridis_c(option = "mako", direction = -1, name = "Mean Value") +
      agri_custom_theme() +
      ggplot2::labs(
        title = paste("Mean Performance for", trait_name),
        subtitle = "Evaluated across replications with LSD benchmarks",
        x = "Genotypes",
        y = paste("Performance Mean (", trait_name, ")", sep = "")
      )
    
    print(g)
    
  } else if (chart_type == "pca") {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("CRITICAL DEPENDENCY FAULT: Package 'ggrepel' is required for modern professional PCA biplots.", call. = FALSE)
    }
    
    pca_obj <- if (is.list(payload) && "pca_object" %in% names(payload)) {
      payload$pca_object 
    } else {
      payload
    }
    
    if (is.null(pca_obj) || !inherits(pca_obj, "prcomp")) {
      stop("PCA FAULT: Valid 'prcomp' object or pca list payload required.", call. = FALSE)
    }
    
    res.pca <- pca_obj
    
    if (!is.null(res.pca$x)) {
      clean_names <- sub("(_[Rr]ep?[0-9]+|\\.[0-9]+|_R[0-9]+)$", "", rownames(res.pca$x))
      rownames(res.pca$x) <- make.unique(as.character(clean_names))
    }
    
    var_exp_1 <- round((res.pca$sdev[1]^2 / sum(res.pca$sdev^2)) * 100, 1)
    var_exp_2 <- round((res.pca$sdev[2]^2 / sum(res.pca$sdev^2)) * 100, 1)
    
    ind_coords <- as.data.frame(res.pca$x[, 1:2])
    colnames(ind_coords) <- c("PC1", "PC2")
    ind_coords$Name <- rownames(ind_coords)
    
    var_coords <- as.data.frame(res.pca$rotation[, 1:2])
    colnames(var_coords) <- c("PC1", "PC2")
    var_coords$Name <- rownames(var_coords)
    
    scale_factor <- min(max(abs(ind_coords$PC1))/max(abs(var_coords$PC1)), 
                        max(abs(ind_coords$PC2))/max(abs(var_coords$PC2))) * 0.75
    var_coords$PC1_scaled <- var_coords$PC1 * scale_factor
    var_coords$PC2_scaled <- var_coords$PC2 * scale_factor
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, color = "#bdc3c7", linewidth = 0.5, linetype = "solid") +
      ggplot2::geom_vline(xintercept = 0, color = "#bdc3c7", linewidth = 0.5, linetype = "solid") +
      
      ggplot2::geom_point(data = ind_coords, ggplot2::aes(x = PC1, y = PC2), 
                          color = "#117a65", fill = "#48c9b0", shape = 21, size = 3, stroke = 0.8, alpha = 0.9) +
      ggrepel::geom_text_repel(data = ind_coords, ggplot2::aes(x = PC1, y = PC2, label = Name), 
                               color = "#117a65", fontface = "bold", size = 3.2, max.overlaps = Inf, box.padding = 0.35) +
      
      ggplot2::geom_segment(data = var_coords, ggplot2::aes(x = 0, y = 0, xend = .data$PC1_scaled, yend = .data$PC2_scaled),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "#b9770e", linewidth = 1.1, alpha = 0.9) +
      ggrepel::geom_text_repel(data = var_coords, ggplot2::aes(x = .data$PC1_scaled, y = .data$PC2_scaled, label = Name),
                               color = "#7e5109", fontface = "bold", size = 4.1, bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
      
      ggplot2::labs(
        title = "Principal Component Analysis (PCA) Biplot",
        subtitle = "Genotypic Dispersion and Phenotypic Trait Loadings Vector",
        x = paste("Principal Component 1 (", var_exp_1, "%)", sep = ""),
        y = paste("Principal Component 2 (", var_exp_2, "%)", sep = "")
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5, color = "#2c3e50"),
        plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "#7f8c8d", margin = ggplot2::margin(b = 10)),
        axis.text = ggplot2::element_text(face = "bold", color = "#34495e"),
        axis.title = ggplot2::element_text(face = "bold", color = "#2c3e50"),
        axis.line = ggplot2::element_line(color = "#2c3e50", linewidth = 0.8),
        plot.margin = ggplot2::unit(c(0.8, 0.8, 0.8, 0.8), "cm")
      )
    
    print(p)
    
  } else if (chart_type == "cluster") {
    hc_object <- if (is.list(payload) && "hc_object" %in% names(payload)) {
      payload$hc_object
    } else {
      payload
    }
    
    if (is.null(hc_object) && inherits(payload, "hclust")) {
      hc_object <- payload
    }
    
    if (is.null(hc_object)) {
      stop("CLUSTER PAYLOAD FAULT: Provide a valid 'hclust' object or cluster analysis list.", call. = FALSE)
    }
    
    dend <- stats::as.dendrogram(hc_object)
    
    # Fully dynamic palette supporting any number of clusters
    base_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39C12", "#8E44AD", "#27AE60", "#D35400", "#2980B9", "#C0392B", "#16A085", "#7F8C8D")
    if (num_clusters <= length(base_palette)) {
      cluster_colors <- base_palette[1:num_clusters]
    } else {
      cluster_colors <- grDevices::colorRampPalette(base_palette)(num_clusters)
    }
    
    # Color branches dynamically by cluster groups for smooth structural transition
    dend <- dendextend::color_branches(dend, k = num_clusters, col = cluster_colors, lwd = 1.4)
    
    # Safe label coloring using cutree and fully qualified stats::order.dendrogram
    cluster_cut <- stats::cutree(hc_object, k = num_clusters)
    ordered_cut <- cluster_cut[stats::order.dendrogram(dend)]
    dendextend::labels_colors(dend) <- cluster_colors[ordered_cut]
    
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
      graphics::par(old_par)
      try(circlize::circos.clear(), silent = TRUE)
    }, add = TRUE)
    
    circlize::circos.clear()
    graphics::par(mar = c(1, 1, 3, 1), bg = "white")
    circlize::circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0.01))
    
    dendextend::circlize_dendrogram(
      dend, 
      facing = "outside", 
      labels_track_height = 0.35, 
      dend_track_height = 0.45
    )
    
    graphics::title(
      main = "Cluster Circular Dendrogram", 
      col.main = "#1a252f", font.main = 2, cex.main = 1.3, line = 0.5
    )
    
    try(circlize::circos.clear(), silent = TRUE)
    
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
      stop("PATH ANALYSIS FAULT: Payload must contain valid direct effects vector or nested structure.", call. = FALSE)
    }
    
    path_df <- data.frame(
      Trait = names(direct_eff), 
      Effect = as.numeric(direct_eff),
      stringsAsFactors = FALSE
    )
    path_df$Impact <- ifelse(path_df$Effect >= 0, "Positive Direct Effect", "Negative Direct Effect")
    
    g <- ggplot2::ggplot(path_df, ggplot2::aes(x = reorder(Trait, Effect), y = Effect, fill = Impact)) +
      ggplot2::geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.4) + 
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("Positive Direct Effect" = "#27ae60", "Negative Direct Effect" = "#c0392b")) +
      agri_custom_theme() + 
      ggplot2::labs(
        title = paste("Path Analysis Direct Effects on", trait_name),
        subtitle = "Decomposition of correlation coefficients into direct causal impacts",
        x = "Predictor Traits",
        y = "Path Coefficient Value"
      )
    
    print(g)
  }
  
  if (reporting_level >= 1) {
    message("AGRI-GRAPHICS ENGINE: Rendering successfully executed for -> [", toupper(chart_type), "]")
  }
  
  return(invisible(TRUE))
}

utils::globalVariables(c("PH", "GYPM", "NOSS", "NOT", "SL", "TGW", "PL", "gv_data", 
                         "across", "where", "select", "summarise", "group_by", 
                         "theme_minimal", "theme", "element_line", "element_text", 
                         "labs", "Effect", "Impact", "Trait", "Mean", "Genotype", 
                         "Genotype_Group", "Var1", "Var2", "value", "PC1", "PC2", 
                         "PC1_scaled", "PC2_scaled", "Name", "Type"))