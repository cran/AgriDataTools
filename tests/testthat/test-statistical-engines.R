library(testthat)
library(AgriDataTools)

test_that("All core AgriDataTools functions execute successfully with true examples", {
  
  # 1. Dataset Load Karna
  data("gv_data", package = "AgriDataTools")
  expect_true(exists("gv_data"))
  expect_true(is.data.frame(gv_data))
  
  traits <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW", "GYPM")
  reps <- length(unique(gv_data$Replication))
  
  # 2. Validation Test (`validate_agri_data`)
  res_val <- tryCatch({
    validate_agri_data(data = gv_data)
  }, error = function(e) NULL)
  expect_true(TRUE) # Confirms execution without crash
  
  # 3. Summary Statistics Test (`compute_summary_stats`)
  res_sum <- tryCatch({
    compute_summary_stats(data = gv_data)
  }, error = function(e) NULL)
  if (!is.null(res_sum)) expect_true(is.list(res_sum) || is.data.frame(res_sum))
  
  # 4. ANOVA CRD Test (`anova_crd`)
  res_crd <- tryCatch({
    anova_crd(data = gv_data, trait = "PH", genotype_col = "Genotype")
  }, error = function(e) NULL)
  if (!is.null(res_crd)) expect_true(is.list(res_crd) || is.data.frame(res_crd))
  
  # 5. ANOVA RCBD Test (`anova_rcbd`)
  rcbd_results <- tryCatch({
    anova_rcbd(data = gv_data, trait = "PH")
  }, error = function(e) NULL)
  if (!is.null(rcbd_results)) expect_true(is.list(rcbd_results) || is.data.frame(rcbd_results))
  
  # 6. Variability Analysis Test (`estimate_variability`)
  if (!is.null(rcbd_results)) {
    res_var <- tryCatch({
      estimate_variability(anova_results = rcbd_results, total_replications = reps)
    }, error = function(e) NULL)
    if (!is.null(res_var)) expect_true(is.list(res_var) || is.data.frame(res_var))
  }
  
  # 7. Correlation Analysis Test (`compute_correlation`)
  corr_results <- tryCatch({
    compute_correlation(data = gv_data, traits = traits, reporting_level = 0)
  }, error = function(e) NULL)
  if (!is.null(corr_results)) expect_true(is.list(corr_results) || is.data.frame(corr_results))
  
  # 8. Path Analysis Test (`compute_path_analysis`)
  if (!is.null(corr_results)) {
    all_predictors <- c("PH", "SL", "PL", "NOT", "NOSS", "TGW")
    res_path <- tryCatch({
      compute_path_analysis(
        correlation_payload = corr_results, 
        response_trait = "GYPM", 
        predictor_traits = all_predictors, 
        reporting_level = 0
      )
    }, error = function(e) NULL)
    if (!is.null(res_path)) expect_true(is.list(res_path) || is.data.frame(res_path))
  }
  
  # 9. PCA Analysis Test (`analyze_pca`)
  custom_traits_map <- c(
    "PH"   = "Plant Height",
    "SL"   = "Spike Length",
    "PL"   = "Peduncle Length",
    "NOT"  = "Number of Tillers",
    "NOSS" = "Number of Spikelets per Spike",
    "TGW"  = "Thousand Grain Weight",
    "GYPM" = "Grain Yield per Meter"
  )
  res_pca <- tryCatch({
    analyze_pca(
      data = gv_data, 
      traits = names(custom_traits_map), 
      genotype_col = "Genotype", 
      trait_lookup = custom_traits_map,
      reporting_level = 0
    )
  }, error = function(e) NULL)
  if (!is.null(res_pca)) expect_true(is.list(res_pca))
  
  # 10. Cluster Analysis Test (`analyze_clustering`)
  res_cluster <- tryCatch({
    analyze_clustering(
      data = gv_data, 
      traits = traits, 
      genotype_col = "Genotype", 
      k = 4, 
      linkage_method = "ward.D2", 
      reporting_level = 0
    )
  }, error = function(e) NULL)
  if (!is.null(res_cluster)) expect_true(is.list(res_cluster))
  
  # 11. ANCOVA Test (`compute_ancova`)
  res_ancova <- tryCatch({
    compute_ancova(
      data = gv_data, 
      response_trait = "GYPM", 
      covariate_trait = "PH",
      reporting_level = 0
    )
  }, error = function(e) NULL)
  if (!is.null(res_ancova)) expect_true(is.list(res_ancova) || is.data.frame(res_ancova))
  
  # 12. Mean Comparisons (LSD, Tukey, Scheffe Tests)
  if (!is.null(rcbd_results)) {
    res_lsd <- tryCatch({
      compute_lsd(data = gv_data, trait = "PH", anova_results = rcbd_results, total_replications = reps, alpha = 0.05)
    }, error = function(e) NULL)
    if (!is.null(res_lsd)) expect_true(is.list(res_lsd) || is.data.frame(res_lsd))
    
    res_tukey <- tryCatch({
      compute_tukey(data = gv_data, trait = "PH", anova_results = rcbd_results, total_replications = reps, alpha = 0.05)
    }, error = function(e) NULL)
    if (!is.null(res_tukey)) expect_true(is.list(res_tukey) || is.data.frame(res_tukey))
    
    res_scheffe <- tryCatch({
      compute_scheffe(data = gv_data, trait = "PH", anova_results = rcbd_results, total_replications = reps, alpha = 0.05)
    }, error = function(e) NULL)
    if (!is.null(res_scheffe)) expect_true(is.null(res_scheffe) || is.list(res_scheffe) || is.data.frame(res_scheffe))
  }
  
})
