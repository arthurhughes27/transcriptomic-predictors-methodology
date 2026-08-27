# analysis/pipeline_comparisons/prevac_ad26mva/02_compare_engineering.R
#
# Compare feature-engineering choices for predicting PREVAC Ad26/MVA
# (+ placebo) day-180 antibody titer from day-7 gene expression, against
# the reference pipeline (z-scored gene-level transform, mean gene-set
# aggregation, elastic net model; see R/pipeline_defaults.R).
#
# GENESET-ONLY: this is one of "the main analyses" (see this folder's
# README.md) - only gene-set aggregation methods (median, max, 1st PC,
# GSVA, ssGSEA) applied on top of a z-scored gene-level transform are
# compared here, against the mean-aggregation reference. Gene-level
# (no-aggregation) engineering - z-score vs. no transform - is compared
# separately, as a SUPPLEMENTARY analysis, in
# analysis/supplementary/prevac_ad26mva/01_compare_engineering_genewise.R.
#
# treatment is not passed to this comparison - it's not used by any
# engineering option, and treatment_predictor = FALSE (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# regardless.
#
# Run analysis/pipeline_comparisons/prevac_ad26mva/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))
source(fs::path("R", "metrics_io.R"))

analysis_data <- readRDS(fs::path("output", "results", "prevac_ad26mva_analysis_data.rds"))
single <- analysis_data$single
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_ad26mva")

reference_params <- reference_pipeline_params(genesets)

option_choices <- list(
  "Gene-set: median" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "median"),
  "Gene-set: max"    = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "max"),
  "Gene-set: 1st PC" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "pc1"),
  "Gene-set: GSVA"   = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "gsva", gsva_min_size = 2),
  "Gene-set: ssGSEA" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "ssgsea", ssgsea_min_size = 2)
)

res <- run_or_load_comparison(
  dataset = "prevac_ad26mva", label = "engineering",
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "engineering", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "prevac_ad26mva", category = "engineering")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: gene-set aggregation (PREVAC Ad26/MVA)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_engineering_prevac_ad26mva.pdf")

gc()
rm(list = ls())
