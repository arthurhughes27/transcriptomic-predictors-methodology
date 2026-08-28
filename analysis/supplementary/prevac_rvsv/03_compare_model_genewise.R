# analysis/supplementary/prevac_rvsv/03_compare_model_genewise.R
#
# SUPPLEMENTARY comparison: predictive-model choices for predicting PREVAC
# rVSV (+ placebo) day-180 antibody titer from day-7 gene expression,
# against the gene-wise reference pipeline (z-scored gene-level transform,
# variance top-7,500 pre-filter, elastic net model; see
# R/pipeline_defaults.R::raw_gene_reference_params()).
#
# This is NOT sourced by
# analysis/pipeline_comparisons/prevac_rvsv/00_master.R and does not run as
# part of analysis/master_analysis.R - run it (and its siblings
# 01_compare_engineering_genewise.R/02_compare_selection_genewise.R) via
# analysis/supplementary/prevac_rvsv/00_master.R, or
# analysis/supplementary/master_supplementary.R for all three datasets.
#
# This is NOT treated as a/the reference approach for this dataset - see
# 04_compare_model.R (this folder's main, geneset-only model comparison)
# for that. Uses the SINGLE (day-7-only) dataset, exactly as
# 04_compare_model.R does - the fold-change transform is an engineering
# choice, not something model choice needs to hold fixed.
#
# Run analysis/pipeline_comparisons/prevac_rvsv/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))
source(fs::path("R", "metrics_io.R"))

analysis_data <- readRDS(fs::path("output", "results", "prevac_rvsv_analysis_data.rds"))
single <- analysis_data$single

figure_path <- fs::path("output", "figures", "supplementary", "prevac_rvsv")

reference_params <- raw_gene_reference_params()

option_choices <- list(
  "Linear regression"         = list(method = "lm", inner_folds = 10, metric = "r2", scale = TRUE, compute_importance = TRUE),
  "Lasso"                     = list(method = "lasso", inner_folds = 10, metric = "r2", scale = TRUE, compute_importance = TRUE),
  "Ridge"                     = list(method = "ridge", inner_folds = 10, metric = "r2", scale = TRUE, compute_importance = TRUE),
  "Random forest"             = list(method = "ranger", inner_folds = 10, metric = "r2", scale = TRUE, compute_importance = TRUE),
  "Support vector regression" = list(method = "svr", inner_folds = 10, metric = "r2", scale = TRUE, compute_importance = TRUE)
)

res <- run_or_load_comparison(
  dataset = "prevac_rvsv", label = "model_genewise",
  metrics_subdir = "supplementary",
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "model", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "prevac_rvsv", category = "model_genewise", subdir = "supplementary")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise model choice (PREVAC rVSV)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_model_genewise_prevac_rvsv.pdf")

gc()
rm(list = ls())
