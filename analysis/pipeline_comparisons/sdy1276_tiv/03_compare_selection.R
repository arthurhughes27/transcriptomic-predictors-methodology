# analysis/pipeline_comparisons/sdy1276_tiv/03_compare_selection.R
#
# Compare feature-selection choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against the reference pipeline
# (no feature selection; see R/pipeline_defaults.R).
#
# SDY1276 is a single-arm study with only one post-vaccination expression
# timepoint used as predictors, so there is no natural two-group contrast
# (e.g. treatment vs. placebo, or pre- vs. post-vaccination) to power the
# selection methods that require one: TLS screening with RISE
# (`rise.screen()` compares a "treated" vs. "control" group) and
# differential-expression testing with dearseq (which tests a
# `variables2test` group contrast). Those methods are more naturally compared
# on the placebo-controlled Ebola/PREVAC datasets, in a future script.
#
# Here we compare selection methods that instead score each engineered
# feature directly against the continuous response: variance filtering,
# absolute correlation (Spearman and Pearson), and univariate
# regression-based screening ("relative gain").
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))

analysis_data <- readRDS(fs::path("data", "derived", "sdy1276_tiv_analysis_data.rds"))
single <- analysis_data$single

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")

reference_params <- reference_pipeline_params()

option_choices <- list(
  "Variance (top 100)" = list(method = "variance", top_n = 100),
  "Variance (top 500)" = list(method = "variance", top_n = 500),
  "Correlation - Spearman (top 100)"   = list(method = "spearman", top_n = 100),
  "Correlation - Spearman (|r| > 0.3)" = list(method = "spearman", threshold = 0.3),
  "Correlation - Pearson (top 100)"    = list(method = "pearson", top_n = 100),
  "Correlation - Pearson (|r| > 0.3)"  = list(method = "pearson", threshold = 0.3),
  "Univariate regression screening (threshold = 0)" = list(
    method = "relative_gain", threshold = 0,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse", relative_gain_seed = 12345
  ),
  "Univariate regression screening (threshold = 0.1)" = list(
    method = "relative_gain", threshold = 0.1,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse", relative_gain_seed = 12345
  )
)

res <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "selection", option_choices = option_choices,
  reference_params = reference_params
)

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: feature selection (SDY1276, TIV)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_selection_sdy1276_tiv.pdf")

gc()
rm(list = ls())
