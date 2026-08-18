# analysis/pipeline_comparisons/prevac_ad26mva/04_compare_model.R
#
# Compare predictive-model choices for predicting PREVAC Ad26/MVA
# (+ placebo) day-180 antibody titer from day-7 gene expression, against the
# reference pipeline (elastic net; see R/pipeline_defaults.R).
#
# treatment is not passed to this comparison - no model option uses it, and
# treatment_predictor = FALSE (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# regardless.
#
# NOTE: the "Linear regression" option below assumes `predictomics` exposes
# an unregularised OLS model under method = "lm" - see
# sdy1276_tiv/04_compare_model.R for the same caveat.
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
  "Linear regression"         = list(method = "lm", inner_folds = 10, metric = "r2"),
  "Lasso"                     = list(method = "lasso", inner_folds = 10, metric = "r2"),
  "Ridge"                     = list(method = "ridge", inner_folds = 10, metric = "r2"),
  "Random forest"             = list(method = "ranger", inner_folds = 10, metric = "r2"),
  "Support vector regression" = list(method = "svr", inner_folds = 10, metric = "r2")
)

res <- run_or_load_comparison(
  dataset = "prevac_ad26mva", label = "model",
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "model", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "prevac_ad26mva", category = "model")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: model choice (PREVAC Ad26/MVA)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_model_prevac_ad26mva.pdf")

gc()
rm(list = ls())
