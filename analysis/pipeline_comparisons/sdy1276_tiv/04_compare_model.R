# analysis/pipeline_comparisons/sdy1276_tiv/04_compare_model.R
#
# Compare predictive-model choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against the reference pipeline
# (elastic net; see R/pipeline_defaults.R).
#
# NOTE: the "Linear regression" option below assumes `predictomics` exposes
# an unregularised OLS model under method = "lm". This mirrors the naming
# convention of the other model methods used elsewhere in this repository
# ("lasso", "ridge", "ranger", "svr", "glmnet"), but has not been verified
# against the package source, since neither R nor the underlying gene
# expression data are available in the environment this script was written
# in. Please check this method string against the `predictomics` source
# before running, and correct it if needed.
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
  "Linear regression"         = list(method = "lm", inner_folds = 5, metric = "r2"),
  "Lasso"                     = list(method = "lasso", inner_folds = 5, metric = "r2"),
  "Ridge"                     = list(method = "ridge", inner_folds = 5, metric = "r2"),
  "Random forest"             = list(method = "ranger", inner_folds = 5, metric = "r2"),
  "Support vector regression" = list(method = "svr", inner_folds = 5, metric = "r2")
)

res <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "model", option_choices = option_choices,
  reference_params = reference_params
)

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: model choice (SDY1276, TIV)")

save_pipeline_comparison_plot(p, figure_path, "comparison_model_sdy1276_tiv.pdf")
