# analysis/pipeline_comparisons/prevac_rvsv/05_find_best_model.R
#
# Find the best-performing analytical pipeline for PREVAC rVSV (+ placebo)
# via the targeted, 3-round greedy coordinate-ascent search implemented in
# R/best_pipeline_search.R::find_best_pipeline() (see that file's header for
# the full method and its trade-offs), then keep the winning pipeline's fit
# - already produced as a side effect of round 3's compare_pipelines() call,
# no separate refit needed - for downstream interpretation (selected
# features, model coefficients, etc.).
#
# Unlike SDY1276, PREVAC has a placebo arm, so RISE and dearseq are
# available as selection candidates here in their "classic" (treatment vs.
# placebo) mode - see R/best_pipeline_search.R's header. treatment is
# passed only for that screening contrast; treatment_predictor = FALSE (the
# default in R/run_comparison.R::run_pipeline_comparison()) keeps it out of
# the model itself throughout the search.
#
# Run analysis/pipeline_comparisons/prevac_rvsv/01_prepare_data.R first.
#
# The greedy search itself (3 compare_pipelines() calls) is re-run only if
# no saved result exists yet at results_dir/best_pipeline_prevac_rvsv.rds
# - once found, `best`/`best_fit` are just reloaded, so the plotting code
# below can be iterated on without re-searching every time.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)
library(patchwork)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "metrics_labels.R"))
source(fs::path("R", "best_pipeline_search.R"))
source(fs::path("R", "panel_helpers.R"))

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

best_path <- fs::path(results_dir, "best_pipeline_prevac_rvsv.rds")
best_fit_path <- fs::path(results_dir, "best_model_fit_prevac_rvsv.rds")

if (fs::file_exists(best_path) && fs::file_exists(best_fit_path)) {

  best <- readRDS(best_path)
  best_fit <- readRDS(best_fit_path)

} else {

  analysis_data <- readRDS(fs::path("output", "results", "prevac_rvsv_analysis_data.rds"))
  single <- analysis_data$single
  genesets <- analysis_data$genesets

  best <- find_best_pipeline(
    X = single$X, Y = single$Y, covariates = single$covariates,
    treatment = single$treatment,
    genesets = genesets,
    dearseq_mode = "classic"
  )

  saveRDS(best, file = best_path)

  # The winning model-round pipeline's fit is already a full predictomics
  # predict_cv() object (compare_pipelines() fits every candidate, including
  # the winner, via predict_cv() internally) - reused here rather than
  # re-fitting an identical model from scratch.
  best_fit <- best$round_results$model$fits[[best$winners$model$pipeline]]
  saveRDS(best_fit, file = best_fit_path)
}

cat("Best pipeline search: PREVAC rVSV (+ placebo)\n")
print(best$summary)

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_rvsv")
fs::dir_create(figure_path)

p_summary <- plot_best_model_summary(best_fit, best, title = "PREVAC (rVSV): best pipeline")
print(p_summary)
ggsave(
  fs::path(figure_path, "best_model_summary_prevac_rvsv.pdf"),
  p_summary, width = 13.5, height = 8, dpi = 300
)

rm(list = ls())
