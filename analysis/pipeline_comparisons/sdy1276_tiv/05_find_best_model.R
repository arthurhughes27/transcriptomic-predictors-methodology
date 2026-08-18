# analysis/pipeline_comparisons/sdy1276_tiv/05_find_best_model.R
#
# Find the best-performing analytical pipeline for SDY1276 (TIV) via the
# targeted, 3-round greedy coordinate-ascent search implemented in
# R/best_pipeline_search.R::find_best_pipeline() (see that file's header for
# the full method and its trade-offs), then keep the winning pipeline's fit
# - already produced as a side effect of round 3's compare_pipelines() call,
# no separate refit needed - for downstream interpretation (selected
# features, model coefficients, etc.).
#
# SDY1276 has no placebo arm, so RISE and dearseq (which would need either
# the paired baseline-vs-post contrast this single-data-only search
# deliberately avoids, or a treatment-vs-placebo contrast this dataset
# doesn't have) are not offered as selection candidates here - see
# R/best_pipeline_search.R's header for the full reasoning.
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R first.
#
# The greedy search itself (3 compare_pipelines() calls) is re-run only if
# no saved result exists yet at results_dir/best_pipeline_sdy1276_tiv.rds
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

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

best_path <- fs::path(results_dir, "best_pipeline_sdy1276_tiv.rds")
best_fit_path <- fs::path(results_dir, "best_model_fit_sdy1276_tiv.rds")

if (fs::file_exists(best_path) && fs::file_exists(best_fit_path)) {

  best <- readRDS(best_path)
  best_fit <- readRDS(best_fit_path)

} else {

  analysis_data <- readRDS(fs::path("output", "results", "sdy1276_tiv_analysis_data.rds"))
  single <- analysis_data$single
  genesets <- analysis_data$genesets

  best <- find_best_pipeline(
    X = single$X, Y = single$Y, covariates = single$covariates,
    genesets = genesets,
    dearseq_mode = NULL
  )

  saveRDS(best, file = best_path)

  # The winning model-round pipeline's fit is already a full predictomics
  # predict_cv() object (compare_pipelines() fits every candidate, including
  # the winner, via predict_cv() internally) - reused here rather than
  # re-fitting an identical model from scratch.
  best_fit <- best$round_results$model$fits[[best$winners$model$pipeline]]
  saveRDS(best_fit, file = best_fit_path)
}

cat("Best pipeline search: SDY1276 (TIV)\n")
print(best$summary)

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")
fs::dir_create(figure_path)

p_summary <- plot_best_model_summary(best_fit, best, title = "SDY1276 (TIV): best pipeline")
print(p_summary)
ggsave(
  fs::path(figure_path, "best_model_summary_sdy1276_tiv.pdf"),
  p_summary, width = 13, height = 6, dpi = 300
)

rm(list = ls())
