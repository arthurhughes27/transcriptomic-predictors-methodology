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

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "best_pipeline_search.R"))

analysis_data <- readRDS(fs::path("data", "derived", "sdy1276_tiv_analysis_data.rds"))
single <- analysis_data$single
genesets <- analysis_data$genesets

best <- find_best_pipeline(
  X = single$X, Y = single$Y, covariates = single$covariates,
  genesets = genesets,
  dearseq_mode = NULL
)

cat("Best pipeline search: SDY1276 (TIV)\n")
print(best$summary)

derived_data_dir <- fs::path("data", "derived")
fs::dir_create(derived_data_dir)

saveRDS(best, file = fs::path(derived_data_dir, "best_pipeline_sdy1276_tiv.rds"))

# The winning model-round pipeline's fit is already a full predictomics
# predict_cv() object (compare_pipelines() fits every candidate, including
# the winner, via predict_cv() internally) - reused here rather than
# re-fitting an identical model from scratch.
best_fit <- best$round_results$model$fits[[best$winners$model$pipeline]]
saveRDS(best_fit, file = fs::path(derived_data_dir, "best_model_fit_sdy1276_tiv.rds"))

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")
fs::dir_create(figure_path)

p_fit <- tryCatch(plot(best_fit), error = function(e) {
  message("[best model] plot(best_fit) failed: ", conditionMessage(e))
  NULL
})
if (!is.null(p_fit)) {
  print(p_fit)
  ggsave(fs::path(figure_path, "best_model_fit_sdy1276_tiv.pdf"), p_fit, width = 6, height = 4.5, dpi = 300)
}

if (!is.null(best$selection_params)) {
  p_stability <- tryCatch(
    plot_selection_stability(best_fit, top_n = 20, type = "explicit", plot_type = "frequency"),
    error = function(e) {
      message("[best model] plot_selection_stability() failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(p_stability)) {
    print(p_stability)
    ggsave(
      fs::path(figure_path, "best_model_selection_stability_sdy1276_tiv.pdf"),
      p_stability, width = 7, height = 5, dpi = 300
    )
  }
}

rm(list = ls())
