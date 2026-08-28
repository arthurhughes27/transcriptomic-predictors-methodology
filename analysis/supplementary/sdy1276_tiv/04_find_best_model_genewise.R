# analysis/supplementary/sdy1276_tiv/04_find_best_model_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/pipeline_comparisons/sdy1276_tiv/05_find_best_model.R: find the
# best-performing GENE-WISE (no gene-set aggregation) analytical pipeline
# for SDY1276 (TIV) via the same targeted, 3-round greedy coordinate-ascent
# search, restricted to gene-level options and the gene-wise reference
# pipeline throughout (see
# R/best_pipeline_search.R::find_best_pipeline_genewise()'s header).
#
# SDY1276 has no placebo arm, so RISE and dearseq are not offered as
# selection candidates here either (dearseq_mode = NULL) - same reasoning
# as the main search's own 05_find_best_model.R.
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R first.
#
# The greedy search itself (3 compare_pipelines() calls) is re-run only if
# no saved result exists yet at
# output/results/supplementary/best_pipeline_genewise_sdy1276_tiv.rds - once
# found, `best`/`best_fit` are just reloaded, so the plotting code below can
# be iterated on without re-searching every time.

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

results_dir <- fs::path("output", "results", "supplementary")
fs::dir_create(results_dir)

best_path <- fs::path(results_dir, "best_pipeline_genewise_sdy1276_tiv.rds")
best_fit_path <- fs::path(results_dir, "best_model_fit_genewise_sdy1276_tiv.rds")

if (fs::file_exists(best_path) && fs::file_exists(best_fit_path)) {

  best <- readRDS(best_path)
  best_fit <- readRDS(best_fit_path)

} else {

  analysis_data <- readRDS(fs::path("output", "results", "sdy1276_tiv_analysis_data.rds"))
  single <- analysis_data$single

  best <- find_best_pipeline_genewise(
    X = single$X, Y = single$Y, covariates = single$covariates,
    dearseq_mode = NULL
  )

  saveRDS(best, file = best_path)

  best_fit <- best$round_results$model$fits[[best$winners$model$pipeline]]
  saveRDS(best_fit, file = best_fit_path)
}

cat("Best gene-wise pipeline search: SDY1276 (TIV)\n")
print(best$summary)

figure_path <- fs::path("output", "figures", "supplementary", "sdy1276_tiv")
fs::dir_create(figure_path)

p_summary <- plot_best_model_summary(
  best_fit, best, title = "SDY1276 (TIV): best gene-wise pipeline",
  describe_fn = describe_best_pipeline_genewise
)
print(p_summary)
ggsave(
  fs::path(figure_path, "best_model_summary_genewise_sdy1276_tiv.pdf"),
  p_summary, width = 13.5, height = 8, dpi = 300
)

rm(list = ls())
