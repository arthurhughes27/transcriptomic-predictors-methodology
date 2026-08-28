# analysis/supplementary/prevac_ad26mva/04_find_best_model_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/pipeline_comparisons/prevac_ad26mva/05_find_best_model.R: find
# the best-performing GENE-WISE (no gene-set aggregation) analytical
# pipeline for PREVAC Ad26/MVA via the same targeted, 3-round greedy
# coordinate-ascent search, restricted to gene-level options and the
# gene-wise reference pipeline throughout (see
# R/best_pipeline_search.R::find_best_pipeline_genewise()'s header).
#
# PREVAC has a placebo arm, so RISE and dearseq (gene-level) are available
# as selection candidates here in "classic" mode (dearseq_mode = "classic",
# contrasting Ad26/MVA vs. placebo) - same reasoning as the main search's
# own 05_find_best_model.R.
#
# Run analysis/pipeline_comparisons/prevac_ad26mva/01_prepare_data.R first.
#
# The greedy search itself (3 compare_pipelines() calls) is re-run only if
# no saved result exists yet at
# output/results/supplementary/best_pipeline_genewise_prevac_ad26mva.rds -
# once found, `best`/`best_fit` are just reloaded, so the plotting code
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

results_dir <- fs::path("output", "results", "supplementary")
fs::dir_create(results_dir)

best_path <- fs::path(results_dir, "best_pipeline_genewise_prevac_ad26mva.rds")
best_fit_path <- fs::path(results_dir, "best_model_fit_genewise_prevac_ad26mva.rds")

if (fs::file_exists(best_path) && fs::file_exists(best_fit_path)) {

  best <- readRDS(best_path)
  best_fit <- readRDS(best_fit_path)

} else {

  analysis_data <- readRDS(fs::path("output", "results", "prevac_ad26mva_analysis_data.rds"))
  single <- analysis_data$single

  best <- find_best_pipeline_genewise(
    X = single$X, Y = single$Y, covariates = single$covariates,
    treatment = single$treatment,
    dearseq_mode = "classic"
  )

  saveRDS(best, file = best_path)

  best_fit <- best$round_results$model$fits[[best$winners$model$pipeline]]
  saveRDS(best_fit, file = best_fit_path)
}

cat("Best gene-wise pipeline search: PREVAC Ad26/MVA\n")
print(best$summary)

figure_path <- fs::path("output", "figures", "supplementary", "prevac_ad26mva")
fs::dir_create(figure_path)

p_summary <- plot_best_model_summary(
  best_fit, best, title = "PREVAC Ad26/MVA: best gene-wise pipeline",
  describe_fn = describe_best_pipeline_genewise
)
print(p_summary)
ggsave(
  fs::path(figure_path, "best_model_summary_genewise_prevac_ad26mva.pdf"),
  p_summary, width = 13.5, height = 5, dpi = 300
)

rm(list = ls())
