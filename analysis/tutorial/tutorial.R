# analysis/tutorial/tutorial.R
#
# A self-contained demonstration of the two main entry points of the
# `predictomics` package - predict_cv() and compare_pipelines() - on
# simulated data (predictomics::simulate_predictomics_data()), for the
# thesis's package-tutorial section, in the same style as the SurrogateRank
# tutorial.
#
# Deliberately independent of this repository's real, privacy-protected
# data (data/df_merged_all.rds etc.) and of the pipeline-comparison
# infrastructure under R/ (pipeline_defaults.R, run_comparison.R,
# metrics_io.R) built for those analyses: this script only illustrates
# package usage, so it stays fully self-contained and runnable by anyone
# with just `predictomics` installed - it is not part of the chapter's
# substantive analysis and is not sourced from analysis/master_analysis.R.
#
# Run this script from the repository root. It writes:
#   - output/figures/tutorial/*.pdf : one PDF per plot, for \includegraphics{}
#   - output/text/tutorial/*.txt    : one plain-text file per console
#     output, for pasting into a \begin{minted}{text} ... \end{minted} block
#     (or \verbatiminput{}/\lstinputlisting{} if preferred)
#
# See analysis/tutorial/README.md for how each output maps onto the chapter.

library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "plotting.R"))
source(fs::path("R", "text_output.R"))

figures_dir <- fs::path("output", "figures", "tutorial")
text_dir    <- fs::path("output", "text", "tutorial")


# =============================================================================
# Part 1: predict_cv() - a single cross-validated pipeline
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Simulate data
#
# simulate_predictomics_data() generates gene expression data (X), a
# continuous response (Y), a binary treatment indicator, covariates (age,
# sex), and a partition of the genes into named genesets. Two of the ten
# genesets are "signal" genesets whose mean expression drives Y; the rest
# are pure noise, so we know in advance which features a well-behaved
# pipeline should recover.
# -----------------------------------------------------------------------------
sim <- simulate_predictomics_data(
  n                 = 60,  # samples
  p                 = 200, # genes
  n_genesets        = 10,
  geneset_size      = 20,
  n_signal_genesets = 2,
  seed              = 1
)

save_console_output(capture.output(str(sim)), text_dir, "01_simulate_str.txt")

# -----------------------------------------------------------------------------
# 1.2 Fit a single pipeline with predict_cv()
#
# We aggregate genes into geneset-level mean expression (z-scored beforehand)
# and fit a lasso model, which performs embedded feature selection. Both the
# engineering and model steps are refit within every CV fold by default
# (outside_cv = FALSE), so the reported performance is not optimistically
# biased. model_params$scale = TRUE and compute_importance = TRUE additionally
# let us interpret which genesets the model actually relied on.
# -----------------------------------------------------------------------------
engineering_params <- list(
  method        = "engineer",
  col_transform = "z",
  genesets      = sim$genesets,
  agg_method    = "mean"
)

model_params <- list(
  method             = "lasso",
  scale              = TRUE,
  compute_importance = TRUE
)

fit <- predict_cv(
  Y                  = sim$Y,
  X                  = sim$X,
  engineering_params = engineering_params,
  model_params       = model_params,
  covariates         = sim$covariates,
  cv_type            = "kfold",
  folds              = 10,
  seed               = 12345,
  verbose            = FALSE
)

# -----------------------------------------------------------------------------
# 1.3 Inspect the result
# -----------------------------------------------------------------------------
save_console_output(capture.output(print(fit)), text_dir, "02_predict_cv_print.txt")

p_fit <- plot(fit)
save_pipeline_comparison_plot(p_fit, figures_dir, "tutorial_predict_cv_scatter.pdf",
                              width = 7, height = 5)

# Embedded selection stability: how consistently the lasso picked out each
# geneset across the 10 outer folds.
p_stability <- plot_selection_stability(fit, type = "embedded")
save_pipeline_comparison_plot(p_stability$frequency, figures_dir,
                              "tutorial_selection_stability.pdf", width = 7, height = 5)

# Feature importance: mean standardised |coefficient| across folds, which is
# meaningful regardless of how many genesets ended up being selected.
p_importance <- plot_feature_importance(fit)
save_pipeline_comparison_plot(p_importance, figures_dir,
                              "tutorial_feature_importance.pdf", width = 7, height = 5)

# The genesets that actually generated the signal, for comparison against
# the plots above.
save_console_output(capture.output(print(sim$signal_genesets)), text_dir,
                    "03_signal_genesets.txt")


# =============================================================================
# Part 2: compare_pipelines() - comparing alternative pipelines
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Compare three modelling choices against the lasso reference
#
# compare_pipelines() refits a reference pipeline plus one pipeline per entry
# of option_choices, all under an identical CV split, and returns a results
# table plus a comparison plot. Here we vary option_type = "model", holding
# the geneset-level engineering step from Part 1 fixed for every pipeline.
# Called directly (rather than via R/run_comparison.R::run_pipeline_comparison())
# since this dataset is small enough that the managed `future` parallel
# backend built for the real pipeline-comparison analyses isn't needed here.
# -----------------------------------------------------------------------------
cmp <- compare_pipelines(
  Y                  = sim$Y,
  X                  = sim$X,
  option_type        = "model",
  option_choices     = list(
    "Random forest" = list(method = "ranger"),
    "SVR"           = list(method = "svr", scale = TRUE)
  ),
  reference_params   = list(
    engineering_params = engineering_params,
    selection_params   = NULL,
    model_params        = model_params
  ),
  covariates         = sim$covariates,
  cv_type            = "kfold",
  folds              = 10,
  seed               = 12345,
  metric             = "sRMSE",
  verbose            = FALSE
)

# -----------------------------------------------------------------------------
# 2.2 Inspect the result
# -----------------------------------------------------------------------------
save_console_output(capture.output(print(cmp)), text_dir, "04_compare_pipelines_print.txt")

# A preview of the full results table, one row per pipeline (baseline,
# reference, and each option).
save_console_output(capture.output(print(cmp$results, row.names = FALSE)),
                    text_dir, "05_compare_pipelines_results.txt")

p_cmp <- plot(cmp)
save_pipeline_comparison_plot(p_cmp, figures_dir, "tutorial_compare_pipelines.pdf",
                              width = 8, height = 5)

message("Tutorial complete. Figures written to ", figures_dir,
       "; console outputs written to ", text_dir, ".")

rm(list = ls())
