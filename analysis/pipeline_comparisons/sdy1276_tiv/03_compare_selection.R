# analysis/pipeline_comparisons/sdy1276_tiv/03_compare_selection.R
#
# Compare feature-selection choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against the reference pipeline
# (no feature selection; see R/pipeline_defaults.R).
#
# Two comparisons are run, reported as two separate figures:
#
#   1. Selection methods that score each (mean-aggregated) engineered
#      feature directly against the continuous response: variance filtering,
#      absolute correlation (Spearman and Pearson), and univariate
#      regression-based screening ("relative gain"). SDY1276 is a single-arm
#      study with only one post-vaccination expression timepoint used as
#      predictors, so there is no treatment-vs-placebo or pre-vs-post
#      contrast to power a *classic* RISE/dearseq comparison here (see
#      `run_selection()`'s "classic" `dearseq_mode`, and `rise_paired =
#      FALSE`); those are more naturally compared on the placebo-controlled
#      Ebola/PREVAC datasets, in a future script.
#   2. RISE and dearseq's *paired* modes, which instead contrast each
#      participant's baseline (day 0) vs. day-1 gene expression to screen
#      for vaccine-responsive genes - this SDY1276 dataset does have that
#      contrast available (see 01_prepare_data.R's paired dataset). Both
#      return individual gene names (not gene-set scores), and are
#      implemented as a two-step process rather than a single
#      `compare_pipelines()` call - see the comment above that section for
#      why.
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
paired <- analysis_data$paired
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")

reference_params <- reference_pipeline_params(genesets)

## ---- 1. Selection methods scored against the continuous response ----------

option_choices <- list(
  "Variance (top 25)" = list(method = "variance", top_n = 25),
  "Variance (top 100)" = list(method = "variance", top_n = 100),
  "Correlation - Spearman (top 25)"   = list(method = "spearman", top_n = 25),
  "Correlation - Spearman (|r| > 0.5)" = list(method = "spearman", threshold = 0.5),
  "Correlation - Pearson (top 25)"    = list(method = "pearson", top_n = 100),
  "Correlation - Pearson (|r| > 0.5)"  = list(method = "pearson", threshold = 0.5),
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

## ---- 2. Paired RISE / dearseq (baseline vs. day-1 expression contrast) ----
#
# Both `rise_paired = TRUE` and `dearseq_mode = "paired"` screen genes on the
# PAIRED (baseline + day-1) dataset, contrasting each participant's
# timepoint == 1 (day 1) vs. timepoint == 0 (day 0) expression. That screen
# returns a set of individual gene names.
#
# For RISE, `predictomics::predict_cv()`/`compare_pipelines()` handle the
# "discard the pre-vaccination rows, then model on the selected genes" step
# automatically once `rise_paired = TRUE` and `individual_id`/`timepoint` are
# supplied. dearseq's paired mode does not have the same automatic row-drop -
# per the package's own documentation, it "operates upstream of engineering"
# and only filters *columns* (genes), leaving row handling to the caller.
#
# Rather than rely on two different automatic behaviours (and since neither
# selection method's output - individual gene names - lines up with the
# reference pipeline's gene-SET-aggregated columns; see
# `raw_gene_reference_params()`), both are run here as an explicit two-step
# process that works identically for both methods and mirrors exactly what
# you described: (1) select genes on the paired data, via
# `predictomics::run_selection()` directly; (2) discard the pre-vaccination
# (day 0) rows entirely and evaluate the selected genes on the single (day-1)
# dataset, via `compare_pipelines(option_type = "predictors")` against the
# same raw-gene reference used for the "Gene-set: none" option in
# 02_compare_engineering.R.
#
# NOTE: a handful of the genes RISE/dearseq select from the paired dataset
# may be absent from `single$gene_names` (each dataset's complete-case
# filtering in 01_prepare_data.R is applied independently, so a gene present
# in both timepoints of the paired data could still have been dropped from
# the single/day-1-only dataset, or vice versa). These are silently excluded
# via `intersect()` below.

rise_paired_fit <- predictomics::run_selection(
  X_train       = paired$X,
  Y_train       = paired$Y,
  covariates    = paired$covariates,
  individual_id = paired$participant_id,
  timepoint     = paired$timepoint,
  params = list(
    method            = "rise",
    rise_paired       = TRUE,
    top_n             = 25,
    rise_power_want_s = 0.8,
    rise_p_correction = "BH"
  )
)

dearseq_paired_fit <- predictomics::run_selection(
  X_train       = paired$X,
  Y_train       = NULL,
  covariates    = paired$covariates,
  individual_id = paired$participant_id,
  timepoint     = paired$timepoint,
  params = list(
    method        = "dearseq",
    dearseq_mode  = "paired",
    dearseq_level = "gene",
    threshold     = 0.05
  )
)

rise_paired_genes    <- intersect(rise_paired_fit$selected_features, single$gene_names)
dearseq_paired_genes <- intersect(dearseq_paired_fit$selected_features, single$gene_names)

cat(sprintf(
  "Paired RISE selected %d gene(s); paired dearseq selected %d gene(s) (%d/%d after intersecting with the single-timepoint dataset's genes).\n",
  length(rise_paired_fit$selected_features), length(dearseq_paired_fit$selected_features),
  length(rise_paired_genes), length(dearseq_paired_genes)
))

predictor_options <- list(
  "RISE (paired)"                = single$X[, rise_paired_genes, drop = FALSE],
  "Dearseq (paired, gene-level)" = single$X[, dearseq_paired_genes, drop = FALSE]
)

res_paired <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "predictors", option_choices = predictor_options,
  reference_params = raw_gene_reference_params()
)

p_paired <- plot(res_paired, metric = "R2") +
  ggtitle("Pipeline comparison: paired feature selection (SDY1276, TIV)")

print(p_paired)

save_pipeline_comparison_plot(p_paired, figure_path, "comparison_selection_paired_sdy1276_tiv.pdf")

gc()
rm(list = ls())
