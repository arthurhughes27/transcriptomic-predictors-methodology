# analysis/pipeline_comparisons/sdy1276_tiv/03_compare_selection.R
#
# Compare feature-selection choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression.
#
# Two comparisons are run, reported as two separate figures, one per
# engineering scale (RISE and dearseq behave differently with respect to
# gene-set aggregation - see below - so a single figure spanning both scales
# isn't meaningful):
#
#   1. Geneset-level: reference is the usual mean-BTM-aggregated pipeline
#      (`reference_pipeline_params()`). Selection methods score each
#      aggregated feature: variance filtering, absolute correlation
#      (Spearman/Pearson), univariate regression-based screening ("relative
#      gain"), and dearseq's *paired* mode at the geneset level
#      (`dearseq_level = "geneset"`), which is compatible with geneset
#      engineering. RISE is not included here - paired RISE screening always
#      operates on the raw gene-level matrix, so `predictomics::predict_cv()`
#      now rejects `rise_paired = TRUE` combined with
#      `engineering_params$genesets` outright.
#   2. Gene-wise: reference is z-score only, no gene-set aggregation
#      (`raw_gene_reference_params()`). The same selection methods are
#      compared, with top_n thresholds rescaled from the geneset-level
#      values above to the ~20,000-gene scale (see comment below), plus
#      paired RISE and dearseq's paired mode at the gene level
#      (`dearseq_level = "gene"`).
#
# Both comparisons contrast each participant's baseline (day 0) vs. day-1
# gene expression to power RISE/dearseq's *paired* modes (SDY1276 has no
# treatment-vs-placebo arm for their *classic* modes - those are more
# naturally compared on the placebo-controlled Ebola/PREVAC datasets, in a
# future script).
#
# Both comparisons pass the PAIRED dataset (both timepoints) as the top-level
# X/Y/individual_id/timepoint to a single `compare_pipelines(option_type =
# "selection")` call. Per predictomics' paired row-discard parity handling:
# any pipeline that discards pre-treatment rows internally (rise_paired,
# dearseq_mode = "paired") screens on both arms and then models on
# post-treatment rows only, while every other pipeline in the *same* call
# (including the reference and baseline) is automatically restricted to
# post-treatment (day-1) rows first - keeping every pipeline compared on an
# identical, independent (one row per participant) sample, in one call and
# one figure, regardless of whether a given option happens to be paired.
#
# NOTE: because both comparisons now source their post-treatment-only data
# from restricting the PAIRED dataset (rather than from the independently
# complete-case-filtered `single` dataset used in 02_compare_engineering.R
# and 04_compare_model.R), the set of genes available here may differ
# slightly: paired's complete-case filtering (in 01_prepare_data.R) drops a
# gene if it has missing values at *either* timepoint, which can be stricter
# than single's own (day-1-only) filtering.
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))
source(fs::path("R", "metrics_io.R"))

analysis_data <- readRDS(fs::path("data", "derived", "sdy1276_tiv_analysis_data.rds"))
paired <- analysis_data$paired
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")

## ---- 1. Geneset-level feature selection ------------------------------------

reference_params_geneset <- reference_pipeline_params(genesets)

option_choices_geneset <- list(
  "Variance (top 25)" = list(method = "variance", top_n = 25),
  "Variance (top 100)" = list(method = "variance", top_n = 100),
  "Correlation - Spearman (top 25)"   = list(method = "spearman", top_n = 25),
  "Correlation - Spearman (|r| > 0.5)" = list(method = "spearman", threshold = 0.5),
  "Correlation - Pearson (top 25)"    = list(method = "pearson", top_n = 100),
  "Correlation - Pearson (|r| > 0.5)"  = list(method = "pearson", threshold = 0.5),
  "Univariate regression screening (threshold = 0)" = list(
    method = "relative_gain", threshold = 0,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse"
  ),
  "Dearseq (alpha = 0.05)" = list(
    method = "dearseq", dearseq_mode = "paired", dearseq_level = "geneset",
    genesets = genesets, threshold = 0.05
  )
)

res_geneset <- run_pipeline_comparison(
  X = paired$X, Y = paired$Y, covariates = paired$covariates,
  individual_id = paired$participant_id, timepoint = paired$timepoint,
  option_type = "selection", option_choices = option_choices_geneset,
  reference_params = reference_params_geneset
)

save_comparison_metrics(res_geneset, dataset = "sdy1276_tiv", category = "selection_geneset")

p_geneset <- plot(res_geneset, metric = "R2") +
  ggtitle("Pipeline comparison: geneset-level feature selection (SDY1276, TIV)")

print(p_geneset)

save_pipeline_comparison_plot(p_geneset, figure_path, "comparison_selection_geneset_sdy1276_tiv.pdf")

gc()

## ---- 2. Gene-wise (raw) feature selection ----------------------------------
#
# top_n thresholds are rescaled ~100x from their geneset-level counterparts
# above, matching the ~100x jump from a few hundred BTM gene sets to
# ~20,000 raw genes (the same scale reflected in
# raw_gene_reference_params()'s own top_n = 7,500 variance pre-filter).
# Threshold-based options (correlation |r|, relative-gain, dearseq p-value)
# are on a fixed (dimensionless / p-value) scale and so are left unchanged.

reference_params_genewise <- raw_gene_reference_params()

option_choices_genewise <- list(
  "Variance (top 500)" = list(method = "variance", top_n = 500),
  "Correlation - Spearman (top 500)"  = list(method = "spearman", top_n = 500),
  "Correlation - Spearman (|r| > 0.5)"  = list(method = "spearman", threshold = 0.5),
  "Correlation - Pearson (top 500)"  = list(method = "pearson", top_n = 500),
  "Correlation - Pearson (|r| > 0.5)"   = list(method = "pearson", threshold = 0.5),
  "Univariate regression screening (threshold = 0)" = list(
    method = "relative_gain", threshold = 0,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse" 
  ),
  "RISE (top 500)" = list(
    method = "rise", rise_paired = TRUE, top_n = 500,
    rise_power_want_s = 0.8, rise_p_correction = "BH"
  ),
  "Dearseq (alpha = 0.05)" = list(
    method = "dearseq", dearseq_mode = "paired", dearseq_level = "gene",
    threshold = 0.05
  )
)

res_genewise <- run_pipeline_comparison(
  X = paired$X, Y = paired$Y, covariates = paired$covariates,
  individual_id = paired$participant_id, timepoint = paired$timepoint,
  option_type = "selection", option_choices = option_choices_genewise,
  reference_params = reference_params_genewise
)

save_comparison_metrics(res_genewise, dataset = "sdy1276_tiv", category = "selection_genewise")

p_genewise <- plot(res_genewise, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise feature selection (SDY1276, TIV)")

print(p_genewise)

save_pipeline_comparison_plot(p_genewise, figure_path, "comparison_selection_genewise_sdy1276_tiv.pdf")

gc()
rm(list = ls())
