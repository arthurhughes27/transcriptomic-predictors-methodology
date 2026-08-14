# analysis/pipeline_comparisons/sdy1276_tiv/03_compare_selection.R
#
# Compare feature-selection choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against the reference pipeline
# (mean gene-set aggregation; see R/pipeline_defaults.R::reference_pipeline_params()).
# This is the ONLY feature-selection comparison treated as a "reference"
# comparison for this dataset - see
# 03b_compare_selection_genewise.R for the supplementary, gene-wise (raw
# gene, no aggregation) comparison, which uses a different reference
# pipeline and is not used to characterise "the" reference approach for
# this dataset (e.g. in analysis/pipeline_comparisons/visualize_comparisons.R's
# baseline/reference context figure).
#
# Selection methods score each aggregated feature: variance filtering,
# absolute correlation (Spearman/Pearson), univariate regression-based
# screening ("relative gain"), and dearseq's *paired* mode at the geneset
# level (`dearseq_level = "geneset"`), which is compatible with geneset
# engineering. RISE is not included here - paired RISE screening always
# operates on the raw gene-level matrix, so `predictomics::predict_cv()`
# rejects `rise_paired = TRUE` combined with `engineering_params$genesets`
# outright (see 03b_compare_selection_genewise.R, where RISE is compared
# against the raw-gene reference instead).
#
# This comparison passes the PAIRED dataset (both timepoints) as the
# top-level X/Y/individual_id/timepoint to a single
# `compare_pipelines(option_type = "selection")` call, since dearseq's
# paired mode needs both arms to screen genes. Per predictomics' paired
# row-discard parity handling: any pipeline that discards pre-treatment
# rows internally (dearseq_mode = "paired") screens on both arms and then
# models on post-treatment rows only, while every other pipeline in the
# same call (including the reference and baseline) is automatically
# restricted to post-treatment (day-1) rows first - keeping every pipeline
# compared on an identical, independent (one row per participant) sample.
#
# NOTE: because this comparison sources its post-treatment-only data from
# restricting the PAIRED dataset (rather than from the independently
# complete-case-filtered `single` dataset used in 02_compare_engineering.R
# and 04_compare_model.R), the set of genes available here may differ
# slightly: paired's complete-case filtering (in 01_prepare_data.R) drops a
# gene if it has missing values at *either* timepoint, which can be
# stricter than single's own (day-1-only) filtering.
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

res_geneset <- run_or_load_comparison(
  dataset = "sdy1276_tiv", label = "selection_geneset",
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
rm(list = ls())
