# analysis/pipeline_comparisons/sdy1276_tiv/03b_compare_selection_genewise.R
#
# SUPPLEMENTARY comparison: feature selection at the gene-wise (raw gene,
# no gene-set aggregation) level, for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against
# R/pipeline_defaults.R::raw_gene_reference_params() (z-score only, no
# aggregation).
#
# This is NOT treated as a/the reference approach for this dataset - only
# the geneset-mean pipeline (03_compare_selection.R, and
# 02_compare_engineering.R/04_compare_model.R) is. This comparison exists
# to characterise feature-selection choices at the gene-wise scale in their
# own right (a different, incompatible reference pipeline is required for
# RISE and gene-level dearseq - see below), not to compete with the main
# geneset-level comparison for "the" dataset-level baseline.
#
# Adds RISE (paired mode) and dearseq's *paired* mode at the gene level
# (`dearseq_level = "gene"`), on top of the same variance/correlation/
# relative-gain methods compared at the geneset level in
# 03_compare_selection.R. top_n thresholds are rescaled ~100x from that
# script's geneset-level values, matching the ~100x jump from a few hundred
# BTM gene sets to ~20,000 raw genes (the same scale reflected in
# raw_gene_reference_params()'s own top_n = 7,500 variance pre-filter).
# Threshold-based options (correlation |r|, relative-gain, dearseq p-value)
# are on a fixed (dimensionless/p-value) scale and so are left unchanged.
#
# As in 03_compare_selection.R, this passes the PAIRED dataset as the
# top-level X/Y/individual_id/timepoint, relying on predictomics' paired
# row-discard parity handling so RISE/dearseq's paired options and the
# non-paired options/reference all end up compared on the same
# post-treatment-only sample within this one call.
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

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")

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
