# analysis/pipeline_comparisons/prevac_rvsv/03b_compare_selection_genewise.R
#
# SUPPLEMENTARY comparison: feature selection at the gene-wise (raw gene,
# no gene-set aggregation) level, for predicting PREVAC rVSV (+ placebo)
# day-180 antibody titer from day-7 gene expression, against
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
# Adds RISE and dearseq's "classic" mode at the gene level
# (`dearseq_level = "gene"`), on top of the same variance/correlation/
# relative-gain methods compared at the geneset level in
# 03_compare_selection.R. Both use the placebo arm as their treatment
# contrast, same as 03_compare_selection.R's dearseq option - see that
# script's header for the "classic" mode/treatment_predictor details.
#
# Run analysis/pipeline_comparisons/prevac_rvsv/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))
source(fs::path("R", "metrics_io.R"))

analysis_data <- readRDS(fs::path("data", "derived", "prevac_rvsv_analysis_data.rds"))
single <- analysis_data$single

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_rvsv")

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
    method = "rise", top_n = 500,
    rise_power_want_s = 0.8, rise_p_correction = "BH"
  ),
  "Dearseq (gene, alpha = 0.05)" = list(
    method = "dearseq", dearseq_mode = "classic", dearseq_level = "gene",
    threshold = 0.05
  )
)

res_genewise <- run_or_load_comparison(
  dataset = "prevac_rvsv", label = "selection_genewise",
  X = single$X, Y = single$Y, covariates = single$covariates,
  treatment = single$treatment,
  option_type = "selection", option_choices = option_choices_genewise,
  reference_params = reference_params_genewise
)

save_comparison_metrics(res_genewise, dataset = "prevac_rvsv", category = "selection_genewise")

p_genewise <- plot(res_genewise, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise feature selection (PREVAC rVSV)")

print(p_genewise)

save_pipeline_comparison_plot(p_genewise, figure_path, "comparison_selection_genewise_prevac_rvsv.pdf")

gc()
rm(list = ls())
