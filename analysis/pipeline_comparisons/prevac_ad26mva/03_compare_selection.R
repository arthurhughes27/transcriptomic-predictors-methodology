# analysis/pipeline_comparisons/prevac_ad26mva/03_compare_selection.R
#
# Compare feature-selection choices for predicting PREVAC Ad26/MVA
# (+ placebo) day-180 antibody titer from day-7 gene expression.
#
# Two comparisons are run, reported as two separate figures, split by
# engineering scale (as in sdy1276_tiv/03_compare_selection.R): RISE always
# screens the raw gene-level matrix and is never compatible with gene-set
# aggregation; dearseq's gene-level mode (`dearseq_level = "gene"`) is
# likewise never compatible with gene-set aggregation. So a
# geneset-aggregated reference can only host dearseq's geneset-level mode,
# while a raw-gene reference can host both RISE and dearseq's gene-level
# mode.
#
#   1. Geneset-level: reference is the usual mean-BTM-aggregated pipeline
#      (`reference_pipeline_params()`). Selection methods score each
#      aggregated feature: variance filtering, absolute correlation
#      (Spearman/Pearson), univariate regression-based screening ("relative
#      gain"), and dearseq's "classic" mode at the geneset level
#      (`dearseq_level = "geneset"`).
#   2. Gene-wise: reference is z-score only, no gene-set aggregation
#      (`raw_gene_reference_params()`). The same selection methods are
#      compared, plus RISE and dearseq's "classic" mode at the gene level.
#
# Unlike SDY1276 (TIV), PREVAC has a placebo arm, so RISE and dearseq are
# used here in their "classic" mode rather than "paired": RISE contrasts
# `treatment == 1` (Ad26/MVA) vs. `treatment == 0` (placebo); dearseq's
# default `dearseq_mode = "classic"` does the same. Neither needs
# `individual_id`/`timepoint` or predictomics' paired row-discard handling
# (see sdy1276_tiv/03_compare_selection.R for that machinery) - "classic"
# mode operates on the single (day-7-only) dataset directly, with every row
# retained for modelling. `treatment` is passed only for this screening
# contrast; `treatment_predictor = FALSE` (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# itself.
#
# Run analysis/pipeline_comparisons/prevac_ad26mva/01_prepare_data.R first.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)

source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "run_comparison.R"))
source(fs::path("R", "plotting.R"))
source(fs::path("R", "metrics_io.R"))

analysis_data <- readRDS(fs::path("data", "derived", "prevac_ad26mva_analysis_data.rds"))
single <- analysis_data$single
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_ad26mva")

## ---- 1. Geneset-level feature selection ------------------------------------

reference_params_geneset <- reference_pipeline_params(genesets)

option_choices_geneset <- list(
  "Variance (top 25)" = list(method = "variance", top_n = 25),
  "Variance (top 100)" = list(method = "variance", top_n = 100),
  "Correlation - Spearman (top 25)"   = list(method = "spearman", top_n = 25),
  "Correlation - Spearman (|r| > 0.5)" = list(method = "spearman", threshold = 0.5),
  "Correlation - Pearson (top 25)"    = list(method = "pearson", top_n = 25),
  "Correlation - Pearson (|r| > 0.5)"  = list(method = "pearson", threshold = 0.5),
  "Univariate regression screening (threshold = 0)" = list(
    method = "relative_gain", threshold = 0,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse", relative_gain_seed = 12345
  ),
  "Univariate regression screening (threshold = 0.1)" = list(
    method = "relative_gain", threshold = 0.1,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse", relative_gain_seed = 12345
  ),
  "Dearseq (geneset, alpha = 0.05)" = list(
    method = "dearseq", dearseq_mode = "classic", dearseq_level = "geneset",
    genesets = genesets, threshold = 0.05
  )
)

res_geneset <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  treatment = single$treatment,
  option_type = "selection", option_choices = option_choices_geneset,
  reference_params = reference_params_geneset
)

save_comparison_metrics(res_geneset, dataset = "prevac_ad26mva", category = "selection_geneset")

p_geneset <- plot(res_geneset, metric = "R2") +
  ggtitle("Pipeline comparison: geneset-level feature selection (PREVAC Ad26/MVA)")

print(p_geneset)

save_pipeline_comparison_plot(p_geneset, figure_path, "comparison_selection_geneset_prevac_ad26mva.pdf")

gc()

## ---- 2. Gene-wise (raw) feature selection ----------------------------------

reference_params_genewise <- raw_gene_reference_params()

option_choices_genewise <- list(
  "Variance (top 500)" = list(method = "variance", top_n = 500),
  "Correlation - Spearman (top 500)"  = list(method = "spearman", top_n = 500),
  "Correlation - Spearman (|r| > 0.5)"  = list(method = "spearman", threshold = 0.5),
  "Correlation - Pearson (top 500)"  = list(method = "pearson", top_n = 500),
  "Correlation - Pearson (|r| > 0.5)"   = list(method = "pearson", threshold = 0.5),
  "Univariate regression screening (threshold = 0)" = list(
    method = "relative_gain", threshold = 0,
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse", relative_gain_seed = 12345
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

res_genewise <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  treatment = single$treatment,
  option_type = "selection", option_choices = option_choices_genewise,
  reference_params = reference_params_genewise
)

save_comparison_metrics(res_genewise, dataset = "prevac_ad26mva", category = "selection_genewise")

p_genewise <- plot(res_genewise, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise feature selection (PREVAC Ad26/MVA)")

print(p_genewise)

save_pipeline_comparison_plot(p_genewise, figure_path, "comparison_selection_genewise_prevac_ad26mva.pdf")

gc()
rm(list = ls())
