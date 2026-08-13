# analysis/pipeline_comparisons/prevac_rvsv/03_compare_selection.R
#
# Compare feature-selection choices for predicting PREVAC rVSV (+ placebo)
# day-180 antibody titer from day-7 gene expression, against the reference
# pipeline (mean gene-set aggregation; see
# R/pipeline_defaults.R::reference_pipeline_params()). This is the ONLY
# feature-selection comparison treated as a "reference" comparison for this
# dataset - see 03b_compare_selection_genewise.R for the supplementary,
# gene-wise (raw gene, no aggregation) comparison, which uses a different
# reference pipeline and is not used to characterise "the" reference
# approach for this dataset (e.g. in
# analysis/pipeline_comparisons/visualize_comparisons.R's baseline/reference
# context figure).
#
# Selection methods score each aggregated feature: variance filtering,
# absolute correlation (Spearman/Pearson), univariate regression-based
# screening ("relative gain"), and dearseq's "classic" mode at the geneset
# level (`dearseq_level = "geneset"`).
#
# PREVAC has a placebo arm, so dearseq is used here in its "classic" mode:
# `dearseq_mode = "classic"` (the default) contrasts `treatment == 1`
# (rVSV) vs. `treatment == 0` (placebo). This needs no
# `individual_id`/`timepoint` or predictomics' paired row-discard handling
# (see sdy1276_tiv/03_compare_selection.R for that machinery) - "classic"
# mode operates on the single (day-7-only) dataset directly, with every row
# retained for modelling. `treatment` is passed only for this screening
# contrast; `treatment_predictor = FALSE` (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# itself.
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
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_rvsv")

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
    relative_gain_inner_folds = 5, relative_gain_metric = "rmse"
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

save_comparison_metrics(res_geneset, dataset = "prevac_rvsv", category = "selection_geneset")

p_geneset <- plot(res_geneset, metric = "R2") +
  ggtitle("Pipeline comparison: geneset-level feature selection (PREVAC rVSV)")

print(p_geneset)

save_pipeline_comparison_plot(p_geneset, figure_path, "comparison_selection_geneset_prevac_rvsv.pdf")

gc()
rm(list = ls())
