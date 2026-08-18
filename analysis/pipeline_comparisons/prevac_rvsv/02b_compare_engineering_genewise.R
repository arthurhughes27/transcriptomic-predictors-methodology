# analysis/pipeline_comparisons/prevac_rvsv/02b_compare_engineering_genewise.R
#
# SUPPLEMENTARY comparison: gene-level (no gene-set aggregation) engineering
# choices for predicting PREVAC rVSV (+ placebo) day-180 antibody titer
# from day-7 gene expression - no transform vs. z-score - against the same
# reference pipeline used throughout this folder (z-scored gene-level
# transform, mean gene-set aggregation, elastic net model; see
# R/pipeline_defaults.R::reference_pipeline_params()).
#
# This is NOT treated as a/the reference approach for this dataset - see
# 02_compare_engineering.R (this folder's main, geneset-only engineering
# comparison) for that. This comparison exists to characterise gene-level
# transform choices in their own right, alongside gene-set aggregation as
# an alternative, rather than to compete with the geneset-only comparison
# for "the" dataset-level baseline (see this folder's README.md for the
# geneset/gene-wise separation).
#
# Unlike sdy1276_tiv/02b_compare_engineering_genewise.R, individual-level
# fold-change isn't compared here: 01_prepare_data.R doesn't build a paired
# (baseline + day-7) dataset for this study (see that script's header for
# why), so this uses the single (day-7-only) dataset throughout.
#
# treatment is not passed to this comparison - it's not used by any
# engineering option, and treatment_predictor = FALSE (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# regardless.
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

analysis_data <- readRDS(fs::path("output", "results", "prevac_rvsv_analysis_data.rds"))
single <- analysis_data$single
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "prevac_rvsv")

reference_params <- reference_pipeline_params(genesets)

option_choices <- list(
  "Gene-level: none"    = list(method = "engineer", col_transform = "none", genesets = NULL, agg_method = "mean"),
  "Gene-level: z-score" = list(method = "engineer", col_transform = "z",    genesets = NULL, agg_method = "mean")
)

res <- run_or_load_comparison(
  dataset = "prevac_rvsv", label = "engineering_genewise",
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "engineering", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "prevac_rvsv", category = "engineering_genewise")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise transformation (PREVAC rVSV)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_engineering_genewise_prevac_rvsv.pdf")

gc()
rm(list = ls())
