# analysis/pipeline_comparisons/sdy1276_tiv/02_compare_engineering.R
#
# Compare feature-engineering choices for predicting SDY1276 (TIV) day-28
# antibody titer from day-1 gene expression, against the reference pipeline
# (z-scored gene-level transform, mean gene-set aggregation, elastic net
# model; see R/pipeline_defaults.R).
#
# Two comparisons are run, matching the two axes of feature engineering
# considered in Chapter 5, Table 5.2, and reported as two separate figures:
#
#   1. Gene-set-level aggregation (gene-wise/no aggregation, mean, median,
#      max, 1st PC, GSVA, ssGSEA of Blood Transcriptional Modules), each
#      applied on top of a z-scored gene-level transform, using the single
#      (day-1) expression timepoint. Mean aggregation is the reference and
#      so is also included here as an explicit option, to show it on the
#      same comparison plot as the alternatives.
#   2. Gene-level transformation only: z-score vs. individual-level
#      fold-change from baseline. The fold-change option needs *both*
#      baseline (day 0) and day-1 expression per participant, so this
#      comparison uses the paired dataset built in 01_prepare_data.R and is
#      necessarily a separate `compare_pipelines()` call/figure from (1).
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
single <- analysis_data$single
paired <- analysis_data$paired
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")

reference_params <- reference_pipeline_params(genesets)

## ---- 1. Gene-set aggregation (single timepoint) ---------------------------

aggregation_options <- list(
  "Gene-wise"   = list(method = "engineer", col_transform = "z", genesets = NULL, agg_method = "mean"),
  "Gene-set: median" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "median"),
  "Gene-set: max"    = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "max"),
  "Gene-set: 1st PC" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "pc1"),
  "Gene-set: GSVA"   = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "gsva", gsva_min_size = 2),
  "Gene-set: ssGSEA" = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "ssgsea", ssgsea_min_size = 2)
)

res_aggregation <- run_pipeline_comparison(
  X = single$X, Y = single$Y, covariates = single$covariates,
  option_type = "engineering", option_choices = aggregation_options,
  reference_params = reference_params
)

save_comparison_metrics(
  res_aggregation, dataset = "sdy1276_tiv", category = "engineering",
  label = "engineering_aggregation"
)

gc()

p_aggregation <- plot(res_aggregation, metric = "R2") +
  ggtitle("Pipeline comparison: gene-set aggregation (SDY1276, TIV)")

print(p_aggregation)

save_pipeline_comparison_plot(
  p_aggregation, figure_path, "comparison_engineering_aggregation_sdy1276_tiv.pdf"
)

## ---- 2. Gene-level transformation (paired baseline + day-1 timepoints) ----

gene_level_options <- list(
  "Gene-level: z-score" = list(
    method = "engineer", col_transform = "z", gene_level_fc = FALSE,
    genesets = NULL, agg_method = "mean"
  ),
  "Gene-level: fold-change" = list(
    method = "engineer", col_transform = "none", gene_level_fc = TRUE,
    genesets = NULL, agg_method = "mean"
  )
)

# The reference pipeline is re-evaluated here on the paired dataset (rather
# than reusing res_aggregation's reference fit above), so that all three
# gene-level options are compared on an identical dataset and CV split.
res_gene_level <- run_pipeline_comparison(
  X = paired$X, Y = paired$Y, covariates = paired$covariates,
  timepoint = paired$timepoint, individual_id = paired$participant_id,
  option_type = "engineering", option_choices = gene_level_options,
  reference_params = reference_params
)

save_comparison_metrics(
  res_gene_level, dataset = "sdy1276_tiv", category = "engineering",
  label = "engineering_gene_level"
)

p_gene_level <- plot(res_gene_level, metric = "R2") +
  ggtitle("Pipeline comparison: gene-level transformation (SDY1276, TIV)")

print(p_gene_level)

save_pipeline_comparison_plot(
  p_gene_level, figure_path, "comparison_engineering_gene_level_sdy1276_tiv.pdf"
)

gc()
rm(list = ls())
