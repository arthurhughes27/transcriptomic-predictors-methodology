# analysis/supplementary/sdy1276_tiv/01_compare_engineering_genewise.R
#
# NOTE: this used to live at analysis/pipeline_comparisons/sdy1276_tiv/
# 02b_compare_engineering_genewise.R - moved here (and the "02b" cut to "01")
# so this SUPPLEMENTARY, gene-wise comparison is entirely separate from the
# main analysis: it is not sourced by
# analysis/pipeline_comparisons/sdy1276_tiv/00_master.R and does not run as part of
# analysis/master_analysis.R. Run it (and its sibling
# 02_compare_selection_genewise.R) via analysis/supplementary/sdy1276_tiv/00_master.R,
# or analysis/supplementary/master_supplementary.R for all three datasets.
#
# SUPPLEMENTARY comparison: gene-level (no gene-set aggregation) engineering
# choices for predicting SDY1276 (TIV) day-28 antibody titer from day-1
# gene expression - no transform, z-score, and individual-level fold-change
# from baseline - against the same reference pipeline used throughout this
# folder (z-scored gene-level transform, mean gene-set aggregation, elastic
# net model; see R/pipeline_defaults.R::reference_pipeline_params()).
#
# This is NOT treated as a/the reference approach for this dataset - see
# 02_compare_engineering.R (this folder's main, geneset-only engineering
# comparison) for that. This comparison exists to characterise gene-level
# transform choices in their own right, alongside gene-set aggregation as
# an alternative, rather than to compete with the geneset-only comparison
# for "the" dataset-level baseline (see this folder's README.md for the
# geneset/gene-wise separation).
#
# The individual-level fold-change transformation needs *both* baseline
# (day 0) and day-1 expression per participant, so - unlike
# 02_compare_engineering.R - this comparison uses the PAIRED dataset built
# in 01_prepare_data.R. Per predictomics' paired row-discard parity
# handling, the fold-change option (which discards pre-treatment rows
# internally) models on post-treatment (day-1) rows only, while every other
# pipeline in the same call (including "no transform"/"z-score" and the
# reference) is automatically restricted to post-treatment rows first -
# keeping every option compared on an identical, independent (one row per
# participant) sample.
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

analysis_data <- readRDS(fs::path("output", "results", "sdy1276_tiv_analysis_data.rds"))
paired <- analysis_data$paired
genesets <- analysis_data$genesets

figure_path <- fs::path("output", "figures", "supplementary", "sdy1276_tiv")

reference_params <- reference_pipeline_params(genesets)

option_choices <- list(
  "Gene-level: none"        = list(method = "engineer", col_transform = "none", gene_level_fc = FALSE, genesets = NULL, agg_method = "mean"),
  "Gene-level: z-score"     = list(method = "engineer", col_transform = "z",    gene_level_fc = FALSE, genesets = NULL, agg_method = "mean"),
  "Gene-level: fold-change" = list(method = "engineer", col_transform = "none", gene_level_fc = TRUE,  genesets = NULL, agg_method = "mean")
)

res <- run_or_load_comparison(
  dataset = "sdy1276_tiv", label = "engineering_genewise",
  metrics_subdir = "supplementary",
  X = paired$X, Y = paired$Y, covariates = paired$covariates,
  individual_id = paired$participant_id, timepoint = paired$timepoint,
  option_type = "engineering", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "sdy1276_tiv", category = "engineering_genewise", subdir = "supplementary")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise transformation (SDY1276, TIV)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_engineering_genewise_sdy1276_tiv.pdf")

gc()
rm(list = ls())
