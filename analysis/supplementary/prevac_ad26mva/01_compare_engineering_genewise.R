# analysis/supplementary/prevac_ad26mva/01_compare_engineering_genewise.R
#
# SUPPLEMENTARY comparison: gene-level (no gene-set aggregation) engineering
# choices for predicting PREVAC Ad26/MVA (+ placebo) day-180 antibody titer
# from day-7 gene expression - individual-level fold-change from baseline
# vs. the gene-wise reference pipeline (z-scored gene-level transform,
# variance top-7,500 pre-filter, elastic net model; see
# R/pipeline_defaults.R::raw_gene_reference_params()).
#
# Neither "Gene-level: z-score" nor "Gene-level: none" (raw, untransformed
# expression) is offered as an explicit option here: with
# model_params$scale = TRUE always on in this repo, caret centres and
# scales every predictor before fitting regardless of col_transform, so a
# column-wise linear transform like z-scoring or "none" becomes numerically
# identical to the reference by the time the model sees it - there is no
# distinct "no transform" pipeline left to compare. Fold-change survives as
# a genuinely different (nonlinear, individual-level baseline-relative)
# transform, so it remains the sole engineering option here.
#
# This is NOT sourced by
# analysis/pipeline_comparisons/prevac_ad26mva/00_master.R and does not run
# as part of analysis/master_analysis.R - run it (and its siblings
# 02_compare_selection_genewise.R/03_compare_model_genewise.R) via
# analysis/supplementary/prevac_ad26mva/00_master.R, or
# analysis/supplementary/master_supplementary.R for all three datasets.
#
# This is NOT treated as a/the reference approach for this dataset - see
# 02_compare_engineering.R (this folder's main, geneset-only engineering
# comparison) for that. This comparison exists to characterise gene-level
# transform choices in their own right, rather than to compete with the
# geneset-only comparison for "the" dataset-level baseline.
#
# The individual-level fold-change transformation needs *both* baseline
# (P+0D) and day-7 (P+7D) expression per participant, so - like
# analysis/supplementary/sdy1276_tiv/01_compare_engineering_genewise.R -
# this comparison uses the PAIRED dataset built in 01_prepare_data.R. Per
# predictomics' paired row-discard parity handling, the fold-change option
# (which discards pre-vaccination rows internally) models on
# post-vaccination (day-7) rows only, while the reference is automatically
# restricted to post-vaccination rows first - keeping the option and
# reference compared on an identical, independent (one row per participant)
# sample.
#
# treatment is not passed to this comparison - it's not used by any
# engineering option, and treatment_predictor = FALSE (the default in
# R/run_comparison.R::run_pipeline_comparison()) keeps it out of the model
# regardless.
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

analysis_data <- readRDS(fs::path("output", "results", "prevac_ad26mva_analysis_data.rds"))
paired <- analysis_data$paired

figure_path <- fs::path("output", "figures", "supplementary", "prevac_ad26mva")

reference_params <- raw_gene_reference_params()

option_choices <- list(
  "Gene-level: fold-change" = list(method = "engineer", col_transform = "none", gene_level_fc = TRUE, genesets = NULL, agg_method = "mean")
)

res <- run_or_load_comparison(
  dataset = "prevac_ad26mva", label = "engineering_genewise",
  metrics_subdir = "supplementary",
  X = paired$X, Y = paired$Y, covariates = paired$covariates,
  individual_id = paired$participant_id, timepoint = paired$timepoint,
  option_type = "engineering", option_choices = option_choices,
  reference_params = reference_params
)

save_comparison_metrics(res, dataset = "prevac_ad26mva", category = "engineering_genewise", subdir = "supplementary")

p <- plot(res, metric = "R2") +
  ggtitle("Pipeline comparison: gene-wise transformation (PREVAC Ad26/MVA)")

print(p)

save_pipeline_comparison_plot(p, figure_path, "comparison_engineering_genewise_prevac_ad26mva.pdf")

gc()
rm(list = ls())
