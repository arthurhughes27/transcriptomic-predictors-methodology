# analysis/pipeline_comparisons/prevac_ad26mva/01_prepare_data.R
#
# Build and cache the analysis-ready dataset used by the pipeline-comparison
# scripts in this folder (02-04).
#
# Dataset: PREVAC Ad26/MVA vs. placebo, predicting the anti-Ebola-glycoprotein
# binding antibody titer at day 180 (`ab_p_180`) from whole-blood gene
# expression at day 7 post-vaccination, adjusting for age, sex, and race
# (Chapter 5, Section 5.4).
#
# The placebo arm is included in every prediction task here (extra examples
# of low/no response), but *not* as a model predictor: a binary `treatment`
# indicator (1 = Ad26/MVA, 0 = placebo) is built for use as the contrast in
# "classic" RISE/dearseq feature selection (see 03_compare_selection.R),
# while `treatment_predictor = FALSE` is used throughout (02-04) so
# treatment itself never enters the model as a feature.
#
# The placebo arm already gives "classic" RISE/dearseq a treatment contrast
# to screen on, so paired mode (baseline vs. post) isn't needed for the
# MAIN (geneset-level) comparisons in this folder - only the single
# post-vaccination timepoint is used there.
#
# A paired (baseline P+0D + post P+7D) dataset IS also built below, purely
# to support the individual-level gene-level fold-change transformation
# offered as a SUPPLEMENTARY, gene-wise engineering option (see
# analysis/supplementary/prevac_ad26mva/01_compare_engineering_genewise.R) -
# mirroring sdy1276_tiv/01_prepare_data.R's own paired dataset, built for
# the same reason.
#
# NOTE: PREVAC's `race` covariate is recorded as a constant "Unknown" for
# every participant (see
# analysis/preprocessing/preprocessing_clinical_harmonisation.R). It's
# still included here for consistency with the covariate set used elsewhere
# in this chapter (age, sex, race), but carries no information for this
# dataset - worth checking it doesn't cause a degenerate design matrix when
# actually run.
#
# Run this script once before running 02_compare_engineering.R,
# 03_compare_selection.R, or 04_compare_model.R.

library(dplyr)
library(fs)

source(fs::path("R", "data_io.R"))
source(fs::path("R", "gene_columns.R"))
source(fs::path("R", "build_dataset.R"))
source(fs::path("R", "pipeline_defaults.R"))

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

df_merged_all <- load_merged_data()
genesets <- load_genesets()

treatment_arm <- "prevac-Ad26MVA"
study_vaccine_groups <- c(treatment_arm, "prevac-placebo")

prevac_ad26mva_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = study_vaccine_groups,
  timepoint = "P+7D",
  response_col = "ab_p_180",
  covariate_names = default_covariates,
  treatment_arm = treatment_arm,
  log_response = TRUE
)

prevac_ad26mva_paired <- build_paired_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = study_vaccine_groups,
  baseline_timepoint = "P+0D",
  post_timepoint = "P+7D",
  response_col = "ab_p_180",
  covariate_names = default_covariates,
  log_response = TRUE
)

cat(sprintf(
  "PREVAC Ad26/MVA (+ placebo): %d participants (%d Ad26/MVA, %d placebo), %d gene features (single); %d participants (paired timepoints).\n",
  length(prevac_ad26mva_single$participant_id),
  sum(prevac_ad26mva_single$treatment == 1),
  sum(prevac_ad26mva_single$treatment == 0),
  length(prevac_ad26mva_single$gene_names),
  length(unique(prevac_ad26mva_paired$participant_id))
))

saveRDS(
  list(
    single = prevac_ad26mva_single,
    paired = prevac_ad26mva_paired,
    genesets = genesets,
    covariate_names = default_covariates
  ),
  file = fs::path(results_dir, "prevac_ad26mva_analysis_data.rds")
)

rm(list = ls())
