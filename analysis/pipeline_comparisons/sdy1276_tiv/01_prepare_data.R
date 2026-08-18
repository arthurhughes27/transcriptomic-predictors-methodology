# analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R
#
# Build and cache the analysis-ready datasets used by the pipeline-comparison
# scripts in this folder (02-04). Centralising data preparation here means
# all three comparisons (feature engineering, feature selection, model
# choice) analyse an identical dataset, and the filtering logic used to build
# it is written and reviewed in one place only.
#
# Dataset: SDY1276 (TIV / Influenza (IN)), predicting the mean neutralising
# antibody titer across the three vaccine strains at day 28 (`ab_p_28`, see
# analysis/preprocessing/preprocessing_is2_immuneresponse.R for how this
# mean is computed) from whole-blood gene expression at day 1
# post-vaccination, adjusting for age, sex, and race (Chapter 5, Section 5.4).
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

study_vaccine_groups <- "SDY1276-Influenza (IN)"

# Single-timepoint dataset (day-1 expression only): used for the feature
# selection and model-choice comparisons, and for every feature-engineering
# option except the individual-level fold-change transform.
tiv_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = study_vaccine_groups,
  timepoint = "P+1D",
  response_col = "ab_p_28",
  covariate_names = default_covariates,
  log_response = TRUE
)

# Paired baseline (day 0) + post-vaccination (day 1) dataset: needed only to
# compute the individual-level fold-change gene-level transformation, which
# requires both timepoints per participant.
tiv_paired <- build_paired_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = study_vaccine_groups,
  baseline_timepoint = "P+0D",
  post_timepoint = "P+1D",
  response_col = "ab_p_28",
  covariate_names = default_covariates,
  log_response = TRUE
)

cat(sprintf(
  "SDY1276 TIV: %d participants (single timepoint), %d participants (paired timepoints), %d gene features.\n",
  length(tiv_single$participant_id),
  length(unique(tiv_paired$participant_id)),
  length(tiv_single$gene_names)
))

saveRDS(
  list(
    single = tiv_single,
    paired = tiv_paired,
    genesets = genesets,
    covariate_names = default_covariates
  ),
  file = fs::path(results_dir, "sdy1276_tiv_analysis_data.rds")
)

rm(list = ls())
