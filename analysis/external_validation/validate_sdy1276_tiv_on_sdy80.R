# analysis/external_validation/validate_sdy1276_tiv_on_sdy80.R
#
# Validate SDY1276 (TIV)'s best pipeline (see
# analysis/pipeline_comparisons/sdy1276_tiv/05_find_best_model.R) on SDY80,
# an independent influenza (IN) vaccine study with the same paired
# baseline-vs-post structure and no placebo arm.
#
# Per R/external_validation.R's header: this does NOT apply the already-fit
# SDY1276 model directly to SDY80 (expression scales differ too much across
# studies for that to be meaningful), and does NOT re-run feature selection
# on SDY80 either. Instead: the winning engineering and model choices are
# carried over unchanged, the exact set of features SDY1276's winning
# selection step chose (recovered via
# get_discovery_selected_features()) is fixed and applied directly to SDY80,
# and cross-validated performance is re-measured from scratch on SDY80 using
# that fixed panel. If SDY1276's winning pipeline had no selection step at
# all, every engineered feature is carried over unchanged instead.
#
# ASSUMPTIONS ABOUT SDY80 (data/ is gitignored, so these should be checked
# against the actual merged data the first time this script is run):
#   - study_vaccine == "SDY80-Influenza (IN)" identifies SDY80's TIV
#     recipients, mirroring SDY1276's own "SDY1276-Influenza (IN)" - both
#     studies go through the same generic IS2 preprocessing pipeline (see
#     analysis/preprocessing/preprocessing_is2_clinical.R), where
#     study_vaccine = "<study_accession>-<vaccine_name>" and vaccine_name is
#     derived from pathogen + vaccine_type ("Influenza (IN)" for inactivated
#     influenza), so this should hold for any IS2 study administering TIV.
#   - Gene expression at "P+1D" (day 1 post-vaccination) and the antibody
#     response `ab_p_28` (day 28) are both available for SDY80, mirroring
#     SDY1276's own design (see sdy1276_tiv/01_prepare_data.R) - both derive
#     from the same generic IS2 immune-response processing (see
#     analysis/preprocessing/preprocessing_is2_immuneresponse.R), which does
#     not hardcode timepoints/outcomes per study.
#
# Unlike sdy1276_tiv/01_prepare_data.R, no paired (baseline + post) dataset
# is built here - fixing the feature panel upfront removes the only reason
# paired-mode RISE/dearseq (needing baseline+post pairing) was ever used in
# the discovery search, so only the single post-vaccination timepoint is
# needed for SDY80, exactly as for the PREVAC datasets' own searches.
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/00_master.R (through at
# least 05_find_best_model.R) first, so that
# output/results/best_pipeline_sdy1276_tiv.rds exists.

library(dplyr)
library(fs)
library(predictomics)
library(ggplot2)
library(patchwork)

source(fs::path("R", "data_io.R"))
source(fs::path("R", "gene_columns.R"))
source(fs::path("R", "build_dataset.R"))
source(fs::path("R", "pipeline_defaults.R"))
source(fs::path("R", "metrics_labels.R"))
source(fs::path("R", "external_validation.R"))
source(fs::path("R", "panel_helpers.R"))

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

discovery <- readRDS(fs::path(results_dir, "sdy1276_tiv_analysis_data.rds"))
best <- readRDS(fs::path(results_dir, "best_pipeline_sdy1276_tiv.rds"))

# --- Recover the fixed feature panel from SDY1276 (the discovery dataset) --
# SDY1276 is single-arm (no placebo), so treatment is NULL here - matching
# how find_best_pipeline() was itself called for this dataset.
selected_features <- get_discovery_selected_features(
  best,
  X = discovery$single$X, Y = discovery$single$Y,
  covariates = discovery$single$covariates, treatment = NULL
)

# --- Build the SDY80 validation dataset -------------------------------------
df_merged_all <- load_merged_data()

sdy80_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = "SDY80-Influenza (IN)",
  timepoint = "P+1D",
  response_col = "ab_p_28",
  covariate_names = default_covariates,
  log_response = TRUE
)

X_validation <- sdy80_single$X

# Built AFTER the validation dataset so gene-set restriction (for
# aggregating engineering) can be checked against SDY80's own available
# gene panel, not just the discovery-selected feature list - see
# restrict_engineering_for_validation()'s docs for why this matters even
# when there's no upstream selection step to intersect against.
restricted <- restrict_engineering_for_validation(
  best, selected_features, discovery$genesets,
  validation_gene_names = colnames(X_validation)
)

if (!is.null(restricted$fixed_gene_features)) {
  available <- intersect(restricted$fixed_gene_features, colnames(X_validation))
  missing   <- setdiff(restricted$fixed_gene_features, colnames(X_validation))
  if (length(missing) > 0) {
    message(sprintf(
      "[validation] %d of %d discovery-selected genes are not available in SDY80 and will be dropped: %s",
      length(missing), length(restricted$fixed_gene_features), paste(missing, collapse = ", ")
    ))
  }
  X_validation <- X_validation[, available, drop = FALSE]
}

cat(sprintf(
  "SDY1276 (TIV) -> SDY80 validation: %d participants, %d features.\n",
  length(sdy80_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = sdy80_single$Y, covariates = sdy80_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(results_dir, "validation_fit_sdy1276_tiv_on_sdy80.rds"))

figure_path <- fs::path("output", "figures", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  model_params = best$model_params,
  title = "SDY1276 (TIV) best pipeline: validated on SDY80",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_sdy1276_tiv_on_sdy80.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
