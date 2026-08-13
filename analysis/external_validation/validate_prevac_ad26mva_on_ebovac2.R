# analysis/external_validation/validate_prevac_ad26mva_on_ebovac2.R
#
# Validate PREVAC Ad26/MVA's best pipeline (see
# analysis/pipeline_comparisons/prevac_ad26mva/05_find_best_model.R) on
# EBOVAC2, an independent Ad26/MVA Ebola vaccine study.
#
# Per R/external_validation.R's header: this does NOT apply the already-fit
# PREVAC model directly to EBOVAC2, and does NOT re-run feature selection on
# EBOVAC2 either. Instead: the winning engineering and model choices are
# carried over unchanged, the exact set of features PREVAC Ad26/MVA's
# winning selection step chose (recovered via
# get_discovery_selected_features(), needing `treatment` there only to
# reproduce that upfront selection exactly as the discovery search made it -
# e.g. if it picked "classic"-mode RISE/dearseq, screening on the
# Ad26/MVA-vs-placebo contrast) is fixed and applied directly to EBOVAC2, and
# cross-validated performance is re-measured from scratch on EBOVAC2 using
# that fixed panel. If PREVAC Ad26/MVA's winning pipeline had no selection
# step at all, every engineered feature is carried over unchanged instead.
#
# EBOVAC2 has no placebo arm to fit here (only the active Ad26/MVA arm,
# "ebovac2-Ad26MVA", is used - see analysis/preprocessing/
# preprocessing_clinical_harmonisation.R), and fixing the feature panel
# upfront removes the only reason a placebo contrast was needed in the first
# place, so - unlike prevac_ad26mva/01_prepare_data.R - no `treatment`
# vector is built or needed for the validation fit itself.
#
# ASSUMPTIONS ABOUT EBOVAC2 (data/ is gitignored, so these should be checked
# against the actual merged data the first time this script is run):
#   - Gene expression at "P+7D" (day 7 post-vaccination) is available for
#     EBOVAC2, mirroring PREVAC Ad26/MVA's own day-7 design (see
#     prevac_ad26mva/01_prepare_data.R) - EBOVAC2's harmonised timepoints
#     include "P+7D" (see analysis/preprocessing/preprocessing_ebovac2.R).
#   - EBOVAC2 lacks a day-180 antibody measurement (ab_p_180 is NA for every
#     EBOVAC2 row - see analysis/preprocessing/
#     preprocessing_clinical_harmonisation.R) but does have a day-364/365
#     measurement (ab_p_365), its latest available timepoint - used here as
#     the validation response in place of PREVAC's ab_p_180.
#
# Run analysis/pipeline_comparisons/prevac_ad26mva/00_master.R (through at
# least 05_find_best_model.R) first, so that
# data/derived/best_pipeline_prevac_ad26mva.rds exists.

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

derived_data_dir <- fs::path("data", "derived")
fs::dir_create(derived_data_dir)

discovery <- readRDS(fs::path(derived_data_dir, "prevac_ad26mva_analysis_data.rds"))
best <- readRDS(fs::path(derived_data_dir, "best_pipeline_prevac_ad26mva.rds"))

# --- Recover the fixed feature panel from PREVAC Ad26/MVA (discovery) ------
selected_features <- get_discovery_selected_features(
  best,
  X = discovery$single$X, Y = discovery$single$Y,
  covariates = discovery$single$covariates, treatment = discovery$single$treatment
)

restricted <- restrict_engineering_for_validation(best, selected_features, discovery$genesets)

# --- Build the EBOVAC2 validation dataset -----------------------------------
df_merged_all <- load_merged_data()

ebovac2_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = "ebovac2-Ad26MVA",
  timepoint = "P+7D",
  response_col = "ab_p_365",
  covariate_names = default_covariates,
  log_response = TRUE
)

X_validation <- ebovac2_single$X

if (!is.null(restricted$fixed_gene_features)) {
  available <- intersect(restricted$fixed_gene_features, colnames(X_validation))
  missing   <- setdiff(restricted$fixed_gene_features, colnames(X_validation))
  if (length(missing) > 0) {
    message(sprintf(
      "[validation] %d of %d discovery-selected genes are not available in EBOVAC2 and will be dropped: %s",
      length(missing), length(restricted$fixed_gene_features), paste(missing, collapse = ", ")
    ))
  }
  X_validation <- X_validation[, available, drop = FALSE]
}

cat(sprintf(
  "PREVAC Ad26/MVA -> EBOVAC2 validation: %d participants, %d features.\n",
  length(ebovac2_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = ebovac2_single$Y, covariates = ebovac2_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(derived_data_dir, "validation_fit_prevac_ad26mva_on_ebovac2.rds"))

figure_path <- fs::path("output", "figures", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  title = "PREVAC Ad26/MVA best pipeline: validated on EBOVAC2",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_prevac_ad26mva_on_ebovac2.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
