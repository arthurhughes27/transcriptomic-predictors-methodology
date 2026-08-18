# analysis/external_validation/validate_prevac_rvsv_on_hamburg.R
#
# Validate PREVAC rVSV's best pipeline (see
# analysis/pipeline_comparisons/prevac_rvsv/05_find_best_model.R) on
# Hamburg, an independent rVSV Ebola vaccine study.
#
# Per R/external_validation.R's header: this does NOT apply the already-fit
# PREVAC model directly to Hamburg, and does NOT re-run feature selection on
# Hamburg either. Instead: the winning engineering and model choices are
# carried over unchanged, the exact set of features PREVAC rVSV's winning
# selection step chose (recovered via get_discovery_selected_features(),
# needing `treatment` there only to reproduce that upfront selection exactly
# as the discovery search made it - e.g. if it picked "classic"-mode
# RISE/dearseq, screening on the rVSV-vs-placebo contrast) is fixed and
# applied directly to Hamburg, and cross-validated performance is
# re-measured from scratch on Hamburg using that fixed panel. If PREVAC
# rVSV's winning pipeline had no selection step at all, every engineered
# feature is carried over unchanged instead.
#
# Hamburg has no placebo arm at all (every participant is "hamburg-rVSV" -
# see analysis/preprocessing/preprocessing_clinical_harmonisation.R), and
# fixing the feature panel upfront removes the only reason a placebo
# contrast was needed in the first place, so - unlike
# prevac_rvsv/01_prepare_data.R - no `treatment` vector is built or needed
# for the validation fit itself.
#
# ASSUMPTIONS ABOUT HAMBURG (data/ is gitignored, so these should be checked
# against the actual merged data the first time this script is run):
#   - Gene expression at "P+7D" (day 7 post-vaccination) is available for
#     Hamburg, mirroring PREVAC rVSV's own day-7 design (see
#     prevac_rvsv/01_prepare_data.R) - Hamburg's harmonised timepoints
#     include "d7" -> "P+7D" (see analysis/preprocessing/
#     preprocessing_clinical_harmonisation.R).
#   - Hamburg lacks a day-364/365 antibody measurement (ab_p_365 is NA for
#     every Hamburg row) but does have a day-180 measurement (ab_p_180) -
#     used here as the validation response, matching PREVAC's own ab_p_180.
#
# Run analysis/pipeline_comparisons/prevac_rvsv/00_master.R (through at
# least 05_find_best_model.R) first, so that
# output/results/best_pipeline_prevac_rvsv.rds exists.

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

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

discovery <- readRDS(fs::path(results_dir, "prevac_rvsv_analysis_data.rds"))
best <- readRDS(fs::path(results_dir, "best_pipeline_prevac_rvsv.rds"))

# --- Recover the fixed feature panel from PREVAC rVSV (discovery) ----------
selected_features <- get_discovery_selected_features(
  best,
  X = discovery$single$X, Y = discovery$single$Y,
  covariates = discovery$single$covariates, treatment = discovery$single$treatment
)

restricted <- restrict_engineering_for_validation(best, selected_features, discovery$genesets)

# --- Build the Hamburg validation dataset -----------------------------------
df_merged_all <- load_merged_data()

hamburg_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = "hamburg-rVSV",
  timepoint = "P+7D",
  response_col = "ab_p_180",
  covariate_names = default_covariates,
  log_response = F
)

X_validation <- hamburg_single$X

if (!is.null(restricted$fixed_gene_features)) {
  available <- intersect(restricted$fixed_gene_features, colnames(X_validation))
  missing   <- setdiff(restricted$fixed_gene_features, colnames(X_validation))
  if (length(missing) > 0) {
    message(sprintf(
      "[validation] %d of %d discovery-selected genes are not available in Hamburg and will be dropped: %s",
      length(missing), length(restricted$fixed_gene_features), paste(missing, collapse = ", ")
    ))
  }
  X_validation <- X_validation[, available, drop = FALSE]
}

cat(sprintf(
  "PREVAC rVSV -> Hamburg validation: %d participants, %d features.\n",
  length(hamburg_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = hamburg_single$Y, covariates = hamburg_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(results_dir, "validation_fit_prevac_rvsv_on_hamburg.rds"))

figure_path <- fs::path("output", "figures", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  title = "PREVAC rVSV best pipeline: validated on Hamburg",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_prevac_rvsv_on_hamburg.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
