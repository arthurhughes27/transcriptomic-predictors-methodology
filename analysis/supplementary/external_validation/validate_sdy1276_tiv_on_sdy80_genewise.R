# analysis/supplementary/external_validation/validate_sdy1276_tiv_on_sdy80_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/external_validation/validate_sdy1276_tiv_on_sdy80.R: validate
# SDY1276 (TIV)'s best GENE-WISE pipeline (see
# analysis/supplementary/sdy1276_tiv/04_find_best_model_genewise.R) on
# SDY80, using the exact same discovery -> validation methodology described
# in R/external_validation.R's header (fixed feature panel carried over
# from the discovery search, engineering/model choices unchanged, no
# re-selection on SDY80) - just with the gene-wise best-pipeline result and
# gene-wise discovery dataset in place of the geneset-level ones.
#
# ASSUMPTIONS ABOUT SDY80: as in
# analysis/external_validation/validate_sdy1276_tiv_on_sdy80.R - see that
# script's header.
#
# Run analysis/supplementary/sdy1276_tiv/00_master.R (through at least
# 04_find_best_model_genewise.R) first, so that
# output/results/supplementary/best_pipeline_genewise_sdy1276_tiv.rds
# exists.

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
source(fs::path("R", "best_pipeline_search.R"))

results_dir <- fs::path("output", "results")
supplementary_results_dir <- fs::path(results_dir, "supplementary")
fs::dir_create(supplementary_results_dir)

discovery <- readRDS(fs::path(results_dir, "sdy1276_tiv_analysis_data.rds"))
best <- readRDS(fs::path(supplementary_results_dir, "best_pipeline_genewise_sdy1276_tiv.rds"))

# --- Recover the fixed feature panel from SDY1276 (the discovery dataset) --
# SDY1276 is single-arm (no placebo), so treatment is NULL here - matching
# how find_best_pipeline_genewise() was itself called for this dataset.
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

# Engineering never aggregates for the gene-wise best pipeline, so
# restrict_engineering_for_validation() only ever restricts the fixed gene
# feature panel to genes actually present in SDY80 - see that function's
# docs.
restricted <- restrict_engineering_for_validation(
  best, selected_features, genesets = NULL,
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
  "SDY1276 (TIV) gene-wise best pipeline -> SDY80 validation: %d participants, %d features.\n",
  length(sdy80_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = sdy80_single$Y, covariates = sdy80_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(supplementary_results_dir, "validation_fit_genewise_sdy1276_tiv_on_sdy80.rds"))

figure_path <- fs::path("output", "figures", "supplementary", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  model_params = best$model_params,
  title = "SDY1276 (TIV) best gene-wise pipeline: validated on SDY80",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_genewise_sdy1276_tiv_on_sdy80.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
