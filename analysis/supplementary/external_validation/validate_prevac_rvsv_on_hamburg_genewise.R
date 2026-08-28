# analysis/supplementary/external_validation/validate_prevac_rvsv_on_hamburg_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/external_validation/validate_prevac_rvsv_on_hamburg.R: validate
# PREVAC rVSV's best GENE-WISE pipeline (see
# analysis/supplementary/prevac_rvsv/04_find_best_model_genewise.R) on
# Hamburg, using the exact same discovery -> validation methodology
# described in R/external_validation.R's header - just with the gene-wise
# best-pipeline result and gene-wise discovery dataset in place of the
# geneset-level ones.
#
# ASSUMPTIONS ABOUT HAMBURG: as in
# analysis/external_validation/validate_prevac_rvsv_on_hamburg.R - see that
# script's header (including `log_response = FALSE`, matched here).
#
# Run analysis/supplementary/prevac_rvsv/00_master.R (through at least
# 04_find_best_model_genewise.R) first, so that
# output/results/supplementary/best_pipeline_genewise_prevac_rvsv.rds
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

discovery <- readRDS(fs::path(results_dir, "prevac_rvsv_analysis_data.rds"))
best <- readRDS(fs::path(supplementary_results_dir, "best_pipeline_genewise_prevac_rvsv.rds"))

# --- Recover the fixed feature panel from PREVAC rVSV (discovery) ----------
selected_features <- get_discovery_selected_features(
  best,
  X = discovery$single$X, Y = discovery$single$Y,
  covariates = discovery$single$covariates, treatment = discovery$single$treatment
)

# --- Build the Hamburg validation dataset -----------------------------------
df_merged_all <- load_merged_data()

hamburg_single <- build_prediction_dataset(
  df_merged = df_merged_all,
  study_vaccine_groups = "hamburg-rVSV",
  timepoint = "P+7D",
  response_col = "ab_p_180",
  covariate_names = default_covariates,
  log_response = FALSE
)

X_validation <- hamburg_single$X

# Engineering never aggregates for the gene-wise best pipeline, so
# restrict_engineering_for_validation() only ever restricts the fixed gene
# feature panel to genes actually present in Hamburg - see that function's
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
      "[validation] %d of %d discovery-selected genes are not available in Hamburg and will be dropped: %s",
      length(missing), length(restricted$fixed_gene_features), paste(missing, collapse = ", ")
    ))
  }
  X_validation <- X_validation[, available, drop = FALSE]
}

cat(sprintf(
  "PREVAC rVSV gene-wise best pipeline -> Hamburg validation: %d participants, %d features.\n",
  length(hamburg_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = hamburg_single$Y, covariates = hamburg_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(supplementary_results_dir, "validation_fit_genewise_prevac_rvsv_on_hamburg.rds"))

figure_path <- fs::path("output", "figures", "supplementary", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  model_params = best$model_params,
  title = "PREVAC rVSV best gene-wise pipeline: validated on Hamburg",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_genewise_prevac_rvsv_on_hamburg.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
