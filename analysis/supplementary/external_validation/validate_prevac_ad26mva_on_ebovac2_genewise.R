# analysis/supplementary/external_validation/validate_prevac_ad26mva_on_ebovac2_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/external_validation/validate_prevac_ad26mva_on_ebovac2.R:
# validate PREVAC Ad26/MVA's best GENE-WISE pipeline (see
# analysis/supplementary/prevac_ad26mva/04_find_best_model_genewise.R) on
# EBOVAC2, using the exact same discovery -> validation methodology
# described in R/external_validation.R's header - just with the gene-wise
# best-pipeline result and gene-wise discovery dataset in place of the
# geneset-level ones.
#
# ASSUMPTIONS ABOUT EBOVAC2: as in
# analysis/external_validation/validate_prevac_ad26mva_on_ebovac2.R - see
# that script's header.
#
# Run analysis/supplementary/prevac_ad26mva/00_master.R (through at least
# 04_find_best_model_genewise.R) first, so that
# output/results/supplementary/best_pipeline_genewise_prevac_ad26mva.rds
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

discovery <- readRDS(fs::path(results_dir, "prevac_ad26mva_analysis_data.rds"))
best <- readRDS(fs::path(supplementary_results_dir, "best_pipeline_genewise_prevac_ad26mva.rds"))

# --- Recover the fixed feature panel from PREVAC Ad26/MVA (discovery) ------
selected_features <- get_discovery_selected_features(
  best,
  X = discovery$single$X, Y = discovery$single$Y,
  covariates = discovery$single$covariates, treatment = discovery$single$treatment
)

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

# Engineering never aggregates for the gene-wise best pipeline, so
# restrict_engineering_for_validation() only ever restricts the fixed gene
# feature panel to genes actually present in EBOVAC2 - see that function's
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
      "[validation] %d of %d discovery-selected genes are not available in EBOVAC2 and will be dropped: %s",
      length(missing), length(restricted$fixed_gene_features), paste(missing, collapse = ", ")
    ))
  }
  X_validation <- X_validation[, available, drop = FALSE]
}

cat(sprintf(
  "PREVAC Ad26/MVA gene-wise best pipeline -> EBOVAC2 validation: %d participants, %d features.\n",
  length(ebovac2_single$participant_id), ncol(X_validation)
))

fit <- run_validation_fit(
  X = X_validation, Y = ebovac2_single$Y, covariates = ebovac2_single$covariates,
  engineering_params = restricted$engineering_params, model_params = best$model_params
)

saveRDS(fit, file = fs::path(supplementary_results_dir, "validation_fit_genewise_prevac_ad26mva_on_ebovac2.rds"))

figure_path <- fs::path("output", "figures", "supplementary", "external_validation")
fs::dir_create(figure_path)

p_summary <- plot_validation_summary(
  fit,
  model_params = best$model_params,
  title = "PREVAC Ad26/MVA best gene-wise pipeline: validated on EBOVAC2",
  subtitle = validation_subtitle(best, ncol(X_validation))
)
print(p_summary)
ggsave(
  fs::path(figure_path, "validation_summary_genewise_prevac_ad26mva_on_ebovac2.pdf"),
  p_summary, width = 12, height = 5, dpi = 300
)

rm(list = ls())
