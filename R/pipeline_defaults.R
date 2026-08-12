# R/pipeline_defaults.R
#
# Shared defaults for the pipeline-comparison scripts under
# analysis/pipeline_comparisons/, so that the "reference" pipeline and the
# baseline covariate set are defined identically in every comparison
# (engineering, selection, model choice), rather than being redefined (and
# risking drift) in each script.
#
# The reference pipeline: a z-scored gene-level transform, mean aggregation
# of Blood Transcriptional Modules, a relaxed variance pre-filter (top 7,500
# gene sets) to keep model fitting tractable, and an elastic net model.

#' Baseline covariates used throughout the pipeline-comparison analyses
#' (Chapter 5, Section 5.4): age, sex, and race.
# default_covariates <- c("age", "sex", "race", "ab_p_0")
default_covariates <- c("age", "sex", "race")

#' The reference (minimal-complexity) analytical pipeline.
#'
#' `genesets` must be supplied explicitly (rather than picked up from a
#' global variable) so this function has no hidden dependency on the calling
#' script's environment.
#'
#' @param genesets Named list of gene sets (see R/data_io.R::load_genesets()),
#'   used for the reference's mean gene-set aggregation.
#' @param model_inner_folds Number of inner CV folds for model hyperparameter
#'   tuning (see Chapter 5, Section 5.3.1).
#' @param model_metric Metric used to select model hyperparameters in the
#'   inner CV loop.
reference_pipeline_params <- function(genesets, model_inner_folds = 5, model_metric = "r2") {
  list(
    engineering_params = list(
      method = "engineer",
      col_transform = "z",
      gene_level_fc = FALSE,
      genesets = genesets,
      agg_method = "mean"
    ),
    selection_params = list(method = "variance", top_n = 7500),
    model_params = list(method = "glmnet", inner_folds = model_inner_folds, metric = model_metric)
  )
}

#' A "raw gene" pipeline: z-scored gene-level transform, no gene-set
#' aggregation, the same relaxed variance pre-filter as
#' `reference_pipeline_params()` (to keep model fitting tractable on ~20,000
#' raw genes), and an elastic net model.
#'
#' Used as the comparator for feature-selection methods that operate on (and
#' return) individual gene names rather than aggregated gene-set scores -
#' e.g. the paired RISE/dearseq options in 03_compare_selection.R - since
#' those selected gene names would not match the reference pipeline's
#' aggregated (gene-set-named) columns.
#'
#' @param model_inner_folds,model_metric As in reference_pipeline_params().
#' @param variance_top_n Number of genes retained by the variance pre-filter.
#'   Matches reference_pipeline_params()'s default (7,500) unless overridden.
raw_gene_reference_params <- function(model_inner_folds = 5, model_metric = "r2",
                                       variance_top_n = 7500) {
  list(
    engineering_params = list(
      method = "engineer",
      col_transform = "z",
      gene_level_fc = FALSE,
      genesets = NULL,
      agg_method = "mean"
    ),
    selection_params = list(method = "variance", top_n = variance_top_n),
    model_params = list(method = "glmnet", inner_folds = model_inner_folds, metric = model_metric)
  )
}
