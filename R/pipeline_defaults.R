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
