# R/run_comparison.R
#
# Thin wrapper around `predictomics::compare_pipelines()` that manages the
# `future` parallel backend and applies consistent outer cross-validation
# defaults across the pipeline-comparison scripts.
#
# The `future` plan is always reset to sequential on exit (via `on.exit`),
# including if `compare_pipelines()` errors, which avoids leaving background
# worker processes running after a failed script.

#' Run a `predictomics::compare_pipelines()` comparison with a managed
#' parallel backend.
#'
#' @param X,Y,covariates As for `predictomics::compare_pipelines()`.
#' @param option_type One of "engineering", "selection", "model".
#' @param option_choices Named list of option configurations to compare
#'   against the reference pipeline.
#' @param reference_params The reference pipeline definition (see
#'   R/pipeline_defaults.R).
#' @param treatment,treatment_predictor,timepoint,individual_id Optional
#'   arguments passed through to `compare_pipelines()`.
#' @param cv_type,folds,seed Outer cross-validation settings.
#' @param n_workers Number of parallel workers for `future::multisession`.
#' @param verbose Passed through to `compare_pipelines()`.
run_pipeline_comparison <- function(X, Y, covariates,
                                     option_type, option_choices, reference_params,
                                     treatment = NULL, treatment_predictor = FALSE,
                                     timepoint = NULL, individual_id = NULL,
                                     cv_type = "kfold", folds = 5, seed = 12345,
                                     n_workers = 7, verbose = TRUE) {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  predictomics::compare_pipelines(
    X = X,
    Y = Y,
    option_type = option_type,
    option_choices = option_choices,
    reference_params = reference_params,
    treatment = treatment,
    treatment_predictor = treatment_predictor,
    verbose = verbose,
    covariates = covariates,
    cv_type = cv_type,
    folds = folds,
    seed = seed,
    outside_cv = FALSE,
    timepoint = timepoint,
    individual_id = individual_id
  )
}
