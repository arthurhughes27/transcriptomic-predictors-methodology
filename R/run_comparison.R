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
                                     n_workers = 6, verbose = TRUE) {
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
    individual_id = individual_id,
    diagnostics = "summary"
  )
}

#' Run `run_pipeline_comparison()`, but skip it entirely if a previous run
#' was already cached to `cache_path` - so a comparison script (02, 02b, 03,
#' 03b, 04) can be re-run (e.g. to re-generate its figure after a plotting
#' tweak) without re-fitting every pipeline in it.
#'
#' @param cache_path Path (e.g. from `R/metrics_io.R::comparison_cache_path()`)
#'   to save/load the full `predictomics_comparison` object as an `.rds`
#'   file. If it already exists, it's loaded and returned as-is - NOT
#'   re-validated against the arguments below, so deleting the file (or
#'   passing a different `cache_path`) is how to force a re-run after
#'   changing an option, the reference pipeline, or the data.
#' @param ... Forwarded to `run_pipeline_comparison()`.
#'
#' @return A `predictomics_comparison` object, freshly fit or loaded from
#'   cache.
run_or_load_comparison <- function(cache_path, ...) {
  if (fs::file_exists(cache_path)) {
    message("[run_comparison] Loading cached comparison from ", cache_path,
            " (delete this file to force a re-run).")
    return(readRDS(cache_path))
  }

  res <- run_pipeline_comparison(...)

  fs::dir_create(fs::path_dir(cache_path))
  saveRDS(res, file = cache_path)

  res
}
