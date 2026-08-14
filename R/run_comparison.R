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

#' Run `run_pipeline_comparison()`, but skip it entirely if this exact
#' comparison's results are already saved - so a comparison script (02, 02b,
#' 03, 03b, 04) can be re-run (e.g. to re-generate its figure after a
#' plotting tweak) without re-fitting every pipeline in it.
#'
#' Reuses `R/metrics_io.R::save_comparison_metrics()`'s own save file
#' (`metrics_dir()/<dataset>__<label>.rds`) as the cache, rather than a
#' second, separate cache - every 02/02b/03/03b/04 script already calls
#' `save_comparison_metrics()` right after this, so a metrics file existing
#' already means "this comparison has been run before" (including runs from
#' before this caching wrapper existed).
#'
#' The trade-off: only `results` (not the full `predictomics_comparison`
#' object - `fits`, `call`, etc. aren't saved by `save_comparison_metrics()`)
#' survives the round trip through the metrics file. That's sufficient for
#' every 02/02b/03/03b/04 script, which only ever calls `plot(res, metric = "R2")`
#' - `plot.predictomics_comparison()` reads only `x$results` and
#' `x$option_type` (`x$metric` only matters as a default when `metric` isn't
#' passed explicitly, which these scripts always do) - so a minimal
#' `predictomics_comparison`-classed stand-in with just those two fields is
#' enough to reproduce the same figure.
#'
#' @param dataset,label As for `save_comparison_metrics()` - identify the
#'   `metrics_dir()` file this comparison's results are saved to/loaded
#'   from.
#' @param option_type As for `run_pipeline_comparison()` - also needed here
#'   (independent of `...`) to label a cache-loaded stand-in object's
#'   `option_type` field for `plot()`'s title.
#' @param ... Forwarded to `run_pipeline_comparison()` (which also expects
#'   `option_type` among its named arguments - safe to pass once and have it
#'   matched to both this function's explicit `option_type` and `...`, since
#'   R resolves named arguments before `...`).
#'
#' @return A `predictomics_comparison` object (or, when loaded from cache, a
#'   minimal stand-in with the same `results`/`option_type` fields `plot()`
#'   needs).
run_or_load_comparison <- function(dataset, label, option_type, ...) {
  metrics_path <- fs::path(metrics_dir(), paste0(dataset, "__", label, ".rds"))

  if (fs::file_exists(metrics_path)) {
    message("[run_comparison] Loading cached comparison results from ", metrics_path,
            " (delete this file to force a re-run).")
    results <- readRDS(metrics_path)
    return(structure(list(results = results, option_type = option_type),
                      class = "predictomics_comparison"))
  }

  run_pipeline_comparison(option_type = option_type, ...)
}
