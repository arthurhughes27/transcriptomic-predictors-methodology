# R/metrics_io.R
#
# Save and re-load per-comparison results tables
# (`predictomics_comparison$results`, as returned by
# `predictomics::compare_pipelines()`) across the pipeline-comparison
# scripts, tagged with `dataset` and `category`, so they can later be pooled
# into a single data frame spanning every dataset and category (see
# analysis/pipeline_comparisons/collect_metrics.R) for cross-dataset
# visualisation.
#
# "category" identifies which methodological axis a given comparison
# varies: "engineering", "engineering_genewise", "selection_geneset",
# "selection_genewise", or "model". Geneset-level and gene-wise
# engineering/selection are each kept as two distinct categories (rather
# than one "engineering"/"selection" category) because the "main" analyses
# (engineering, selection_geneset, model) are geneset-only - the gene-wise
# ones are supplementary and, for selection specifically, compare against a
# different reference pipeline (mean-gene-set-aggregated vs. raw-gene; see
# R/pipeline_defaults.R) - so pooling them without that distinction would
# silently mix a main and a supplementary comparison, or two different
# baselines.

#' Directory where per-comparison metrics tables are cached.
metrics_dir <- function() {
  fs::path("data", "derived", "metrics")
}

#' Save one `compare_pipelines()` comparison's results table, tagged with
#' `dataset` and `category`, to a dedicated file in `metrics_dir()`.
#'
#' Call this once per `run_pipeline_comparison()`/`compare_pipelines()` call
#' (i.e. once per figure) - once within each of `02_compare_engineering.R`,
#' `02b_compare_engineering_genewise.R`, `04_compare_model.R`, and twice
#' within `03_compare_selection.R`/`03b_compare_selection_genewise.R`'s
#' pairing (once for the geneset-level comparison, once for the gene-wise
#' one - each its own script/category, unlike engineering's split).
#'
#' `label` (defaulting to `category`) identifies the *file* a comparison is
#' saved to, while `category` identifies the *grouping* used when pooling
#' across datasets (see `load_all_comparison_metrics()`); for every current
#' script, one `category` maps to exactly one `label`/file per dataset, so
#' this distinction rarely matters in practice - it exists for the (now
#' hypothetical) case of a category needing more than one
#' `compare_pipelines()` call/figure.
#'
#' @param res A `predictomics_comparison` object.
#' @param dataset Character scalar identifying the dataset/vaccine arm,
#'   e.g. "sdy1276_tiv", "prevac_rvsv", "prevac_ad26mva".
#' @param category Character scalar identifying the methodological axis
#'   varied in this comparison: one of "engineering", "engineering_genewise",
#'   "selection_geneset", "selection_genewise", "model".
#' @param label Character scalar identifying this specific comparison/file,
#'   for cases where more than one comparison shares the same
#'   `dataset`/`category` (see Details). Defaults to `category`.
#'
#' @return The tagged results data frame, invisibly.
save_comparison_metrics <- function(res, dataset, category, label = category) {
  fs::dir_create(metrics_dir())

  results <- res$results
  results$dataset <- dataset
  results$category <- category
  results$comparison <- label

  saveRDS(
    results,
    file = fs::path(metrics_dir(), paste0(dataset, "__", label, ".rds"))
  )

  invisible(results)
}

#' Load and row-bind every saved comparison results table in
#' `metrics_dir()` into a single long-format data frame.
#'
#' @return A data frame with columns `pipeline`, `role`
#'   ("baseline"/"reference"/"option"), `RMSE`, `sRMSE`, `R2`, `SpearmanR`,
#'   `dataset`, `category`, and `comparison` - one row per pipeline per
#'   saved comparison. Note that the "Reference"/"Baseline" `role` rows are
#'   local to whichever specific `comparison` produced them - when a
#'   `category` spans more than one `comparison` (see
#'   `save_comparison_metrics()`), each `comparison`'s options must be
#'   compared against its *own* Reference/Baseline rows (matched on
#'   `dataset`, `category`, *and* `comparison`), not pooled across
#'   `comparison`s within the category.
load_all_comparison_metrics <- function() {
  files <- fs::dir_ls(metrics_dir(), glob = "*.rds")

  if (length(files) == 0L) {
    stop(
      "[metrics_io] No saved comparison metrics found in ", metrics_dir(),
      ". Run each dataset folder's 02/03/04 comparison scripts first.",
      call. = FALSE
    )
  }

  dplyr::bind_rows(lapply(files, readRDS))
}
