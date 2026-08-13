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
# varies: "engineering", "selection_geneset", "selection_genewise", or
# "model". Geneset-level and gene-wise feature selection are kept as two
# distinct categories (rather than one "selection" category) because they
# compare against two different reference pipelines (mean-gene-set-
# aggregated vs. raw-gene; see R/pipeline_defaults.R), so pooling them
# without that distinction would silently mix two different baselines.

#' Directory where per-comparison metrics tables are cached.
metrics_dir <- function() {
  fs::path("data", "derived", "metrics")
}

#' Save one `compare_pipelines()` comparison's results table, tagged with
#' `dataset` and `category`, to a dedicated file in `metrics_dir()`.
#'
#' Call this once per `run_pipeline_comparison()`/`compare_pipelines()` call
#' (i.e. once per figure) - so twice within a `03_compare_selection.R`
#' script (once for the geneset-level comparison, once for the gene-wise
#' one), once within `04_compare_model.R`, and once or twice within
#' `02_compare_engineering.R` (twice for sdy1276_tiv/, which splits gene-set
#' aggregation and gene-level transformation into two separate calls/figures
#' that are still both the "engineering" axis; once for the two PREVAC
#' folders, which compare both in a single call).
#'
#' Because sdy1276_tiv/02_compare_engineering.R's two calls share the same
#' `category` ("engineering") but must not overwrite each other's saved
#' file, `label` (defaulting to `category`) identifies the *file*, while
#' `category` identifies the *grouping* used when pooling across datasets
#' (see `load_all_comparison_metrics()`) - so both of that script's saves
#' still end up tagged `category = "engineering"` in the pooled data,
#' despite being saved to two different files.
#'
#' @param res A `predictomics_comparison` object.
#' @param dataset Character scalar identifying the dataset/vaccine arm,
#'   e.g. "sdy1276_tiv", "prevac_rvsv", "prevac_ad26mva".
#' @param category Character scalar identifying the methodological axis
#'   varied in this comparison: one of "engineering", "selection_geneset",
#'   "selection_genewise", "model".
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
