# R/metrics_analysis.R
#
# Compute each option's performance relative to its own comparison's
# reference pipeline (delta_r2 = R2(option) - R2(reference)), so that
# methodological choices can be compared fairly across datasets with very
# different inherent predictability (see the chat discussion this
# implements: differences, not ratios - R2 can be <= 0, so a ratio is
# undefined/unstable near a small or negative reference, and differences -
# unlike ratios - aggregate to a well-defined mean when averaged across
# datasets).
#
# Requires R/metrics_labels.R to be sourced first (for
# canonicalize_option_label()).

#' Compute delta-R2 (vs. each comparison's own reference) for every option
#' in a pooled metrics data frame.
#'
#' @param all_metrics The combined data frame produced by
#'   `analysis/pipeline_comparisons/collect_metrics.R`
#'   (`R/metrics_io.R::load_all_comparison_metrics()`'s output): one row per
#'   pipeline per saved comparison, with columns `pipeline`, `role`, `R2`
#'   (among other metrics), `dataset`, `category`, `comparison`.
#'
#' @details
#' For each `(dataset, category, comparison)` group (i.e. each individual
#' `compare_pipelines()` call that was saved), every `role == "option"` row's
#' `R2` is compared against that *same group's* `role == "reference"` row's
#' `R2` - never against a reference from a different comparison, dataset, or
#' category. `role == "baseline"` rows are dropped entirely: they answer a
#' different question (does gene expression help at all, vs. covariates
#' alone) from the one these plots address (which choice, among those that
#' do use gene expression, performs best).
#'
#' Option labels are then canonicalized (`canonicalize_option_label()`) so
#' the same underlying choice, labelled slightly differently across dataset
#' folders (e.g. "paired" vs "classic" RISE/dearseq, or SDY1276's "Gene-wise"
#' vs PREVAC's "Gene-level: z-score"), is treated as one option. Where this
#' canonicalization causes more than one row to share a
#' `(dataset, category, canonical_option)` combination - which happens for
#' SDY1276's engineering comparison, whose two separate `compare_pipelines()`
#' calls both contribute a "z-score, no aggregation" option (see
#' `R/metrics_labels.R`'s comments) - those rows' `delta_r2` values are
#' averaged, so each dataset contributes at most one point per option to the
#' downstream plots.
#'
#' Each comparison's own `role == "reference"` row is also included as a
#' row in the output, labelled via `reference_option_label()` (e.g.
#' "Elastic net (reference)" for `category = "model"`) and flagged
#' `is_reference = TRUE`; its `delta_r2` is 0 by construction (it's being
#' compared to itself). This is deliberate, not a check you need to remove:
#' the reference pipeline is a real, valid choice for its own category and
#' belongs in the ranking/heatmap alongside the options it's the baseline
#' for. Filter on `is_reference` in a caller that doesn't want it (e.g. the
#' forest/dot plot, where a "distance from itself" point adds nothing).
#'
#' @return A data frame with columns `dataset`, `category`,
#'   `canonical_option`, `delta_r2` (mean, if more than one raw row
#'   contributed), `is_reference`, and `n_raw` (how many raw rows were
#'   averaged - 1 unless canonicalization merged rows).
compute_relative_metrics <- function(all_metrics) {

  references <- all_metrics %>%
    dplyr::filter(.data$role == "reference") %>%
    dplyr::select(dataset, category, comparison, reference_r2 = R2)

  options_with_ref <- all_metrics %>%
    dplyr::filter(.data$role == "option") %>%
    dplyr::inner_join(references, by = c("dataset", "category", "comparison")) %>%
    dplyr::transmute(
      dataset = .data$dataset,
      category = .data$category,
      canonical_option = canonicalize_option_label(.data$pipeline, .data$category),
      delta_r2 = .data$R2 - .data$reference_r2,
      is_reference = FALSE
    )

  reference_rows <- all_metrics %>%
    dplyr::filter(.data$role == "reference") %>%
    dplyr::transmute(
      dataset = .data$dataset,
      category = .data$category,
      canonical_option = paste0(reference_option_label(.data$category), " (reference)"),
      delta_r2 = 0,
      is_reference = TRUE
    )

  dplyr::bind_rows(options_with_ref, reference_rows) %>%
    dplyr::group_by(.data$dataset, .data$category, .data$canonical_option, .data$is_reference) %>%
    dplyr::summarise(
      delta_r2 = mean(.data$delta_r2),
      n_raw = dplyr::n(),
      .groups = "drop"
    )
}
