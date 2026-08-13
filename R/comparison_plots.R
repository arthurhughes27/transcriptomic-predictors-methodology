# R/comparison_plots.R
#
# Cross-dataset visualisations built on the pooled pipeline-comparison
# metrics (see R/metrics_analysis.R and
# analysis/pipeline_comparisons/collect_metrics.R):
#
#   * plot_relative_effects(): a "mini forest plot" of delta_r2, one row per
#     option, one point per dataset, plus a diamond marking the mean
#     delta_r2 across datasets. Drawn once per category (never faceted
#     across categories - one category's choices are compared against a
#     different reference pipeline than another's, so mixing categories on
#     one axis would be misleading). Options consistently to the right of
#     zero, across all dataset colours, are choices that robustly help;
#     options with dataset points scattered on both sides of zero are
#     dataset-dependent.
#   * plot_relative_heatmap(): as above but options x datasets, cell fill =
#     delta_r2, cell label = rank within that dataset (1 = best option in
#     that dataset/category). A compact companion view of the same numbers.
#   * plot_reference_context(): NOT about delta_r2 - a grouped bar chart of
#     each dataset's raw (absolute) Baseline and Reference R2, for context
#     on how much inherent signal each dataset has (see that function's own
#     docs).
#
# plot_relative_effects()/plot_relative_heatmap() take the *already-computed*
# long-format output of R/metrics_analysis.R::compute_relative_metrics(),
# pre-filtered to a single category, so they have no dependency on how
# delta_r2 was derived.

#' A mini forest plot of delta_r2 by option, one point per dataset plus a
#' cross-dataset mean.
#'
#' @param df_category A data frame for ONE category, with columns
#'   `dataset`, `canonical_option`, `delta_r2` (as returned by
#'   `compute_relative_metrics()`, filtered to one `category`).
#' @param title Plot title.
plot_relative_effects <- function(df_category, title) {

  option_order <- df_category %>%
    dplyr::group_by(.data$canonical_option) %>%
    dplyr::summarise(mean_delta = mean(.data$delta_r2), .groups = "drop") %>%
    dplyr::arrange(.data$mean_delta) %>%
    dplyr::pull(.data$canonical_option)

  df_category <- df_category %>%
    dplyr::mutate(
      canonical_option = factor(.data$canonical_option, levels = option_order),
      dataset_label = dataset_display_name(.data$dataset)
    )

  summary_df <- df_category %>%
    dplyr::group_by(.data$canonical_option) %>%
    dplyr::summarise(mean_delta = mean(.data$delta_r2), .groups = "drop")

  ggplot2::ggplot(df_category, ggplot2::aes(x = .data$delta_r2, y = .data$canonical_option)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$dataset_label), size = 3, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(x = .data$mean_delta, y = .data$canonical_option),
      shape = 23, size = 3.5, fill = "black", colour = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = expression(Delta * R^2 ~ "vs. reference pipeline"),
      y = NULL,
      colour = "Dataset",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

#' A heatmap of delta_r2 by option x dataset, with within-dataset rank
#' annotated in each cell.
#'
#' @param df_category As for `plot_relative_effects()`.
#' @param title Plot title.
plot_relative_heatmap <- function(df_category, title) {

  option_order <- df_category %>%
    dplyr::group_by(.data$canonical_option) %>%
    dplyr::summarise(mean_delta = mean(.data$delta_r2), .groups = "drop") %>%
    dplyr::arrange(.data$mean_delta) %>%
    dplyr::pull(.data$canonical_option)

  df_category <- df_category %>%
    dplyr::mutate(dataset_label = dataset_display_name(.data$dataset)) %>%
    dplyr::group_by(.data$dataset_label) %>%
    dplyr::mutate(rank = rank(-.data$delta_r2, ties.method = "min")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(canonical_option = factor(.data$canonical_option, levels = option_order))

  fill_limit <- max(abs(df_category$delta_r2), na.rm = TRUE)

  ggplot2::ggplot(
    df_category,
    ggplot2::aes(x = .data$dataset_label, y = .data$canonical_option, fill = .data$delta_r2)
  ) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$rank), size = 3.5) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b", mid = "white", high = "#8DA0CB", midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      fill = expression(Delta * R^2),
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
                   plot.title = ggplot2::element_text(size = 12))
}

#' A grouped bar chart of raw (absolute, not delta) R2 for the Baseline
#' (covariates only) and Reference (geneset-mean-aggregated) pipelines,
#' one pair of bars per dataset.
#'
#' Unlike `plot_relative_effects()`/`plot_relative_heatmap()`, this is
#' deliberately about ABSOLUTE performance, not performance relative to a
#' reference - it exists to show how much inherent signal (achievable R2
#' with a sensible default pipeline) each dataset has, as context for
#' interpreting the delta_r2 plots. Bars (rather than points) are used
#' because, for an absolute magnitude with a meaningful zero, position from
#' that zero is the more accurate encoding than the point-based encoding
#' used for delta_r2 (which is about deviation from a reference, not
#' magnitude).
#'
#' @param df A data frame with columns `dataset`, `role_label` ("Baseline"
#'   or "Reference"), `R2`. Expected to have exactly one row per
#'   `dataset` x `role_label` combination - the caller is responsible for
#'   supplying a single, unambiguous R2 per cell (see
#'   analysis/pipeline_comparisons/visualize_comparisons.R for which
#'   category's Baseline/Reference rows are used, and why).
#' @param title Plot title.
plot_reference_context <- function(df, title) {

  df <- df %>%
    dplyr::mutate(
      dataset_label = dataset_display_name(.data$dataset),
      role_label = factor(.data$role_label, levels = c("Baseline", "Reference"))
    )

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$dataset_label, y = .data$R2, fill = .data$role_label)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c(Baseline = "grey75", Reference = "#8DA0CB")) +
    ggplot2::labs(
      x = NULL, y = expression("Cross-validated " * R^2),
      fill = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
