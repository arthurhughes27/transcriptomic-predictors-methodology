# R/comparison_plots.R
#
# Two cross-dataset visualisations of delta_r2 (see R/metrics_analysis.R),
# each drawn once per category (never faceted across categories - one
# category's choices are compared against a different reference pipeline
# than another's, so mixing categories on one axis would be misleading):
#
#   * plot_relative_effects(): a "mini forest plot" - one row per option,
#     one point per dataset, plus a diamond marking the mean delta_r2
#     across datasets. Options consistently to the right of zero, across
#     all dataset colours, are choices that robustly help; options with
#     dataset points scattered on both sides of zero are dataset-dependent.
#   * plot_relative_heatmap(): options x datasets, cell fill = delta_r2,
#     cell label = rank within that dataset (1 = best option in that
#     dataset/category). A compact companion view of the same numbers.
#
# Both take the *already-computed* long-format output of
# R/metrics_analysis.R::compute_relative_metrics(), pre-filtered to a single
# category, so they have no dependency on how delta_r2 was derived.

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
      ggplot2::aes(colour = .data$dataset_label), size = 3, alpha = 0.85
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
      title = title,
      caption = "Diamond = mean across datasets. Reference pipeline's own delta_r2 is 0 by construction, not shown."
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "grey40"))
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
      low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      fill = expression(Delta * R^2),
      title = title,
      caption = "Cell label = rank within that dataset/category (1 = best). Blank cells = option not compared for that dataset."
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "grey40"),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
}
