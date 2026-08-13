# analysis/pipeline_comparisons/visualize_comparisons.R
#
# Cross-dataset visualisation of pipeline-comparison results: for each
# methodological category (engineering, selection_geneset,
# selection_genewise, model), compare options fairly across datasets by
# looking at delta_r2 (each option's R2 minus its own comparison's
# reference R2) rather than raw R2 - see R/metrics_analysis.R for why a
# difference, computed within each comparison before pooling, is the fair
# way to do this (it cancels out how inherently predictable a given
# dataset's response is, which raw R2 would not).
#
# Produces two separate plots per category (NOT faceted into one figure per
# plot type - each category is its own file, since different categories are
# evaluated against different reference pipelines and so aren't on a
# directly comparable delta_r2 scale to be shown side-by-side):
#   * a "mini forest plot" (R/comparison_plots.R::plot_relative_effects()) -
#     excludes each category's own reference row (compute_relative_metrics()'s
#     is_reference flag): a point at delta_r2 = 0 showing the reference's
#     distance from itself adds nothing here.
#   * a rank/delta_r2 heatmap (R/comparison_plots.R::plot_relative_heatmap()) -
#     includes the reference row (e.g. "Elastic net (reference)" for
#     category = "model"), participating in the within-dataset ranking like
#     any other option, since it's a real, valid choice for that category.
#
# That's 4 categories x 2 plot types = 8 figures, saved to
# output/figures/pipeline_comparisons/cross_dataset/.
#
# Run analysis/pipeline_comparisons/collect_metrics.R first (and, before
# that, every dataset folder's 01-04 scripts - see each folder's README.md).

library(dplyr)
library(fs)
library(ggplot2)

source(fs::path("R", "metrics_labels.R"))
source(fs::path("R", "metrics_analysis.R"))
source(fs::path("R", "comparison_plots.R"))
source(fs::path("R", "plotting.R"))

all_metrics <- readRDS(fs::path("data", "derived", "pipeline_comparison_metrics.rds"))

relative_metrics <- compute_relative_metrics(all_metrics)

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "cross_dataset")

categories <- unique(relative_metrics$category)

for (cat in categories) {

  df_category <- relative_metrics %>% dplyr::filter(.data$category == cat)
  cat_label <- category_display_name(cat)
  n_options <- dplyr::n_distinct(df_category$canonical_option)

  df_effects <- df_category %>% dplyr::filter(!.data$is_reference)

  p_effects <- plot_relative_effects(
    df_effects,
    title = paste0("Relative performance by choice: ", cat_label)
  )
  print(p_effects)
  save_pipeline_comparison_plot(
    p_effects, figure_path, paste0("delta_r2_dotplot_", cat, ".pdf"),
    width = 9, height = max(3, 0.45 * (n_options - 1) + 1.5)
  )

  p_heatmap <- plot_relative_heatmap(
    df_category,
    title = paste0("Relative performance by choice: ", cat_label)
  )
  print(p_heatmap)
  save_pipeline_comparison_plot(
    p_heatmap, figure_path, paste0("delta_r2_heatmap_", cat, ".pdf"),
    width = 8, height = max(3, 0.45 * n_options + 1.5)
  )
}

rm(list = ls())
