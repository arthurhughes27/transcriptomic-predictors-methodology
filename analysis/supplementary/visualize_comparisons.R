# analysis/supplementary/visualize_comparisons.R
#
# SUPPLEMENTARY ONLY: cross-dataset visualisation of the two gene-wise
# supplementary categories (engineering_genewise, selection_genewise) -
# compare options fairly across datasets by looking at delta_r2 (each
# option's R2 minus its own comparison's reference R2) rather than raw R2 -
# see R/metrics_analysis.R for why a difference, computed within each
# comparison before pooling, is the fair way to do this. Structurally the
# same approach as analysis/pipeline_comparisons/visualize_comparisons.R
# (the main analysis' equivalent script), just reading from the completely
# separate supplementary metrics pool and category set - never mixed with
# the main analysis' three categories.
#
# Produces two separate plots per category (NOT faceted into one figure per
# plot type - each category is its own file, since the two categories are
# evaluated against different reference pipelines (mean-gene-set-aggregated
# vs. raw-gene; see R/pipeline_defaults.R) and so aren't on a directly
# comparable delta_r2 scale to be shown side-by-side):
#   * a "mini forest plot" (R/comparison_plots.R::plot_relative_effects())
#   * a rank/delta_r2 heatmap (R/comparison_plots.R::plot_relative_heatmap())
#
# That's 2 categories x 2 plot types = 4 figures, saved to
# output/figures/supplementary/cross_dataset/. Unlike the main analysis'
# visualize_comparisons.R, there is no reference-context figure here: the
# gene-wise reference pipelines are explicitly NOT treated as a/the
# reference approach for a dataset (see each dataset's
# 01_compare_engineering_genewise.R/02_compare_selection_genewise.R
# headers), so there's no "inherent signal" context to show for them.
#
# Run analysis/supplementary/collect_metrics.R first (and, before that,
# every dataset folder's analysis/supplementary/<dataset>/00_master.R).

library(dplyr)
library(fs)
library(ggplot2)

source(fs::path("R", "metrics_labels.R"))
source(fs::path("R", "metrics_analysis.R"))
source(fs::path("R", "comparison_plots.R"))
source(fs::path("R", "plotting.R"))

all_metrics <- readRDS(fs::path("output", "results", "supplementary", "pipeline_comparison_metrics_supplementary.rds"))

relative_metrics <- compute_relative_metrics(all_metrics)

figure_path <- fs::path("output", "figures", "supplementary", "cross_dataset")

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
    width = 10, height = max(3, 0.45 * (n_options - 1) + 1.5)
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
