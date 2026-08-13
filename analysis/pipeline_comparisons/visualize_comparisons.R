# analysis/pipeline_comparisons/visualize_comparisons.R
#
# Cross-dataset visualisation of pipeline-comparison results: for each
# methodological category (engineering, engineering_genewise,
# selection_geneset, selection_genewise, model), compare options fairly
# across datasets by
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
# That's 5 categories x 2 plot types = 10 figures, plus one further figure
# (see below), saved to output/figures/pipeline_comparisons/cross_dataset/.
#
# A ninth figure, produced last, puts the above into context: a grouped bar
# chart of each dataset's raw (absolute) Baseline and Reference R2
# (R/comparison_plots.R::plot_reference_context()) - how much inherent
# signal (achievable R2 with a sensible default pipeline) each dataset has,
# independent of any methodological choice. This is the one place raw R2,
# rather than delta_r2, is shown. It uses ONLY the `category == "model"`
# Baseline/Reference rows (not "engineering"'s or "selection_geneset"'s -
# see the comment above that section for why), and only the geneset-mean
# reference pipeline - the gene-wise/raw-gene reference
# (`category == "selection_genewise"`) is a supplementary comparison (see
# each dataset folder's 03b_compare_selection_genewise.R) and is not
# treated as a/the reference approach for a dataset.
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

## ---- Reference/baseline context (absolute R2, not delta_r2) ---------------
#
# ONLY category == "model" is used: its Reference/Baseline rows are fit on
# the `single` (day-1/day-7-only) dataset in every dataset folder, with no
# ambiguity about which underlying sample was used. "engineering" now also
# fits the reference pipeline on `single` only in every folder (since
# 02_compare_engineering.R was split into a geneset-only main comparison and
# a separate 02b_compare_engineering_genewise.R supplementary one - see that
# split's rationale in each folder's README.md), but "selection_geneset"
# still has a sample-ambiguity issue for sdy1276_tiv specifically: its
# 03_compare_selection.R uses the `paired`-restricted sample (needed for
# dearseq's paired mode), so its reference fit there need not exactly match
# a `single`-only fit of the same pipeline - if those two ever disagree,
# that's a real discrepancy to go investigate at the data-construction level
# (most likely: `single` and `paired`'s independent complete-case gene
# filtering in 01_prepare_data.R admit slightly different gene sets and/or
# participants), not something to paper over here by averaging. "model" is
# used for this figure precisely because it's the one category guaranteed
# unambiguous across all three dataset folders.

reference_context <- all_metrics %>%
  dplyr::filter(.data$category == "model", .data$role %in% c("baseline", "reference")) %>%
  dplyr::mutate(
    role_label = ifelse(.data$role == "baseline", "Baseline", "Reference")
  ) %>%
  dplyr::select(dataset, role_label, R2)

p_reference_context <- plot_reference_context(
  reference_context,
  title = "Cross-validation results: baseline and reference pipelines"
)
print(p_reference_context)
save_pipeline_comparison_plot(
  p_reference_context, figure_path, "reference_context.pdf",
  width = 8, height = 4.5
)

rm(list = ls())
