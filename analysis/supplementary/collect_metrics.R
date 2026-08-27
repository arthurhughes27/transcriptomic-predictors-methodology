# analysis/supplementary/collect_metrics.R
#
# Pool the per-comparison metrics tables saved by each dataset folder's
# 01_compare_engineering_genewise.R / 02_compare_selection_genewise.R
# scripts (via R/metrics_io.R::save_comparison_metrics(..., subdir =
# "supplementary")) into a single, long-format data frame spanning every
# dataset and (supplementary) category, and cache it for downstream
# cross-dataset visualisation.
#
# SUPPLEMENTARY ONLY: this reads metrics_dir("supplementary") (i.e.
# output/results/supplementary/metrics), completely separate from the main
# analysis' own output/results/metrics pooled by
# analysis/pipeline_comparisons/collect_metrics.R - the two are never
# mixed.
#
# Columns: pipeline (option/Reference/Baseline label), role
# ("baseline"/"reference"/"option"), RMSE, sRMSE, R2, SpearmanR, dataset
# (e.g. "sdy1276_tiv", "prevac_rvsv", "prevac_ad26mva"), category
# ("engineering_genewise", "selection_genewise"), and comparison.
#
# Run every dataset folder's analysis/supplementary/<dataset>/00_master.R
# first, and re-run this script whenever any of those are re-run.

library(dplyr)
library(fs)

source(fs::path("R", "metrics_io.R"))

all_metrics <- load_all_comparison_metrics(subdir = "supplementary")

cat("Supplementary (gene-wise) pipeline-comparison metrics collected:\n")
print(
  all_metrics %>%
    dplyr::count(dataset, category, comparison) %>%
    dplyr::arrange(dataset, category, comparison)
)

results_dir <- fs::path("output", "results", "supplementary")
fs::dir_create(results_dir)

saveRDS(all_metrics, file = fs::path(results_dir, "pipeline_comparison_metrics_supplementary.rds"))

rm(list = ls())
