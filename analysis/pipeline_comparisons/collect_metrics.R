# analysis/pipeline_comparisons/collect_metrics.R
#
# Pool the per-comparison metrics tables saved by each dataset folder's
# 02_compare_engineering.R / 03_compare_selection.R / 04_compare_model.R
# scripts (via R/metrics_io.R::save_comparison_metrics()) into a single,
# long-format data frame spanning every dataset and category, and cache it
# for downstream cross-dataset visualisation.
#
# MAIN ANALYSIS ONLY: this reads metrics_dir() (i.e.
# output/results/metrics), which only ever contains the three "main"
# categories below - the gene-wise supplementary comparisons
# (engineering_genewise, selection_genewise; see analysis/supplementary/)
# are saved to a completely separate output/results/supplementary/metrics
# store (subdir = "supplementary") and pooled by
# analysis/supplementary/collect_metrics.R instead, never mixed in here.
#
# Columns: pipeline (option/Reference/Baseline label), role
# ("baseline"/"reference"/"option"), RMSE, sRMSE, R2, SpearmanR, dataset
# (e.g. "sdy1276_tiv", "prevac_rvsv", "prevac_ad26mva"), category
# ("engineering", "selection_geneset", "model"), and comparison (identifies
# which specific compare_pipelines() call a row came from - see
# R/metrics_io.R for why this can differ from category).
#
# This script does not itself compute anything relative to a reference (see
# R/metrics_io.R::load_all_comparison_metrics()'s docs on why that must be
# done per-comparison, not pooled) - it only assembles the raw, long-format
# table that a later visualisation script will need for that.
#
# Run every dataset folder's 01-04 scripts first (see each folder's
# README.md), and re-run this script whenever any of those are re-run.

library(dplyr)
library(fs)

source(fs::path("R", "metrics_io.R"))

all_metrics <- load_all_comparison_metrics()

cat("Pipeline-comparison metrics collected:\n")
print(
  all_metrics %>%
    dplyr::count(dataset, category, comparison) %>%
    dplyr::arrange(dataset, category, comparison)
)

results_dir <- fs::path("output", "results")
fs::dir_create(results_dir)

saveRDS(all_metrics, file = fs::path(results_dir, "pipeline_comparison_metrics.rds"))

rm(list = ls())
