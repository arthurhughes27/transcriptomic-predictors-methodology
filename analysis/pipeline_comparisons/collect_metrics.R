# analysis/pipeline_comparisons/collect_metrics.R
#
# Pool the per-comparison metrics tables saved by each dataset folder's
# 02_compare_engineering.R / 02b_compare_engineering_genewise.R /
# 03_compare_selection.R / 03b_compare_selection_genewise.R /
# 04_compare_model.R scripts (via R/metrics_io.R::save_comparison_metrics())
# into a single, long-format data frame spanning every dataset and category,
# and cache it for downstream cross-dataset visualisation.
#
# Columns: pipeline (option/Reference/Baseline label), role
# ("baseline"/"reference"/"option"), RMSE, sRMSE, R2, SpearmanR, dataset
# (e.g. "sdy1276_tiv", "prevac_rvsv", "prevac_ad26mva"), category
# ("engineering", "engineering_genewise", "selection_geneset",
# "selection_genewise", "model"), and comparison (identifies which specific
# compare_pipelines() call a row came from - see R/metrics_io.R for why this
# can differ from category).
#
# This script does not itself compute anything relative to a reference (see
# R/metrics_io.R::load_all_comparison_metrics()'s docs on why that must be
# done per-comparison, not pooled) - it only assembles the raw, long-format
# table that a later visualisation script will need for that.
#
# Run every dataset folder's 01-04(b) scripts first (see each folder's
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

derived_data_dir <- fs::path("data", "derived")
fs::dir_create(derived_data_dir)

saveRDS(all_metrics, file = fs::path(derived_data_dir, "pipeline_comparison_metrics.rds"))

rm(list = ls())
