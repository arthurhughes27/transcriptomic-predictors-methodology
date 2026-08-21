# analysis/pipeline_comparisons/prevac_rvsv/00_master.R
#
# Master script to run the PREVAC rVSV (+ placebo) pipeline-comparison
# scripts in order. See README.md for details on each step.
#
# Each source() call below re-specifies its own full path, rather than
# reusing a shared folder variable, because every script in this folder
# ends with rm(list = ls()) - since source() evaluates in the calling
# (global) environment by default, that would also wipe a shared variable
# defined here.

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "01_prepare_data.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "02_compare_engineering.R"))

gc()

# source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "02b_compare_engineering_genewise.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "03_compare_selection.R"))

gc()

# source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "03b_compare_selection_genewise.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "04_compare_model.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "05_find_best_model.R"))

gc()
