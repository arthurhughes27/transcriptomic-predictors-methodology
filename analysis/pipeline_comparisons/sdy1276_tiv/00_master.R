# analysis/pipeline_comparisons/sdy1276_tiv/00_master.R
#
# Master script to run the SDY1276 (TIV) pipeline-comparison scripts in
# order. See README.md for details on each step.
#
# The gene-wise SUPPLEMENTARY comparisons for sdy1276_tiv (gene-level
# engineering, gene-wise selection) are NOT run here - they live entirely
# separately at analysis/supplementary/sdy1276_tiv/00_master.R (see that folder's
# header for why).
#
# Each source() call below re-specifies its own full path, rather than
# reusing a shared folder variable, because every script in this folder
# ends with rm(list = ls()) - since source() evaluates in the calling
# (global) environment by default, that would also wipe a shared variable
# defined here.

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "01_prepare_data.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "02_compare_engineering.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "03_compare_selection.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "04_compare_model.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "05_find_best_model.R"))

gc()

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "06_plot_cv_by_sex.R"))
