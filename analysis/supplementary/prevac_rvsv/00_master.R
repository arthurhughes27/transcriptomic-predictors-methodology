# analysis/supplementary/prevac_rvsv/00_master.R
#
# Master script for prevac_rvsv's SUPPLEMENTARY (gene-wise) comparisons only - runs
# 01_compare_engineering_genewise.R, 02_compare_selection_genewise.R,
# 03_compare_model_genewise.R, then the gene-wise best-pipeline search
# (04_find_best_model_genewise.R). Not sourced by
# analysis/pipeline_comparisons/prevac_rvsv/00_master.R or
# analysis/master_analysis.R - run this explicitly (or via
# analysis/supplementary/master_supplementary.R for all three datasets).
#
# Depends on the analysis-ready dataset built by
# analysis/pipeline_comparisons/prevac_rvsv/01_prepare_data.R (data preparation
# itself is shared with the main analysis, not duplicated here - only the
# comparisons and their outputs are kept separate) - run that first.
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "supplementary", "prevac_rvsv", "01_compare_engineering_genewise.R"))

gc()

source(fs::path("analysis", "supplementary", "prevac_rvsv", "02_compare_selection_genewise.R"))

gc()

source(fs::path("analysis", "supplementary", "prevac_rvsv", "03_compare_model_genewise.R"))

gc()

source(fs::path("analysis", "supplementary", "prevac_rvsv", "04_find_best_model_genewise.R"))

gc()
