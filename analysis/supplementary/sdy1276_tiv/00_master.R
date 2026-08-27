# analysis/supplementary/sdy1276_tiv/00_master.R
#
# Master script for sdy1276_tiv's SUPPLEMENTARY (gene-wise) comparisons only - runs
# 01_compare_engineering_genewise.R then 02_compare_selection_genewise.R.
# Not sourced by analysis/pipeline_comparisons/sdy1276_tiv/00_master.R or
# analysis/master_analysis.R - run this explicitly (or via
# analysis/supplementary/master_supplementary.R for all three datasets).
#
# Depends on the analysis-ready dataset built by
# analysis/pipeline_comparisons/sdy1276_tiv/01_prepare_data.R (data preparation
# itself is shared with the main analysis, not duplicated here - only the
# comparisons and their outputs are kept separate) - run that first.
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "supplementary", "sdy1276_tiv", "01_compare_engineering_genewise.R"))

gc()

source(fs::path("analysis", "supplementary", "sdy1276_tiv", "02_compare_selection_genewise.R"))

gc()
