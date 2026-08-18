# analysis/descriptive/master_descriptive.R
#
# Master script to run all descriptive analyses in order:
#   - study_descriptions.R: gene-expression sample availability by
#     study/group and timepoint (figure).
#   - dataset_characteristics_table.R: manuscript-ready LaTeX table
#     summarising the exploration and external-validation cohorts (see
#     that script's header).
#
# Run analysis/preprocessing/master_preprocessing.R first (both scripts
# read from data/df_merged_all.rds).
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "descriptive", "study_descriptions.R"))

source(fs::path("analysis", "descriptive", "dataset_characteristics_table.R"))
