# analysis/descriptive/master_descriptive.R
#
# Master script to run all descriptive analyses in order:
#   - study_descriptions.R: gene-expression sample availability by
#     study/group and timepoint (figure).
#   - dataset_characteristics_table.R: manuscript-ready LaTeX table
#     summarising the exploration and external-validation cohorts (see
#     that script's header).
#   - model_hyperparameters_table.R: manuscript-ready LaTeX table of each
#     regression method's hyperparameters and tuning grids (see that
#     script's header) - static (a property of the code, not any dataset),
#     so it doesn't need data/df_merged_all.rds.
#
# Run analysis/preprocessing/master_preprocessing.R first (study_descriptions.R
# and dataset_characteristics_table.R both read from data/df_merged_all.rds).
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "descriptive", "study_descriptions.R"))

source(fs::path("analysis", "descriptive", "dataset_characteristics_table.R"))

source(fs::path("analysis", "descriptive", "model_hyperparameters_table.R"))
