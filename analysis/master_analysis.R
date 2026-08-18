# analysis/master_analysis.R
#
# Top-level master script for the whole chapter: runs preprocessing, then
# every analysis, in order:
#   1. Preprocessing (analysis/preprocessing/master_preprocessing.R) -
#      builds data/df_merged_all.rds and the processed gene-set files
#      everything downstream reads.
#   2. Pipeline comparisons (analysis/pipeline_comparisons/
#      master_pipeline_comparisons.R) - within-study feature-engineering/
#      selection/model-choice comparisons and best-pipeline search, for all
#      three vaccines, plus the cross-dataset comparison figures.
#   3. External validation (analysis/external_validation/
#      master_external_validation.R) - validates each vaccine's best
#      pipeline (from step 2) on an independent dataset for the same
#      vaccine.
#   4. Descriptive (analysis/descriptive/master_descriptive.R) - dataset
#      description figures/tables.
#
# As in every sub-master script, every source() call below re-specifies its
# own full path rather than reusing a shared variable, since the sourced
# scripts (transitively) end with rm(list = ls()), which would otherwise
# wipe such a variable (source() evaluates in the calling/global
# environment by default).

## ---- 1. Preprocessing ------------------------------------------------------

source(fs::path("analysis", "preprocessing", "master_preprocessing.R"))

## ---- 2. Pipeline comparisons (per-vaccine + cross-dataset) -----------------

source(fs::path("analysis", "pipeline_comparisons", "master_pipeline_comparisons.R"))

## ---- 3. External validation -------------------------------------------------

source(fs::path("analysis", "external_validation", "master_external_validation.R"))

## ---- 4. Descriptive ---------------------------------------------------------

source(fs::path("analysis", "descriptive", "master_descriptive.R"))
