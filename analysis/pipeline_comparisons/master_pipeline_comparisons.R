# analysis/pipeline_comparisons/master_pipeline_comparisons.R
#
# Master script for the within-study pipeline-comparison analyses: runs all
# three dataset folders' own master scripts (each of which runs that
# dataset's 01-05 comparison/best-pipeline-search scripts in order - see
# that folder's README.md), then pools every dataset's saved metrics and
# produces the cross-dataset comparison figures.
#
# As in each dataset folder's 00_master.R, every source() call below
# re-specifies its own full path rather than reusing a shared variable,
# since the sourced scripts (transitively) end with rm(list = ls()), which
# would otherwise wipe such a variable (source() evaluates in the calling/
# global environment by default).

## ---- 1. Per-dataset comparisons and best-pipeline search ------------------

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "00_master.R"))

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "00_master.R"))

source(fs::path("analysis", "pipeline_comparisons", "prevac_ad26mva", "00_master.R"))

## ---- 2. Pool metrics across datasets, then visualise -----------------------

source(fs::path("analysis", "pipeline_comparisons", "collect_metrics.R"))

source(fs::path("analysis", "pipeline_comparisons", "visualize_comparisons.R"))
