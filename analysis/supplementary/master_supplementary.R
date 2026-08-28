# analysis/supplementary/master_supplementary.R
#
# Master script for ALL gene-wise supplementary comparisons: runs all three
# dataset folders' own supplementary masters (engineering/selection/model
# comparisons, then the gene-wise best-pipeline search), pools their metrics
# and produces the cross-dataset supplementary figures, then validates each
# vaccine's best gene-wise pipeline on its external dataset (see
# analysis/supplementary/external_validation/).
#
# NOT part of the main analysis: this is not sourced by
# analysis/master_analysis.R, and its outputs (output/results/supplementary/,
# output/figures/supplementary/) are entirely separate from the main
# analysis' own (output/results/, output/figures/pipeline_comparisons/) -
# see this folder's per-dataset scripts and R/metrics_io.R's header for why.
# Run this explicitly whenever the supplementary (gene-wise) comparisons
# are wanted.
#
# Depends on each dataset's analysis-ready data already being built by
# analysis/pipeline_comparisons/<dataset>/01_prepare_data.R (data
# preparation itself is shared with the main analysis, not duplicated here)
# - run analysis/pipeline_comparisons/master_pipeline_comparisons.R (or at
# least each dataset's 01_prepare_data.R) first.
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since the
# sourced scripts (transitively) end with rm(list = ls()), which would
# otherwise wipe such a variable (source() evaluates in the calling/global
# environment by default).

## ---- 1. Per-dataset gene-wise comparisons ----------------------------------

source(fs::path("analysis", "supplementary", "sdy1276_tiv", "00_master.R"))

source(fs::path("analysis", "supplementary", "prevac_rvsv", "00_master.R"))

source(fs::path("analysis", "supplementary", "prevac_ad26mva", "00_master.R"))

## ---- 2. Pool metrics across datasets, then visualise -----------------------

source(fs::path("analysis", "supplementary", "collect_metrics.R"))

source(fs::path("analysis", "supplementary", "visualize_comparisons.R"))

## ---- 3. External validation of each vaccine's best gene-wise pipeline -----

source(fs::path("analysis", "supplementary", "external_validation", "master_external_validation_genewise.R"))
