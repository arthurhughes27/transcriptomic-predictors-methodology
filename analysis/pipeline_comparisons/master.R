# analysis/pipeline_comparisons/master.R
#
# Master script to run all three within-study pipeline-comparison master
# scripts, in order: SDY1276 (TIV), PREVAC rVSV, PREVAC Ad26/MVA.
#
# As in each dataset folder's 00_master.R, every source() call below
# re-specifies its own full path rather than reusing a shared variable,
# since the sourced scripts (transitively) end with rm(list = ls()), which
# would otherwise wipe such a variable (source() evaluates in the calling/
# global environment by default).

source(fs::path("analysis", "pipeline_comparisons", "sdy1276_tiv", "00_master.R"))

source(fs::path("analysis", "pipeline_comparisons", "prevac_rvsv", "00_master.R"))

source(fs::path("analysis", "pipeline_comparisons", "prevac_ad26mva", "00_master.R"))
