# analysis/external_validation/master_external_validation.R
#
# Master script to run all three external-validation scripts in order (one
# per vaccine - see each script's own header for details):
#   SDY1276 (TIV)       -> SDY80
#   PREVAC Ad26/MVA     -> EBOVAC2
#   PREVAC rVSV         -> Hamburg
#
# Each script re-derives the fixed feature panel from its vaccine's saved
# best_pipeline_*.rds (see analysis/pipeline_comparisons/), so run
# analysis/pipeline_comparisons/master_pipeline_comparisons.R (through at
# least each dataset folder's 05_find_best_model.R) first.
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "external_validation", "validate_sdy1276_tiv_on_sdy80.R"))

source(fs::path("analysis", "external_validation", "validate_prevac_ad26mva_on_ebovac2.R"))

source(fs::path("analysis", "external_validation", "validate_prevac_rvsv_on_hamburg.R"))
