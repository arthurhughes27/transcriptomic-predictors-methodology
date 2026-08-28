# analysis/supplementary/external_validation/master_external_validation_genewise.R
#
# SUPPLEMENTARY, gene-wise counterpart to
# analysis/external_validation/master_external_validation.R: runs all three
# gene-wise external-validation scripts in order (one per vaccine - see each
# script's own header for details):
#   SDY1276 (TIV)       -> SDY80
#   PREVAC Ad26/MVA     -> EBOVAC2
#   PREVAC rVSV         -> Hamburg
#
# Each script re-derives the fixed feature panel from its vaccine's saved
# output/results/supplementary/best_pipeline_genewise_*.rds, so run
# analysis/supplementary/master_supplementary.R's per-dataset step (through
# at least each dataset folder's 04_find_best_model_genewise.R) first.
#
# As elsewhere in this repository, every source() call below re-specifies
# its own full path rather than reusing a shared variable, since each
# sourced script ends with rm(list = ls()), which would otherwise wipe such
# a variable (source() evaluates in the calling/global environment by
# default).

source(fs::path("analysis", "supplementary", "external_validation", "validate_sdy1276_tiv_on_sdy80_genewise.R"))

source(fs::path("analysis", "supplementary", "external_validation", "validate_prevac_ad26mva_on_ebovac2_genewise.R"))

source(fs::path("analysis", "supplementary", "external_validation", "validate_prevac_rvsv_on_hamburg_genewise.R"))
