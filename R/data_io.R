# R/data_io.R
#
# Helpers for loading the processed data objects shared across analysis
# scripts. Centralising these paths/loaders means the harmonised `data/`
# layout produced by analysis/preprocessing/ is defined in one place, rather
# than being re-typed (and potentially drifting) in every downstream script.

#' Path to the harmonised, merged clinical + gene-expression dataset produced
#' by analysis/preprocessing/preprocessing_merging_harmonisation.R
merged_data_path <- function() {
  fs::path("data", "df_merged_all.rds")
}

#' Path to the processed Blood Transcriptional Module gene sets produced by
#' analysis/preprocessing/preprocessing_BTM.R
genesets_path <- function() {
  fs::path("data", "BTM_processed.rds")
}

#' Load the harmonised, merged clinical + gene-expression dataset.
load_merged_data <- function(path = merged_data_path()) {
  readRDS(path)
}

#' Load processed gene sets and return them as a plain named list
#' (geneset name -> character vector of gene symbols), ready for use as the
#' `genesets` argument to `predictomics` functions.
load_genesets <- function(path = genesets_path()) {
  gs <- readRDS(path)
  genesets <- gs[["genesets"]]
  names(genesets) <- gs[["geneset.names.descriptions"]]
  genesets
}
