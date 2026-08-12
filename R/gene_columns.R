# R/gene_columns.R
#
# Identify gene-expression columns in the merged clinical + expression data.
#
# The scripts under analysis/application/ historically selected gene columns
# positionally, e.g. `dplyr::select(a1cf:zzz3)`, relying on every gene symbol
# sorting alphabetically between those two bounds. This is fragile: any
# clinical/metadata column that happens to sort inside that range would
# silently be treated as a gene column. Here we instead explicitly enumerate
# the known non-gene ("metadata") columns and treat everything else in the
# merged data frame as a gene column, which is robust to column reordering
# and to the addition of new metadata columns upstream.

#' Non-gene metadata columns present in `df_merged_all.rds`.
#'
#' See analysis/preprocessing/preprocessing_clinical_harmonisation.R (for the
#' clinical/antibody columns) and
#' analysis/preprocessing/preprocessing_merging_harmonisation.R (for how they
#' are merged with gene expression) for the definitive column list.
metadata_columns <- function(df) {
  known <- c(
    "participant_id", "study_accession", "group", "study_vaccine",
    "group_long", "age", "sex", "race", "ethnicity", "time"
  )
  antibody_cols <- grep("^ab_(p|b)_", names(df), value = TRUE)
  intersect(c(known, antibody_cols), names(df))
}

#' Names of gene-expression columns in a merged clinical + expression data
#' frame (i.e. every column that is not a known metadata column).
get_gene_columns <- function(df) {
  setdiff(names(df), metadata_columns(df))
}
