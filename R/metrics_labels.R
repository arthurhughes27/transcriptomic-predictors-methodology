# R/metrics_labels.R
#
# Harmonise pipeline-option labels across the three dataset folders, so that
# "the same choice" (e.g. RISE screening, whether run in its "paired" mode
# for SDY1276 or its "classic" mode for PREVAC) lands on one shared row when
# pooled across datasets for cross-dataset visualisation.
#
# The three dataset folders were written independently and their option
# labels have drifted slightly (typos, inconsistent wording, and
# dataset-specific qualifiers like "(paired)" vs "(classic)" or "(geneset,
# ...)" vs "(gene, ...)"). Most labels are already identical across datasets
# and need no mapping; `.option_label_overrides` below covers the ones that
# aren't, as of the current 02/03/04 scripts in each dataset folder.
#
# IMPORTANT: this table is a manually-curated snapshot, not derived
# automatically from the option lists themselves. If a dataset folder's
# option labels change (new option, renamed option), this table needs a
# matching update, or that option will silently appear as its own
# (unmerged) row instead of being pooled with its counterpart in other
# datasets. canonicalize_option_label() does NOT warn when this drifts, so
# it's worth spot-checking after editing any 02/03/04 script.
#
# Keys are "<category>|<raw pipeline label>"; values are the canonical
# label used in cross-dataset plots. Entries are grouped, and commented,
# by why the raw labels differ:

.option_label_overrides <- c(

  # --- engineering / engineering_genewise: since 02_compare_engineering.R
  # (geneset-only) and 02b_compare_engineering_genewise.R (gene-level-only)
  # use identically-worded option labels across all three dataset folders
  # (a deliberate choice when they were split into separate, geneset-vs-
  # gene-wise scripts/categories - see each folder's 02b script), no
  # override is needed for either category.

  # --- selection_geneset / selection_genewise: dearseq's "(geneset, ...)"/
  # "(gene, ...)" qualifier is redundant once split by category (geneset-
  # level dearseq only ever appears under selection_geneset, gene-level
  # dearseq only under selection_genewise), and the "paired" vs "classic"
  # qualifier SDY1276 vs PREVAC otherwise used is redundant with dataset
  # (SDY1276 is always paired, PREVAC always classic) - both stripped so
  # "dearseq screening" is treated as one canonical choice regardless of
  # which contrast mechanism a given dataset had available for it.
  "selection_geneset|Dearseq (geneset, alpha = 0.05)"  = "Dearseq (alpha = 0.05)",
  "selection_genewise|Dearseq (gene, alpha = 0.05)"    = "Dearseq (alpha = 0.05)"

  # Everything else (variance/correlation/relative-gain options, RISE,
  # model options, unmapped gene-set aggregation methods) is already
  # labelled identically across all three dataset folders and needs no
  # override - canonicalize_option_label() passes those through unchanged.
)

#' Map a dataset-specific pipeline option label to its cross-dataset
#' canonical label.
#'
#' @param pipeline Character vector of raw `pipeline` labels (as saved by
#'   `R/metrics_io.R::save_comparison_metrics()`).
#' @param category Character vector (recycled/paired with `pipeline`) of
#'   `category` values.
#'
#' @return Character vector of canonical labels, the same length as
#'   `pipeline`. Labels with no entry in `.option_label_overrides` are
#'   returned unchanged.
canonicalize_option_label <- function(pipeline, category) {
  key <- paste(category, pipeline, sep = "|")
  mapped <- .option_label_overrides[key]
  ifelse(is.na(mapped), pipeline, unname(mapped))
}

#' Human-readable display name for a `dataset` tag.
dataset_display_name <- function(dataset) {
  lookup <- c(
    sdy1276_tiv    = "SDY1276 (TIV)",
    prevac_rvsv    = "PREVAC (rVSV)",
    prevac_ad26mva = "PREVAC (Ad26/MVA)"
  )
  unname(ifelse(dataset %in% names(lookup), lookup[dataset], dataset))
}

#' Human-readable display name for a `category` tag.
category_display_name <- function(category) {
  lookup <- c(
    engineering          = "Feature engineering",
    engineering_genewise = "Feature engineering (gene-wise)",
    selection_geneset    = "Feature selection (geneset-level)",
    selection_genewise   = "Feature selection (gene-wise)",
    model                = "Model choice"
  )
  unname(ifelse(category %in% names(lookup), lookup[category], category))
}

#' Human-readable name for what a category's reference pipeline actually
#' does, used to label the reference row added to the cross-dataset
#' heatmaps (see `R/metrics_analysis.R::compute_relative_metrics()`) - e.g.
#' "Elastic net" for the model-choice heatmap, shown there as "Elastic net
#' (reference)".
#'
#' Unlike the other lookups in this file, this one can't be derived from the
#' saved metrics tables themselves (they carry option labels, not pipeline
#' parameter values), so it's a manually-curated match to
#' `R/pipeline_defaults.R`'s `reference_pipeline_params()` (used by
#' `engineering`/`engineering_genewise`/`selection_geneset`/`model` - the
#' gene-level engineering comparison is deliberately evaluated against the
#' SAME mean-aggregation reference as the geneset comparison, since "does
#' skipping aggregation help" is the question that reference answers) and
#' `raw_gene_reference_params()` (used by `selection_genewise`).
reference_option_label <- function(category) {
  lookup <- c(
    engineering           = "Mean aggregation",
    engineering_genewise  = "Mean aggregation",
    selection_geneset     = "Variance (top 7,500)",
    selection_genewise    = "Variance (top 7,500)",
    model                 = "Elastic net"
  )
  unname(ifelse(category %in% names(lookup), lookup[category], "Reference"))
}
