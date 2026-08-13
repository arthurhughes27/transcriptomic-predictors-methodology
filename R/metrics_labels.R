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

  # --- engineering: SDY1276 splits gene-level transform vs aggregation
  # across two compare_pipelines() calls (see
  # sdy1276_tiv/02_compare_engineering.R); "Gene-wise" (from the
  # aggregation call: z-score, no aggregation) and "Gene-level: z-score"
  # (from the gene-level call, same config, evaluated on a differently
  # -restricted sample - see that script's header) both denote the same
  # underlying choice as PREVAC's "Gene-level: z-score" (its only
  # single-call comparison). All three canonicalize to one label; when a
  # dataset contributes more than one raw label to it (SDY1276 does),
  # compute_relative_metrics() averages them - see that function's docs.
  "engineering|Gene-wise"                    = "Gene-level: z-score (no aggregation)",
  "engineering|Gene-level: z-score"          = "Gene-level: z-score (no aggregation)",

  # PREVAC's two folders independently named the "no transform, no
  # aggregation" option differently ("no transformation" vs "none");
  # SDY1276 doesn't have an equivalent option in its engineering
  # comparisons at all (so this canonical label ends up PREVAC-only,
  # covering 2 of 3 datasets - expected, not a bug).
  "engineering|Gene-level: no transformation" = "Gene-level: none (no aggregation)",
  "engineering|Gene-level: none"              = "Gene-level: none (no aggregation)",

  # SDY1276-only (needs paired baseline+post data, which PREVAC's folders
  # don't build - see prevac_*/01_prepare_data.R). Renamed only for a
  # consistent "(no aggregation)" suffix alongside its siblings above.
  "engineering|Gene-level: fold-change" = "Gene-level: fold-change (no aggregation)",

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
    engineering         = "Feature engineering",
    selection_geneset   = "Feature selection (geneset-level)",
    selection_genewise  = "Feature selection (gene-wise)",
    model               = "Model choice"
  )
  unname(ifelse(category %in% names(lookup), lookup[category], category))
}
