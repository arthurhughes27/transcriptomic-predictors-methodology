# R/gene_prefilter.R
#
# One-off, GLOBAL (not per-CV-fold) unsupervised variance pre-filter,
# applied upfront to the raw gene matrix passed into the SUPPLEMENTARY,
# gene-wise selection comparison and best-pipeline search
# (analysis/supplementary/*/02_compare_selection_genewise.R,
# R/best_pipeline_search.R::find_best_pipeline_genewise()'s selection
# round, both called with a pre-filtered X).
#
# Restricting to the top `top_n` highest-variance genes BEFORE any pipeline
# is fit means every compared selection method (correlation, RISE,
# dearseq, and especially relative_gain's per-gene nested CV - by far the
# slowest of the six, since it fits `relative_gain_inner_folds` separate
# lm() models per candidate gene, per outer CV fold) only ever scans
# `top_n` genes, not the full ~20,000 - directly cutting their cost by
# roughly the same factor as the reduction in gene count.
#
# This is safe to compute ONCE, outside of cross-validation, because it is
# entirely unsupervised (a function of X alone, never Y) - unlike a
# supervised filter (relative_gain/correlation/RISE/dearseq themselves,
# which all look at Y and so MUST be recomputed independently per training
# fold to avoid leaking test-fold information; see Ambroise & McLachlan,
# PNAS 2002, for the general principle that unsupervised filtering, unlike
# supervised filtering, does not need to be nested inside CV to stay
# leakage-free). It is the same category of preprocessing as this repo's
# per-fold "variance" selection step - just applied once, globally, rather
# than being recomputed inside every fold.
#
# By default this uses the SAME top_n as
# R/pipeline_defaults.R::raw_gene_reference_params()'s own per-fold
# variance selection step, so once applied upstream, the gene-wise
# reference pipeline's own variance filter becomes a (harmless) no-op -
# ncol(X) already equals top_n, so predictomics' run_selection() skips
# scoring entirely and keeps every column (see that function's docs) - and
# the reference effectively reduces to "this fixed panel of top_n genes, no
# further selection", while every compared OPTION still applies its own
# additional filter within that same panel.

#' Restrict `X` to its `top_n` highest-variance columns, computed once on
#' the full matrix (not recomputed per CV fold).
#'
#' @param X A numeric feature matrix, genes/features in columns.
#' @param top_n Number of highest-variance columns to keep. If `ncol(X) <=
#'   top_n`, `X` is returned unchanged (nothing to filter).
#'
#' @return `X`, subset to its `top_n` highest-variance columns (original
#'   column order preserved).
prefilter_by_variance <- function(X, top_n = 5000) {
  if (ncol(X) <= top_n) {
    return(X)
  }
  variances <- apply(X, 2, stats::var)
  keep      <- colnames(X)[order(variances, decreasing = TRUE)[seq_len(top_n)]]
  X[, keep, drop = FALSE]
}
