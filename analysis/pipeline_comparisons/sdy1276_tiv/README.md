# Pipeline comparisons: SDY1276 (TIV)

Scripts implementing the Chapter 5 pipeline comparisons (feature engineering,
feature selection, model choice) on the SDY1276 influenza (TIV) dataset. This
is the first dataset covered by the new `analysis/pipeline_comparisons/`
structure; the plan is to add one sibling folder per dataset (e.g. PREVAC,
Ebola) here as the analysis is extended, and eventually retire the older
scripts under `analysis/application/`.

## Running order

1. `01_prepare_data.R` — builds and caches the analysis-ready dataset(s) to
   `data/derived/sdy1276_tiv_analysis_data.rds`. Run this first, and re-run it
   if the upstream processed data (`data/df_merged_all.rds`,
   `data/BTM_processed.rds`) changes.
2. `02_compare_engineering.R`, `03_compare_selection.R`, `04_compare_model.R`
   — each independently loads the cached dataset and runs one or more
   `predictomics::compare_pipelines()` comparisons (`04_compare_model.R`
   runs one; `02_compare_engineering.R` and `03_compare_selection.R` each
   run two — see below), saving a figure per comparison to
   `output/figures/pipeline_comparisons/sdy1276_tiv/`.

## Design notes / decisions made without the ability to run R

This chapter was written and reviewed in an environment without R installed
and without access to the (privacy-protected) underlying data, so none of
these scripts have been executed. A few choices are worth double-checking
before relying on the output:

- **Reference pipeline.** A z-scored gene-level transform, mean aggregation of
  Blood Transcriptional Modules, a relaxed variance pre-filter (top 7,500
  gene sets, to keep model fitting tractable), and an elastic net model
  (`R/pipeline_defaults.R::reference_pipeline_params()`), matching the
  reference used in the older `analysis/application/` scripts.
  `reference_pipeline_params()` takes `genesets` as an explicit argument
  (rather than reading it from a global variable), so every call site must
  pass it, e.g. `reference_pipeline_params(genesets)`.
- **Covariates.** Follows Chapter 5, Section 5.4 literally: age, sex, race.
  The older scripts additionally included baseline antibody titer (`ab_p_0`)
  as a covariate; this is omitted here since the chapter text does not list
  it. If you intended `ab_p_0` to be included, adjust
  `R/pipeline_defaults.R::default_covariates`.
- **Feature selection: classic RISE/dearseq omitted; paired RISE/dearseq
  included, split by engineering scale.** RISE's and dearseq's *classic*
  modes need a treatment-vs-control contrast, which SDY1276 TIV doesn't have
  (single arm); those are a planned addition for the placebo-controlled
  Ebola/PREVAC datasets. Their *paired* modes instead contrast each
  participant's baseline (day 0) vs. day-1 expression, which this dataset
  does have. `03_compare_selection.R` now runs two `option_type =
  "selection"` comparisons directly against the paired dataset: (1)
  geneset-level, against the usual mean-aggregated reference, including
  dearseq's paired mode at `dearseq_level = "geneset"` (RISE is excluded
  here - `predictomics::predict_cv()` rejects `rise_paired = TRUE` combined
  with geneset engineering outright, since paired RISE screening always
  operates on the raw gene-level matrix); (2) gene-wise, against
  `R/pipeline_defaults.R::raw_gene_reference_params()` (z-score only, no
  aggregation), including both paired RISE and dearseq's paired mode at
  `dearseq_level = "gene"`, with `top_n` thresholds rescaled ~100x from (1)
  to the ~20,000-raw-gene scale. Both comparisons rely on predictomics'
  paired row-discard parity handling (`compare_pipelines()`'s
  `.discards_pretreatment_rows()`/row-parity logic): any option that
  discards pre-treatment rows internally (`rise_paired`, `dearseq_mode =
  "paired"`) screens on both arms then models on post-treatment rows only,
  while every other pipeline in the same call - including the reference and
  baseline - is automatically restricted to post-treatment (day-1) rows, so
  mixed paired/non-paired comparisons are safe in a single call. (This
  replaces an earlier, more convoluted two-step workaround that was needed
  before predictomics generalized this row-parity handling beyond
  `option_type = "engineering"`.)
- **Gene-level fold-change.** Requires paired baseline (day 0) + post (day 1)
  expression per participant, unlike the rest of the pipeline (day-1
  expression only). It is therefore evaluated as a separate
  `compare_pipelines()` call/figure within `02_compare_engineering.R`, rather
  than being merged into the main aggregation-method comparison.
- **`method = "lm"` for linear regression** in `04_compare_model.R` is an
  assumed method string (by analogy with `"lasso"`/`"ridge"`/`"ranger"`/
  `"svr"`/`"glmnet"`) — please verify against the `predictomics` source.

## Reusable helpers

Shared logic lives in `R/` at the repository root (sourced via
`source(fs::path("R", "<file>.R"))`):

- `data_io.R` — load the merged data and processed gene sets.
- `gene_columns.R` — robustly identify gene-expression columns.
- `build_dataset.R` — filter/shape the merged data into `predictomics`-ready
  `(X, Y, covariates, ...)` objects (single-timepoint and paired-timepoint
  variants).
- `pipeline_defaults.R` — the shared reference pipeline and default
  covariate set.
- `run_comparison.R` — wraps `predictomics::compare_pipelines()` with managed
  `future` parallelisation.
- `plotting.R` — consistent figure saving.
