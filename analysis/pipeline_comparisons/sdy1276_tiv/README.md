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
   — each independently loads the cached dataset and runs one
   `predictomics::compare_pipelines()` comparison (or, for feature
   engineering, two — see below), saving a figure to
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
- **Feature selection: RISE and dearseq omitted.** Both require a two-group
  contrast (treated vs. control for RISE; a `variables2test` group contrast
  for dearseq), which SDY1276 TIV doesn't have (single arm, single
  post-vaccination expression timepoint used as predictors). These are
  compared in `03_compare_selection.R`'s header comment as a planned addition
  for the placebo-controlled Ebola/PREVAC datasets.
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
