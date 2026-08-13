# Pipeline comparisons: PREVAC Ad26/MVA (+ placebo)

Scripts implementing the Chapter 5 pipeline comparisons (feature
engineering, feature selection, model choice) on the PREVAC Ad26/MVA Ebola
vaccine arm, with the placebo arm included as extra examples of low/no
response. Sibling folder to `../sdy1276_tiv/` and `../prevac_rvsv/` within
`analysis/pipeline_comparisons/`.

## Running order

1. `01_prepare_data.R` — builds and caches the analysis-ready dataset to
   `data/derived/prevac_ad26mva_analysis_data.rds`. Run this first, and
   re-run it if the upstream processed data (`data/df_merged_all.rds`,
   `data/BTM_processed.rds`) changes.
2. `02_compare_engineering.R`, `03_compare_selection.R` (two comparisons),
   `04_compare_model.R` — each independently loads the cached dataset and
   runs one or more `predictomics::compare_pipelines()` comparisons, saving
   a figure per comparison to
   `output/figures/pipeline_comparisons/prevac_ad26mva/`.

## How this differs from `sdy1276_tiv/`

This folder is structurally identical to `../prevac_rvsv/` (same design,
different vaccine arm) - see that folder's README for the full comparison
against `sdy1276_tiv/`'s single-arm, paired-mode design. In short:

- The placebo arm is included in every prediction task
  (`study_vaccine_groups = c("prevac-Ad26MVA", "prevac-placebo")` in
  `01_prepare_data.R`), with a binary `treatment` indicator (1 = Ad26/MVA,
  0 = placebo) built for use as the "classic" RISE/dearseq screening
  contrast in `03_compare_selection.R`. `treatment_predictor = FALSE`
  throughout keeps treatment out of the model itself.
- Only post-vaccination (day-7) data is used - no paired dataset is built,
  so `02_compare_engineering.R` compares the gene-level and gene-set
  aggregation axes together in one call/figure, and
  `03_compare_selection.R` uses RISE/dearseq's "classic" (treatment-vs-
  placebo) mode rather than "paired", needing none of predictomics' paired
  row-discard-parity machinery.
- `03_compare_selection.R` is still split into geneset-level and gene-wise
  comparisons, since RISE (always raw gene-level) and dearseq's gene-level
  mode are never compatible with gene-set-aggregated engineering.

## Design notes / decisions made without the ability to run R

This was written and reviewed without R installed and without access to the
(privacy-protected) underlying data, so none of these scripts have been
executed. Worth double-checking:

- **Response and predictor timepoint.** `ab_p_180` (log-transformed) as
  response, day-7 (`P+7D`) gene expression as predictors - matches Chapter
  5, Section 5.4's stated Ebolavirus design (same as `../prevac_rvsv/`) and
  the covariates used elsewhere in this chapter (age, sex, race; see
  `R/pipeline_defaults.R::default_covariates`).
- **`race` is a constant "Unknown" for PREVAC** (see
  `analysis/preprocessing/preprocessing_clinical_harmonisation.R`) - kept
  in the covariate set for consistency with the rest of the chapter, but
  it's degenerate here and may need dropping if it causes a design-matrix
  problem when actually run.
- **`top_n` values.** Reused the current scale from `sdy1276_tiv/` (25/100
  at the geneset level, 500 at the gene level) as a starting point; PREVAC's
  sample size and gene/geneset counts after complete-case filtering may
  differ, so these may need retuning. The Ad26/MVA arm in particular may
  have a different sample size than rVSV - worth checking `01_prepare_data.R`'s
  printed participant counts for both arms before assuming shared settings
  are appropriate for both.
- **`method = "lm"` for linear regression** in `04_compare_model.R` is an
  assumed method string - see `sdy1276_tiv/04_compare_model.R`'s note.

## Reusable helpers

Shared logic lives in `R/` at the repository root (see
`../sdy1276_tiv/README.md` for the full list), including
`build_dataset.R::build_prediction_dataset()`'s `treatment_arm` argument,
used by this folder's `01_prepare_data.R`.
