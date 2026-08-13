# Pipeline comparisons: PREVAC rVSV (+ placebo)

Scripts implementing the Chapter 5 pipeline comparisons (feature
engineering, feature selection, model choice) on the PREVAC rVSV Ebola
vaccine arm, with the placebo arm included as extra examples of low/no
response. Sibling folder to `../sdy1276_tiv/` and `../prevac_ad26mva/`
within `analysis/pipeline_comparisons/`.

## Running order

1. `01_prepare_data.R` — builds and caches the analysis-ready dataset to
   `data/derived/prevac_rvsv_analysis_data.rds`. Run this first, and re-run
   it if the upstream processed data (`data/df_merged_all.rds`,
   `data/BTM_processed.rds`) changes.
2. `02_compare_engineering.R`, `03_compare_selection.R` (two comparisons),
   `04_compare_model.R` — each independently loads the cached dataset and
   runs one or more `predictomics::compare_pipelines()` comparisons, saving
   a figure per comparison to
   `output/figures/pipeline_comparisons/prevac_rvsv/`.

## How this differs from `sdy1276_tiv/`

SDY1276 (TIV) is a single-arm study, so its `03_compare_selection.R` had to
rely on RISE/dearseq's *paired* modes (contrasting each participant's
baseline vs. post-vaccination expression) to get a two-group screening
contrast at all, and needed predictomics' paired row-discard-parity
handling to mix those with non-paired selection methods safely. PREVAC has
a placebo arm, which changes the design here:

- **Placebo included, not as a predictor.** Every prediction task here
  includes the placebo arm (`study_vaccine_groups = c("prevac-rVSV",
  "prevac-placebo")` in `01_prepare_data.R`), giving more examples of
  low/no response. A binary `treatment` indicator (1 = rVSV, 0 = placebo) is
  built alongside `X`/`Y`/`covariates`, but `treatment_predictor = FALSE`
  throughout (the default in `R/run_comparison.R::run_pipeline_comparison()`)
  keeps treatment out of the model itself; `treatment` is only passed to
  `03_compare_selection.R`'s comparisons, where it's used as the "classic"
  RISE/dearseq screening contrast.
- **Only post-vaccination (day-7) data is used.** PREVAC does have paired
  baseline/day-7 expression, but since the placebo arm already supplies a
  treatment contrast, "classic" RISE/dearseq modes are used instead of
  "paired" ones - so `01_prepare_data.R` only builds a single-timepoint
  dataset (no `build_paired_dataset()` call, no `individual_id`/`timepoint`
  anywhere in these scripts). This also lets `02_compare_engineering.R`
  compare the gene-level ("none"/"z-score") and gene-set aggregation axes
  together in one call/figure, rather than needing to split across two like
  `sdy1276_tiv/02_compare_engineering.R` does for its fold-change option.
- **`03_compare_selection.R`'s two comparisons (geneset-level and
  gene-wise) are still split**, for the same reason as in `sdy1276_tiv/`:
  RISE always screens the raw gene-level matrix (never compatible with
  gene-set aggregation), and dearseq's gene-level mode
  (`dearseq_level = "gene"`) is likewise never compatible with gene-set
  aggregation - only dearseq's geneset-level mode
  (`dearseq_level = "geneset"`) is. Since neither comparison uses a paired
  selection mode here, predictomics' row-discard-parity machinery isn't
  needed at all - both are ordinary `option_type = "selection"` calls on the
  single dataset.

## Design notes / decisions made without the ability to run R

This was written and reviewed without R installed and without access to the
(privacy-protected) underlying data, so none of these scripts have been
executed. Worth double-checking:

- **Response and predictor timepoint.** `ab_p_180` (log-transformed) as
  response, day-7 (`P+7D`) gene expression as predictors - matches Chapter
  5, Section 5.4's stated Ebolavirus design and the covariates used
  elsewhere in this chapter (age, sex, race; see
  `R/pipeline_defaults.R::default_covariates`).
- **`race` is a constant "Unknown" for PREVAC** (see
  `analysis/preprocessing/preprocessing_clinical_harmonisation.R`) - kept
  in the covariate set for consistency with the rest of the chapter, but
  it's degenerate here and may need dropping if it causes a design-matrix
  problem when actually run.
- **`top_n` values.** Reused the current scale from `sdy1276_tiv/` (25/100
  at the geneset level, 500 at the gene level) as a starting point; PREVAC's
  sample size and gene/geneset counts after complete-case filtering may
  differ, so these may need retuning.
- **`method = "lm"` for linear regression** in `04_compare_model.R` is an
  assumed method string - see `sdy1276_tiv/04_compare_model.R`'s note.

## Reusable helpers

Shared logic lives in `R/` at the repository root (see
`../sdy1276_tiv/README.md` for the full list). New here:
`build_dataset.R::build_prediction_dataset()`'s `treatment_arm` argument,
which builds the binary `treatment` indicator used by this folder's
`03_compare_selection.R`.
