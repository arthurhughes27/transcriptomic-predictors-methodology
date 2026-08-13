# transcriptomic-predictors-methodology

Repository for files to produce outputs for chapter 5 of my PhD thesis: "Methodological challenges in identifying transcriptomic predictors of immunogenicity"

## Overview

This repository holds the analysis code used to evaluate how methodological choices
(feature engineering, feature selection, and modelling strategies) affect the
predictive performance of transcriptomic (gene expression) predictors of vaccine
immunogenicity. Analyses are built around the `predictomics` R package, which
implements a modular, configurable machine-learning pipeline (engineering →
selection → modelling, with cross-validation) for relating whole-blood gene
expression to post-vaccination antibody response.

Raw and processed data are not stored in this repository for privacy reasons; the
scripts here document the harmonisation and analysis logic and assume the expected
files exist locally under `data-raw/` and `data/` (both gitignored).

## Studies

Data are harmonised across several vaccine studies/cohorts:

- **ebovac2** — Ad26/MVA Ebola vaccine study
- **hamburg** — rVSV Ebola vaccine study
- **prevac** — PREVAC trial (rVSV, Ad26/MVA, and placebo arms)
- **is2 / SDY1276** — Influenza (TIV) vaccine study, sourced from ImmuneSpace (IS2)

## Repository structure

```
analysis/
  preprocessing/          # Per-study cleaning and harmonisation, merging into unified datasets
  descriptive/            # Descriptive summaries/figures of study design and data availability
  application/            # Original application of the predictomics pipeline (being superseded)
    reference_models/         # "Reference" pipeline fit per study/vaccine group
    pipeline_comparison/      # Systematic comparison of engineering/selection/model choices
  pipeline_comparisons/   # Current pipeline-comparison analyses, one folder per dataset
    sdy1276_tiv/               # SDY1276 (TIV) engineering/selection/model comparisons
    prevac_rvsv/               # PREVAC rVSV (+ placebo) engineering/selection/model comparisons
    prevac_ad26mva/            # PREVAC Ad26/MVA (+ placebo) engineering/selection/model comparisons
R/                        # Reusable helper functions shared across analysis scripts
```

### Preprocessing (`analysis/preprocessing/`)

Each study has its own preprocessing script (`preprocessing_<study>.R`) that reads
raw study data and standardises it. `preprocessing_is2_*.R` scripts handle the
multi-file IS2/ImmuneSpace pipeline (clinical, expression, immune response, then
merging). Harmonisation scripts then align naming/coding conventions across
studies:

- `preprocessing_clinical_harmonisation.R` — harmonises clinical/antibody data across studies into `df_clinical_all.rds`
- `preprocessing_merging_harmonisation.R` — merges harmonised clinical data with per-study gene expression into `df_merged_all.rds`
- `preprocessing_BTM.R` / `preprocessing_BG3M.R` — prepare Blood Transcriptional Module (and BG3M) gene sets for use in feature engineering/selection
- `preprocessing_master.R` — runs all preprocessing scripts in the correct order

### Descriptive (`analysis/descriptive/`)

- `study_descriptions.R` — summarises antibody and gene expression sampling availability by study/vaccine group and timepoint

### Application (`analysis/application/`)

- `predictomics_test.R` / `predictomics_test_is2.R` — exploratory scripts exercising the `predictomics` pipeline (feature screening with `SurrogateRank`/RISE, `predict_cv`, `dearseq`/`dgsa_seq`)
- `reference_models/` — fits a single "reference" pipeline configuration per vaccine platform (TIV, Ad26/MVA, rVSV); `reference_model_master.R` runs all three
- `pipeline_comparison/` — uses `predictomics::compare_pipelines()` to systematically compare alternative engineering, selection, and modelling choices against the reference pipeline, across studies/vaccine groups
- `ebola_risemeta.R` — RISE-based meta-analysis across the Ebola vaccine studies

### Pipeline comparisons (`analysis/pipeline_comparisons/`)

Reimplementation of the pipeline-comparison analyses, restructured for
modularity and reuse; intended to eventually replace
`analysis/application/pipeline_comparison/` and
`analysis/application/reference_models/`. Each dataset/vaccine-arm has its
own subfolder (`sdy1276_tiv/`, `prevac_rvsv/`, `prevac_ad26mva/`). Within a
folder, a `01_prepare_data.R` script builds and caches an analysis-ready
dataset, and one script per methodological choice (`02_compare_engineering.R`,
`03_compare_selection.R`, `04_compare_model.R`) loads that cached dataset and
runs one or more `predictomics::compare_pipelines()` comparisons against a
shared reference pipeline. Shared logic (data loading, dataset construction,
reference pipeline definition, comparison execution, plot saving) lives in
`R/` at the repository root rather than being duplicated across scripts.

SDY1276 (TIV) is single-arm, so its selection comparisons rely on RISE/
dearseq's *paired* modes (baseline-vs-post-vaccination expression contrast).
The two PREVAC folders (rVSV and Ad26/MVA) instead include each vaccine
arm's placebo group in every prediction task - `treatment_predictor = FALSE`
throughout keeps treatment out of the model itself, but it powers RISE/
dearseq's *classic* (treatment-vs-placebo) mode in `03_compare_selection.R`,
using only post-vaccination gene expression (no paired dataset needed). See
each folder's README.md for details and dataset-specific assumptions.

Every `02_compare_engineering.R`/`03_compare_selection.R`/`04_compare_model.R`
script also saves its comparison's results table (via
`R/metrics_io.R::save_comparison_metrics()`) to `data/derived/metrics/`,
tagged with `dataset` and `category` (`engineering`, `selection_geneset`,
`selection_genewise`, or `model` - geneset-level and gene-wise feature
selection are kept as separate categories, since they compare against two
different reference pipelines).
`analysis/pipeline_comparisons/collect_metrics.R` pools all of these into a
single long-format data frame (`data/derived/pipeline_comparison_metrics.rds`),
for cross-dataset visualisation of which methodological choices help,
independent of any one dataset's baseline predictability.

## Data conventions

Processed data frames follow shared column naming conventions:

- `participant_id` — unique participant identifier (harmonised from study-specific `pid`)
- `study_accession` — study identifier (e.g. `ebovac2`, `hamburg`, `prevac`, IS2 `SDY` accessions)
- `group` / `group_long` — vaccine arm (short/long form), and `study_vaccine` = `study_accession-group`
- `age`, `sex`, `race`, `ethnicity` — harmonised demographic covariates
- `time` — gene-expression sampling timepoint relative to prime (`P+`) or boost (`B+`) vaccination, in hours (`H`) or days (`D`), e.g. `P+0D`, `P+7D`, `B+3H`
- `ab_p_<day>` / `ab_b_<day>` — antibody titres at a given day post-prime/-boost (e.g. `ab_p_28`, `ab_b_21`); not all studies measure all timepoints, so many are `NA`
- Gene expression columns — one column per gene, named by lowercase gene symbol, spanning alphabetically from `a1cf`/`a1bg` to `zzz3`; selected in code via `dplyr::select(a1cf:zzz3)` or similar ranges

Gene set objects (BTM/BG3M) are lists containing `genesets` (gene symbol members),
`geneset.names`, `geneset.descriptions`, and `geneset.names.descriptions`.

## Requirements

Analyses depend on the `predictomics` package (not on CRAN) plus `tidyverse`,
`SurrogateRank`, `dearseq`, `GSA`, `fs`, `future`, and related packages used for
plotting and parallelisation.
