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
  preprocessing/     # Per-study cleaning and harmonisation, merging into unified datasets
  descriptive/        # Descriptive summaries/figures of study design and data availability
  application/        # Application of the predictomics pipeline
    reference_models/       # "Reference" pipeline fit per study/vaccine group
    pipeline_comparison/    # Systematic comparison of engineering/selection/model choices
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
