# analysis/descriptive/dataset_characteristics_table.R
#
# Produce a manuscript-ready LaTeX table summarising basic characteristics
# of all six datasets used in this chapter: the three exploration
# (within-study best-pipeline search) cohorts and the three independent
# external-validation cohorts for the same three vaccines (see
# analysis/pipeline_comparisons/ and analysis/external_validation/,
# respectively).
#
# For each cohort, reports: cohort name, vaccine (+ placebo, where
# applicable), sample size (broken down into active/placebo arms where a
# placebo group exists), the gene-expression timepoint used as predictors,
# the baseline covariates from {age, sex, race} that are both available and
# non-constant in that cohort (computed from the data, not assumed - e.g.
# PREVAC's `race` is known to be constant "Unknown" and so is expected to
# drop out here), and the immunogenicity response (assay/form and
# timepoint).
#
# Sample sizes and covariate availability are computed directly from the
# same R/build_dataset.R::build_prediction_dataset() builder every other
# analysis script in this repository uses, restricted to rows with a
# non-missing response (i.e. the sample actually usable for prediction),
# so this table's N always matches what analysis/pipeline_comparisons/ and
# analysis/external_validation/ actually fit on.
#
# Run this after analysis/preprocessing/ (i.e. whenever data/df_merged_all.rds
# exists).

library(dplyr)
library(fs)

source(fs::path("R", "data_io.R"))
source(fs::path("R", "gene_columns.R"))
source(fs::path("R", "build_dataset.R"))
source(fs::path("R", "pipeline_defaults.R"))

df_merged_all <- load_merged_data()

# --- Cohort specifications ---------------------------------------------------
#
# GE timepoints and response columns/labels mirror exactly what each
# cohort's own 01_prepare_data.R / external-validation script uses (see
# analysis/pipeline_comparisons/*/01_prepare_data.R and
# analysis/external_validation/validate_*.R for the sourcing of each
# assay/timepoint choice - repeated here rather than re-derived, since this
# script's only job is to describe the data, not re-fit anything).

cohort_specs <- list(
  list(
    section = "Exploration cohorts",
    cohort = "SDY1276", vaccine = "Influenza (TIV)",
    study_vaccine_groups = "SDY1276-Influenza (IN)", treatment_arm = NULL,
    timepoint = "P+1D", timepoint_label = "Day 1",
    response_col = "ab_p_28",
    response_label = "Neutralizing Ab titer (log), day 28"
  ),
  list(
    section = "Exploration cohorts",
    cohort = "PREVAC", vaccine = "Ebola (Ad26/MVA) + placebo",
    study_vaccine_groups = c("prevac-Ad26MVA", "prevac-placebo"),
    treatment_arm = "prevac-Ad26MVA",
    timepoint = "P+7D", timepoint_label = "Day 7",
    response_col = "ab_p_180",
    response_label = "Anti-Ebola-GP binding Ab titer (log), day 180"
  ),
  list(
    section = "Exploration cohorts",
    cohort = "PREVAC", vaccine = "Ebola (rVSV) + placebo",
    study_vaccine_groups = c("prevac-rVSV", "prevac-placebo"),
    treatment_arm = "prevac-rVSV",
    timepoint = "P+7D", timepoint_label = "Day 7",
    response_col = "ab_p_180",
    response_label = "Anti-Ebola-GP binding Ab titer (log), day 180"
  ),
  list(
    section = "External validation cohorts",
    cohort = "SDY80", vaccine = "Influenza (TIV)",
    study_vaccine_groups = "SDY80-Influenza (IN)", treatment_arm = NULL,
    timepoint = "P+1D", timepoint_label = "Day 1",
    response_col = "ab_p_28",
    response_label = "Neutralizing Ab titer (log), day 28"
  ),
  list(
    section = "External validation cohorts",
    cohort = "EBOVAC2", vaccine = "Ebola (Ad26/MVA)",
    study_vaccine_groups = "ebovac2-Ad26MVA", treatment_arm = NULL,
    timepoint = "P+7D", timepoint_label = "Day 7",
    response_col = "ab_p_365",
    response_label = "Anti-Ebola-GP binding Ab titer (log), day 365"
  ),
  list(
    section = "External validation cohorts",
    cohort = "Hamburg", vaccine = "Ebola (rVSV)",
    study_vaccine_groups = "hamburg-rVSV", treatment_arm = NULL,
    timepoint = "P+7D", timepoint_label = "Day 7",
    response_col = "ab_p_180",
    response_label = "Anti-Ebola-GP binding Ab titer (log), day 180"
  )
)

# --- Per-cohort N and covariate availability ---------------------------------

covariate_display_names <- c(age = "Age", sex = "Sex", race = "Race")

describe_cohort <- function(spec) {

  built <- build_prediction_dataset(
    df_merged = df_merged_all,
    study_vaccine_groups = spec$study_vaccine_groups,
    timepoint = spec$timepoint,
    response_col = spec$response_col,
    covariate_names = default_covariates,
    treatment_arm = spec$treatment_arm,
    log_response = TRUE
  )

  n_total <- length(built$participant_id)

  n_label <- if (!is.null(spec$treatment_arm)) {
    n_active <- sum(built$treatment == 1)
    n_control <- sum(built$treatment == 0)
    sprintf("%d (%d active / %d placebo)", n_total, n_active, n_control)
  } else {
    as.character(n_total)
  }

  # A covariate counts as an available, informative baseline predictor only
  # if it's both present (build_prediction_dataset() already restricts to
  # covariates that exist in the merged data) AND non-constant (>1 distinct
  # non-missing value) in this specific cohort - e.g. PREVAC's `race` is
  # recorded as a constant "Unknown" for every participant (see
  # analysis/preprocessing/preprocessing_clinical_harmonisation.R) and so is
  # expected to drop out here, without needing to special-case it by name.
  covariate_names_available <- vapply(names(covariate_display_names), function(nm) {
    nm %in% names(built$covariates) &&
      length(unique(stats::na.omit(built$covariates[[nm]]))) > 1
  }, logical(1))

  covariates_label <- if (any(covariate_names_available)) {
    paste(covariate_display_names[names(covariate_display_names)[covariate_names_available]], collapse = ", ")
  } else {
    "None"
  }

  data.frame(
    section = spec$section,
    cohort = spec$cohort,
    vaccine = spec$vaccine,
    n = n_label,
    timepoint = spec$timepoint_label,
    covariates = covariates_label,
    response = spec$response_label,
    stringsAsFactors = FALSE
  )
}

table_data <- do.call(rbind, lapply(cohort_specs, describe_cohort))

cat("Dataset characteristics:\n")
print(table_data)

# --- Render as a manuscript-ready LaTeX table --------------------------------
#
# Hand-built (rather than via knitr::kable()/xtable()) to keep this script
# dependency-free and give full control over the compact, booktabs-style
# layout expected in a manuscript. Requires the `booktabs` LaTeX package
# (\toprule/\midrule/\bottomrule) in the including document's preamble.

escape_latex <- function(x) {
  # fixed = TRUE throughout: each replacement is a plain string built with
  # paste0() before gsub() ever sees it, so there's no regex-replacement
  # backreference syntax to get subtly wrong.
  for (special_char in c("&", "%", "$", "#", "_", "{", "}")) {
    x <- gsub(special_char, paste0("\\", special_char), x, fixed = TRUE)
  }
  x
}

format_data_row <- function(row) {
  paste0(
    paste(escape_latex(c(row$cohort, row$vaccine, row$n, row$timepoint, row$covariates, row$response)), collapse = " & "),
    " \\\\"
  )
}

format_section_row <- function(section_label) {
  sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", escape_latex(section_label))
}

body_lines <- character(0)
for (section_label in unique(table_data$section)) {
  body_lines <- c(body_lines, format_section_row(section_label), "\\addlinespace[2pt]")
  section_rows <- table_data[table_data$section == section_label, , drop = FALSE]
  body_lines <- c(body_lines, vapply(seq_len(nrow(section_rows)), function(i) format_data_row(section_rows[i, ]), character(1)))
  body_lines <- c(body_lines, "\\addlinespace[4pt]")
}
body_lines <- body_lines[-length(body_lines)]  # drop trailing spacer after the last section

latex_lines <- c(
  "% Requires \\usepackage{booktabs} in the preamble.",
  "\\begin{table}[t]",
  "  \\centering",
  "  \\caption{Characteristics of the exploration and external-validation cohorts.}",
  "  \\label{tab:dataset-characteristics}",
  "  \\resizebox{\\textwidth}{!}{%",
  "  \\begin{tabular}{llllll}",
  "    \\toprule",
  "    Cohort & Vaccine & N & Predictor timepoint & Baseline covariates & Response \\\\",
  "    \\midrule",
  paste0("    ", body_lines),
  "    \\bottomrule",
  "  \\end{tabular}%",
  "  }",
  "\\end{table}"
)

table_path <- fs::path("output", "tables")
fs::dir_create(table_path)

writeLines(latex_lines, fs::path(table_path, "dataset_characteristics.tex"))

cat("\nLaTeX table written to ", fs::path(table_path, "dataset_characteristics.tex"), "\n", sep = "")

rm(list = ls())
