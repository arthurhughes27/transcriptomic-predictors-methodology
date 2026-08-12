# R/build_dataset.R
#
# Build analysis-ready (X, Y, covariates, ...) objects from the merged
# clinical + gene-expression data, for use with `predictomics::predict_cv()`
# and `predictomics::compare_pipelines()`.
#
# Two builders are provided:
#   * build_prediction_dataset(): a single post-vaccination expression
#     timepoint used to predict a later immunogenicity outcome. This is the
#     standard design used throughout the pipeline-comparison scripts.
#   * build_paired_dataset(): baseline + post-vaccination expression stacked
#     per participant. This is only needed when comparing the individual-level
#     fold-change gene-level transformation, which requires both timepoints.
#
# Requires: dplyr, rlang (loaded via dplyr); R/gene_columns.R must be sourced
# first.

#' Build a single-timepoint predictive analysis dataset.
#'
#' Filters the merged data to the requested study/vaccine group(s) and
#' expression timepoint, defines the response (optionally log-transformed),
#' and drops any gene column with missing values (predictomics does not
#' currently handle missing predictor values).
#'
#' @param df_merged Merged clinical + expression data (see load_merged_data()).
#' @param study_vaccine_groups Character vector of `study_vaccine` values to keep.
#' @param timepoint Character scalar/vector giving the `time` value(s) to keep.
#' @param response_col Name of the column holding the immunogenicity outcome.
#' @param covariate_names Baseline covariates to retain alongside gene expression.
#' @param log_response If TRUE (default), the response is log-transformed.
#' @param require_complete_cases If TRUE (default), gene columns with any
#'   missing values are dropped.
#'
#' @return A list with elements X (predictor matrix), Y (response vector),
#'   covariates (data frame), participant_id, and gene_names.
build_prediction_dataset <- function(df_merged,
                                      study_vaccine_groups,
                                      timepoint,
                                      response_col,
                                      covariate_names = c("age", "sex", "race"),
                                      log_response = TRUE,
                                      require_complete_cases = TRUE) {

  gene_names <- get_gene_columns(df_merged)

  df_filtered <- df_merged %>%
    dplyr::filter(.data$study_vaccine %in% study_vaccine_groups,
                  .data$time %in% timepoint) %>%
    dplyr::mutate(response = .data[[response_col]]) %>%
    dplyr::filter(!is.na(.data$response)) %>%
    dplyr::select(participant_id, response,
                  dplyr::any_of(covariate_names),
                  dplyr::any_of(gene_names))

  if (log_response) {
    df_filtered <- dplyr::mutate(df_filtered, response = log(.data$response))
  }

  if (require_complete_cases) {
    df_filtered <- df_filtered %>%
      dplyr::select(participant_id, response, dplyr::any_of(covariate_names),
                    dplyr::where(~ !any(is.na(.))))
  }

  gene_names_present <- intersect(gene_names, names(df_filtered))

  list(
    X = df_filtered %>% dplyr::select(dplyr::any_of(gene_names_present)) %>% as.matrix(),
    Y = df_filtered %>% dplyr::pull(response),
    covariates = df_filtered %>% dplyr::select(dplyr::any_of(covariate_names)),
    participant_id = df_filtered %>% dplyr::pull(participant_id),
    gene_names = gene_names_present
  )
}

#' Build a paired baseline + post-vaccination dataset.
#'
#' Required for the individual-level fold-change gene-level transformation,
#' which needs both timepoints per participant to compute a within-individual
#' fold change. Only participants with an expression sample at *both*
#' timepoints are retained.
#'
#' Unlike build_prediction_dataset(), the immunogenicity response is defined
#' once per participant and repeated across both of that participant's
#' expression rows; `predictomics` uses the returned `timepoint` and
#' `participant_id` vectors to compute per-individual fold changes internally
#' before collapsing to one engineered row per participant.
#'
#' @return A list as for build_prediction_dataset(), plus `timepoint` (0 =
#'   baseline, 1 = post-vaccination).
build_paired_dataset <- function(df_merged,
                                  study_vaccine_groups,
                                  baseline_timepoint,
                                  post_timepoint,
                                  response_col,
                                  covariate_names = c("age", "sex", "race"),
                                  log_response = TRUE,
                                  require_complete_cases = TRUE) {

  gene_names <- get_gene_columns(df_merged)

  # Keep only participants with an expression sample at BOTH timepoints
  participants_paired <- df_merged %>%
    dplyr::filter(.data$study_vaccine %in% study_vaccine_groups,
                  .data$time %in% c(baseline_timepoint, post_timepoint)) %>%
    dplyr::distinct(participant_id, time) %>%
    dplyr::count(participant_id) %>%
    dplyr::filter(n == 2) %>%
    dplyr::pull(participant_id)

  df_filtered <- df_merged %>%
    dplyr::filter(.data$study_vaccine %in% study_vaccine_groups,
                  .data$time %in% c(baseline_timepoint, post_timepoint),
                  .data$participant_id %in% participants_paired) %>%
    dplyr::mutate(response = .data[[response_col]]) %>%
    dplyr::filter(!is.na(.data$response)) %>%
    dplyr::arrange(participant_id, .data$time) %>%
    dplyr::select(participant_id, time, response,
                  dplyr::any_of(covariate_names),
                  dplyr::any_of(gene_names))

  if (log_response) {
    df_filtered <- dplyr::mutate(df_filtered, response = log(.data$response))
  }

  if (require_complete_cases) {
    df_filtered <- df_filtered %>%
      dplyr::select(participant_id, time, response, dplyr::any_of(covariate_names),
                    dplyr::where(~ !any(is.na(.))))
  }

  gene_names_present <- intersect(gene_names, names(df_filtered))

  list(
    X = df_filtered %>% dplyr::select(dplyr::any_of(gene_names_present)) %>% as.matrix(),
    Y = df_filtered %>% dplyr::pull(response),
    covariates = df_filtered %>% dplyr::select(dplyr::any_of(covariate_names)),
    participant_id = df_filtered %>% dplyr::pull(participant_id),
    timepoint = as.integer(df_filtered$time == post_timepoint),
    gene_names = gene_names_present
  )
}
