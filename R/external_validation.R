# R/external_validation.R
#
# Shared helpers for validating each vaccine's best pipeline (see
# R/best_pipeline_search.R) on an independent external dataset:
#   PREVAC Ad26/MVA (+ placebo) -> EBOVAC2 ("ebovac2-Ad26MVA")
#   PREVAC rVSV (+ placebo)     -> Hamburg ("hamburg-rVSV")
#   SDY1276 (TIV)               -> SDY80   ("SDY80-Influenza (IN)")
#
# Unlike the within-study best-pipeline search, validation does NOT apply an
# already-fitted model directly to new data (expression scales differ too
# much across studies for that to be meaningful), and it does NOT re-run
# feature selection on the validation data either: the engineering and model
# choices are carried over unchanged, but the SET OF FEATURES selected by
# the winning pipeline on its DISCOVERY dataset is fixed and simply applied
# to the validation dataset, with no further selection step. Cross-validated
# performance is then re-measured from scratch on the validation dataset
# using that fixed feature panel plus the discovery's engineering/model
# choices.
#
# Also unlike the within-study search, no placebo/treatment contrast or
# baseline+post pairing is used here at all: none of EBOVAC2, Hamburg, or
# SDY80 has a placebo arm, and fixing the feature panel upfront removes the
# only reason RISE/dearseq's "classic" (treatment-contrast) or "paired"
# (baseline-vs-post) modes were needed in the first place. Validation
# therefore always uses a single post-vaccination-timepoint dataset (built
# via R/build_dataset.R::build_prediction_dataset(), the same builder
# sdy1276_tiv/01_prepare_data.R uses for its `single` dataset), and always
# fits with `selection_params = NULL`.

#' Recover the fixed set of features selected by a `find_best_pipeline()`
#' winner on its DISCOVERY dataset (i.e. the dataset it was searched on -
#' PREVAC Ad26/MVA, PREVAC rVSV, or SDY1276), for later re-use on the
#' validation dataset.
#'
#' Re-runs `predict_cv()` once, with `outside_cv = TRUE`, using the winning
#' engineering/selection/model params on the FULL discovery dataset. This is
#' the same "apply engineering/selection once to the full dataset" logic
#' `compare_pipelines()`/`find_best_pipeline()` never used during the search
#' itself (which is why the leakage warning `outside_cv = TRUE` emits does
#' not apply to how the result is used here): the predictions and CV
#' performance from this call are discarded entirely - only the resulting
#' fixed feature panel (`dearseq_selection$selected_features` for dearseq,
#' `outside_cv_selection$selected_features` for every other filter method)
#' is kept, precisely because that panel is deliberately meant to be
#' re-applied to an external dataset rather than re-discovered on it.
#'
#' @param best A `find_best_pipeline()` result for the discovery dataset.
#' @param X,Y,covariates,treatment The discovery dataset's SINGLE
#'   (post-vaccination-only) X/Y/covariates(/treatment) - i.e. exactly what
#'   was passed to `find_best_pipeline()` for this dataset.
#'
#' @return A character vector of feature names (gene names if the winning
#'   engineering does not aggregate; gene-set names if it does), in the same
#'   column-name space as the winning pipeline's engineered features. `NULL`
#'   if the winning pipeline had no selection step at all (`best$selection_params`
#'   is `NULL`), meaning every engineered feature should be carried over
#'   unchanged.
get_discovery_selected_features <- function(best, X, Y, covariates, treatment = NULL) {

  if (is.null(best$selection_params)) {
    return(NULL)
  }

  fit <- suppressWarnings(
    predictomics::predict_cv(
      Y = Y, X = X, covariates = covariates,
      treatment = treatment, treatment_predictor = FALSE,
      engineering_params = best$engineering_params,
      selection_params   = best$selection_params,
      model_params        = best$model_params,
      outside_cv = TRUE,
      cv_type = "kfold", folds = 5, seed = 12345,
      verbose = FALSE
    )
  )

  if (!is.null(fit$dearseq_selection)) {
    fit$dearseq_selection$selected_features
  } else {
    fit$outside_cv_selection$selected_features
  }
}

#' Build the engineering/selection specification used to fit the validation
#' dataset: the winning pipeline's engineering choice, restricted (if it
#' aggregates into gene sets) to only the discovery-selected gene sets, with
#' selection itself switched off (`selection_params = NULL`) since the
#' feature panel is already fixed.
#'
#' @param best A `find_best_pipeline()` result for the discovery dataset.
#' @param selected_features As returned by `get_discovery_selected_features()`
#'   (`NULL` if the winning pipeline had no selection step).
#' @param genesets Named list of gene sets (see `R/data_io.R::load_genesets()`).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{engineering_params}{`best$engineering_params`, with `genesets`
#'       restricted to the selected gene sets when engineering aggregates
#'       and a selection step fixed a subset.}
#'     \item{fixed_gene_features}{Only set when engineering does NOT
#'       aggregate and a selection step fixed a subset: the gene names to
#'       subset the validation X matrix to BEFORE fitting (since, unlike the
#'       gene-set case, there's no `engineering_params` field to restrict
#'       instead). `NULL` otherwise.}
#'   }
restrict_engineering_for_validation <- function(best, selected_features, genesets) {

  engineering_params <- best$engineering_params
  aggregates <- !is.null(engineering_params$genesets)

  if (is.null(selected_features)) {
    return(list(engineering_params = engineering_params, fixed_gene_features = NULL))
  }

  if (aggregates) {
    engineering_params$genesets <- genesets[intersect(names(genesets), selected_features)]
    list(engineering_params = engineering_params, fixed_gene_features = NULL)
  } else {
    list(engineering_params = engineering_params, fixed_gene_features = selected_features)
  }
}

#' Run the validation fit: `predict_cv()` with a managed `future` parallel
#' backend (mirroring `R/run_comparison.R::run_pipeline_comparison()`),
#' `selection_params = NULL` throughout (see this file's header).
run_validation_fit <- function(X, Y, covariates, engineering_params, model_params,
                                cv_type = "kfold", folds = 5, seed = 12345,
                                n_workers = 6, verbose = TRUE) {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  predictomics::predict_cv(
    Y = Y, X = X, covariates = covariates,
    engineering_params = engineering_params,
    selection_params   = NULL,
    model_params        = model_params,
    outside_cv = FALSE,
    cv_type = cv_type, folds = folds, seed = seed,
    verbose = verbose
  )
}

#' One-line description of a validation fit's specification: the discovery
#' pipeline's engineering/model choice (unchanged), and how many features
#' its fixed, carried-over selection panel contains (or "None" if the
#' winning pipeline had no selection step to carry over).
#'
#' @param best The DISCOVERY `find_best_pipeline()` result.
#' @param n_features_used Number of features actually used to fit the
#'   validation model (i.e. `ncol()` of the matrix passed to `predict_cv()`),
#'   after any restriction to genes/gene sets actually present in the
#'   validation dataset - may be smaller than the discovery selection's
#'   panel size if some features aren't available there.
validation_subtitle <- function(best, n_features_used) {

  engineering_label <- if (best$winners$engineering$role == "reference") {
    reference_option_label("engineering")
  } else {
    best$winners$engineering$pipeline
  }

  model_label <- if (best$winners$model$role == "reference") {
    reference_option_label("model")
  } else {
    best$winners$model$pipeline
  }

  selection_label <- if (is.null(best$selection_params)) {
    "None"
  } else {
    paste0(n_features_used, " features fixed from discovery selection")
  }

  paste0(
    "Engineering: ", engineering_label,
    " | Selection: ", selection_label,
    " | Model: ", model_label
  )
}

#' Combine a validation fit's CV-prediction plot and selection-frequency
#' stability plot into a single side-by-side figure, labelled A) and B) -
#' the validation analogue of
#' `R/best_pipeline_search.R::plot_best_model_summary()`.
#'
#' Stability is always read via `type = "embedded"`: validation fits are
#' always run with `selection_params = NULL` (see this file's header), so
#' there is never a filter-selection step to diagnose via `type = "explicit"`
#' - only the model's own embedded selection (e.g. elastic net's per-fold
#' non-zero coefficients) is available.
#'
#' @param fit The validation `predict_cv()` fit.
#' @param title Overall figure title (e.g. "PREVAC Ad26/MVA -> EBOVAC2").
#' @param subtitle As built by `validation_subtitle()`.
#' @param selection_top_n `top_n` passed to `plot_selection_stability()`.
plot_validation_summary <- function(fit, title, subtitle, selection_top_n = 20) {

  p_fit <- tryCatch(plot(fit), error = function(e) {
    message("[validation] plot(fit) failed: ", conditionMessage(e))
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "CV prediction plot unavailable") +
      ggplot2::theme_void()
  })

  p_stability <- tryCatch(
    plot_selection_stability(fit, top_n = selection_top_n, type = "embedded", plot_type = "frequency"),
    error = function(e) {
      message("[validation] plot_selection_stability() failed: ", conditionMessage(e))
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "Selection stability plot unavailable") +
        ggplot2::theme_void()
    }
  )

  # As in plot_best_model_summary(): no guides = "collect", so the CV plot
  # keeps its own legend on its own right-hand side rather than it being
  # pulled to the bottom/side of the combined figure.
  (p_fit + p_stability) +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1)) +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      tag_levels = "A",
      tag_suffix = ")",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 13)
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 16, face = "bold")
    )
}
