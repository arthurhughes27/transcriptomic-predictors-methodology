# R/panel_helpers.R
#
# Shared logic for the B) panel of the combined best-model/validation-fit
# summary figures (R/best_pipeline_search.R::plot_best_model_summary(),
# R/external_validation.R::plot_validation_summary()): whichever of
# explicit selection-frequency, embedded selection-frequency, or (when
# NEITHER was performed) predictomics::plot_feature_importance() is
# actually meaningful for the fit at hand.

#' Build the "which feature panel makes sense here" B) panel for a
#' `predict_cv()` fit.
#'
#' Three cases, in order of preference:
#'   1. An explicit filter selection step ran (`selection_params` not
#'      `NULL`, e.g. variance/correlation/RISE/dearseq) -
#'      `plot_selection_stability(type = "explicit")`.
#'   2. No filter step, but the model itself performs embedded selection
#'      (`model_params$method \%in\% c("glmnet", "lasso")`) -
#'      `plot_selection_stability(type = "embedded")`.
#'   3. Neither: no filter step AND the model has no embedded selection of
#'      its own (e.g. "lm", "ridge", "ranger", "svr") - there is no
#'      selection frequency to show at all, so
#'      `predictomics::plot_feature_importance()` is shown instead. This
#'      needs the fit to have been produced with
#'      `model_params$compute_importance = TRUE` (see
#'      R/best_pipeline_search.R::model_search_menu()/
#'      default_search_model_params(), where this is always set) - if it
#'      wasn't (e.g. a fit cached from before that flag existed), the
#'      resulting error is caught and a placeholder shown instead, same as
#'      any other panel-build failure here.
#'
#' @param fit A `predict_cv()`/`compare_pipelines()` fit result.
#' @param selection_params The `selection_params` used to produce `fit`
#'   (`NULL` if none - e.g. always `NULL` for external-validation fits,
#'   which fix their feature panel upfront rather than re-selecting; see
#'   R/external_validation.R's header).
#' @param model_params The `model_params` used to produce `fit`.
#' @param top_n `top_n` passed through to whichever plotting function is
#'   used.
#'
#' @return A `ggplot` object (or a blank placeholder annotated with an
#'   explanatory message, if the appropriate plot could not be built).
build_selection_or_importance_panel <- function(fit, selection_params, model_params, top_n = 20) {

  placeholder <- function(label) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = label) +
      ggplot2::theme_void()
  }

  has_explicit_selection <- !is.null(selection_params)
  has_embedded_selection <- isTRUE(model_params$method %in% c("glmnet", "lasso"))

  if (!has_explicit_selection && !has_embedded_selection) {
    return(tryCatch(
      predictomics::plot_feature_importance(fit, top_n = top_n),
      error = function(e) {
        message("[panel] plot_feature_importance() failed: ", conditionMessage(e))
        placeholder("Feature importance plot unavailable.")
      }
    ))
  }

  stability_type <- if (has_explicit_selection) "explicit" else "embedded"

  tryCatch(
    plot_selection_stability(fit, top_n = top_n, type = stability_type, plot_type = "frequency"),
    error = function(e) {
      message("[panel] plot_selection_stability(type = '", stability_type,
              "') failed: ", conditionMessage(e))
      # As before: an explicit filter step's per-fold diagnostics can be
      # unavailable for reasons unrelated to whether one was configured
      # (e.g. a method that doesn't report them) - fall back to embedded
      # diagnostics rather than giving up on the panel entirely, since the
      # model may still perform its own per-fold selection regardless.
      if (stability_type == "explicit") {
        tryCatch(
          plot_selection_stability(fit, top_n = top_n, type = "embedded", plot_type = "frequency"),
          error = function(e2) {
            message("[panel] plot_selection_stability(type = 'embedded') fallback also failed: ",
                    conditionMessage(e2))
            placeholder("No feature selection performed.")
          }
        )
      } else {
        placeholder("No feature selection performed.")
      }
    }
  )
}
