# R/best_pipeline_search.R
#
# A "smart", targeted search for the best-performing analytical pipeline
# per dataset, rather than an exhaustive grid over every
# engineering x selection x model combination.
#
# GENESET-ONLY: this search is one of "the main analyses" (geneset
# engineering + geneset-level selection + model choice), matching the
# separation used throughout the repo (analysis/pipeline_comparisons/'s
# 02/03 = main, geneset-only; analysis/supplementary/'s 01/02 = supplementary,
# gene-wise-only). Every engineering option offered here aggregates into gene
# sets, so round 2 only ever needs the geneset-level selection menu - there
# is no gene-wise counterpart to this search. Gene-wise engineering/selection
# are compared only in the analysis/supplementary/*/01_compare_engineering_genewise.R
# and 02_compare_selection_genewise.R scripts, never as part of finding "the
# best" pipeline.
#
# Implements a three-round GREEDY COORDINATE-ASCENT search through the
# pipeline stages, in their natural execution order (engineering ->
# selection -> model):
#   1. Compare engineering options (gene-set aggregation methods only)
#      against the usual mean-aggregation reference
#      (R/pipeline_defaults.R::reference_pipeline_params()) - pick the
#      single best (by R2).
#   2. Compare selection options (geneset-level methods only), with
#      engineering FIXED to round 1's winner (not the default reference) -
#      pick the single best.
#   3. Compare model options, with engineering AND selection FIXED to
#      rounds 1-2's winners - pick the single best.
#
# This costs 3 compare_pipelines() calls per dataset (one per round) instead
# of an exhaustive grid (~7 engineering x ~7 selection x ~5 model options
# per dataset). The trade-off: it's greedy, not exhaustive - each round
# conditions only on the PREVIOUS rounds' winners, so it can miss a
# combination where a locally-suboptimal earlier choice would have paired
# better with a later one (an interaction a one-round-at-a-time design
# can't see). That's the accepted cost of a targeted search over an
# exhaustive one.
#
# IMPORTANT SCOPE LIMITATION: this search runs entirely on the SINGLE
# (post-vaccination-only) dataset, for every dataset - never on the paired
# (baseline + post) dataset, even where one exists. This keeps every
# round's reference/options compared on an identical, consistent sample
# throughout the whole search (see
# analysis/pipeline_comparisons/visualize_comparisons.R's reference-context
# section for why that consistency matters - the same "don't mix samples"
# principle applies here). The consequence: RISE and dearseq's *paired*
# modes (which need individual_id/timepoint from the paired dataset) are
# never offered as selection candidates for SDY1276 (TIV), which has no
# placebo arm to power their *classic* mode instead - SDY1276's selection
# menus are therefore limited to variance/correlation/relative-gain.
# PREVAC's two folders do have a placebo arm, so RISE and dearseq's
# *classic* mode (contrasting treatment vs. placebo, needing only
# `treatment`, not pairing) are available there at both scales - pass
# `dearseq_mode = "classic"` to find_best_pipeline() to enable them.
#
# This module defines fresh option menus for all three search rounds
# (rather than re-reading option lists or saved results from
# 02_compare_engineering.R/03_compare_selection.R), so it has no dependency
# on those scripts' exact current option labels - engineering's search
# menu below intentionally mirrors Table 5.2's full option set, independent
# of whatever 02_compare_engineering.R happens to currently compare.

#' Elastic net, the default model used while searching the engineering and
#' selection axes (rounds 1-2), before model choice is itself searched
#' (round 3).
default_search_model_params <- function(inner_folds = 10, metric = "r2") {
  list(method = "glmnet", inner_folds = inner_folds, metric = metric, scale = TRUE,
       compute_importance = TRUE)
}

#' Engineering options compared in round 1: the five gene-set aggregation
#' methods (Table 5.2), same menu for every dataset. Mean aggregation is
#' deliberately excluded as an explicit option - it's the reference itself
#' (`reference_pipeline_params()`), already evaluated as the "Reference" row
#' by `compare_pipelines()`. Gene-level (no-aggregation) options are
#' intentionally NOT offered here - see this file's header: gene-wise
#' engineering is compared only in the supplementary
#' analysis/supplementary/*/01_compare_engineering_genewise.R scripts, never as a candidate for "the
#' best" pipeline.
engineering_search_menu <- function(genesets) {
  list(
    "Gene-set: median"    = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "median"),
    "Gene-set: max"       = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "max"),
    "Gene-set: 1st PC"    = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "pc1"),
    "Gene-set: GSVA"      = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "gsva", gsva_min_size = 2),
    "Gene-set: ssGSEA"    = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "ssgsea", ssgsea_min_size = 2)
  )
}

#' Selection options compared in round 2 (geneset-level only - see this
#' file's header).
#'
#' @param genesets Named list of gene sets, forwarded to the dearseq option.
#' @param dearseq_mode "classic" to include dearseq at the geneset level
#'   (requires `treatment`; PREVAC datasets), or NULL to omit it entirely
#'   (SDY1276 - see this file's header for why "paired" mode isn't offered
#'   in this search).
selection_geneset_search_menu <- function(genesets, dearseq_mode = NULL) {
  menu <- list(
    "Variance (top 25)" = list(method = "variance", top_n = 25),
    "Variance (top 100)" = list(method = "variance", top_n = 100),
    "Correlation - Spearman (top 25)"    = list(method = "spearman", top_n = 25),
    "Correlation - Spearman (|r| > 0.5)" = list(method = "spearman", threshold = 0.5),
    "Correlation - Pearson (top 25)"     = list(method = "pearson", top_n = 25),
    "Correlation - Pearson (|r| > 0.5)"  = list(method = "pearson", threshold = 0.5),
    "Univariate regression screening (threshold = 0)" = list(
      method = "relative_gain", threshold = 0,
      relative_gain_inner_folds = 10, relative_gain_metric = "rmse"
    )
  )
  if (!is.null(dearseq_mode)) {
    menu[["Dearseq (geneset, alpha = 0.05)"]] <- list(
      method = "dearseq", dearseq_mode = dearseq_mode, dearseq_level = "geneset",
      genesets = genesets, threshold = 0.05
    )
  }
  menu
}

#' Model options compared in round 3. Same menu for every dataset.
#'
#' `compute_importance = TRUE` throughout: whichever model wins, its saved
#' fit then supports `predictomics::plot_feature_importance()` - needed for
#' the combined best-model summary figure's B) panel
#' (`plot_best_model_summary()`/`R/panel_helpers.R`) whenever the winner has
#' no explicit filter step AND no embedded selection of its own (i.e.
#' anything other than "lasso"/"glmnet" here).
model_search_menu <- function(inner_folds = 10, metric = "r2") {
  list(
    "Linear regression"         = list(method = "lm", inner_folds = inner_folds, metric = metric, scale = TRUE, compute_importance = TRUE),
    "Lasso"                     = list(method = "lasso", inner_folds = inner_folds, metric = metric, scale = TRUE, compute_importance = TRUE),
    "Ridge"                     = list(method = "ridge", inner_folds = inner_folds, metric = metric, scale = TRUE, compute_importance = TRUE),
    "Random forest"             = list(method = "ranger", inner_folds = inner_folds, metric = metric, scale = TRUE, compute_importance = TRUE),
    "Support vector regression" = list(method = "svr", inner_folds = inner_folds, metric = metric, scale = TRUE, compute_importance = TRUE)
  )
}

#' Pick the best-performing row (by R2) among a compare_pipelines() result's
#' reference and option rows (excluding "baseline", which isn't a candidate
#' pipeline choice for any of the three search rounds).
#'
#' @return A one-row data frame with columns `pipeline`, `role`, `R2`.
.pick_round_winner <- function(res) {
  candidates <- res$results[res$results$role %in% c("reference", "option"), ]
  candidates[which.max(candidates$R2), c("pipeline", "role", "R2")]
}

#' Run the three-round greedy coordinate-ascent search described in this
#' file's header, for one dataset.
#'
#' @param X,Y,covariates,treatment As for `run_pipeline_comparison()`.
#'   Should be the SINGLE (post-vaccination-only) dataset's
#'   X/Y/covariates(/treatment) for the dataset being searched - never the
#'   paired dataset (see this file's header).
#' @param genesets Named list of gene sets (see `R/data_io.R::load_genesets()`).
#' @param dearseq_mode "classic" to make RISE and dearseq available as
#'   selection candidates (requires `treatment`; PREVAC datasets), or NULL
#'   to omit them entirely (SDY1276).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{engineering_params, selection_params, model_params}{The winning
#'       spec for each stage. `selection_params` may be `NULL`, meaning no
#'       selection step at all won round 2 (round 2's own reference is the
#'       "no selection" pipeline, since round 1's winning engineering always
#'       aggregates into a few hundred gene sets, not ~20,000 raw genes - no
#'       pre-filter fallback is needed here).}
#'     \item{winners}{A list of the three rounds' winner rows (from
#'       `.pick_round_winner()`), for reference.}
#'     \item{summary}{A data frame with one row per round: round label,
#'       winner, role ("reference" or "option"), R2.}
#'     \item{round_results}{The three raw `predictomics_comparison` objects
#'       (named `engineering`, `selection`, `model`), so the winning fit for
#'       any round can be recovered via
#'       `round_results[[stage]]$fits[[winners[[stage]]$pipeline]]`.}
#'   }
find_best_pipeline <- function(X, Y, covariates, treatment = NULL, genesets,
                                dearseq_mode = NULL) {

  model_params_default <- default_search_model_params()

  ## ---- Round 1: engineering ------------------------------------------------

  reference_params1 <- reference_pipeline_params(genesets)

  res1 <- run_pipeline_comparison(
    X = X, Y = Y, covariates = covariates, treatment = treatment,
    option_type = "engineering",
    option_choices = engineering_search_menu(genesets),
    reference_params = reference_params1
  )

  winner1 <- .pick_round_winner(res1)
  engineering_params <- if (winner1$role == "reference") {
    reference_params1$engineering_params
  } else {
    engineering_search_menu(genesets)[[winner1$pipeline]]
  }

  ## ---- Round 2: selection (geneset-level only), conditioned on round 1 -----
  #
  # Every round-1 engineering option (and the reference) aggregates into
  # gene sets - see this file's header - so round 2 only ever needs the
  # geneset-level selection menu. Gene-wise selection (RISE, gene-level
  # dearseq) is never a candidate here; it's compared only in the
  # supplementary analysis/supplementary/*/02_compare_selection_genewise.R scripts.

  reference_params2 <- list(
    engineering_params = engineering_params,
    selection_params   = NULL,
    model_params        = model_params_default
  )

  selection_menu <- selection_geneset_search_menu(genesets, dearseq_mode = dearseq_mode)

  res2 <- run_pipeline_comparison(
    X = X, Y = Y, covariates = covariates, treatment = treatment,
    option_type = "selection",
    option_choices = selection_menu,
    reference_params = reference_params2
  )

  winner2 <- .pick_round_winner(res2)
  selection_params <- if (winner2$role == "reference") {
    reference_params2$selection_params
  } else {
    selection_menu[[winner2$pipeline]]
  }

  ## ---- Round 3: model, conditioned on rounds 1-2's winners ------------------

  reference_params3 <- list(
    engineering_params = engineering_params,
    selection_params   = selection_params,
    model_params        = model_params_default
  )

  # diagnostics = "full": unlike rounds 1-2, round 3's winning fit is the
  # one actually saved and reused for downstream interpretation (best_fit -
  # see 05_find_best_model.R and plot_best_model_summary()), so it needs to
  # keep model_params$compute_importance = TRUE's per-feature scores rather
  # than have them stripped by the default "summary" diagnostics (see
  # R/run_comparison.R::run_pipeline_comparison()'s docs) - otherwise
  # predictomics::plot_feature_importance() has nothing to plot for a
  # winner with no explicit or embedded selection.
  res3 <- run_pipeline_comparison(
    X = X, Y = Y, covariates = covariates, treatment = treatment,
    option_type = "model",
    option_choices = model_search_menu(),
    reference_params = reference_params3,
    diagnostics = "full"
  )

  winner3 <- .pick_round_winner(res3)
  model_params <- if (winner3$role == "reference") {
    reference_params3$model_params
  } else {
    model_search_menu()[[winner3$pipeline]]
  }

  summary <- data.frame(
    round  = c("1. engineering", "2. selection", "3. model"),
    winner = c(winner1$pipeline, winner2$pipeline, winner3$pipeline),
    role   = c(winner1$role, winner2$role, winner3$role),
    R2     = c(winner1$R2, winner2$R2, winner3$R2),
    stringsAsFactors = FALSE
  )

  list(
    engineering_params = engineering_params,
    selection_params   = selection_params,
    model_params        = model_params,
    winners             = list(engineering = winner1, selection = winner2, model = winner3),
    summary             = summary,
    round_results        = list(engineering = res1, selection = res2, model = res3)
  )
}

#' Human-readable one-line description of a `find_best_pipeline()` result's
#' winning specification, e.g. "Engineering: Gene-set: 1st PC | Selection:
#' Dearseq (alpha = 0.05) | Model: Ridge". Used as the subtitle on the
#' combined best-model summary figure (`plot_best_model_summary()`) so the
#' figure is self-describing without cross-referencing `best$summary`.
#'
#' A round's winner is labelled with the *reference* option's human name
#' (via `R/metrics_labels.R::reference_option_label()`) whenever that round's
#' winner was the reference row rather than one of the searched options -
#' `find_best_pipeline()`'s own option-menu labels are otherwise used
#' directly, since they're already human-readable.
#'
#' @param best A list as returned by `find_best_pipeline()`.
#' @return A single character string.
describe_best_pipeline <- function(best) {

  engineering_label <- if (best$winners$engineering$role == "reference") {
    reference_option_label("engineering")
  } else {
    best$winners$engineering$pipeline
  }

  selection_label <- if (is.null(best$selection_params)) {
    "None"
  } else if (best$winners$selection$role == "reference") {
    reference_option_label("selection_geneset")
  } else {
    best$winners$selection$pipeline
  }

  model_label <- if (best$winners$model$role == "reference") {
    reference_option_label("model")
  } else {
    best$winners$model$pipeline
  }

  paste0(
    "Engineering: ", engineering_label,
    " | Selection: ", selection_label,
    " | Model: ", model_label
  )
}

#' Combine the winning pipeline's CV-prediction plot and selection-frequency
#' stability plot into a single side-by-side figure, labelled A) and B), with
#' a shared title and a subtitle detailing the winning specification.
#'
#' @param best_fit The winning pipeline's fitted `predictomics` object (a
#'   `predict_cv()`/`compare_pipelines()` fit result).
#' @param best The `find_best_pipeline()` result `best_fit` came from - used
#'   only to build the subtitle via `describe_best_pipeline()`.
#' @param title Overall figure title (e.g. dataset display name).
#' @param selection_top_n `top_n` passed to `plot_selection_stability()`.
#'
#' @return A `patchwork` object (two aligned panels tagged "A)"/"B)", plus
#'   the shared title/subtitle), or, if a panel fails to build, that panel is
#'   replaced with a blank placeholder rather than aborting the whole figure.
plot_best_model_summary <- function(best_fit, best, title, selection_top_n = 35) {

  p_fit <- tryCatch(plot(best_fit), error = function(e) {
    message("[best model] plot(best_fit) failed: ", conditionMessage(e))
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "CV prediction plot unavailable") +
      ggplot2::theme_void()
  })

  # See R/panel_helpers.R::build_selection_or_importance_panel() for the
  # explicit-selection / embedded-selection / feature-importance choice.
  p_stability <- build_selection_or_importance_panel(
    best_fit, best$selection_params, best$model_params, top_n = selection_top_n
  )

  # No guides = "collect" here: collecting would pool p_fit's legend with
  # p_stability's panel, pushing it away from the CV plot it belongs to and
  # over towards the stability panel. Keeping each panel's own guides means
  # p_fit's legend renders directly to the right of p_fit, not the combined
  # figure.
  (p_fit + p_stability) +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1)) +
    patchwork::plot_annotation(
      title = title,
      subtitle = describe_best_pipeline(best),
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
