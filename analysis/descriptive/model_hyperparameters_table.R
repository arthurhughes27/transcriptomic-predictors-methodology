# analysis/descriptive/model_hyperparameters_table.R
#
# Produce a manuscript-ready LaTeX table summarising the hyperparameters and
# tuning grids for each of the six regression methods used in this chapter's
# model-choice comparisons (see R/best_pipeline_search.R::model_search_menu(),
# analysis/pipeline_comparisons/*/04_compare_model.R) and reference pipeline
# (R/pipeline_defaults.R::reference_pipeline_params()/raw_gene_reference_params()).
#
# Every grid below is verified directly against source, not recalled from
# memory of caret/glmnet defaults:
#   - predictomics::run_model() (github.com/arthurhughes27/predictomics,
#     R/run_model.R) - how model_params is turned into a caret::train() call,
#     including which methods get an EXPLICIT tune_grid built by
#     predictomics itself (lasso, ridge, svr) vs. which fall through to
#     caret's own default grid function (glmnet with alpha tuned, ranger).
#   - caret's own per-method grid functions (github.com/topepo/caret,
#     inst/model_sources/files/{glmnet,ranger,svmLinear}.R) and
#     caret::var_seq() (R/misc.R), for the methods predictomics does NOT
#     override.
#
# Two things worth flagging, found only by reading run_model() directly
# rather than assuming standard caret/glmnet defaults:
#   1. model_params$metric (e.g. "r2", used elsewhere in this repo for
#      OUTER performance reporting) is never read by run_model() - inner-CV
#      hyperparameter selection always minimises RMSE
#      (caret::train(..., metric = "RMSE") is hardcoded for every tuned
#      method). This table's footnote makes that explicit rather than
#      implying the grid is searched by R2.
#   2. glmnet's default alpha/lambda grid (used only for "glmnet"/elastic
#      net, since "lasso"/"ridge" get predictomics' own explicit grid
#      instead) is PARTLY data-dependent: alpha is fixed
#      (caret::getModelInfo("glmnet")$grid: seq(0.1, 1, length = tuneLength)),
#      but lambda is refit from the training fold's own glmnet
#      regularisation path each time, not a fixed a priori sequence. Ranger's
#      default mtry grid (caret::var_seq()) is similarly a function of p
#      (the number of predictors reaching the model for a given pipeline),
#      not a fixed sequence. Both are reported here as the exact PROCEDURE
#      rather than literal numbers, since the actual numbers vary by
#      pipeline/fold.
#
# This script's content is static (a property of the code, not of any
# dataset), so - unlike dataset_characteristics_table.R - it does not read
# any data/ files; it can be run standalone at any time.

library(fs)

table_data <- data.frame(
  method = c(
    "Linear regression (OLS)",
    "Elastic net",
    "Lasso",
    "Ridge",
    "Random forest",
    "Support vector regression (linear kernel)"
  ),
  hyperparameters = c(
    "None",
    "$\\alpha$ (mixing parameter), $\\lambda$ (regularisation strength)",
    "$\\lambda$ ($\\alpha = 1$ fixed)",
    "$\\lambda$ ($\\alpha = 0$ fixed)",
    "\\texttt{mtry}, \\texttt{splitrule} (\\texttt{min.node.size} fixed, not tuned)",
    "$C$ (cost)"
  ),
  grid = c(
    "--- (single OLS fit; no inner-CV tuning loop)",
    "$\\alpha \\in \\{0.10,\\, 0.55,\\, 1.00\\}$ (\\texttt{caret} default grid, tune length 3); for each $\\alpha$, 3 $\\lambda$ values taken from the \\texttt{glmnet} regularisation path fit on that training fold (\\texttt{nlambda} = 5, endpoints dropped) - data-dependent, not a fixed sequence",
    "50 values, log$_{10}$-spaced, $\\lambda \\in [10^{-3}, 10^{1}]$ (i.e.\\ 0.001 to 10)",
    "50 values, log$_{10}$-spaced, $\\lambda \\in [10^{-3}, 10^{1}]$ (i.e.\\ 0.001 to 10) - same grid as Lasso",
    "\\texttt{mtry}: 3 values via \\texttt{caret::var\\_seq()} - evenly-spaced integers between 2 and $p$ for $p < 500$, log$_2$-spaced for $p \\geq 500$ ($p$ = number of predictors reaching the model for that pipeline); \\texttt{splitrule} $\\in$ \\{variance, extratrees\\}; \\texttt{min.node.size} fixed at 5",
    "20 values, log$_{10}$-spaced, $C \\in [10^{-2}, 10^{2}]$ (i.e.\\ 0.01 to 100)"
  ),
  stringsAsFactors = FALSE
)

cat("Model hyperparameters and tuning grids:\n")
print(table_data)

# --- Render as a manuscript-ready LaTeX table --------------------------------
#
# Hand-built (as in dataset_characteristics_table.R) to keep this script
# dependency-free. Unlike that table, the hyperparameters/grid columns
# contain deliberate LaTeX markup (math mode, \texttt{}), so they are NOT
# passed through any escaping step here - only the plain-text Method column
# would need escaping, and none of its values contain LaTeX special
# characters.

format_data_row <- function(row) {
  paste0(
    paste(c(row$method, row$hyperparameters, row$grid), collapse = " & "),
    " \\\\"
  )
}

body_lines <- vapply(seq_len(nrow(table_data)), function(i) format_data_row(table_data[i, ]), character(1))

latex_lines <- c(
  "% Requires \\usepackage{booktabs} in the preamble.",
  "\\begin{table}[t]",
  "  \\centering",
  "  \\caption{Hyperparameters and tuning grids for each regression method.}",
  "  \\label{tab:model-hyperparameters}",
  "  \\resizebox{\\textwidth}{!}{%",
  "  \\begin{tabular}{lll}",
  "    \\toprule",
  "    Method & Hyperparameters & Hyperparameter grid \\\\",
  "    \\midrule",
  paste0("    ", body_lines),
  "    \\bottomrule",
  "  \\end{tabular}%",
  "  }",
  "",
  "  \\vspace{4pt}",
  "  \\footnotesize",
  "  \\raggedright",
  "  All tuned methods select hyperparameters via 10-fold inner cross-validation, minimising RMSE (the R\\textsuperscript{2}-based metric recorded elsewhere in this chapter is used only for outer, between-pipeline comparison, not for inner-loop hyperparameter selection). Every predictor matrix is centred and scaled immediately before fitting, uniformly across methods.",
  "\\end{table}"
)

table_path <- fs::path("output", "tables")
fs::dir_create(table_path)

writeLines(latex_lines, fs::path(table_path, "model_hyperparameters.tex"))

cat("\nLaTeX table written to ", fs::path(table_path, "model_hyperparameters.tex"), "\n", sep = "")

rm(list = ls())
