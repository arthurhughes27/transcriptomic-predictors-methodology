# Tutorial: `predictomics` package demonstration

A self-contained demonstration of `predictomics`'s two main entry points,
`predict_cv()` and `compare_pipelines()`, on simulated data - source
material for the thesis's package-tutorial section, in the same style as
the `SurrogateRank` tutorial (Section "The `SurrogateRank` R package").

Unlike every other folder under `analysis/`, this one does not touch this
repository's real, privacy-protected data or the pipeline-comparison
infrastructure in `R/` (`pipeline_defaults.R`, `run_comparison.R`,
`metrics_io.R`) - it only illustrates package usage on simulated data, so
it is fully self-contained and runnable by anyone with just `predictomics`
installed. It is not part of the chapter's substantive analysis and is
deliberately not sourced from `analysis/master_analysis.R`.

## Running it

From the repository root:

```r
source(fs::path("analysis", "tutorial", "tutorial.R"))
```

This writes two kinds of output (neither committed - see `.gitignore`,
which excludes `output/` entirely):

- `output/figures/tutorial/*.pdf` - one PDF per plot produced during the
  tutorial, ready for `\includegraphics{}`.
- `output/text/tutorial/*.txt` - one plain-text file per console output
  (`str()`, `print()`, table previews), ready to paste into a
  `\begin{minted}{text} ... \end{minted}` block, or to bring in directly via
  `\verbatiminput{}`/`\lstinputlisting{}` if you'd rather not paste by hand.

## Mapping outputs onto the chapter

The script is split into two parts, each producing the same three kinds of
material the `SurrogateRank` tutorial uses - narrative code, a console-text
table/summary, and a figure:

| # | What | File | Suggested LaTeX use |
|---|------|------|----------------------|
| 1.1 | `str(sim)` | `output/text/tutorial/01_simulate_str.txt` | `\begin{minted}{text}` block, as in the `SurrogateRank` `str(full_data)` example |
| 1.2 | `predict_cv()` call itself | (paste the R code block directly) | `\begin{minted}{r}` block |
| 1.3 | `print(fit)` | `output/text/tutorial/02_predict_cv_print.txt` | `\begin{minted}{text}` block, or reformat the key lines (RMSE/sRMSE/R2/SpearmanR) into a table |
| 1.3 | Observed-vs-predicted scatter | `output/figures/tutorial/tutorial_predict_cv_scatter.pdf` | `\includegraphics` |
| 1.3 | Embedded selection stability (lasso) | `output/figures/tutorial/tutorial_selection_stability.pdf` | `\includegraphics` |
| 1.3 | Feature importance | `output/figures/tutorial/tutorial_feature_importance.pdf` | `\includegraphics` |
| 1.3 | Ground-truth signal genesets | `output/text/tutorial/03_signal_genesets.txt` | inline text, or a footnote confirming the pipeline recovers the simulated signal |
| 2.1 | `compare_pipelines()` call itself | (paste the R code block directly) | `\begin{minted}{r}` block |
| 2.2 | `print(cmp)` | `output/text/tutorial/04_compare_pipelines_print.txt` | `\begin{minted}{text}` block |
| 2.2 | Full results table | `output/text/tutorial/05_compare_pipelines_results.txt` | source for a hand-built `\begin{table}`, or `\begin{minted}{text}` directly |
| 2.2 | Comparison bar chart | `output/figures/tutorial/tutorial_compare_pipelines.pdf` | `\includegraphics` |

## Reusable helpers

- `R/plotting.R::save_pipeline_comparison_plot()` - same figure-saving
  helper used throughout `analysis/pipeline_comparisons/`, reused here
  rather than duplicated.
- `R/text_output.R::save_console_output()` - the same pattern for captured
  console text, added alongside this tutorial since no equivalent existed
  yet (the rest of the repository saves metrics tables and figures, not
  console-text snippets).

## Editing/extending

`tutorial.R` is deliberately linear and narrated with section comments so
that blocks of it can be copied into `\begin{minted}{r}` environments with
minimal editing - the code as written is what should appear in the chapter,
not a wrapper around it. If you want to demonstrate a different pipeline
axis (e.g. `option_type = "engineering"` or `"selection"` instead of
`"model"` in Part 2, or a paired/`gene_level_fc` design via
`simulate_predictomics_data(design = "paired")`), duplicate Part 2's
structure with the new `option_choices` rather than editing it in place, so
the existing figures/outputs stay reproducible from the committed script.
