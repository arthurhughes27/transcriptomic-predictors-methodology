# R/text_output.R
#
# Small helper to save captured console output (str()/print() text, table
# previews) consistently, analogous to R/plotting.R::save_pipeline_comparison_plot()
# for figures. Used by analysis/tutorial/tutorial.R to produce plain-text
# files ready to paste into the thesis's minted{text} blocks.

#' Save captured console output to `folder/filename`, creating the folder if
#' it does not already exist.
#'
#' @param text A character vector of output lines, typically from
#'   `capture.output()`.
#' @param folder Directory to save into.
#' @param filename File name (including extension, e.g. `"foo.txt"`).
save_console_output <- function(text, folder, filename) {
  fs::dir_create(folder)
  writeLines(text, fs::path(folder, filename))
}
