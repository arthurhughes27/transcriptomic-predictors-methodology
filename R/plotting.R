# R/plotting.R
#
# Small helper to save a `compare_pipelines()` plot consistently across the
# pipeline-comparison scripts (output folder creation, figure size, DPI).

#' Save a pipeline-comparison plot to `folder/filename`, creating the folder
#' if it does not already exist.
save_pipeline_comparison_plot <- function(plot, folder, filename,
                                           width = 9, height = 4.5, dpi = 300) {
  fs::dir_create(folder)
  ggplot2::ggsave(
    filename = fs::path(folder, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}
