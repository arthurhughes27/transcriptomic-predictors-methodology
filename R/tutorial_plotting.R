# R/tutorial_plotting.R
#
# Tutorial-only figure styling: every figure in analysis/tutorial/tutorial.R
# is given a visual "this came directly from code" cue matching the
# chapter's minted code-block style (\setminted{bgcolor=codebg, frame=lines,
# ...}, with \definecolor{codebg}{HTML}{F0F2FF}) - the same F0F2FF
# background fill, plus black horizontal rules along the top and bottom
# edges only (minted's frame=lines draws top/bottom rules, not a full box -
# mirrored here rather than a 4-sided border).
#
# Deliberately kept separate from R/plotting.R::save_pipeline_comparison_plot()
# (used by every other figure in the repo), which is left unstyled - this
# framing is specific to the tutorial's figures.

#' Wrap a ggplot in the tutorial's code-block-style frame: an F0F2FF
#' background behind the whole plot (panel, legend, and margins alike) and a
#' thin black rule spanning the full width at the top and bottom edges.
#'
#' Implemented as a 3-row `gtable` (thin top border row / the plot itself /
#' thin bottom border row) rather than a `ggplot2::theme()` border, since
#' `element_rect()` only supports a uniform 4-sided border, not top/bottom
#' only.
#'
#' @param plot A `ggplot` object.
#' @param fill Background colour, matching `codebg` in the thesis's LaTeX
#'   preamble.
#' @param border_colour Colour of the top/bottom rules.
#' @param line_width `lwd` (in points) of the top/bottom rules.
#' @param frame_height `grid::unit()` height of each border row (just needs
#'   to comfortably contain `line_width`'s rule).
#'
#' @return A `gtable` object - pass this to `save_tutorial_plot()` (or
#'   `ggplot2::ggsave()` directly; `ggsave()` accepts any grid grob, not
#'   just `ggplot` objects).
add_tutorial_frame <- function(plot, fill = "#F0F2FF", border_colour = "black",
                                line_width = 1.2, frame_height = grid::unit(6, "pt")) {

  plot <- plot + ggplot2::theme(
    plot.background    = ggplot2::element_rect(fill = fill, colour = NA),
    panel.background   = ggplot2::element_rect(fill = fill, colour = NA),
    legend.background  = ggplot2::element_rect(fill = fill, colour = NA),
    legend.key         = ggplot2::element_rect(fill = fill, colour = NA)
  )

  # Filled rectangle (so the border row itself is also codebg-coloured, not
  # left blank/white) plus a centred horizontal rule spanning the full row
  # width, combined into one grob so both are added to the gtable together.
  border_row <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = fill, col = NA)),
    grid::linesGrob(
      x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(0.5, 0.5), "npc"),
      gp = grid::gpar(col = border_colour, lwd = line_width)
    )
  )

  gt <- gtable::gtable(
    widths  = grid::unit(1, "npc"),
    heights = grid::unit.c(frame_height, grid::unit(1, "null"), frame_height)
  )
  gt <- gtable::gtable_add_grob(gt, ggplot2::ggplotGrob(plot), t = 2, l = 1, name = "tutorial-plot")
  gt <- gtable::gtable_add_grob(gt, border_row, t = 1, l = 1, name = "tutorial-top-border")
  gt <- gtable::gtable_add_grob(gt, border_row, t = 3, l = 1, name = "tutorial-bottom-border")

  gt
}

#' Save a tutorial figure with `add_tutorial_frame()`'s styling applied,
#' creating `folder` if it does not already exist - the tutorial's
#' counterpart to `R/plotting.R::save_pipeline_comparison_plot()`.
save_tutorial_plot <- function(plot, folder, filename, width = 9, height = 4.5, dpi = 300) {
  fs::dir_create(folder)
  ggplot2::ggsave(
    filename = fs::path(folder, filename),
    plot = add_tutorial_frame(plot),
    width = width,
    height = height,
    dpi = dpi
  )
}
