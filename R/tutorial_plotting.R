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
#' background filling the entire figure, and a thin black rule spanning the
#' full width exactly at the top and bottom edges of that same background.
#'
#' Implemented as a single-cell `gtable` stack, drawn bottom-to-top: (1) an
#' opaque background rectangle filling the whole cell, (2) the plot itself,
#' with its own `plot.background`/`panel.background` made fully transparent
#' so (1) shows through everywhere the plot doesn't otherwise paint, and (3)
#' the two border rules, positioned at the cell's own top/bottom edges (not
#' a separate row) so they sit exactly at the background's own edges.
#'
#' This two-part design (transparent plot over an explicit background,
#' rather than colouring `plot.background`/`panel.background` directly and
#' relying on them to fill the whole figure) matters for a plot using
#' `coord_fixed()`/`coord_equal()` (e.g. a predicted-vs-observed scatter):
#' such a plot's *panel* is forced to a square smaller than the figure's
#' full extent, so a fill on `panel.background` alone would only colour
#' that square, not the figure - drawing the background as a separate,
#' full-cell rectangle underneath sidesteps that entirely, regardless of
#' what aspect ratio the plot itself uses.
#'
#' @param plot A `ggplot` object.
#' @param fill Background colour, matching `codebg` in the thesis's LaTeX
#'   preamble.
#' @param border_colour Colour of the top/bottom rules.
#' @param line_width `lwd` (in points) of the top/bottom rules.
#'
#' @return A `gtable` object - pass this to `save_tutorial_plot()` (or
#'   `ggplot2::ggsave()` directly; `ggsave()` accepts any grid grob, not
#'   just `ggplot` objects).
add_tutorial_frame <- function(plot, fill = "#F0F2FF", border_colour = "black",
                                line_width = 1.2) {

  # Transparent, not codebg-filled: the background is painted separately,
  # underneath, by the explicit rectGrob below - see this function's docs
  # for why (coord_fixed()/coord_equal() plots).
  plot <- plot + ggplot2::theme(
    plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
    legend.background  = ggplot2::element_rect(fill = NA, colour = NA),
    legend.key         = ggplot2::element_rect(fill = NA, colour = NA)
  )

  background <- grid::rectGrob(gp = grid::gpar(fill = fill, col = NA))

  # Inset by half the stroke width so the full rule is drawn inside the
  # figure (a line centred exactly on the edge would have half its width
  # clipped by the device boundary) while still reading as flush with the
  # background's own edge.
  top_line <- grid::linesGrob(
    x = grid::unit(c(0, 1), "npc"),
    y = grid::unit(1, "npc") - grid::unit(line_width / 2, "pt"),
    gp = grid::gpar(col = border_colour, lwd = line_width)
  )
  bottom_line <- grid::linesGrob(
    x = grid::unit(c(0, 1), "npc"),
    y = grid::unit(0, "npc") + grid::unit(line_width / 2, "pt"),
    gp = grid::gpar(col = border_colour, lwd = line_width)
  )

  gt <- gtable::gtable(
    widths  = grid::unit(1, "npc"),
    heights = grid::unit(1, "npc")
  )
  gt <- gtable::gtable_add_grob(gt, background, t = 1, l = 1, z = 1, name = "tutorial-background")
  gt <- gtable::gtable_add_grob(gt, ggplot2::ggplotGrob(plot), t = 1, l = 1, z = 2, name = "tutorial-plot")
  gt <- gtable::gtable_add_grob(gt, top_line, t = 1, l = 1, z = 3, name = "tutorial-top-border")
  gt <- gtable::gtable_add_grob(gt, bottom_line, t = 1, l = 1, z = 3, name = "tutorial-bottom-border")

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
