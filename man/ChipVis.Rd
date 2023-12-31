% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChipVis.R
\name{ChipVis}
\alias{ChipVis}
\title{ChipVis Function}
\usage{
ChipVis(
  object = NULL,
  plot.type = c("both", "profile", "heatmap"),
  sample.order = NULL,
  group.order = NULL,
  group2.order = NULL,
  rowsplit.order = NULL,
  rowgroup.order = NULL,
  facet_profile_params = list(),
  facet_heatmap_params = list(),
  theme_params = list(),
  vline_params = list(),
  axisx_params = list(),
  geom_line_params = list(),
  rel_height = c(0.3, 0.7),
  multi_legend_barheight = 1,
  multi.heatmap.col = NULL,
  geom_rect_params = list(),
  draw_rect = TRUE,
  rect_ratio = 0.1,
  merge_facet = FALSE,
  add_CI95_layer = TRUE,
  geom_ribbon_params = list(),
  facet_profile = TRUE,
  add_nested_line = FALSE
)
}
\arguments{
\item{object}{An object containing the ChIP-seq data which from retriveData function.
Must have "profile" and "heatmap" components.}

\item{plot.type}{Type of plot to generate. Can be "both" (default), "profile", or "heatmap".}

\item{sample.order}{Order of the samples for plotting. Default is NULL.}

\item{group.order}{Order of the groups for plotting. Default is NULL.}

\item{group2.order}{Order of the secondary groups for plotting. Default is NULL.}

\item{rowsplit.order}{Order of the heatmap row split. Default is NULL.}

\item{rowgroup.order}{Order of the heatmap row split groups. Default is NULL.}

\item{facet_profile_params}{A list parameters for customizing the profile facet layer.
Details see ggh4x::facet_nested().}

\item{facet_heatmap_params}{A list parameters for customizing the heatmap facet layer.
Details see ggh4x::facet_nested().}

\item{theme_params}{A list parameters for customizing the theme of the plot.
Details see ggplot2::theme().}

\item{vline_params}{A list parameters for customizing the vertical line layer.
Details see ggplot2::geom_vline().}

\item{axisx_params}{A list parameters for customizing the x-axis scale.
Details see ggplot2::scale_x_continues().}

\item{geom_line_params}{A list parameters for customizing the line layer.
Details see ggplot2::geom_line().}

\item{rel_height}{Relative heights of the profile and heatmap plots when plot.type
is "both". Default is c(0.3,0.7).}

\item{multi_legend_barheight}{Height of the colorbar legend for multiple
colors in the heatmap. Default is 1.}

\item{multi.heatmap.col}{Colors for the multiple colors in the heatmap. Default is NULL.}

\item{geom_rect_params}{A list parameters for customizing the rectangle layer.
Details see ggplot2::geom_rect().}

\item{draw_rect}{Whether to draw a rectangle layer. Default is TRUE.}

\item{rect_ratio}{Ratio of the rectangle width to the maximum density range.
Default is 0.1.}

\item{merge_facet}{Whether to merge the facet layers for different binding sites.
Default is FALSE.}

\item{add_CI95_layer}{Whether to add a confidence interval (CI) 95\% layer to the
profile plot. Default is TRUE.}

\item{geom_ribbon_params}{A list parameters for customizing the CI 95\% layer.
Details see ggplot2::geom_ribbon().}

\item{facet_profile}{Whether to include the profile facet layer. Default is TRUE.}

\item{add_nested_line}{It is used when plot.type = "both" and you want to add
nested line for the upper profile plot. Default is FALSE.}
}
\value{
A ggplot object representing the ChIP-seq data visualization.
}
\description{
Visualize ChIP-seq data using both profile and heatmap plots. This function use
output data from "retriveData" function and use ggplot2 package to produce profile
and heatmap plots. Most parameters can be accepted in a list form with ggplot
arguments.
}
\author{
Jun Zhang
}
