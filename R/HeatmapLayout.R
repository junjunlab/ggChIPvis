#' Heatmap Layout Function
#'
#' @author Jun Zhang
#'
#' @description
#' This function creates a layout for a heatmap plot with using "gtabel" R
#' package. And this is used in "ChipHeatmap" function.
#'
#' @param mat A matrix for the heatmap data which from
#' normalizeToMatrix/getTagMatrix/parseDeeptools.
#' @param heatmap.size A numeric vector specifying the size of the heatmap in
#' relative units (default: c(0.7, 0.5)).
#' @param heatmap.xaxis.h The height of the heatmap x-axis labels in relative
#' units (default: 0.1).
#' @param heatmap.blegend.h The height of the heatmap bottom legend in relative
#' units (default: 0.035).
#' @param heatmap.rlegend.w The width of the heatmap right legend in relative
#' units (default: 0.035).
#' @param heatmap.profile.h The height of the heatmap profile area in relative
#' units (default: 0.2).
#' @param heatmap.label.h The height of the heatmap label area in relative
#' units (default: 0.075).
#' @param heatmap.lanno.w The width of the heatmap left annotation area in
#' relative units (default: 0.05).
#' @param heatmap.lanno.label.w The width of the heatmap left annotation labels
#' in relative units (default: 0.1).
#' @param draw.profile Logical indicating whether to draw the heatmap profile
#' area (default: TRUE).
#' @param legend.pos The position of the heatmap legend ("right" or "bottom",
#' default: "right").
#' @param plot.size A numeric vector specifying the size of the overall plot
#' (default: c(1, 1)).
#' @param newpage Logical indicating whether to start a new page for the plot
#' (default: TRUE).
#' @param draw.layout Logical indicating whether to draw the heatmap layout
#' (default: FALSE).
#'
#' @return A gtable object representing the heatmap layout.
#'
#' @import grid
#' @import gtable
#' @importFrom RColorBrewer brewer.pal
HeatmapLayout <- function(mat = NULL,
                          heatmap.size = c(0.7,0.5),
                          heatmap.xaxis.h = 0.1,
                          heatmap.blegend.h = 0.035,
                          heatmap.rlegend.w = 0.035,
                          heatmap.profile.h = 0.2,
                          heatmap.label.h = 0.075,
                          heatmap.lanno.w = 0.05,
                          heatmap.lanno.label.w = 0.1,
                          draw.profile = TRUE,
                          legend.pos = c("right","bottom"),
                          plot.size = c(1,1),
                          newpage = TRUE,
                          draw.layout = FALSE){
  legend.pos <- match.arg(legend.pos,c("right","bottom"))

  # ============================================================================
  # design heatmap
  # ============================================================================
  if(!is.null(mat)){
    xscale = c(1,ncol(mat))
    yscale = c(0,nrow(mat))

  }else{
    xscale = c(0,1)
    yscale = c(0,1)
  }

  heatmap.size <- heatmap.size
  heatmap.xaxis.h <- heatmap.xaxis.h

  if(legend.pos == "right"){
    heatmap.blegend.h <- 0
  }else{
    heatmap.rlegend.w <- 0
  }

  heatmap.profile.h <- ifelse(draw.profile == TRUE,heatmap.profile.h,0)

  heatmap.label.h <- heatmap.label.h
  heatmap.lanno.w <- heatmap.lanno.w
  heatmap.lanno.label.w <- heatmap.lanno.label.w

  # ============================================================================
  # layout
  # ============================================================================
  main.p <- gtable::gtable(widths = unit(c(heatmap.lanno.label.w,heatmap.lanno.w,
                                           heatmap.size[1],heatmap.rlegend.w),"npc"),
                           heights = unit(c(heatmap.label.h,heatmap.profile.h,
                                            heatmap.size[2],
                                            heatmap.xaxis.h,heatmap.blegend.h),"npc"),
                           name = "main",
                           vp = viewport(width = plot.size[1],height = plot.size[2]))


  # gtable_show_layout(main.p)

  return(main.p)
}
