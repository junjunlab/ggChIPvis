#' MultiHeatmap
#'
#' @author Jun Zhang
#'
#' @description
#' This function creates a multi-panel heatmap plot based on multiple matrices
#' using ChipHeatmap function.
#'
#' @param mat.list A list of matrices to be plotted which from
#' normalizeToMatrix/getTagMatrix/parseDeeptools.
#' @param row.split A vector specifying the groups to which each row belongs.
#' This vector should be of the same length as the number of rows in the matrices.
#' Default is NULL.
#' @param sample.label A vector of labels for the samples. This vector should
#'  be of the same length as mat.list. Default is NULL.
#' @param heatmap_rank_method A character string specifying the method to rank
#' the rows in the heatmap. The available options are "weighting", "sum",
#' "mean", and "median". Default is "weighting". Note "weighting" can not
#' be used for data which is from "getTagMatrix" function.
#' @param quantile.threshold A numeric value specifying the quantile threshold
#' for determining the range of the color scale. Default is 0.99.
#' @param add.panel.rect Logical value indicating whether to add a rectangle
#' around the entire heatmap plot. Default is FALSE.
#' @param panel.rect.gp An object of class "gpar" specifying graphical
#' parameters for the panel rectangle. Default is gpar().
#' @param keep.row.order Logical value indicating whether to keep the row order
#' of the first plot for the rest of the plots. Default is TRUE.
#' @param scale.range Logical value indicating whether to scale the color range
#' of each heatmap based on the quantile threshold. Default is FALSE.
#' @param heatmap.col A list of color schemes for each heatmap. If NULL, the
#' default color scheme is used. Default is NULL.
#' @param keep.one.left.annobar Logical value indicating whether to keep only
#' one left annotation panel for all heatmaps. Default is FALSE.
#' @param keep.one.legend Logical value indicating whether to keep only one
#' legend panel for all heatmaps. Default is FALSE.
#' @param ChipHeatmap.params A list of additional parameters to be passed to
#' the ChipHeatmap function. Default is list().
#' @param respect Logical value indicating whether to respect the layout
#' parameters specified by widths and heights when using multiple heatmaps.
#' Default is FALSE.
#' @param plot.size A numeric vector specifying the size of the overall plot.
#' The first element represents the width and the second element represents
#' the height. Default is c(0.9,0.9).
#' @param draw.profile Logical value indicating whether to draw a profile plot
#' below each heatmap. Default is FALSE.
#' @param keep.one.line.legend Logical value indicating whether to keep only
#' one line legend panel for all profile plots. Default is FALSE.
#' @param keep.one.profile.yaxis Logical value indicating whether to keep only
#' one y-axis label for all profile plots. Default is FALSE.
#' @param scale.y.range Logical value indicating whether to scale the y-axis
#' range of each profile plot based on the quantile threshold. Default is FALSE.
#' @param sample.group1 A vector specifying the grouping of samples for the
#' first annotation panel. Default is NULL.
#' @param sample.group2 A vector specifying the grouping of samples for the
#' second annotation panel. Default is NULL.
#' @param anno.panel.height Numeric value specifying the relative height of
#' the annotation panel. Default is 0.1.
#' @param draw.anno.fun.params.g1 A list of additional parameters to be passed
#' to the annotation function for the first group of samples.
#' Default is list(), details see draw.anno.fun().
#' @param draw.anno.fun.params.g2 A list of additional parameters to be passed
#' to the annotation function for the second group of samples.
#' Default is list(), details see draw.anno.fun().
#'
#' @return None
#'
#' @export
multiHeatmap <- function(mat.list = NULL,
                         row.split = NULL,
                         sample.label = NULL,
                         heatmap_rank_method = c("weighting","sum","mean","median"),
                         quantile.threshold = 0.99,
                         add.panel.rect = FALSE,
                         panel.rect.gp = gpar(),
                         keep.row.order = TRUE,
                         scale.range = FALSE,
                         heatmap.col = NULL,
                         keep.one.left.annobar = FALSE,
                         keep.one.legend = FALSE,
                         ChipHeatmap.params = list(),
                         respect = FALSE,
                         plot.size = c(0.9,0.9),
                         draw.profile = FALSE,
                         keep.one.line.legend = FALSE,
                         keep.one.profile.yaxis = FALSE,
                         scale.y.range = FALSE,
                         sample.group1 = NULL,
                         sample.group2 = NULL,
                         anno.panel.height = 0.1,
                         draw.anno.fun.params.g1 = list(),
                         draw.anno.fun.params.g2 = list()){
  heatmap_rank_method <- match.arg(heatmap_rank_method,c("weighting","sum","mean","median"))
  # ============================================================================
  # panels
  # ============================================================================
  ncol <- length(mat.list)

  # sample.group1 = NULL
  # sample.group2 = NULL
  if(is.null(sample.group1)){
    nrow <- 1
    heights = unit(1,"npc")
    ht.pos.row <- 1
  }else{
    if(is.null(sample.group2)){
      nrow <- 2
      heights = unit(c(anno.panel.height,1 - anno.panel.height),"npc")
      ht.pos.row <- 2
    }else{
      nrow <- 3
      heights = unit(c(anno.panel.height,anno.panel.height,1 - 2*anno.panel.height),"npc")
      ht.pos.row <- 3
    }
  }


  layout <- grid.layout(nrow = nrow,ncol = ncol,
                        widths = unit(rep(1/ncol,ncol),"npc"),
                        heights = heights,
                        respect = respect)
  # ============================================================================
  # basic settings
  # ============================================================================
  # whether multiple heatmap row orders be same with the first plot
  # keep.row.order = TRUE
  if(keep.row.order == TRUE){
    # get the first ranked row orders
    ods <- get.row.order(norm.mat = mat.list[[1]],
                         row.split = row.split,
                         heatmap_rank_method = heatmap_rank_method)

    row.order <- ods$genes
  }else{
    row.order <- NULL
  }

  # get the all matrces data range
  # scale.range = FALSE
  if(scale.range == TRUE){
    col.rg <- get.matlist.range(mat.list = mat.list,quantile.threshold = quantile.threshold)
  }else{
    col.rg <- NULL
  }

  # heatmap.col = NULL
  if(!is.null(heatmap.col)){
    if(!is.list(heatmap.col)){
      heatmap.col <- lapply(1:ncol, function(x){heatmap.col})
    }else{
      heatmap.col <- heatmap.col
    }
  }else{
    heatmap.col <- lapply(1:ncol, function(x){c("white","#379237")})
  }

  # keep.one.left.annobar = FALSE
  if(keep.one.left.annobar == TRUE){
    draw.anno.left <- c("TRUE",rep("FALSE",ncol - 1))
  }else{
    draw.anno.left <- rep("TRUE",ncol)
  }

  # keep.one.legend = FALSE
  if(scale.range == TRUE){
    if(keep.one.legend == TRUE){
      draw.legend <- c(rep("FALSE",ncol - 1),"TRUE")
      legend.pos <- rep("right",ncol)
    }
  }else{
    draw.legend <- rep("TRUE",ncol)
    legend.pos <- rep("bottom",ncol)
  }

  # whether draw profile
  if(draw.profile == TRUE){
    # keep.one.line.legend = FALSE
    if(keep.one.line.legend == TRUE){
      draw.line.legend <- c(rep("FALSE",ncol - 1),"TRUE")
    }else{
      draw.line.legend <- rep("TRUE",ncol)
    }

    # keep.one.profile.yaxis = FALSE
    if(keep.one.profile.yaxis == TRUE){
      profile.yaxis <- c("TRUE",rep("FALSE",ncol - 1))
    }else{
      profile.yaxis <- rep("TRUE",ncol)
    }

    # scale.y.range = FALSE
    if(scale.y.range == TRUE){
      profile.y.range <- get.matlist.range(mat.list = mat.list,
                                           type = "ave",
                                           quantile.threshold = quantile.threshold,
                                           row.split = row.split)
    }else{
      profile.y.range <- NULL
    }
  }else{
    draw.line.legend <- rep("TRUE",ncol)
    profile.yaxis <- rep("TRUE",ncol)
    profile.y.range <- NULL
  }

  # ============================================================================
  # panels
  # ============================================================================
  # plot.size = c(0.85,0.85)
  # respect = FALSE
  grid.newpage()
  pushViewport(viewport(x = 0.5,y = 0.5,
                        width = plot.size[1],height = plot.size[2],
                        layout = layout))

  # add.panel.rect = TRUE
  # panel.rect.gp = gpar()
  if(add.panel.rect == TRUE){
    grid.rect(gp = panel.rect.gp)
  }
  # ===========================================
  # draw.anno.fun.params.g1 = list()
  # draw.anno.fun.params.g2 = list()

  # draw anno
  if(nrow == 2){
    do.call(draw.anno.fun,modifyList(list(group.list = sample.group1,
                                          row.pos = 1),
                                     draw.anno.fun.params.g1))
  }else if(nrow == 3){
    do.call(draw.anno.fun,modifyList(list(group.list = sample.group1,
                                          row.pos = 1),
                                     draw.anno.fun.params.g1))

    do.call(draw.anno.fun,modifyList(list(group.list = sample.group2,
                                          row.pos = 2),
                                     draw.anno.fun.params.g2))
  }

  # ===========================================================
  # draw heatmap
  # ===========================================================
  for (j in 1:ncol) {
    # panel grid
    pushViewport(viewport(layout.pos.row = ht.pos.row,layout.pos.col = j))
    # grid.rect(gp = gpar(fill = "grey90"))

    # draw heatmap
    # ChipHeatmap.params = list()
    do.call(ChipHeatmap,
            modifyList(list(norm.mat = mat.list[[j]],
                            sample.label = sample.label[j],
                            heatmap_rank_method = heatmap_rank_method,
                            legend.pos = legend.pos[j],
                            row.split = row.split,
                            row.order = row.order,
                            col.range = col.rg,
                            ht.col = heatmap.col[[j]],
                            quantile.threshold = quantile.threshold,
                            draw.anno.left = draw.anno.left[j],
                            draw.legend = draw.legend[j],
                            draw.profile = draw.profile,
                            draw.line.legend = draw.line.legend[j],
                            profile.yaxis = profile.yaxis[j],
                            profile.y.range = profile.y.range),
                       ChipHeatmap.params))

    popViewport()
  }
}
