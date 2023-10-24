#' grid.yaxis2 function
#' @author Jun Zhang
#'
#' This function creates a y-axis in a graphical grid.
#'
#' @param at numeric vector specifying the locations of tick marks on the axis.
#' @param breaks integer specifying the number of breaks on the axis.
#' @param labels character vector providing custom labels for the tick marks.
#' @param tick.len numeric specifying the length of the tick marks.
#' @param label.space numeric specifying the space between the tick marks and labels.
#' @param side character specifying the side to place the axis ("left" or "right").
#' @param ticks.gp gpar() seetings for the axis ticks.
#' @param text.gp gpar() seetings for the axis label.
#' @param hjust label hjust.
#' @param vjust label vjust.
#' @param rot numeric specifying the rotation angle for the tick mark labels.
#' @param draw if draw == "FALSE" , then return a gTree object.
#' @param vp viewport seetings.
#'
#' @import grid
#'
#' @export
grid.yaxis2 <- function(at = NULL,
                        breaks = 5,
                        labels = NULL,
                        tick.len = 0.5,
                        label.space = 0.25,
                        side = c("left","right"),
                        rot = 0,
                        ticks.gp = NULL,
                        text.gp = NULL,
                        hjust = NULL,
                        vjust = NULL,
                        draw = TRUE,
                        vp = NULL){
  # labels and ticks
  if(is.null(at) || is.null(labels)){
    if(is.null(vp)){
      at <- grid.pretty(current.viewport()$yscale,n = breaks)
    }else{
      at <- grid.pretty(vp$yscale,n = breaks)
    }
    labels <- as.character(at)
  }else{
    at <- at
    labels <- as.character(labels)
  }

  # axis position
  side <- match.arg(side,c("left","right"))
  if(side == "left"){
    tck.x0 = unit(0, "npc")
    tck.x1 = unit(-tick.len, "lines")
    text.x = unit(-tick.len - label.space,"lines")
    text.just = "right"
  }else{
    tck.x0 = unit(1, "npc")
    tck.x1 = unit(1, "npc") + unit(tick.len, "lines")
    text.x = unit(abs(tick.len) + abs(label.space),"lines") + unit(1, "npc")
    text.just = "left"
  }

  if(draw == TRUE){
    grid.segments(x0 = 0,x1 = 0,y0 = 0,y1 = 1,vp = vp)
    grid.segments(y0 = unit(at, "native"),y1 = unit(at, "native"),
                  x0 = tck.x0,x1 = tck.x1,
                  gp = ticks.gp,vp = vp)

    grid.text(label = labels,
              y = unit(at, "native"),x = text.x,
              rot = rot,
              just = text.just,hjust = hjust,vjust = vjust,
              gp = text.gp,vp = vp)
  }else{
    seg.main <- segmentsGrob(x0 = 0,x1 = 0,y0 = 0,y1 = 1)

    ticks <- segmentsGrob(y0 = unit(at, "native"),y1 = unit(at, "native"),
                          x0 = tck.x0,x1 = tck.x1,
                          gp = ticks.gp)

    text <- textGrob(label = labels,
                     y = unit(at, "native"),x = text.x,
                     rot = rot,
                     just = text.just,hjust = hjust,vjust = vjust,
                     gp = text.gp)

    comb <- gTree(children = gList(seg.main,ticks,text),vp = vp)
    # comb <- gList(seg.main,ticks,text)
    return(comb)
  }
}



#' This function creates a x-axis in a graphical grid.
#'
#' @param at numeric vector specifying the locations of tick marks on the axis.
#' @param breaks integer specifying the number of breaks on the axis.
#' @param labels character vector providing custom labels for the tick marks.
#' @param tick.len numeric specifying the length of the tick marks.
#' @param label.space numeric specifying the space between the tick marks and labels.
#' @param side character specifying the side to place the axis ("bottom" or "top").
#' @param ticks.gp gpar() seetings for the axis ticks.
#' @param text.gp gpar() seetings for the axis label.
#' @param hjust label hjust.
#' @param vjust label vjust.
#' @param rot numeric specifying the rotation angle for the tick mark labels.
#' @param draw if draw == "FALSE" , then return a gTree object.
#' @param vp viewport seetings.
#'
#' @export
grid.xaxis2 <- function(at = NULL,
                        breaks = 5,
                        labels = NULL,
                        tick.len = 0.5,
                        label.space = 0.5,
                        side = c("bottom","top"),
                        rot = 0,
                        ticks.gp = NULL,
                        text.gp = NULL,
                        hjust = NULL,
                        vjust = NULL,
                        draw = TRUE,
                        vp = NULL){
  # labels and ticks
  if(is.null(at) || is.null(labels)){
    if(is.null(vp)){
      at <- grid.pretty(current.viewport()$yscale,n = breaks)
    }else{
      at <- grid.pretty(vp$xscale,n = breaks)
    }
    labels <- as.character(at)
  }else{
    at <- at
    labels <- as.character(labels)
  }

  # axis position
  side <- match.arg(side,c("bottom","top"))
  if(side == "bottom"){
    tck.y0 = unit(0, "npc")
    tck.y1 = unit(-tick.len, "lines")
    text.y = unit(-tick.len - label.space,"lines")
  }else{
    tck.y0 = unit(1, "npc")
    tck.y1 = unit(1, "npc") + unit(tick.len, "lines")
    text.y = unit(abs(tick.len) + abs(label.space),"lines") + unit(1, "npc")
  }

  if(draw == TRUE){
    grid.segments(x0 = 0,x1 = 1,y0 = 0,y1 = 0,vp = vp)
    grid.segments(x0 = unit(at, "native"),x1 = unit(at, "native"),
                  y0 = tck.y0,y1 = tck.y1,
                  gp = ticks.gp,vp = vp)

    grid.text(label = labels,
              x = unit(at, "native"),y = text.y,
              rot = rot,
              hjust = hjust,vjust = vjust,
              gp = text.gp,vp = vp)
  }else{
    seg.main <- segmentsGrob(x0 = 0,x1 = 1,y0 = 0,y1 = 0)

    ticks <- segmentsGrob(x0 = unit(at, "native"),x1 = unit(at, "native"),
                          y0 = tck.y0,y1 = tck.y1,
                          gp = ticks.gp)

    text <- textGrob(label = labels,
                     x = unit(at, "native"),y = text.y,
                     rot = rot,
                     hjust = hjust,vjust = vjust,
                     gp = text.gp)

    comb <- gTree(children = gList(seg.main,ticks,text),vp = vp)
    # comb <- gList(seg.main,ticks,text)
    return(comb)
  }

}




#' Create a Color Key Grid
#'
#' This function creates a color key grid using the grid graphics system.
#'
#' @param x Numeric vector specifying the range for the x-axis.
#' @param color Vector of color values to be used for the color key.
#' @param color.n Number of colors to generate.
#' @param ticks.side Side of the color key to display tick marks
#' ("left", "right", "top", or "bottom").
#' @param nbreaks the colorbar break numbers, default 5.
#' @param ... other args passed by grid.xaxis2()/grid.yaxis2().
#' @param pos Position of the color key ("h" for horizontal, "v" for vertical).
#'
#' @return NULL
#'
#' @importFrom grDevices colorRampPalette
#'
#' @export
grid.colorkey <- function(x = NULL,
                          color = NULL,
                          color.n = 100,
                          ticks.side = c("left","right","top","bottom"),
                          pos = c("h","v"),
                          nbreaks = 5,...){
  pos <- match.arg(pos,c("h","v"))
  ticks.side <- match.arg(ticks.side,c("left","right","top","bottom"))

  # check position
  if(pos == "v"){
    x_scale <- c(0,1)
    y_scale <- range(as.numeric(x))

    xpos <- 0.5
    ypos <- seq(0,1, length = color.n)

    r_width = unit(1, "npc")
    r_height = 1/(color.n - 1)
  }else{
    y_scale <- c(0,1)
    x_scale <- range(as.numeric(x))

    ypos <- 0.5
    xpos <- seq(0,1, length = color.n)

    r_width = 1/(color.n - 1)
    r_height = unit(1, "npc")
  }

  # canvas
  pushViewport(viewport(angle = 0,
                        yscale = y_scale,
                        xscale = x_scale))
  # assign colors
  if(is.null(color)){
    cols <- c("blue","white","red")
  }else{
    cols <- color
  }

  col_p <- colorRampPalette(cols)(color.n)

  # just
  if(pos == "v"){
    just <- c("bottom",rep("centre",color.n-2),"top")

    # loop to create color
    for (i in 1:color.n) {
      grid.rect(x = xpos,
                y = ypos[i],
                height = r_height,
                width = r_width,
                just = just[i],
                gp = gpar(col = col_p[i], fill = col_p[i]))
    }
  }else{
    just <- c("left",rep("centre",color.n-2),"right")

    # loop to create color
    for (i in 1:color.n) {
      grid.rect(x = xpos[i],
                y = ypos,
                height = r_height,
                width = r_width,
                just = just[i],
                gp = gpar(col = col_p[i], fill = col_p[i]))
    }
  }

  grid.rect(gp = gpar(fill = NA))

  # add axis
  if(pos == "h"){
    grid.xaxis2(side = ticks.side,tick.len = 0.25,breaks = nbreaks,...)
  }else{
    grid.yaxis2(side = ticks.side,tick.len = 0.25,breaks = nbreaks,...)
  }
  popViewport()
}


#' Create a Color Key Grid
#'
#' This function creates a color key grid using the grid graphics system.
#'
#' @param x Numeric vector specifying the range for the x-axis.
#' @param color Vector of color values to be used for the color key.
#' @param color.n Number of colors to generate.
#' @param ticks.side Side of the color key to display tick marks
#' ("left", "right", "top", or "bottom").
#' @param nbreaks the colorbar break numbers, default 5.
#' @param pos Position of the color key ("h" for horizontal, "v" for vertical).
#'
#' @return return a gTree object
#'
#' @export
grid.colorkey2 <- function(x = NULL,
                           color = NULL,
                           color.n = 100,
                           ticks.side = c("left","right","top","bottom"),
                           pos = c("h","v"),
                           nbreaks = 5){
  pos <- match.arg(pos,c("h","v"))
  ticks.side <- match.arg(ticks.side,c("left","right","top","bottom"))

  # check position
  if(pos == "v"){
    x_scale <- c(0,1)
    y_scale <- range(as.numeric(x))

    xpos <- 0.5
    ypos <- seq(0,1, length = color.n)

    r_width = unit(1, "npc")
    r_height = 1/(color.n - 1)
  }else{
    y_scale <- c(0,1)
    x_scale <- range(as.numeric(x))

    ypos <- 0.5
    xpos <- seq(0,1, length = color.n)

    r_width = 1/(color.n - 1)
    r_height = unit(1, "npc")
  }

  # canvas
  vp <- viewport(angle = 0,yscale = y_scale,xscale = x_scale)
  # assign colors
  if(is.null(color)){
    cols <- c("blue","white","red")
  }else{
    cols <- color
  }

  col_p <- colorRampPalette(cols)(color.n)

  # just
  just <- c(0,rep(0.5,color.n-2),1)

  if(pos == "v"){
    # to create color
    rect.grob <- rectGrob(x = xpos,
                          y = ypos,
                          height = r_height,
                          width = r_width,
                          vjust = just,
                          gp = gpar(col = col_p, fill = col_p),
                          vp = vp)
  }else{
    # to create color
    rect.grob <- rectGrob(x = xpos,
                          y = ypos,
                          height = r_height,
                          width = r_width,
                          hjust = just,
                          gp = gpar(col = col_p, fill = col_p))
  }

  border.grob <- rectGrob(gp = gpar(fill = NA),vp = vp)

  glist <- gTree(children = gList(rect.grob,border.grob),vp = vp)
}




#' Mapping a numeric range to colors
#'
#' @param cols the colors to be mapped, default c("white","red").
#' @param cols.n the numbers color to be generated, default 100.
#' @param range.val the data range to be mapped.
#'
#' @return a data.frame
#' @export
match.col2val <- function(cols = c("white","red"),
                          cols.n = 100,
                          range.val = NULL){
  col_p <-grDevices::colorRampPalette(cols)(cols.n)
  cut_range <- cut(range.val,cols.n)

  labs <- levels(cut_range)
  names(labs) <- col_p

  # prepare color dataframe
  left <- as.numeric(sapply(strsplit(labs, "\\(|,|\\]"),"[",2))
  right <- as.numeric(sapply(strsplit(labs, "\\(|,|\\]"),"[",3))

  range_col <- data.frame(col = names(labs),left = left,right = right)

  return(range_col)
}




#' Getting ranked row orders from a matrix
#'
#' @param norm.mat normalized matrix with attributes
#' @param row.split row.split
#' @param heatmap_rank_method heatmap_rank_method
#' @param quantile.threshold quantile.threshold
#'
#' @return list orders
#' @export
get.row.order <- function(norm.mat = NULL,
                          row.split = NULL,
                          heatmap_rank_method = NULL,
                          quantile.threshold = 1){
  # ===================================================
  qth <- quantile(norm.mat,c(0,quantile.threshold))

  mat <- apply(norm.mat, c(1,2), function(x){
    if(x > qth[2]){qth[2]}else{x}
  })

  mat <- EnrichedHeatmap::copyAttr(norm.mat,mat)

  # ===================================================

  if(is.null(row.split)){
    # row.split <- rep("genes",nrow(mat))

    # check data source
    if("deeptoolsMat" %in% class(norm.mat)){
      row.split <- rep(attr(norm.mat,"group_labels"),attr(norm.mat,"group_numbers"))
    }else{
      binding_sites <- nrow(norm.mat)
      row.split <- rep(paste0("Binding sites: ",binding_sites),binding_sites)
    }
  }

  splits <- table(row.split)

  # loop get each split order
  lapply(seq_along(splits), function(x){
    # idx <- base::grep(names(splits)[x],row.split)
    idx <- base::grep(paste0("^",names(splits)[x],"$"),row.split)
    tmp_mat <- mat[idx,]

    # assign attributes
    class(tmp_mat) <- class(mat)
    attr(tmp_mat, "upstream_index") <- attr(mat, "upstream_index")
    attr(tmp_mat, "downstream_index") <-  attr(mat, "downstream_index")
    attr(tmp_mat, "target_index") <-  attr(mat, "target_index")

    # ============================================================================
    # reoder the matrix heatmap rank method
    # ============================================================================
    if(heatmap_rank_method == "sum"){
      row_order <- order(rowSums(data.frame(tmp_mat)))
    }else if(heatmap_rank_method == "mean"){
      row_order <- order(rowMeans(data.frame(tmp_mat)))
    }else if(heatmap_rank_method == "median"){
      row_order <- order(apply(data.frame(tmp_mat),1,median))
    }else if(heatmap_rank_method == "weighting"){
      row_order <- order(EnrichedHeatmap::enriched_score(tmp_mat),decreasing = FALSE)
    }

    return(row_order)
  }) -> order.list

  names(order.list) <- names(splits)

  return(order.list)
}




#' Getting whole range from list of matrice
#'
#' @param mat.list mat.list
#' @param type the type to calculate range, "raw" for heatmap, "ave" for profile.
#' @param quantile.threshold quantile.threshold
#' @param row.split row.split
#'
#' @return range of value
#' @export
get.matlist.range <- function(mat.list = NULL,
                              type = c("raw","ave"),
                              quantile.threshold = 1,
                              row.split = NULL){
  type <- match.arg(type,c("raw","ave"))

  # x = 1
  sapply(seq_along(mat.list), function(x){
    # what type to calculate
    if(type == "raw"){
      if(quantile.threshold == 1){
        rg <- range(mat.list[[x]])
      }else{
        rg <- quantile(mat.list[[x]],c(0,quantile.threshold))
      }
    }else{
      rg <- get.col.average(norm.mat = mat.list[[x]],
                            row.split = row.split,
                            quantile.threshold = quantile.threshold)
    }

  }) |> range() -> rg.all

  return(rg.all)
}




#' Calculate mean/median for columns of different rowsplit clusters
#'
#' @param norm.mat norm.mat
#' @param row.split row.split
#' @param fun fun
#' @param quantile.threshold quantile.threshold
#'
#' @return range of value
#' @export
get.col.average <- function(norm.mat = NULL,
                            row.split = NULL,
                            fun = c("mean","median"),
                            quantile.threshold = 1){
  fun <- match.arg(fun,c("mean","median"))
  # ===================================================
  qth <- quantile(norm.mat,c(0,quantile.threshold))

  mat <- apply(norm.mat, c(1,2), function(x){
    if(x > qth[2]){qth[2]}else{x}
  })
  # ===================================================
  if(is.null(row.split)){
    row.split <- rep("genes",nrow(mat))
  }

  splits <- table(row.split)

  # loop get each split order
  # x = 2
  sapply(seq_along(splits), function(x){
    # idx <- base::grep(names(splits)[x],row.split)
    idx <- base::grep(paste0("^",names(splits)[x],"$"),row.split)
    tmp_mat <- mat[idx,]

    ave.val.rg <- range(apply(tmp_mat,2,match.fun(fun)))
    return(ave.val.rg)
  }) |> range() -> val.rg

  return(val.rg)
}




#' Getting panel boundary from group list
#'
#' @param group goup list
#'
#' @return list of boundary
#' @export
get.group.bound <- function(group = NULL){
  group.n <- sapply(group,length)
  start.panel <- cumsum(group.n) - (group.n - 1)
  end.panel <- cumsum(group.n)

  res <- list(start.panel = start.panel,end.panel = end.panel)
  return(res)
}




#' Draw sample annnotations
#'
#' @param group.list group.list
#' @param row.pos row.pos
#' @param group.anno.rect.fill group.anno.rect.fill
#' @param group.anno.rect.col group.anno.rect.col
#' @param group.anno.line.col group.anno.line.col
#' @param group.anno.line.lwd group.anno.line.lwd
#' @param group.anno.text.face group.anno.text.face
#' @param group.anno.text.size group.anno.text.size
#' @param sample.anno.type sample.anno.type
#'
#' @return grid object
#' @export
draw.anno.fun <- function(group.list = NULL,
                          row.pos = 1,
                          group.anno.rect.fill = NULL,
                          group.anno.rect.col = NULL,
                          group.anno.line.col = NULL,
                          group.anno.line.lwd = 2,
                          group.anno.text.face = "bold",
                          group.anno.text.size = 12,
                          sample.anno.type = c("rect","line")){
  sample.anno.type <- match.arg(sample.anno.type,c("rect","line"))

  gp.list <- get.group.bound(group = group.list)

  # group.anno.rect.fill = NULL
  # group.anno.rect.col = NULL
  if(!is.null(group.anno.rect.fill)){
    rect.fill <- group.anno.rect.fill
  }else{
    rect.fill <- rep("grey90",length(group.list))
  }

  if(!is.null(group.anno.rect.col)){
    rect.col <- group.anno.rect.col
  }else{
    rect.col <- rep("black",length(group.list))
  }

  # group.anno.line.col = NULL
  # group.anno.line.lwd = 2
  if(!is.null(group.anno.line.col)){
    line.col <- group.anno.line.col
  }else{
    line.col <- rep("black",length(group.list))
  }

  # draw anno
  for (i in seq_along(group.list)) {
    pushViewport(viewport(layout.pos.row = row.pos,
                          layout.pos.col = (gp.list$start.panel[i]):(gp.list$end.panel[i])))
    # sample.anno.type = c("rect","line)
    if(sample.anno.type == "rect"){
      grid.rect(gp = gpar(fill = rect.fill[i],col = rect.col[i]))
    }else{
      grid.segments(x0 = 0.02,x1 = 0.98,y0 = 0,y1 = 0,
                    gp = gpar(col = line.col[i],lwd = group.anno.line.lwd))
    }

    # group.anno.text.face = "bold"
    # group.anno.text.size = 10
    grid.text(label = names(group.list)[i],
              just = ifelse(sample.anno.type == "line","top","centre"),
              gp = gpar(fontface = group.anno.text.face,
                        fontsize = group.anno.text.size))
    upViewport()
  }
}




#' Geeting attributes from matrix
#'
#' @param mat normalizedMatrix/getTagMatrix/parseDeeptools object
#'
#' @return axis_name breaks vline.x text.x.hjust
#' @export
getAttrs <- function(mat = NULL){
  class.mat <- attr(mat,"class")

  if("normalizedMatrix" %in% class.mat || "deeptoolsMat" %in% class.mat){
    # attributes
    upstream_index = attr(mat, "upstream_index")
    downstream_index = attr(mat, "downstream_index")
    target_index = attr(mat, "target_index")

    n1 = length(upstream_index)
    n2 = length(target_index)
    n3 = length(downstream_index)
    n = n1 + n2 + n3
    nbin <- ncol(mat)

    extend = attr(mat, "extend")

    upstream_flipped = attr(mat, "upstream_flipped")
    if(is.null(upstream_flipped)) upstream_flipped = FALSE

    # xtcks <- seq(-extend[1],extend[2],length.out = ncol(mat_trim))

    axis_name = NULL

    # axis labels
    if(n1 && n2 && n3) {
      axis_name = c(paste0("-", extend[1]), "start", "end", extend[2])
      breaks = c(1,n1 + 1,n1 + n2 + 1,nbin)
      vline.x = c(n1 + 1,n1 + n2 + 1)
      text.x.hjust = c(0,0.5,0.5,1)
    } else if(n1 && !n2 && n3) {
      axis_name = c(paste0("-", extend[1]), "start", extend[2])
      breaks = c(1,n1 + 1,nbin)
      vline.x = c(n1 + 1)
      text.x.hjust = c(0,0.5,1)
    } else if(!n1 && n2 && n3) {
      axis_name = c("start", "end", extend[2])
      breaks = c(n1 + 1,n1 + n2 + 1,nbin)
      vline.x = c(n1 + n2 + 1)
      text.x.hjust = c(0,0.5,1)
    } else if(n1 && n2 && !n3) {
      axis_name = c(paste0("-", extend[1]), "start", "end")
      breaks = c(1,n1 + 1,n1 + n2 + 1)
      vline.x = c(n1 + 1)
      text.x.hjust = c(0,0.5,1)
    } else if(!n1 && n2 && !n3) {
      axis_name = c("start", "end")
      breaks = c(n1 + 1,n1 + n2 + 1)
      vline.x = c()
      text.x.hjust = c(0,1)
    } else if(n1 && !n2 && !n3) {
      axis_name = c(paste0("-", extend[1]), "start")
      breaks = c(1,n1 + 1)
      vline.x = c()
      text.x.hjust = c(0,1)
    } else if(!n1 && !n2 && n3) {
      if(upstream_flipped) {
        axis_name = c("start", extend[2])
        breaks = c(n1 + 1,nbin)
        vline.x = c()
        text.x.hjust = c(0,1)
      } else {
        axis_name = c("end", extend[2])
        breaks = c(n1 + n2 + 1,nbin)
        vline.x = c()
        text.x.hjust = c(0,1)
      }
    }
  }else{
    upstream <- attr(mat, "upstream")
    downstream <- attr(mat, "downstream")
    binning_Flag <- attr(mat,"is.binning")
    type <- attr(mat,"type")

    extend = c(upstream,downstream)

    body_Flag <- FALSE
    if(type == "body"){
      body_Flag <- TRUE
      label <- attr(mat,"label")
    }

    if(binning_Flag){
      nbin <- dim(mat)[2]
    }else{
      nbin <- ncol(mat)
    }

    if(type == "body"){
      upstreamPer <- floor(upstream/1000)*0.1
      downstreamPer <- floor(downstream/1000)*0.1

      breaks=c(1,
               floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
               # floor(nbin*((upstreamPer+0.25)/(1+upstreamPer+downstreamPer))),
               # floor(nbin*((upstreamPer+0.5)/(1+upstreamPer+downstreamPer))),
               # floor(nbin*((upstreamPer+0.75)/(1+upstreamPer+downstreamPer))),
               floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))),
               nbin)

      axis_name=c(paste0("-",upstream),
                  "start", # label[1],
                  # "25%",
                  # "50%",
                  # "75%",
                  "end",  # label[2],
                  paste0(downstream))

      vline.x <- c(floor(nbin*(upstreamPer/(1+upstreamPer+downstreamPer))),
                   floor(nbin*((upstreamPer+1)/(1+upstreamPer+downstreamPer))))

      text.x.hjust <- c(0,0.5,0.5,1)
    }else{
      if(binning_Flag){
        breaks = c(1,
                   # floor(nbin*(downstream*0.5/(downstream+upstream))),
                   floor(nbin*(downstream/(downstream+upstream))),
                   # floor(nbin*((downstream + upstream*0.5)/(downstream+upstream))),
                   nbin)

      }else{
        breaks = c(1,
                   # floor(downstream*0.5),
                   (downstream + 1),
                   # (downstream + 1 + floor(upstream * 0.5)),
                   upstream+downstream+1)
      }

      axis_name = c((-1*downstream),
                    # floor(-1*downstream*0.5),
                    0,
                    # floor(upstream*0.5),
                    upstream)

      vline.x <- c()
      text.x.hjust <- c(0,0.5,1)
    }
  }


  attrs <- list(axis_name = axis_name,
                breaks = breaks,
                vline.x = vline.x,
                text.x.hjust = text.x.hjust)

  return(attrs)
}
