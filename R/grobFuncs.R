#' grid.xaxis2 function
#'
#' @author Jun Zhang
#'
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
#' @param name grob name.
#' @param gp gpar() seetings.
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
                        name = NULL, gp = NULL, vp = NULL) {
  lst <- list(at = at,
              breaks = breaks,
              labels = labels,
              tick.len = tick.len,
              label.space = label.space,
              side = side,
              rot = rot,
              ticks.gp = ticks.gp,
              text.gp = text.gp,
              hjust = hjust,
              vjust = vjust,
              draw = draw,
              name = name, gp = gp, vp = vp,
              cl = "xaxis2")

  if(draw == TRUE){
    grid.draw(do.call(gTree,lst))
  }else{
    do.call(gTree,lst)
  }
}



# makeContent method for grid.xaxis2
#' @export
makeContent.xaxis2 <- function(x) {
  # labels and ticks
  if(is.null(x$at) || is.null(x$labels)){
    if(is.null(x$vp)){
      at <- grid.pretty(current.viewport()$xscale,n = x$breaks)
    }else{
      at <- grid.pretty(x$vp$xscale,n = x$breaks)
    }
    labels <- as.character(at)
  }else{
    at <- x$at
    labels <- as.character(x$labels)
  }

  # axis position
  side <- match.arg(x$side,c("bottom","top"))
  if(side == "bottom"){
    tck.y0 = unit(0, "npc")
    tck.y1 = unit(-x$tick.len, "lines")
    text.y = unit(-x$tick.len - x$label.space,"lines")
  }else{
    tck.y0 = unit(1, "npc")
    tck.y1 = unit(1, "npc") + unit(x$tick.len, "lines")
    text.y = unit(abs(x$tick.len) + abs(x$label.space),"lines") + unit(1, "npc")
  }

  seg.main <- segmentsGrob(x0 = 0,x1 = 1,y0 = 0,y1 = 0,name = "main")

  ticks <- segmentsGrob(x0 = unit(at, "native"),x1 = unit(at, "native"),
                        y0 = tck.y0,y1 = tck.y1,
                        gp = x$ticks.gp,name = "ticks")

  text <- textGrob(label = labels,
                   x = unit(at, "native"),y = text.y,
                   rot = x$rot,
                   hjust = x$hjust,vjust = x$vjust,
                   gp = x$text.gp,name = "labels")

  setChildren(x, gList(seg.main,ticks,text))
}



#' grid.yaxis2 function
#'
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
#' @param name grob name.
#' @param gp gpar() seetings.
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
                        name = NULL, gp = NULL, vp = NULL) {
  lst <- list(at = at,
              breaks = breaks,
              labels = labels,
              tick.len = tick.len,
              label.space = label.space,
              side = side,
              rot = rot,
              ticks.gp = ticks.gp,
              text.gp = text.gp,
              hjust = hjust,
              vjust = vjust,
              draw = draw,
              name = name, gp = gp, vp = vp,
              cl = "yaxis2")

  if(draw == TRUE){
    grid.draw(do.call(gTree,lst))
  }else{
    do.call(gTree,lst)
  }
}


# makeContent method for grid.yaxis2
#' @export
makeContent.yaxis2 <- function(x) {
  # labels and ticks
  if(is.null(x$at) || is.null(x$labels)){
    if(is.null(x$vp)){
      at <- grid.pretty(current.viewport()$yscale,n = x$breaks)
    }else{
      at <- grid.pretty(x$vp$yscale,n = x$breaks)
    }
    labels <- as.character(at)
  }else{
    at <- x$at
    labels <- as.character(x$labels)
  }

  # axis position
  side <- match.arg(x$side,c("left","right"))
  if(side == "left"){
    tck.x0 = unit(0, "npc")
    tck.x1 = unit(-x$tick.len, "lines")
    text.x = unit(-x$tick.len - x$label.space,"lines")
    text.just = "right"
  }else{
    tck.x0 = unit(1, "npc")
    tck.x1 = unit(1, "npc") + unit(x$tick.len, "lines")
    text.x = unit(abs(x$tick.len) + abs(x$label.space),"lines") + unit(1, "npc")
    text.just = "left"
  }

  seg.main <- segmentsGrob(x0 = 0,x1 = 0,y0 = 0,y1 = 1,name = "main")

  ticks <- segmentsGrob(y0 = unit(at, "native"),y1 = unit(at, "native"),
                        x0 = tck.x0,x1 = tck.x1,
                        gp = x$ticks.gp,name = "ticks")

  text <- textGrob(label = labels,
                   y = unit(at, "native"),x = text.x,
                   rot = x$rot,
                   just = text.just,hjust = x$hjust,vjust = x$vjust,
                   gp = x$text.gp,name = "labels")

  setChildren(x, gList(seg.main,ticks,text))
}

# ==============================================================================
# color legend
# ==============================================================================


#' Create a Color Key for Grid Graphics
#'
#' This function creates a color key for use in grid graphics.
#'
#' @param x The values for which colors will be displayed. Default is NULL.
#' @param color The colors to be used in the color key. Default is NULL.
#' @param color.n The number of colors to generate. Default is 100.
#' @param value.limits The limits of the values. Default is NULL.
#' @param ticks.side The side on which to display tick marks. Default is "left".
#' @param pos The position of the color key ('v' for vertical, 'h' for horizontal).
#' Default is "v".
#' @param nbreaks The number of breaks for the color key. Default is 5.
#' @param draw Logical. Should the color key be drawn immediately? Default is TRUE.
#' @param xaxis2.params Parameters for the x-axis. Default is an empty list.
#' @param yaxis2.params Parameters for the y-axis. Default is an empty list.
#' @param name The name of the color key. Default is NULL.
#' @param gp A graphical parameter object. Default is NULL.
#' @param vp A viewport object. Default is NULL.
#'
#' @return If \code{draw} is \code{TRUE}, the color key is drawn. Otherwise,
#' a grob is returned.
#'
#' @seealso \code{\link[grid]{grid.draw}}, \code{\link[grid]{gTree}}
#'
#' @export
grid.colorkey <- function(x = NULL,
                          color = NULL,
                          color.n = 100,
                          value.limits = NULL,
                          ticks.side = c("left","right","top","bottom"),
                          pos = c("v","h"),
                          nbreaks = 5,
                          draw = TRUE,
                          xaxis2.params = list(),
                          yaxis2.params = list(),
                          name = NULL, gp = NULL, vp = NULL){
  pos <- match.arg(pos,c("v","h"))
  ticks.side <- match.arg(ticks.side,c("left","right","top","bottom"))

  lst <- list(x = x,
              color = color,
              color.n = color.n,
              value.limits = value.limits,
              ticks.side = ticks.side,
              pos = pos,
              nbreaks = nbreaks,
              draw = draw,
              xaxis2.params = xaxis2.params,
              yaxis2.params = yaxis2.params,
              name = name, gp = gp, vp = vp,
              cl = "colorkey")

  if(draw == TRUE){
    grid.draw(do.call(gTree,lst))
  }else{
    do.call(gTree,lst)
  }
}


# makeContext method for grid.colorkey
#' @export
makeContext.colorkey <- function(x) {
  # check position
  if(x$pos == "v"){
    x_scale <- c(0,1)
    if(is.null(x$x)){
      y_scale <- c(0,1)
    }else{
      y_scale <- range(as.numeric(x$x))
    }
  }else{
    y_scale <- c(0,1)
    if(is.null(x$x)){
      x_scale <- c(0,1)
    }else{
      x_scale <- range(as.numeric(x$x))
    }
  }

  # viewport
  if (is.null(x$vp)){
    x$vp <- viewport(yscale = y_scale,xscale = x_scale)
  }else{
    x$vp <- x$vp
  }

  x
}


# makeContent method for grid.colorkey
#' @export
makeContent.colorkey <- function(x){
  # check position
  if(x$pos == "v"){
    xpos <- 0.5
    ypos <- seq(0,1, length = x$color.n)

    r_width = unit(1, "npc")
    r_height = 1/(x$color.n - 1)
  }else{
    ypos <- 0.5
    xpos <- seq(0,1, length = x$color.n)

    r_width = 1/(x$color.n - 1)
    r_height = unit(1, "npc")
  }

  # assign colors
  if(is.null(x$color)){
    cols <- c("white","red")
  }else{
    cols <- x$color
  }

  if(is.null(x$value.limits)){
    limits <- range(x$x)
  }else{
    limits <- x$value.limits
  }

  new.x <- seq(range(x$x)[1],range(x$x)[2],length.out = x$color.n)

  # col_p <- colorRampPalette(cols)(x$color.n)
  col_fun <- circlize::colorRamp2(limits, cols)
  col_p <- col_fun(sort(new.x))
  value <- circlize::col2value(col_p, col_fun = col_fun)

  # just
  just <- c(0,rep(0.5,x$color.n-2),1)

  if(x$pos == "v"){
    # to create color
    rect.grob <- rectGrob(x = xpos,
                          y = ypos,
                          height = r_height,
                          width = r_width,
                          vjust = just,
                          gp = gpar(col = col_p, fill = col_p))
  }else{
    # to create color
    rect.grob <- rectGrob(x = xpos,
                          y = ypos,
                          height = r_height,
                          width = r_width,
                          hjust = just,
                          gp = gpar(col = col_p, fill = col_p))
  }

  border.grob <- rectGrob(gp = gpar(fill = NA))

  # add axis
  if(x$pos == "h"){
    # xaxis2.params = list()
    axis2 <- do.call(grid.xaxis2,modifyList(list(side = x$ticks.side,
                                                 tick.len = 0.25,
                                                 draw = FALSE,
                                                 breaks = x$nbreaks),
                                            x$xaxis2.params))
  }else{
    # yaxis2.params = list()
    axis2 <- do.call(grid.yaxis2,modifyList(list(side = x$ticks.side,
                                                 tick.len = 0.25,
                                                 draw = FALSE,
                                                 breaks = x$nbreaks),
                                            x$yaxis2.params))
  }

  # combine grobs
  setChildren(x, gList(rect.grob,border.grob,axis2))
}
