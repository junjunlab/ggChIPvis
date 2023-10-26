#' Match Values to Colors and Retrieve Color Data Frame
#'
#' This function maps values to colors based on specified color values and
#' optional value limits, and returns a data frame with color, left value, and
#' right value columns.
#'
#' @param x A numeric vector of values to be matched with colors. Default is NULL.
#' @param cols A character vector of color values. Default is c("white", "red").
#' @param color.n The number of color values to generate. Default is 100.
#' @param value.limits A numeric vector specifying the value limits for mapping
#' to colors. Default is NULL.
#'
#' @return A data frame with columns: 'col' for colors, 'left' for left values,
#' and 'right' for right values.
#'
#' @seealso \code{\link[circlize]{colorRamp2}}, \code{\link[circlize]{col2value}}
#'
#' @importFrom circlize colorRamp2 col2value
#'
#' @export
match.col2val <- function(x = NULL,
                          cols = c("white","red"),
                          color.n = 100,
                          value.limits = NULL){
  if(is.null(value.limits)){
    limits <- range(x)
  }else{
    limits <- value.limits
  }

  new.x <- seq(range(x)[1],range(x)[2],length.out = color.n)

  col_fun <- circlize::colorRamp2(limits, cols)
  col_p <- col_fun(new.x)
  value <- circlize::col2value(col_p, col_fun = col_fun)

  value1 <- value[1:(color.n - 1)]
  value2 <- value[2:color.n]

  df <- data.frame(col = c(col_p[1:(color.n - 1)],col_p[color.n]),
                   left = c(value1,value[color.n]),
                   right = c(value2,value[color.n] + 1))

  return(df)
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
