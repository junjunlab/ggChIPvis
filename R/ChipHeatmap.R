#' Create a heatmap of chip data
#'
#' @author Jun Zhang
#'
#' @description
#' This function creates a single heatmap of chip data using grid.
#'
#' @param norm.mat A numeric matrix containing normalized chip data which from
#' normalizeToMatrix/getTagMatrix/parseDeeptools.
#' @param sample.label A character string specifying the label for the sample.
#' Default is "sample name".
#' @param label.rect.gp A \code{gpar} object specifying the graphic parameters
#' for the label rectangle. Default is \code{gpar(fill = "grey90")}.
#' @param label.gp A \code{gpar} object specifying the graphic parameters for
#' the label. Default is \code{gpar(fontface = "bold")}.
#' @param row.split A numeric vector specifying the row split index.
#' Default is \code{NULL}.
#' @param row.order A numeric vector specifying the row order index.
#' Default is \code{NULL}.
#' @param ht.col A character vector specifying the colors for heatmap.
#' Default is \code{c("white","red")}.
#' @param ht.col.legend.nbreaks The number of colorbar legend breaks. Default is 4.
#' @param ht.xaxis.rot A numeric value specifying the rotation angle of
#' the x-axis labels. Default is 0.
#' @param col.range A numeric vector specifying the range of colors for
#' the heatmap. Default is \code{NULL}.
#' @param heatmap_rank_method A character string specifying the method for
#' ranking the heatmap. The available options are "weighting", "sum",
#' "mean", and "median". Default is "weighting". Note "weighting" can not
#' be used for data which is from "getTagMatrix" function.
#' @param draw.legend A logical value indicating whether to draw the legend.
#' Default is \code{TRUE}.
#' @param legend.pos A character vector specifying the position of the legend.
#' Default is \code{c("right","bottom")}.
#' @param cluster.gap A numeric value specifying the gap between clusters.
#' Default is 0.01.
#' @param keep.cluster.panel.same A logical value indicating whether to keep
#' the cluster panel the same. Default is \code{FALSE}.
#' @param quantile.threshold A numeric value specifying the quantile threshold.
#' Default is 1.
#' @param draw.anno.left A logical value indicating whether to draw the left
#' annotation. Default is \code{TRUE}.
#' @param draw.anno.left.shape A character vector specifying the shape of the
#' left annotation. Default is \code{c("tri","rect")}.
#' @param anno.left.col A character vector specifying the color of the left
#' annotation. Default is \code{NULL}.
#' @param draw.profile A logical value indicating whether to draw the profile.
#' Default is \code{FALSE}.
#' @param profile.line.col A character vector specifying the color of the
#' profile line. Default is \code{NULL}.
#' @param profile.line.size A numeric value specifying the size of the profile
#' line. Default is 1.
#' @param profile.y.range A numeric vector specifying the y-range of the
#' profile. Default is \code{NULL}.
#' @param profile.yaxis A logical value indicating whether to draw the y-axis
#' of the profile. Default is \code{TRUE}.
#' @param draw.line.legend A logical value indicating whether to draw the line
#' legend. Default is \code{TRUE}.
#' @param profile.legend.size A numeric vector specifying the size of the line
#' legend. Default is \code{c(0.25,0.5)}.
#' @param profile.legend.fontsize A numeric value specifying the font size of
#' the line legend. Default is 8.
#' @param interpolate A logical value indicating whether to interpolate the
#' heatmap. Default is \code{FALSE}.
#' @param HeatmapLayout.params A list specifying the parameters for the heatmap
#' layout. Default is \code{list()}.
#' @param panel.col.space A numeric value specifying the space between panels
#' in the column direction. Default is 2.
#' @param panel.row.space A numeric value specifying the space between panels
#' in the row direction. Default is 2.
#' @param draw.anno.left.diretion The arrow of heatmap left annotation
#' bar direction "down"(default) or "up".
#' @param add.profile.rect Whether add rect border for the profile plot.
#' Default is TRUE.
#'
#' @return A gtable object representing the heatmap.
#'
#' @import grid
#' @import gtable
#' @export
ChipHeatmap <- function(norm.mat = NULL,
                        sample.label = "sample name",
                        label.rect.gp = gpar(fill = "grey90"),
                        label.gp = gpar(fontface = "bold"),
                        row.split = NULL,
                        row.order = NULL,
                        ht.col = c("white","red"),
                        ht.col.legend.nbreaks = 4,
                        ht.xaxis.rot = 0,
                        col.range = NULL,
                        heatmap_rank_method = c("weighting","sum","mean","median"),
                        draw.legend = TRUE,
                        legend.pos = c("right","bottom"),
                        cluster.gap = 0.01,
                        keep.cluster.panel.same = FALSE,
                        quantile.threshold = 1,
                        draw.anno.left = TRUE,
                        draw.anno.left.shape = c("tri","rect"),
                        draw.anno.left.diretion = c("down","up"),
                        anno.left.col = NULL,
                        draw.profile = FALSE,
                        add.profile.rect = TRUE,
                        profile.line.col = NULL,
                        profile.line.size = 1,
                        profile.y.range = NULL,
                        profile.yaxis = TRUE,
                        draw.line.legend = TRUE,
                        profile.legend.size = c(0.25,0.5),
                        profile.legend.fontsize = 8,
                        interpolate = FALSE,
                        HeatmapLayout.params = list(),
                        panel.col.space = 2,
                        panel.row.space = 2){
  # ============================================================================
  # check args
  # ============================================================================
  heatmap_rank_method <- match.arg(heatmap_rank_method,c("weighting","sum","mean","median"))
  draw.anno.left.shape <- match.arg(draw.anno.left.shape,c("tri","rect"))
  legend.pos <- match.arg(legend.pos,c("right","bottom"))
  draw.anno.left.diretion <- match.arg(draw.anno.left.diretion,c("down","up"))
  # ============================================================================
  # get attrs
  # ============================================================================
  attrs <- getAttrs(mat = norm.mat)

  # whether remove extrme values
  # quantile.threshold = 0.99
  qth <- quantile(norm.mat,c(0,quantile.threshold))

  mat <- apply(norm.mat, c(1,2), function(x){
    if(x > qth[2]){qth[2]}else{x}
  })

  mat <- EnrichedHeatmap::copyAttr(norm.mat,mat)

  # ============================================================================
  # calculate each rowsplit cluster panel coordinates
  # ============================================================================
  # panel.gap <- 0.01
  if(is.null(row.split)){
    # row.split <- rep("genes",nrow(mat))

    # check data source
    if("deeptoolsMat" %in% class(norm.mat)){
      row.split <- rep(attr(norm.mat,"group_labels"),attr(norm.mat,"group_numbers"))
    }else{
      binding_sites <- nrow(norm.mat)
      row.split <- rep(paste0("Binding sites: ",binding_sites),binding_sites)
    }

    splits <- table(row.split)
    n.split <- length(splits)

    # calculate each panels y coordinate
    if("deeptoolsMat" %in% class(norm.mat)){
      if(n.split > 1){
        # calculate each panels y coordinate
        if(keep.cluster.panel.same == TRUE){
          htmp <- (1 - (n.split - 1)*cluster.gap)/n.split
          h.panels <- rep(htmp,n.split)
        }else{
          h.panels <- ((1 - (n.split - 1)*cluster.gap)/nrow(mat))*splits
        }

        h.panels.y <- as.numeric(h.panels/2 + c(0,cumsum(h.panels)[1:(n.split - 1)]))
        h.panels.y <- h.panels.y + cluster.gap*(0:(n.split - 1))
      }else{
        h.panels <- 1
        # h.panels.y <- ht.y
      }
    }else{
      h.panels <- 1
      # h.panels.y <- ht.y
    }

  }else{
    splits <- table(row.split)
    n.split <- length(splits)

    # calculate each panels y coordinate
    # keep.cluster.panel.same = FALSE
    if(keep.cluster.panel.same == TRUE){
      htmp <- (1 - (n.split - 1)*cluster.gap)/n.split
      h.panels <- rep(htmp,n.split)
    }else{
      h.panels <- ((1 - (n.split - 1)*cluster.gap)/nrow(mat))*splits
    }

    h.panels.y <- as.numeric(h.panels/2 + c(0,cumsum(h.panels)[1:(n.split - 1)]))
    h.panels.y <- h.panels.y + cluster.gap*(0:(n.split - 1))
  }

  # ============================================================================
  # design heatmap layout
  # ============================================================================
  # HeatmapLayout.params = list()
  main.vp <-
    do.call(HeatmapLayout,modifyList(list(mat = mat,
                                          draw.profile = draw.profile,
                                          legend.pos = legend.pos,
                                          newpage = FALSE),
                                     HeatmapLayout.params))

  # ============================================================================
  # draw triangle legend
  # ============================================================================
  # draw.anno.left = TRUE
  # anno.left.col = NULL
  if(is.null(anno.left.col)){
    tri.col <- rep("black",n.split)
  }else{
    tri.col <- anno.left.col
  }

  # draw.anno.left.shape = "rect"
  # grid.newpage()
  if(draw.anno.left == TRUE){
    # ==========================================================================
    # loop generate left anno bar grobs
    lapply(seq_along(splits), function(i){
      # annobar
      if(n.split == 1){
        vp.left.anno <- viewport()
      }else{
        vp.left.anno <- viewport(y = ifelse(n.split == 1,0.5,1 - h.panels.y[i]),
                                 width = unit(1,"npc"),height = h.panels[i])
      }

      # draw.anno.left.shape = "tri"
      if(draw.anno.left.shape == "tri"){
        # draw.anno.left.diretion = c("down","up)
        if(draw.anno.left.diretion == "down"){
          left.anno.grob <- polygonGrob(x = c(1,0,1),y = c(0,1,1),
                                        gp = gpar(fill = tri.col[i],col = tri.col[i]),
                                        vp = vp.left.anno)
        }else{
          left.anno.grob <- polygonGrob(x = c(0,1,1),y = c(0,1,0),
                                        gp = gpar(fill = tri.col[i],col = tri.col[i]),
                                        vp = vp.left.anno)
        }
      }else{
        left.anno.grob <- rectGrob(gp = gpar(fill = tri.col[i],col = tri.col[i]),
                                   vp = vp.left.anno)
      }

      return(left.anno.grob)
    }) -> left.anno.grob.list

    main.vp <- gtable::gtable_add_grob(main.vp,
                                       grobs = left.anno.grob.list,
                                       t = 3,l = 2,
                                       name = names(splits))

    # grid.draw(main.vp)

    # ==========================================================================

    # loop generate left anno bar label grobs
    lapply(seq_along(splits), function(i){
      # annobar label
      vp.left.anno.label <- viewport(x = 0.5,
                                     y = ifelse(n.split == 1,0.5,1 - h.panels.y[i]),
                                     width = 1,just = "centre",
                                     height = h.panels[i])

      idx <- base::grep(paste0("^",names(splits)[i],"$"),row.split)
      tmp_mat <- mat[idx,]

      left.anno.label.grob <- textGrob(label = paste0(names(splits)[i],
                                                      "\n(n=",nrow(tmp_mat),")"),
                                       rot = 90,vp = vp.left.anno.label)

      return(left.anno.label.grob)
    }) -> left.anno.label.grob.list

    main.vp <- gtable::gtable_add_grob(main.vp,
                                       grobs = left.anno.label.grob.list,
                                       t = 3,l = 1,
                                       name = names(splits),
                                       clip = "off")

    # grid.draw(main.vp)
  }

  # ============================================================================
  # heatmap main viewport
  # ============================================================================
  # i = 1
  # heatmap_rank_method = "sum"
  lapply(seq_along(splits), function(i){
    # idx <- base::grep(names(splits)[i],row.split)
    idx <- base::grep(paste0("^",names(splits)[i],"$"),row.split)
    tmp_mat <- mat[idx,]

    # assign attributes
    class(tmp_mat) <- class(mat)
    attr(tmp_mat, "upstream_index") <- attr(mat, "upstream_index")
    attr(tmp_mat, "downstream_index") <-  attr(mat, "downstream_index")
    attr(tmp_mat, "target_index") <-  attr(mat, "target_index")

    # ============================================================================
    # viewport
    vp.heatmap <- viewport(x = 0.5,
                           y = ifelse(n.split == 1,0.5,1 - h.panels.y[i]),
                           width = 1,
                           height = h.panels[i],
                           xscale = c(1,ncol(tmp_mat)),yscale = c(0,nrow(tmp_mat)))
    # ============================================================================
    # reoder the matrix heatmap rank method
    # ============================================================================
    # row.order = NULL
    # heatmap_rank_method = "weighting"
    if(is.null(row.order)){
      if(heatmap_rank_method == "sum"){
        row_order <- order(rowSums(data.frame(tmp_mat)))
      }else if(heatmap_rank_method == "mean"){
        row_order <- order(rowMeans(data.frame(tmp_mat)))
      }else if(heatmap_rank_method == "median"){
        row_order <- order(apply(data.frame(tmp_mat),1,median))
      }else if(heatmap_rank_method == "weighting"){
        row_order <- order(EnrichedHeatmap::enriched_score(tmp_mat),decreasing = FALSE)
      }
    }else{
      row_order <- row.order
    }

    new.mat <- tmp_mat[rev(row_order),]
    # ============================================================================
    # draw heatmap
    # ============================================================================
    df.coord <- suppressWarnings(reshape2::melt(as.matrix(new.mat)))
    colnames(df.coord)[1:2] <- c("y","x")

    # color mapping range
    # col.range = NULL
    if(is.null(col.range)){
      range.val <- range(mat)
    }else{
      range.val <- col.range
    }

    range_col <- match.col2val(cols = ht.col,cols.n = 100,range.val = range.val)

    # add colors for each value
    # x = 100
    purrr::map_df(1:nrow(range_col), function(x){
      tmp <- range_col[x,]

      tmp.coord <- subset(df.coord,value >= tmp$left & value < tmp$right)

      if(nrow(tmp.coord) > 0){
        tmp.coord$col_f <- tmp$col
        return(tmp.coord)
      }else{
        return(NULL)
      }
    }) -> df.coord.new

    # matrix with colors filled
    col.mat <- suppressMessages(df.coord.new |> dplyr::select(-value) |>
                                  reshape2::dcast(y ~ x) |> dplyr::select(-y))

    rownames(col.mat) <- 1:nrow(col.mat)
    col.mat <- col.mat[as.character(1:nrow(col.mat)),]

    # ============================================================================
    # heatmap
    # interpolate = FALSE
    heatmap.grob <- rasterGrob(image = as.matrix(col.mat),width = 1,height = 1,
                               interpolate = interpolate,
                               vp = vp.heatmap)

    # vertical line
    vline.grob <- segmentsGrob(x0 = unit(attrs$vline.x,"native"),x1 = unit(attrs$vline.x,"native"),
                               y0 = 0,y1 = 1,
                               gp = gpar(lty = "dashed"),
                               vp = vp.heatmap)

    cluster.border.grob <- rectGrob(gp = gpar(fill = NA,col = "black",lwd = 1),vp = vp.heatmap)

    comb <- gTree(children = gList(heatmap.grob,vline.grob,cluster.border.grob))

    return(comb)
  }) -> heatmap.body.grobs.list

  main.vp <- gtable::gtable_add_grob(main.vp,
                                     grobs = heatmap.body.grobs.list,
                                     name = if(n.split == 1){1}else{1:n.split},
                                     t = 3,l = 3)

  # grid.draw(main.vp)

  # ============================================================================

  # draw x axis for heatmap
  # ht.xaxis.rot = 0
  # xaxis.grob <- suppressWarnings(xaxisGrob(at = attrs$breaks,label = attrs$axis_name,
  #                                          vp = viewport(xscale = c(1,ncol(mat))),
  #                                          edits = gEdit(gPath="labels",rot = ht.xaxis.rot)))

  if(ht.xaxis.rot == 0){
    hjust <- attrs$text.x.hjust
  }else{
    hjust <- 0.5
  }

  xaxis.grob <- grid.xaxis2(at = attrs$breaks,
                            labels = attrs$axis_name,
                            hjust = hjust,
                            rot = ht.xaxis.rot,
                            draw = FALSE,vp = viewport(xscale = c(1,ncol(mat))))

  main.vp <- gtable::gtable_add_grob(main.vp,
                                     grobs = xaxis.grob,
                                     t = 3,l = 3,
                                     name = "xaxis",clip = "off")

  # grid.draw(main.vp)
  # ============================================================================
  # draw profile plot
  # ============================================================================
  # grid.newpage()
  # draw.profile = FALSE
  # get range for each split
  # profile.y.range = NULL
  if(is.null(profile.y.range)){
    prg <- get.col.average(norm.mat = mat,row.split = row.split,
                           quantile.threshold = quantile.threshold) |>
      extendrange(f = 0.1)
  }else{
    prg <- profile.y.range
  }


  # draw
  # draw.profile = TRUE
  if(draw.profile == TRUE){
    vp.profile <- viewport(xscale = c(1,ncol(mat)),
                           yscale = c(0,prg[2]),just = "centre")

    # add.profile.rect = TRUE
    if(add.profile.rect == TRUE){
      profile.rect.grob <- rectGrob(vp = vp.profile)
    }else{
      profile.rect.grob <- zeroGrob()
    }

    # profie.yaxis = TRUE
    if(profile.yaxis == TRUE){
      # profile.yaxis.params = list()
      # yaxis.grob <- do.call(grid.yaxis2,modifyList(list(breaks = 3,draw = FALSE,vp = vp.profile),
      #                                              profile.yaxis.params))
      yaxis.grob <- yaxisGrob(vp = vp.profile)
      # yaxis.grob <- grid.yaxis2(breaks = 1,draw = FALSE,vp = vp.profile)
    }else{
      yaxis.grob <- zeroGrob()
    }

    # loop to draw each split cluster
    # profile.line.col = NULL
    # profile.line.size = 1
    if(is.null(profile.line.col)){
      line.col <- RColorBrewer::brewer.pal(9,"Set1")
    }else{
      line.col <- profile.line.col
    }

    # i = 2
    # line.grob.list <- list()
    lapply(seq_along(splits), function(i){
      idx <- base::grep(paste0("^",names(splits)[i],"$"),row.split)
      tmp_mat0 <- mat[idx,]

      line.grob <- linesGrob(x = 1:ncol(mat),
                             y = apply(tmp_mat0,2,mean),
                             gp = gpar(col = line.col[i],lwd = profile.line.size),
                             default.units = "native",
                             vp = vp.profile)

      return(line.grob)
    }) -> line.grob.list

    profile.vline.grob <- segmentsGrob(x0 = unit(attrs$vline.x,"native"),
                                       x1 = unit(attrs$vline.x,"native"),
                                       y0 = 0,y1 = 1,gp = gpar(lty = "dashed"),
                                       vp = vp.profile)

    # draw legend for lines
    # profile.legend.size = c(0.5,1)
    # profile.legend.fontsize = 8
    # draw.line.legend = TRUE
    if(draw.line.legend == TRUE){
      vp.profile.legend <- viewport(x = 0,y = 1,
                                    width = profile.legend.size[1],
                                    height = profile.legend.size[2],
                                    just = c("left","top"))

      legend.grob <- legendGrob(labels = names(splits),ncol = 1,
                                hgap = unit(0.1, "lines"), vgap = unit(0.1, "lines"),
                                do.lines = T,
                                gp = gpar(col = line.col[1:n.split],
                                          fontsize = profile.legend.fontsize),
                                vp = vp.profile.legend)
    }else{
      legend.grob <- zeroGrob()
    }

    profile.list <- list(profile.rect.grob,yaxis.grob,
                         profile.vline.grob,legend.grob)

    main.vp <- gtable::gtable_add_grob(main.vp,
                                       grobs = profile.list,
                                       t = 2,l = 3,
                                       name = c("rect.p","yaxis.p","vline.p","legend.p"),
                                       clip = "off")

    main.vp <- gtable::gtable_add_grob(main.vp,
                                       grobs = line.grob.list,
                                       t = 2,l = 3,
                                       name = names(splits))

    # grid.draw(main.vp)
  }

  # ============================================================================
  # add heatmap labels
  # ============================================================================
  # grid.newpage()
  # label.rect.gp = gpar(fill = "grey90")
  # label.gp = gpar(fontface = "bold.italic")

  vp.sample.label <- viewport(just = "centre")

  sample.label.rect.grob <- rectGrob(gp = label.rect.gp,
                                     vp = vp.sample.label)

  sample.label.label.grob <- textGrob(label = sample.label,gp = label.gp,
                                      vp = vp.sample.label)

  main.vp <- gtable::gtable_add_grob(main.vp,
                                     grobs = list(sample.label.rect.grob,
                                                  sample.label.label.grob),
                                     t = 1,l = 3,
                                     name = c("rect","sample.label"))

  # grid.draw(main.vp)
  # ============================================================================
  # draw color legend
  # ============================================================================
  # legend.pos = "right"

  # legend color range
  if(is.null(col.range)){
    legend.col.rg <- range(mat)
  }else{
    legend.col.rg <- col.range
  }

  # draw colorbar
  # draw.legend = TRUE
  if(draw.legend == TRUE){
    if(legend.pos == "right"){
      vp.rlegend <- viewport(yscale = legend.col.rg)

      colorbar.grob <- grid.colorkey2(x = legend.col.rg,pos = "v",ticks.side = "right",
                                      color = ht.col)

      # axis.grob <- yaxisGrob(main = FALSE,vp = vp.rlegend)
      axis.grob <- grid.yaxis2(breaks = ht.col.legend.nbreaks,side = "right",draw = FALSE,vp = vp.rlegend)

      panel.pos <- c(3,4)
    }else if(legend.pos == "bottom"){
      vp.blegend <- viewport(xscale = legend.col.rg)

      colorbar.grob <- grid.colorkey2(x = legend.col.rg,pos = "h",ticks.side = "bottom",
                                      color = ht.col)

      # axis.grob <- xaxisGrob(vp = vp.blegend)
      axis.grob <- grid.xaxis2(breaks = ht.col.legend.nbreaks,side = "bottom",draw = FALSE,vp = vp.blegend)

      panel.pos <- c(5,3)
    }

    # add to panel
    main.vp <- gtable::gtable_add_grob(main.vp,
                                       grobs = list(colorbar.grob,axis.grob),
                                       t = panel.pos[1],l = panel.pos[2],
                                       name = c("colorbar","colorbar.label"),
                                       clip = "off")

  }

  # ============================================================================
  # panel space
  # ============================================================================
  # panel.col.space = 2
  # panel.row.space = 2
  main.vp <- gtable::gtable_add_col_space(main.vp,width = unit(panel.col.space,"mm"))
  main.vp <- gtable::gtable_add_row_space(main.vp,height = unit(panel.row.space,"mm"))

  # draw final plot
  grid.draw(main.vp)
}
