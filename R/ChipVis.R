globalVariables(c("group", "group2", "value", "x", "y", "yend","split.group"))

#' ChipVis Function
#'
#' @author Jun Zhang
#'
#' @description
#' Visualize ChIP-seq data using both profile and heatmap plots. This function use
#' output data from "retriveData" function and use ggplot2 package to produce profile
#' and heatmap plots. Most parameters can be accepted in a list form with ggplot
#' arguments.
#'
#' @param object An object containing the ChIP-seq data which from retriveData function.
#' Must have "profile" and "heatmap" components.
#' @param plot.type Type of plot to generate. Can be "both" (default), "profile", or "heatmap".
#' @param sample.order Order of the samples for plotting. Default is NULL.
#' @param group.order Order of the groups for plotting. Default is NULL.
#' @param group2.order Order of the secondary groups for plotting. Default is NULL.
#' @param rowsplit.order Order of the heatmap row split. Default is NULL.
#' @param rowgroup.order Order of the heatmap row split groups. Default is NULL.
#' @param facet_profile_params A list parameters for customizing the profile facet layer.
#' Details see ggh4x::facet_nested().
#' @param facet_heatmap_params A list parameters for customizing the heatmap facet layer.
#' Details see ggh4x::facet_nested().
#' @param theme_params A list parameters for customizing the theme of the plot.
#' Details see ggplot2::theme().
#' @param vline_params A list parameters for customizing the vertical line layer.
#' Details see ggplot2::geom_vline().
#' @param axisx_params A list parameters for customizing the x-axis scale.
#' Details see ggplot2::scale_x_continues().
#' @param geom_line_params A list parameters for customizing the line layer.
#' Details see ggplot2::geom_line().
#' @param rel_height Relative heights of the profile and heatmap plots when plot.type
#' is "both". Default is c(0.3,0.7).
#' @param multi_legend_barheight Height of the colorbar legend for multiple
#' colors in the heatmap. Default is 1.
#' @param multi.heatmap.col Colors for the multiple colors in the heatmap. Default is NULL.
#' @param geom_rect_params A list parameters for customizing the rectangle layer.
#' Details see ggplot2::geom_rect().
#' @param draw_rect Whether to draw a rectangle layer. Default is TRUE.
#' @param rect_ratio Ratio of the rectangle width to the maximum density range.
#' Default is 0.1.
#' @param merge_facet Whether to merge the facet layers for different binding sites.
#' Default is FALSE.
#' @param add_CI95_layer Whether to add a confidence interval (CI) 95% layer to the
#' profile plot. Default is TRUE.
#' @param geom_ribbon_params A list parameters for customizing the CI 95% layer.
#' Details see ggplot2::geom_ribbon().
#' @param facet_profile Whether to include the profile facet layer. Default is TRUE.
#' @param add_nested_line It is used when plot.type = "both" and you want to add
#' nested line for the upper profile plot. Default is FALSE.
#'
#' @return A ggplot object representing the ChIP-seq data visualization.
#'
#' @import ggplot2 ggh4x
#'
#' @export
ChipVis <- function(object = NULL,
                    plot.type = c("both","profile","heatmap"),
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
                    rel_height = c(0.3,0.7),
                    multi_legend_barheight = 1,
                    multi.heatmap.col = NULL,
                    geom_rect_params = list(),
                    draw_rect = TRUE,
                    rect_ratio = 0.1,
                    merge_facet = FALSE,
                    add_CI95_layer = TRUE,
                    geom_ribbon_params = list(),
                    facet_profile = TRUE,
                    add_nested_line = FALSE){
  plot.type <- match.arg(plot.type,c("both","profile","heatmap"))
  attrs <- attributes(object)
  # ============================================================================
  # process data
  # ============================================================================
  # get data
  df.profile <- object$profile
  df.heatmap <- object$heatmap

  # sample orders
  # sample.order = NULL
  if(!is.null(sample.order)){
    df.profile$sample <- factor(df.profile$sample,levels = sample.order)
    df.heatmap$sample <- factor(df.heatmap$sample,levels = sample.order)
  }

  # group orders
  # group.order = NULL
  # group2.order = NULL
  if(!is.null(attrs$group_sample)){
    if(!is.null(group.order)){
      df.profile$group <- factor(df.profile$group,levels = group.order)
      df.heatmap$group <- factor(df.heatmap$group,levels = group.order)
    }

    if(!is.null(group2.order)){
      df.profile$group2 <- factor(df.profile$group2,levels = group2.order)
      df.heatmap$group2 <- factor(df.heatmap$group2,levels = group2.order)
    }
  }

  # rowsplit orders
  # rowsplit.order = NULL
  if(!is.null(rowsplit.order)){
    df.profile$split <- factor(df.profile$split,levels = rowsplit.order)
    df.heatmap$split <- factor(df.heatmap$split,levels = rowsplit.order)
  }

  # rowgroup.order orders
  # rowgroup.order = NULL
  if(!is.null(rowgroup.order)){
    df.profile$split.group <- factor(df.profile$split.group,levels = rowgroup.order)
    df.heatmap$split.group <- factor(df.heatmap$split.group,levels = rowgroup.order)
  }
  # ============================================================================
  # layers
  # ============================================================================
  # profile factet layer
  if(merge_facet == TRUE){
    # for different binding sites
    if(is.null(attrs$group_sample)){
      cols = vars(sample,split)
    }else{
      cols = vars(group,sample,split)
    }
    rows = NULL
  }else{
    if(is.null(attrs$group_sample)){
      cols = vars(sample)
    }else{
      if(!is.null(attrs$group2_sample)){
        cols = vars(group2,group,sample)
      }else{
        cols = vars(group,sample)
      }
    }

    if(is.null(attrs$group_split)){
      rows = vars(split)
    }else{
      rows = vars(split.group,split)
    }
  }

  # facet_profile_params = list()
  # facet_profile = TRUE
  if(facet_profile == TRUE){
    facet_layer.profile <- do.call(facet_nested,
                                   modifyList(list(cols = cols,
                                                   solo_line = TRUE),
                                              facet_profile_params))
  }else{
    facet_layer.profile <- NULL
  }

  # heatmap factet layer
  # facet_heatmap_params = list()
  facet_layer.heatmap <- do.call(facet_nested,
                                 modifyList(list(cols = cols,
                                                 rows = rows,
                                                 solo_line = TRUE,
                                                 switch = "y",
                                                 scales = "free",space = "free_y",
                                                 independent = "y"),
                                            facet_heatmap_params))

  # ============================================
  # theme layer
  # theme_params = list()
  # vline_params = list()
  # axisx_params = list()

  plot.theme <-
    theme_bw() +
    do.call(theme,
            modifyList(list(strip.placement = "outside",
                            panel.grid = element_blank(),
                            strip.text = element_text(face = "bold",size = rel(1)),
                            legend.background = element_rect(colour = "black"),
                            axis.text = element_text(color = "black"),
                            axis.text.x = element_text(hjust = attrs$axis_textx_hjust)),
                       theme_params))

  plot.vline <-
    do.call(geom_vline,
            modifyList(list(xintercept = attrs$vline_x,lty = "dashed"),
                       vline_params))

  plot.xscale <-
    do.call(scale_x_continuous,
            modifyList(list(expand = c(0,0),
                            breaks = attrs$axis_breaks,
                            labels = attrs$axis_name),
                       axisx_params))

  # geom line layer
  # geom_line_params = list()

  # which variable for line color
  tmp.df.p <- subset(df.profile,sample == unique(df.profile$sample)[1])
  tmp.sp <- unique(tmp.df.p$split)

  if(length(tmp.sp) == 1){
    line.col.mapping <- "sample"
  }else{
    line.col.mapping <- "split"
  }

  geom.line <- do.call(geom_line,
                       modifyList(list(mapping = aes_string(x = "x",y = "density",
                                                            color = line.col.mapping)),
                                  geom_line_params))

  # rect layer
  # geom_rect_params = list()
  # draw_rect = TRUE
  if(draw_rect == TRUE){
    if("start" %in% attrs$axis_name && "end" %in% attrs$axis_name){
      st <- attrs$axis_breaks[match("start",attrs$axis_name)]
      ed <- attrs$axis_breaks[match("end",attrs$axis_name)]

      # rect_ratio = 0.1
      rect_width = max(extendrange(range(df.profile$density)))*rect_ratio
      rect.df <- data.frame(x = st,xend = ed,y = 0,yend = -rect_width)

      rect.layer <- do.call(geom_rect,
                            modifyList(list(data = rect.df,
                                            mapping = aes(xmin = st,xmax = ed,
                                                          ymin = yend,ymax = y),
                                            fill = "grey80",color = "black"),
                                       geom_rect_params))

      mult = c(0,0.05)
    }else{
      rect.layer <- NULL
      mult = c(0.05,0.05)
    }
  }else{
    rect.layer <- NULL
    mult = c(0.05,0.05)
  }

  # IC polygon layer
  # geom_ribbon_params = list()
  if(add_CI95_layer == TRUE){
    ribbon.layer <- do.call(geom_ribbon,
                            modifyList(list(mapping = aes_string(x = "x",y = "density",
                                                                 ymin = "lower_ci",ymax = "upper_ci",
                                                                 fill = line.col.mapping),
                                            show.legend = FALSE,
                                            color = NA,
                                            alpha = 0.25),
                                       geom_ribbon_params))
  }else{
    ribbon.layer <- NULL
  }

  # ============================================================================
  # plots
  # ============================================================================
  # profile
  ppfofile <-
    ggplot(df.profile) +
    ribbon.layer +
    geom.line +
    theme_bw() +
    facet_layer.profile +
    xlab("") + ylab("Average density") +
    plot.theme +
    plot.vline +
    plot.xscale +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    guides(color = guide_legend(frame.colour = "black")) +
    rect.layer +
    scale_y_continuous(expand = expansion(mult = mult))

  # heatmap
  pmain <-
    ggplot(df.heatmap) +
    facet_layer.heatmap +
    xlab("") + ylab("") +
    coord_cartesian(expand = 0) +
    plot.theme +
    plot.xscale +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  # ============================================================================
  # whether supply multiple colors for heatmap
  # multi.heatmap.col = NULL

  if(!is.null(multi.heatmap.col)){
    if(!is.null(sample.order)){
      samples <- levels(df.heatmap$sample)
    }else{
      samples <- as.character(unique(df.heatmap$sample))
    }

    names(multi.heatmap.col) <- samples

    # ==================================================
    # generate raster layer
    # ==================================================
    # s = "CBX7"
    lapply(samples, function(s){
      tmp <- subset(df.heatmap,sample == s)

      # re-assign aes var
      aes.res <- aes(x = x,y = y,var = value)
      names(aes.res)[3] <- s

      # re-filter levels for y axis
      sp <- if(is.null(levels(tmp$split))){
        unique(tmp$split)
      }else{
        levels(tmp$split)
      }

      # x = 3
      lapply(seq_along(sp), function(x){
        tmpp <- tmp[tmp$split %in% sp[x],]
        lv.prefix <- sapply(strsplit(tmpp$split[1],split = "\n"), "[",1)
        lvs <- data.frame(levels = levels(tmpp$y))
        lvs <- lvs[startsWith(lvs$levels,lv.prefix),]
        tmpp$y <- factor(tmpp$y,levels = lvs)

        # raster layer
        return(geom_raster(data = tmpp,mapping = aes.res))
      }) |> unlist() -> raster.layer.list

      return(raster.layer.list)
    }) |> unlist() -> raster.layer

    # ==================================================
    # re-filter levels for y axis
    # ==================================================
    sp <- if(is.null(levels(df.heatmap$split))){
      unique(df.heatmap$split)
    }else{
      levels(df.heatmap$split)
    }

    # check whether draw different peaks for each sample
    if(merge_facet == FALSE){
      samp <- rep(samples,length(sp))
      spt <- rep(sp,each = length(samples))
    }else{
      sample_split_info <- unique(df.heatmap[,c("sample","split")])
      samp <- sample_split_info$sample
      spt <- sample_split_info$split
    }


    # generate y axis scales layer
    # s = 7
    lapply(seq_along(spt), function(s){
      tmp <- subset(df.heatmap,sample == samp[s])
      tmpp <- tmp[tmp$split %in% spt[s],]

      lv.prefix <- sapply(strsplit(tmpp$split[1],split = "\n"), "[",1)
      lvs <- data.frame(levels = levels(tmpp$y))
      lvs <- lvs[startsWith(lvs$levels,lv.prefix),]

      return(scale_y_discrete(limits = lvs))
    }) -> yscales.layer

    # colorbar layer
    lapply(seq_along(samples), function(s){
      colourbar <- guide_colourbar(barheight = unit(multi_legend_barheight, "cm"),
                                   frame.colour = "black",
                                   ticks.colour = "black",
                                   title = samples[s],
                                   order = s)
    }) -> colorbar.layer

    pheatmap <-
      pmain +
      raster.layer +
      scale_fill_multi(aesthetics = samples,
                       # Providing colours as a list distributes list-elements over different scales
                       colors = multi.heatmap.col,
                       guide = colorbar.layer) +
      plot.vline +
      facetted_pos_scales(y = yscales.layer)

  }else{
    pheatmap <-
      pmain +
      geom_raster(aes(x = x,y = y,fill = value)) +
      plot.vline +
      scale_fill_gradient(low = "grey95",high = "#1450A3",name = "signal") +
      guides(fill = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }

  # ============================================================================
  # which type to plot
  if(plot.type == "profile"){
    p <- ppfofile
  }else if(plot.type == "heatmap"){
    p <- pheatmap
  }else{
    ppfofile.new <- ppfofile +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(b = 0,unit = "cm"))

    # whether add nest line
    if(add_nested_line == TRUE){
      ppfofile.new <- ppfofile.new +
        theme(strip.background = element_blank(),
              ggh4x.facet.nestline = element_line(colour = "black",linewidth = 1))

      if(!is.null(attrs$group_split)){
        pheatmap.new <- pheatmap +
          theme(strip.background = element_blank(),
                ggh4x.facet.nestline = element_line(colour = "black",linewidth = 1))
      }else{
        pheatmap.new <- pheatmap +
          theme(strip.text.x = element_blank(),
                strip.background.x = element_blank())
      }
    }

    # combine
    # rel_height = c(0.3,0.7)
    p <- aplot::plot_list(gglist = list(ppfofile.new,pheatmap.new),
                          ncol = 1,heights = rel_height)
  }

  return(p)
}

