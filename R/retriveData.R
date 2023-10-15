#' Retrieve Data
#'
#' This function retrieves profile and heatmap data (EnrichedHeatmap::normalizeToMatrix
#' and ChIPseeker::getTagMatrix) from a list of matrices. Some codes origin from R package
#' EnrichedHeatmap/ChIPseeker with some modifications.
#'
#' @param mat.list A list of matrices containing the data which output from
#' normalizeToMatrix/getTagMatrix.
#' @param row.split A character vector specifying the group of the peaks.
#' @param sample.names A character vector specifying the names of the samples.
#' @param rm.extreme.value Logical value indicating whether to remove extreme binding
#' signal values. Default FALSE.
#' @param group.sample A character vector specifying the group of each sample.
#' @param group2.sample A character vector specifying the second group of each sample.
#' @param aggregate.fun A character vector specifying the aggregate function(mean/meidan) to
#' use for profile data. Default "mean".
#' @param ChIPseeker Logical value indicating whether ChIPseeker data is used.
#'
#' @return A list containing the profile and heatmap data.
#' @export
retriveData <- function(mat.list = NULL,
                        row.split = NULL,
                        sample.names = NULL,
                        rm.extreme.value = FALSE,
                        group.sample = NULL,
                        group2.sample = NULL,
                        aggregate.fun = c("mean","median"),
                        ChIPseeker = FALSE){
  myfun <- match.arg(aggregate.fun,c("mean","median"))
  # ============================================================================
  # loop extract profile data
  # ============================================================================
  # x = 1
  purrr::map_df(seq_along(mat.list),function(x){
    norm_mat <- mat.list[[x]]
    # extend = attr(norm_mat, "extend")

    # get extend range
    # ChIPseeker = TRUE
    if(ChIPseeker == TRUE){
      extend = c(attr(norm_mat, "upstream"),attr(norm_mat, "downstream"))
    }else{
      extend = attr(norm_mat, "extend")
    }

    # row split for binding sites
    if(is.null(row.split)){
      binding_sites <- nrow(mat.list[[x]])
      row_split <- rep(paste0("Binding sites: ",binding_sites),binding_sites)
    }else{
      row_split <- row.split
    }
    # =====================
    # extract profile data
    # group for rows
    # x = "A"
    purrr::map_df(sort(unique(row_split)),function(split){
      idx <- base::grep(split,row_split)

      fun.apply <- match.fun(myfun)
      mean_mat <- apply(norm_mat[idx,], 2, fun.apply)
      ci_95_lw <- apply(norm_mat[idx,], 2, function(x){Rmisc::CI(x)[3]})
      ci_95_hg <- apply(norm_mat[idx,], 2, function(x){Rmisc::CI(x)[1]})

      # profile to df
      mat_df <- data.frame(density = mean_mat,
                           # x = seq(-extend[1],extend[2],length.out = length(mean_mat)),
                           x = 1:ncol(norm_mat),
                           sample = sample.names[x],
                           split = split,
                           lower_ci = ci_95_lw,
                           upper_ci = ci_95_hg)
      return(mat_df)
    }) -> df.profile

    # add group
    if(!is.null(group.sample)){
      df.profile$group <- group.sample[x]
      if(!is.null(group2.sample)){
        df.profile$group2 <- group2.sample[x]
      }
    }

    return(df.profile)
  }) -> df.profile

  # ============================================================================
  # loop extract heatmap data
  # ============================================================================
  # x = 1
  # split = "group3"
  purrr::map_df(seq_along(mat.list),function(x){
    norm_mat <- mat.list[[x]]
    rownames(norm_mat) <- 1:nrow(norm_mat)
    # extend = attr(norm_mat, "extend")

    # get extend range
    # ChIPseeker = TRUE
    if(ChIPseeker == TRUE){
      extend = c(attr(norm_mat, "upstream"),attr(norm_mat, "downstream"))
      norm_mat <- t(apply(norm_mat, 1, function(x) x/max(x)))

    }else{
      extend = attr(norm_mat, "extend")
    }

    # row split for binding sites
    if(is.null(row.split)){
      binding_sites <- nrow(mat.list[[x]])
      row_split <- rep(paste0("Binding sites: ",binding_sites),binding_sites)
    }else{
      row_split <- row.split
    }

    # =====================
    # extract heatmap data
    # split = "group3"
    purrr::map_df(sort(unique(row_split)),function(split){
      idx <- base::grep(split,row_split)

      tmp_mat <- data.frame(norm_mat[idx,])
      # colnames(tmp_mat) <- seq(-extend[1],extend[2],length.out = ncol(tmp_mat))
      colnames(tmp_mat) <- 1:ncol(tmp_mat)
      rownames(tmp_mat) <- 1:nrow(tmp_mat)
      tmp_mat$y <- paste0(split,"_",rownames(tmp_mat))

      df.long <- reshape2::melt(tmp_mat,id.vars = "y")
      colnames(df.long)[2] <- c("x")
      df.long$split <- paste0(split,"\n(n=",length(idx),")")
      df.long$sample <- sample.names[x]
      df.long$x <- as.numeric(unfactor(df.long$x))

      # binding site order
      if(ChIPseeker == TRUE){
        ii <- order(rowSums(data.frame(norm_mat[idx,])))
        tmp_mat <- tmp_mat[ii,]

        df.long$y <- factor(df.long$y,levels = paste0(split,"_",ii))
      }else{
        row_order <- order(EnrichedHeatmap::enriched_score(norm_mat[idx,]),
                           decreasing = TRUE)

        df.long$y <- factor(df.long$y,levels = paste0(split,"_",rev(row_order)))
      }

      return(df.long)
    }) -> df.heatmap

    # add group
    if(!is.null(group.sample)){
      df.heatmap$group <- group.sample[x]
      if(!is.null(group2.sample)){
        df.heatmap$group2 <- group2.sample[x]
      }
    }

    return(df.heatmap)
  }) -> df.heatmap

  # deal with extreme values
  if(rm.extreme.value == TRUE){
    q99 <- quantile(df.heatmap$value,c(0,0.99))
    df.heatmap.new <- df.heatmap |>
      mutate(value = ifelse(value > q99[2],q99[2],value))
  }else{
    df.heatmap.new <- df.heatmap
  }

  # ============================================================================
  # get attributes
  # ============================================================================
  class.mat <- attr(mat.list[[1]],"class")

  if("normalizedMatrix" %in% class.mat){
    upstream_index = attr(mat.list[[1]], "upstream_index")
    downstream_index = attr(mat.list[[1]], "downstream_index")
    target_index = attr(mat.list[[1]], "target_index")

    extend = attr(mat.list[[1]], "extend")

    n1 = length(upstream_index)
    n2 = length(target_index)
    n3 = length(downstream_index)
    n = n1 + n2 + n3

    upstream_flipped = attr(mat.list[[1]], "upstream_flipped")
    if(is.null(upstream_flipped)) upstream_flipped = FALSE

    # xtcks <- seq(-extend[1],extend[2],length.out = ncol(mat.list[[1]]))

    # axis labels
    if(n1 && n2 && n3) {
      axis_name = c(paste0("-", extend[1]), "start", "end", extend[2])
      breaks = c(1,n1 + 1,n1 + n2 + 1,n)
      vline.x = c(n1 + 1,n1 + n2 + 1)
      text.x.hjust = c(0,0.5,0.5,1)
    } else if(n1 && !n2 && n3) {
      axis_name = c(paste0("-", extend[1]), "start", extend[2])
      breaks = c(1,n1 + 1,n)
      vline.x = c(n1 + 1)
      text.x.hjust = c(0,0.5,0.5,1)
    } else if(!n1 && n2 && n3) {
      axis_name = c("start", "end", extend[2])
      breaks = c(n1 + 1,n1 + n2 + 1,n)
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
        breaks = c(n1 + 1,n)
        vline.x = c()
        text.x.hjust = c(0,1)
      } else {
        axis_name = c("end", extend[2])
        breaks = c(n1 + n2 + 1,n)
        vline.x = c()
        text.x.hjust = c(0,1)
      }
    }
  }else{
    upstream <- attr(mat.list[[1]], "upstream")
    downstream <- attr(mat.list[[1]], "downstream")
    binning_Flag <- attr(mat.list[[1]],"is.binning")
    type <- attr(mat.list[[1]],"type")

    extend = c(upstream,downstream)
    # xtcks <- seq(-extend[1],extend[2],length.out = ncol(mat.list[[1]]))

    # if(type == "body"){
    #   axis_name <- grid.pretty(range = range(xtcks),n = )
    #   breaks <- axis_name
    #   vline.x = c()
    #   text.x.hjust = c(0,1)
    # }
    # attributes(mat.list[[1]])

    body_Flag <- FALSE
    if(type == "body"){
      body_Flag <- TRUE
      label <- attr(mat.list[[1]],"label")
    }

    if(binning_Flag){
      nbin <- dim(mat.list[[1]])[2]
    }else{
      nbin <- ncol(mat.list[[1]])
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

  # ============================================================================
  # return final results
  # ============================================================================
  # return
  res <- list(profile = df.profile,
              heatmap = df.heatmap.new)

  class(res) <- "chipData"

  res <- structure(res,
                   group_sample = group.sample,
                   group2_sample = group2.sample,
                   axis_name = axis_name,
                   axis_breaks = breaks,
                   vline_x = vline.x,
                   axis_textx_hjust = text.x.hjust)

  return(res)
}
