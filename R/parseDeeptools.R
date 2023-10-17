#' Parse deeptools output
#'
#' @author Jun Zhang
#'
#' @description
#' This function parses the deeptools output and returns a list of matrices.
#' Each matrix represents the data for a specific sample and contains additional
#' attributes like upstream, downstream, extend, and target indices.
#'
#' @param deeptools_output Path to the deeptools output file.
#' @return A list of matrices with additional attributes.
#'
#' @details
#' The output of deeptools in gzipped text format is read using `data.table::fread()`.
#' The function extracts meta information from the first column of the data and
#' processes it to get relevant values for upstream, downstream, body, bin size,
#' group labels, group boundaries, sample labels, and sample boundaries.
#' Then, it loops through the sample labels and creates a matrix for each sample,
#' filling any missing values with zeros. Lastly, it assigns additional attributes
#' to each matrix based on the extracted values and returns a list of matrices.
#'
#' @import data.table
#' @importFrom stringr str_extract_all
#'
#' @export
parseDeeptools <- function(deeptools_output = NULL){
  # ============================================================================
  # # prepare meta info
  # ============================================================================
  gz.data <- suppressWarnings(data.table::fread(deeptools_output))

  meta <- colnames(gz.data)[1]
  tmp <- sapply(strsplit(meta,split = '\\{|\\}'), "[",2)

  labels <- c('upstream','downstream','body','bin size',
              'group_labels','group_boundaries','sample_labels','sample_boundaries')

  pattern <- '\\[(.*?)\\]'
  match.val <- str_extract_all(tmp, pattern)[[1]]
  vals <- sapply(strsplit(match.val,split = "\\[|\\]"), "[",2)[c(1,2,3,4,8,9,10,11)]

  upstream <- as.numeric(unlist(strsplit(vals[1],split = ",")))
  downstream <- as.numeric(unlist(strsplit(vals[2],split = ",")))
  body <- as.numeric(unlist(strsplit(vals[3],split = ",")))
  binsize <- as.numeric(unlist(strsplit(vals[4],split = ",")))
  group_labels <- gsub('"',"",unlist(strsplit(vals[5],split = ",")),fixed = T)
  group_boundaries <- as.numeric(unlist(strsplit(vals[6],split = ",")))[-1]
  sample_labels <- gsub('"',"",unlist(strsplit(vals[7],split = ",")),fixed = T)
  sample_boundaries <- as.numeric(unlist(strsplit(vals[8],split = ",")))

  # ============================================================================
  # loop prepare matrix
  # ============================================================================

  # x = 1
  mat.list <- lapply(1:(length(sample_labels)), function(x){
    tmp <- data.frame(gz.data)
    mat <- tmp[,c((7 + sample_boundaries[x]):(6 + sample_boundaries[x + 1]))]

    # fill NA
    mat[is.na(mat)] <- 0

    # check data from "reference-point" or "scale-regions"
    if(body[x] == 0){
      target_index <- integer(0)
    }else{
      target_index <- (upstream[x]/binsize[x] + 1):(upstream[x]/binsize[x] + body[x]/binsize[x])
    }

    # peaks datasets numbers
    if(length(group_boundaries) == 1){
      group_numbers <- group_boundaries
    }else{
      group_numbers <- c(group_boundaries[1],Reduce(f = "-",x = rev(group_boundaries)))
    }

    # assign attributes for matrix
    res <- structure(as.matrix(mat),
                     dim = c(nrow(mat),ncol(mat)),
                     upstream = upstream[x],
                     downstream = downstream[x],
                     extend = c(upstream[x],downstream[x]),
                     upstream_index = 1:(upstream[x]/binsize[x]),
                     target_index = target_index,
                     downstream_index = (upstream[x]/binsize[x] + body[x]/binsize[x] + 1):(ncol(mat)),
                     sample_name = sample_labels[x],
                     group_labels = group_labels,
                     group_numbers = group_numbers)

    # assign class
    class(res) <- c("deeptoolsMat","matrix")

    return(res)
  })

  return(mat.list)
}
