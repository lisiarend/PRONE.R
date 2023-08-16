

#' Function returning some values on the numbers of NA in the data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#'
#' @return list with total amount of values in the data, amount of NA values, and the percentage of NAs
#' @export
#'
get_NA_overview <- function(se, ain="log2"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  na_nr <- sum(is.na(dt))
  tot_nr <- dim(dt)[1] * dim(dt)[2]
  na_perc <- na_nr/tot_nr * 100
  return(data.table::data.table("Total.Values"=c(tot_nr), "NA.Values" = c(na_nr), "NA.Percentage"= c(na_perc)))
}


#' Barplot showing the number of samples per condition
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#'
#' @return ggplot object
#' @export
#'
plot_condition_overview <- function(se, condition = NULL){
  # get condition
  condition <- get_condition_value(se, condition)

  # prepare data
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))
  dt <- md %>% dplyr::count(get(condition))
  colnames(dt) <- c(condition, "n")

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual",]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector[54:73])

  # plot
  p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(condition), y=get("n"), fill=get(condition))) +
    ggplot2::geom_col() +
    ggplot2::labs(x=condition, y="Number of Samples") +
    ggplot2::scale_fill_manual(name=condition, values=col_vector) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=0.5, vjust=0.5))
  return(p)
}


#' Boxplots of intensities of specific markers
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param markers vector of the IDs of the markers to plot
#' @param id_column String specifying the column of the rowData of the SummarizedExperiment object which includes the IDs of the markers
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param shape_by String specifying the column to shape the samples (If NULL or "No", no shaping of samples is done.)
#' @param facet_norm Boolean indicating whether to facet by normalization method (TRUE) or not (FALSE)
#' @param facet_marker Boolean indicating whether to facet by comparison (TRUE) or not (FALSE). Only valid if facet_norm = FALSE.
#'
#' @return ggplot object
#' @export
#'
plot_markers_boxplots <- function(se, ain, markers, id_column = "Protein.IDs", color_by = NULL, shape_by = NULL, facet_norm = TRUE, facet_marker = FALSE){
  # check input parameters
  color_by <- get_color_value(se, color_by)
  shape_by <- get_shape_value(se, shape_by)
  ain <- check_input_assays(se, ain)

  # check id_column
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(!id_column %in% colnames(rd)){
    # if id_column not in rowData
    stop(paste0(id_column, " not in rowData of the SummarizedExperiment object!"))
  }

  # prepare data
  overall_dt <- NULL
  for(assay in c(ain)){
    dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[assay]])
    dt <- cbind(data.table::as.data.table(SummarizedExperiment::rowData(se))[[id_column]], dt)
    melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = id_column)
    melted_dt$Method <- assay
    if(is.null(overall_dt)){
      overall_dt <- melted_dt
    } else {
      overall_dt <- rbind(overall_dt, melted_dt)
    }
  }

  # merge colData
  cd <- data.table::as.data.table(SummarizedExperiment::colData(se))
  overall_dt <- merge(overall_dt, cd, by = "Column")
  # subset intensities by markers
  found_markers <- stringr::str_detect(overall_dt[[id_column]], paste(markers, collapse = "|") )
  overall_dt <- overall_dt[found_markers,]

  # color vector
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector[54:73])

  # facet by normalization method
  plots <- list()
  if(facet_norm){
    for(marker in markers){
      found_markers <- stringr::str_detect(overall_dt[[id_column]], marker )
      tmp <- overall_dt[found_markers,]
      p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), color = get(color_by))) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(x = color_by, color = color_by, shape = shape_by, y = "Intensity") +
        ggplot2::scale_color_manual(values = col_vector) +
        ggplot2::facet_wrap(~Method)
      if(is.null(shape_by)){
        p <- p + ggplot2::geom_point(alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      }
      plots[[marker]] <- p
    }
  # facet by marker
  } else if(facet_marker){
    for(method in ain){
      tmp <- overall_dt[overall_dt$Method == method,]
      p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), color = get(color_by))) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(x = color_by, color = color_by, shape = shape_by, y = "Intensity") +
        ggplot2::scale_color_manual(values = col_vector) +
        ggplot2::facet_wrap(~Protein.IDs)
      if(is.null(shape_by)){
        p <- p + ggplot2::geom_point(alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      }
      plots[[method]] <- p
    }
  # do not facet by anything
  } else {
    for(method in ain){
      for(marker in markers){
        tmp <- overall_dt[overall_dt$Method == method,]
        found_markers <- stringr::str_detect(tmp[[id_column]], marker )
        tmp <- tmp[found_markers,]
        p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), color = get(color_by))) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(x = color_by, color = color_by, shape = shape_by, y = "Intensity") +
          ggplot2::scale_color_manual(values = col_vector)
        if(is.null(shape_by)){
          p <- p + ggplot2::geom_point(alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
        } else {
          p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
        }
        plots[[paste0(method, "_", marker)]] <- p
      }
    }
  }
  return(plots)
}


