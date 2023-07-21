

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



