

#' Function Returning Some Values on the Numbers of NA in the data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
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
  return(data.table::data.table("Total values"=c(tot_nr), "NA values" = c(na_nr), "NA Percentage"= c(na_perc)))
}
