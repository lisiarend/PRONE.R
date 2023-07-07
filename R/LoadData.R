
#' Load proteomics data into a SummarizedExperiment
#'
#' @param data tabular data table with rows = proteins and columns = samples (such as proteinGroups.txt of MaxQuant)
#' @param md experimental design table (requires a column named "Column" for the column names of the sample intensities in data)
#' @param protein_column name of the column in data containing the protein IDs
#' @param gene_column name of the column in data containing the gene names
#' @param ref_samples reference samples if TMT experiment provided (names of samples)
#' @param batch_column name of the column in md defining the batches
#' @param condition_column name of the column in md defining the condition (can still be changed afterwards)
#' @param label_column name of the column in md containing simple sample names (for visualization)
#' @importFrom magrittr %>%
#'
#' @return SummarizedExperiment object
#' @export
#'
load_data <- function(data, md, protein_column = "Protein.IDs", gene_column = "Gene.Names", ref_samples = NULL, batch_column = NULL, condition_column = NULL, label_column = NULL){
  # convert to data.table
  data <- data.table::as.data.table(data)
  md <- data.table::as.data.table(md)
  # assays
  raw <- data[,md$Column, with=FALSE]
  raw[raw == 0] <- NA
  log2 <- log2(raw)
  # rowData
  rowData <- data[, !colnames(data) %in% md$Column, with = FALSE]
  rowData$IDs <- rownames(raw)
  rowData <- rowData %>% dplyr::rename("Protein.IDs" = protein_column, "Gene.Names" = gene_column)
  # colData
  colData <- md
  # metadata
  metadata <- list(condition = condition_column, batch = batch_column, refs = ref_samples, label = label_column)
  # SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = raw, log2 = log2), colData = colData, rowData = rowData, metadata = metadata)
  return(se)
}
