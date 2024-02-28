#' # HarmonizR
#'
#' #' Helper function to prepare input data for the HarmonizR function
#' #'
#' #' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' #' @param ain String which assay should be used as input (default raw)
#' #' @param batch_column  column name of batch (if NULL, batch saved in SummarizedExperiment will be taken)
#' #'
#' #' @return list of two data frames
#' #' @export
#' #'
#' prepare_data_for_harmonizR <- function(se, ain = "raw", batch_column = NULL){
#'   # remove reference samples
#'   se <- remove_reference_samples(se)
#'   # prepare intensity file
#'   dt <- as.data.frame(SummarizedExperiment::assays(se)[[ain]])
#'   row.names(dt) <- SummarizedExperiment::rowData(se)$Protein.IDs
#'   if(ain != "raw"){
#'     dt <- 2^dt
#'   }
#'
#'   # prepare description
#'   md <- data.table::as.data.table(SummarizedExperiment::colData(se))
#'   md <- md[, c("Column", "Batch")]
#'   colnames(md) <- c("ID", "batch")
#'   md$sample <- seq(1, nrow(md))
#'   md <- md[, c("ID", "sample", "batch")]
#'   # check if batch column numeric
#'   batch_list <- md[["batch"]]
#'   if(length(is.na(as.numeric(batch_list)))>0){
#'     # convert batch column to numerics
#'     batch_list <- as.numeric(factor(batch_list))
#'     md$batch <- batch_list
#'   }
#'
#'   return(list(as.data.frame(dt), as.data.frame(md)))
#' }
#'
#' #' HarmonizR proteomics data
#' #'
#' #' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' #' @param ain String which assay should be used as input (default raw)
#' #' @param aout String which assay should be used to save normalized data (default HarmonizR_ComBat)
#' #' @param batch_column  column name of batch (if NULL, batch saved in SummarizedExperiment will be taken)
#' #' @param algorithm String specifying which algorithm should be used (ComBat or limma)
#' #' @param combat_mode Integer speciying which ComBat mode should be executed (1,2,3,4)
#' #'
#' #' @return SummarizedExperiment containing the HarmonizR data as assay
#' #' @export
#' #'
#' #' @details
#' #' combat_mode = 1 -> par.prior = TRUE, mean.only = FALSE
#' #' combat_mode = 2 -> par.prior = TRUE, mean.only = TRUE
#' #' combat_mode = 3 -> par.prior = FALSE, mean.only = FALSE
#' #' combat_mode = 4 -> par.prior = FALSE, mean.only = TRUE
#' #' par.prior -> TRUE indicates parametric adjustments will be used
#' #' mean.only -> TRUE indicates that ComBat only corrects the mean of the batch effect (no scale adjustment)
#' run_harmonizR <- function(se, ain = "raw", aout = "HarmonizR_ComBat", batch_column = NULL, algorithm = "ComBat", combat_mode = 1){
#'   # Check input
#'   stopifnot(algorithm %in% c("ComBat", "limma"))
#'   stopifnot(combat_mode %in% c(1,2,3,4))
#'   if(is.null(batch_column)){
#'     if(is.null(S4Vectors::metadata(se)$batch)){
#'       # print error
#'       stop("No batch column provided!")
#'     } else {
#'       batch_column <- S4Vectors::metadata(se)$batch
#'       message("Batch column of SummarizedExperiment used!")
#'     }
#'   } else {
#'     # check if batch column in SummarizedExperiment
#'     if(! batch_column %in% colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))){
#'       stop(paste0("No column named ", batch_column, " in SummarizedExperiment!"))
#'     }
#'   }
#'
#'   # Prepare input data
#'   data <- prepare_data_for_harmonizR(se, ain = ain, batch_column)
#'   # Run HarmonizR
#'   harmonized_data <- HarmonizR::harmonizR(data_as_input = data[[1]], description_as_input = data[[2]], algorithm = algorithm, ComBat_mode = 1)
#'   # Reconstruct Original Order
#'   harmonized_dt <- harmonized_data[row.names(data[[1]]),]
#'
#'   # Add data to SummarizedExperiment
#'   harmonized_dt <- log2(harmonized_dt)
#'   SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(harmonized_dt)
#'   return(se)
#' }
#'
