## ----- Normalization Methods ----- ##

#' Total Intensity Normalization
#' Intensities of each variable in a sample are divided with the sum of intensities
#' of all variables in the sample and multiplied with the median or mean of sum of intensities
#' of all variables in all samples. Raw data is taken as input (not log2).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default TotalInt_mean or TotalInt_median)
#' @param type String whether to use median or mean to calculate the scaling factor
#' @importFrom SummarizedExperiment assay
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay
#' @export
#'
globalIntNorm <- function(se, ain = "raw", aout="GlobalMedian", type = "median"){
  dt <- SummarizedExperiment::assays(se)[[ain]]
  dt <- as.data.frame(dt)
  if(ain != "raw"){
    dt <- 2^dt
  }
  colSums <- colSums(dt, na.rm = TRUE)
  if(type == "median"){
    colSumsM <- stats::median(colSums)
  } else if (type == "mean"){
    colSumsM <- mean(colSums)
  }
  norm_dt <- matrix(nrow = nrow(dt), ncol = ncol(dt),
                    byrow = TRUE)
  normFunc <- function(colIndex) {
    (dt[rowIndex, colIndex]/colSums[colIndex]) *
      colSumsM
  }
  for (rowIndex in seq_len(nrow(dt))) {
    norm_dt[rowIndex, ] <- vapply(seq_len(ncol(dt)),
                                  normFunc, 0)
  }
  norm_dt <- log2(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Total Intensity Normalization Using Mean For Calculation of Scaling Factors
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default GlobalMean)
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay
#' @export
#'
globalMeanNorm <- function(se, ain = "raw", aout = "GlobalMean"){
  se <- globalIntNorm(se, ain, aout, type = "mean")
  return(se)
}

#' Total Intensity Normalization Using Median For Calculation of Scaling Factors
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default GlobalMedian)
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay
#' @export
#'
globalMedianNorm <- function(se, ain = "raw", aout = "GlobalMedian"){
  se <- globalIntNorm(se, ain, aout, type = "median")
  return(se)
}

#' Median Normalization
#' The intensity of each protein group in a given sample is divided by the median of the
#' intensities of all protein groups in that sample and then multiplied by the mean of
#' median of sum of intensities of all protein groups in all samples.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default mean)
#'
#' @return SummarizedExperiment containing the median normalized data as assay
#' @export
#'
medianNorm <- function(se, ain = "log2", aout = "Median"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # find median of each sample
  sample_med <- apply(dt, 2, stats::median, na.rm=TRUE) # columns
  # find mean of medians
  mean_med <- mean(sample_med, na.rm=TRUE)
  # divide data by median
  norm_dt <- t(t(dt)/sample_med)
  # multiply data by mean of medians
  norm_dt <- norm_dt * mean_med
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Mean Normalization
#' The intensity of each protein group in a given sample is divided by the mean of the
#' intensities of all protein groups in that sample and then multiplied by the mean of
#' mean of sum of intensities of all protein groups in all samples.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default mean)
#'
#' @return SummarizedExperiment containing the mean normalized data as assay
#' @export
#'
meanNorm <- function(se, ain = "log2", aout="Mean"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # find means of each sample
  sample_mean <- apply(dt, 2, mean, na.rm=TRUE)
  # find mean of means of each sample
  mean_mean <- mean(sample_mean, na.rm=TRUE)
  # divide data by mean
  norm_dt <- t(t(dt)/sample_mean)
  # multiply data by mean of means
  norm_dt <- norm_dt * mean_mean
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Internal Reference Scaling Normalization
#' IRS makes different measurements of the same thing all exactly the same and puts
#' all of the intensities on the same scale.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default IRS)
#'
#' @return SummarizedExperiment containing the IRS normalized data as assay
#' @export
#'
irsNorm <- function(se, ain="raw", aout="IRS"){
  # extract necessary info of SummarizedExperiment
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  if(ain != "raw"){
    dt <- 2^dt
  }

  batch <- S4Vectors::metadata(se)$batch
  refs <- S4Vectors::metadata(se)$refs
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))

  # get md of reference samples
  refs_md <- md[md$Column %in% refs,]
  # separate data by batch
  dt_list <- lapply(unique(md[,batch]), function(b){
    md_chunk <- md[md[,batch] == b,]
    dt_chunk <- dt[, md_chunk$Column, with=FALSE]
    return(dt_chunk)
  })
  names(dt_list) <- unique(md[,batch])
  # check if ref_samples have been specified
  if (is.null(refs)){
    # create mock channel with rowsums from each frame
    irs <- lapply(dt_list, function(dt_chunk){
      return(rowSums(dt_chunk, na.rm=TRUE))
    })
    irs <- do.call(cbind, irs)
  } else {
    # take reference sample intensities
    irs <- dt[, refs_md$Column, with=FALSE]
    colnames(irs) <- as.character(refs_md[refs_md$Column %in% refs,][,batch])
  }
  # get the geometric average intensity for each protein
  irs <- tibble::as_tibble(irs)
  irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
  # normalize data
  dt_irs_list <- lapply(names(dt_list), function(b){
    # compute scaling factor vectors
    fac <- irs$average / irs[,b]
    # normalize
    dt_irs_chunk <- dt_list[[b]] * fac[,1]
    return(dt_irs_chunk)
  })
  # reconstruct data after irs normalization
  dt_irs <- do.call(cbind, dt_irs_list)
  dt_irs <- data.table::as.data.table(dt_irs)
  dt_irs <- dt_irs[, colnames(dt), with = FALSE]
  SummarizedExperiment::assay(se, aout, FALSE) <- log2(dt_irs)
  return(se)
}

#' Quantile Normalization of preprocessCore package.
#' Forces distributions of the samples to be the same on the basis of the quantiles of the samples by replacing
#' each protein of a sample with the mean of the corresponding quantile.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default quantile)
#'
#' @return SummarizedExperiment containing the quantile normalized data as assay
#' @export
#'
quantileNorm <- function(se, ain="log2", aout="Quantile"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  norm_dt <- preprocessCore::normalize.quantiles(as.matrix(dt), copy=TRUE)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Variance Stabilization Normalization of limma package
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default vsn)
#'
#' @return SummarizedExperiment containing the vsn normalized data as assay
#' @export
#'
vsnNorm <- function(se, ain="raw", aout="VSN"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  norm_dt <- suppressMessages(limma::normalizeVSN(dt))
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Weighted Trimmed Mean of M Values (TMM) Normalization
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default raw)
#' @param aout String which assay should be used to save normalized data (default TMM)
#'
#' @return SummarizedExperiment containing the TMM normalized data as assay
#' @export
#'
tmmNorm <- function(se, ain="raw", aout="TMM"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  if(ain != "raw"){
    dt <- 2^dt
  }
  tmm <- edgeR::calcNormFactors(stats::na.omit(dt))
  dt_norm <- sweep(dt, 2, tmm, FUN="/")
  dt_norm <- data.table::as.data.table(dt_norm)
  SummarizedExperiment::assay(se, aout, FALSE) <- log2(dt_norm)
  return(se)
}

#' Robust Linear Regression Normalization of NormalyzerDE (uses median values over all samples as reference
#' sample to which all the other samples in the data are normalized to)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default rlr)
#'
#' @return SummarizedExperiment containing the rlr normalized data as assay
#' @export
#'
rlrNorm <- function(se, ain="log2", aout="Rlr"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt <- as.data.frame(dt)
  sample_median <- matrixStats::rowMedians(as.matrix(dt), na.rm = TRUE, useNames = TRUE)
  for (i in 1:ncol(dt)){ # iterate over samples
    sampleA <- dt[,i] # sample to normalize
    lrFit <- MASS::rlm(as.matrix(sampleA) ~
                         sample_median, na.action = stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (sampleA - coefIntercept)/coefSlope
    dt[,i] <- globalFittedRLRCol
  }
  colnames(dt) <- colnames(dt)
  rownames(dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(dt)
  return(se)
}

#' Linear Regression Normalization on MA Transformed Data (similar to Rlr, but data are MA transformed before normalization,
#' (A = median sample, M = difference of that sample to A)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default rlrMA)
#'
#' @return SummarizedExperiment containing the rlrMA normalized data as assay
#'
#' @return SummarizedExperiment containing the rlrMA normalized data as assay
#' @export
#'
rlrMANorm <- function(se, ain="log2",aout="RlrMA"){
  # extract intensities
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt_matrix <- as.matrix(dt)
  # MA transformation
  ref_a <- MatrixGenerics::rowMedians(dt_matrix, na.rm=TRUE, useNames = TRUE) # A = median over all samples --> vector of length (number of proteins)
  for (i in 1:ncol(dt)){ # iterate over samples
    comp_sample <- dt_matrix[, i] # sample to normalize
    m <- comp_sample - ref_a  # M = sample to normalize
    # rlm fit
    fit <- MASS::rlm(m ~ ref_a, na.action = stats::na.exclude)
    fit_values <- stats::predict(fit)
    dt_matrix[, i] <- comp_sample - fit_values
  }
  norm_dt <- data.table::as.data.table(dt_matrix)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Cyclic Linear Regression Normalization on MA Transformed Data (no reference, but MA transformation and normalization
#' of samples done pairwise between two samples, A = average of two samples, M = difference, process iterated through all samples pairs).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default rlrMACyc)
#' @param iterations Number of cyclic iterations to be performed (default 3)
#'
#' @return SummarizedExperiment containing the rlrMACyc normalized data as assay
#'
#' @return SummarizedExperiment containing the rlrMACyc normalized data as assay
#' @export
#'
rlrMACycNorm <- function(se, ain="log2", aout="RlrMACyc", iterations=3){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt_matrix <- as.matrix(dt)
  n <- ncol(dt)
  # iteration over all pairs
  for (k in 1:iterations){
    for (j in 1:n){
      for (i in j:n){
        # MA transformation
        ref_sample <- dt_matrix[, j]
        comp_sample <- dt_matrix[, i]
        m <- comp_sample - ref_sample
        a <- (ref_sample + comp_sample)/2
        # rlm fit
        fit <- MASS::rlm(m ~ a, na.action = stats::na.exclude)
        fit_values <- stats::predict(fit)
        # reorder values
        comp_norm <- comp_sample - fit_values/2
        ref_norm <- ref_sample + fit_values/2
        dt_matrix[, j] <- ref_norm
        dt_matrix[, i] <- comp_norm
      }
    }
  }
  norm_dt <- data.table::as.data.table(dt_matrix)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Cyclic Loess Normalization of limma (two samples of the data are MA transformed and normalized at a time, and all pairs of samples are iterated through)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default cyclicLoess)
#'
#' @return SummarizedExperiment containing the loessCyc normalized data as assay
#' @export
#'
loessCycNorm <- function(se, ain="log2", aout="LoessCyc"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  norm_dt <- limma::normalizeCyclicLoess(as.matrix(dt), method="pairs")
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Fast Loess Normalization of limma (using mean intensities over all the samples as its reference A sample)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default loessF)
#'
#' @return SummarizedExperiment containing the loessF normalized data as assay
#' @export
#'
loessFNorm <- function(se, ain="log2", aout="LoessF"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  norm_dt <- limma::normalizeCyclicLoess(as.matrix(dt), method="fast")
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' EigenMS Normalization
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default eigenMS)
#'
#' @return SummarizedExperiment containing the EigenMS normalized data as assay
#' @export
#'
eigenMSNorm <- function(se, ain="log2", aout="EigenMS"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  prot.info <- cbind(data.frame(SummarizedExperiment::rowData(se)$Protein.IDs), data.frame(SummarizedExperiment::rowData(se)$Protein.IDs))
  colnames(prot.info) <- c("pepIDs", "prID")
  condition <- S4Vectors::metadata(se)$condition
  grps <- factor(data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]], levels = unique(data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]]))
  ints_eig1 <- eig_norm1(m=dt, treatment=grps, prot.info = prot.info)
  ints_norm <- eig_norm2(rv=ints_eig1)
  # present
  prot_present <- ints_eig1$present$pepIDs
  # missing
  prot_missing <- prot.info[!prot.info$pepIDs %in% prot_present,]
  prot_missing$prID <- NULL
  # complete norm_dt
  norm_dt <- data.table::as.data.table(ints_norm$normalized)
  norm_dt$prID <- NULL
  if(nrow(prot_missing)>0){
    # add these proteins with complete NAs
    for(col in names(norm_dt)){ # -1 because of pepIDs columns
      if (col!="pepIDs") prot_missing[,col] <- NA
    }
    norm_dt <- rbind(norm_dt, prot_missing)
    # take ordering of rowData
    ordering <- SummarizedExperiment::rowData(se)$Protein.IDs
    norm_dt <- norm_dt %>% dplyr::arrange(factor(pepIDs, levels = ordering))
    norm_dt$pepIDs <- NULL
  } else {
    norm_dt$pepIDs <- NULL

  }
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Median Absolute Deviation Normalization (substracts the median and divides the data by the median absolute deviation (MAD))
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default mad)
#'
#' @return SummarizedExperiment containing the MAD normalized data as assay
#' @export
#'
medianAbsDevNorm <- function(se, ain="log2", aout="MAD"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt <- as.matrix(dt)
  if (!ain %in% c("raw", "vsn")){
    norm_dt <- NormalyzerDE::performSMADNormalization(dt, noLogTransform = TRUE)
  } else {
    norm_dt <- NormalyzerDE::performSMADNormalization(dt, noLogTransform = FALSE)
  }
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' RobNorm
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default robNorm)
#' @param gamma.0 Numeric representing the exponent of the weighted density. When the sample size
#'                is small, the fitted population of some proteins could be locally trapped such
#'                that the variance of those proteins was very small under a large gamma. To avoid
#'                this, a small gamma is recommended. When sample size smaller than 40, then set
#'                gamma to 0.5 or 0.1.
#'
#' @return SummarizedExperiment containing the robNorm normalized data as assay
#' @export
#'
robNorm <- function(se, ain="log2", aout="RobNorm", gamma.0 = 0.1){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  X.0 <- as.matrix(dt)
  rownames(X.0) <- rownames(dt)
  tryCatch(
    {
      robnorm_res <- RobNorm::RobNorm(X.0, gamma.0=gamma.0, tol=10^(-4), step=200) # default parameter
      norm_dt <- data.table::as.data.table(robnorm_res$norm.data)
      colnames(norm_dt) <- colnames(dt)
      rownames(norm_dt) <- rownames(dt)
      SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
    },
    error = function(e){
      message(paste0("RobNorm: ", e$message))
    }
  )
  return(se)
}

#' limma::removeBatchEffects (limBE)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input (default log2)
#' @param aout String which assay should be used to save normalized data (default limBE)
#'
#' @return SummarizedExperiment containing the limBE normalized data as assay
#' @export
#'
limmaNorm <- function(se, ain = "log2", aout = "limBE"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  batch <- S4Vectors::metadata(se)$batch
  dt_batch <- limma::removeBatchEffect(dt, batch = batch_column)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(dt_batch)
  return(se)
}


## ----- Main Normalization Method Functions ----- ##

#' Function to return available normalization methods' identifier names
#'
#' @return vector of normalization methods
#' @export
#'
get_normalization_methods <- function(){
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "HarmonizR", "limBE")
  return(norm_names)
}

#' Normalize SummarizedExperiment object using different normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_all_normalization_methods())
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
normalize_se_single <- function(se, methods = NULL, gamma.0 = 0.5){
  # vector with available normalization methods
  norm_functions <- norm_functions <- list(globalMeanNorm, globalMedianNorm, medianNorm, meanNorm, irsNorm,
                                           quantileNorm, vsnNorm, loessFNorm, loessCycNorm, rlrNorm,
                                           rlrMANorm, rlrMACycNorm, eigenMSNorm, medianAbsDevNorm, robNorm, tmmNorm, run_harmonizR, limmaNorm)
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "HarmonizR", "limBE")
  names(norm_functions) <- norm_names

  # retrieve normalization methods & check if all methods available
  if(!is.null(methods)){
    not_available_methods <- methods[!methods %in% norm_names]
    if(length(not_available_methods) > 0){
      warning(paste0(paste0(not_available_methods, collapse = ", "), " normalization methods not available!"))
    }
  } else {
    message("All available normalization methods will be performed.")
    methods <- norm_names
  }

  # check if IRS, limBE or HarmonizR can be executed
  batch <-S4Vectors::metadata(se)$batch
  if(sum(c("limBE", "IRS", "HarmonizR") %in% methods) > 0){
    if(is.null(batch)){
      stop("No batch specified! Batch need to be specified for limBE, IRS, and HarmonizR in the SummarizedExperiment under metadata(se)$batch!")
    }
  }
  refs <- S4Vectors::metadata(se)$refs
  if("IRS" %in% methods){
    if(is.null(refs)){
      stop("No reference samples specified! Reference samples need to be specified for IRS in the SummarizedExperiment under metadata(se)$refs!")
    }
  }

  # normalization
  for(method in methods){
    func <- norm_functions[[method]]
    if(method == "RobNorm"){
      se <- func(se, aout = method, gamma.0 = gamma.0)
    } else {
      se <- func(se, aout = method)
    }
    # TODO: error handling
    message(paste0(method, " completed."))
  }
  return(se)
}


#' Normalize SummarizedExperiment object using combinations of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_all_normalization_methods())
#' @param ains Vector of assays of SummarizedExperiment object to apply the normalization methods (e.g. if you want to perform Median normalization on IRS-normalized data)
#' @param combination_pattern String to give name to combination of methods (e.g. IRS_on_Median --> "_on_")
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
 normalize_se_combination <- function(se, methods, ains, combination_pattern = "_on_", gamma.0 = 0.5){

  # vector with available normalization methods
  norm_functions <- norm_functions <- list(globalMeanNorm, globalMedianNorm, medianNorm, meanNorm, irsNorm,
                                           quantileNorm, vsnNorm, loessFNorm, loessCycNorm, rlrNorm,
                                           rlrMANorm, rlrMACycNorm, eigenMSNorm, medianAbsDevNorm, robNorm, tmmNorm)
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM")
  names(norm_functions) <- norm_names

  # retrieve normalization methods & check if all methods available
  if(!is.null(methods)){
    # check if all methods available
    not_available_methods <- methods[!methods %in% norm_names]
    if(length(not_available_methods) > 0){
      warning(paste0(paste0(not_available_methods, collapse = ", "), " normalization methods not available!"))
    }
  } else {
    message("All available normalization methods will be performed.")
    methods <- norm_names
  }

  # perform combination of normalization
  for(ain in ains){
    # check if ain already in se --> if not: perform now
    if(! ain %in% names(SummarizedExperiment::assays(se))){
      message(paste0(ain, " normalization not yet performed. Single ", ain, " normalization performed now."))
      se <- normalize_se_single(se, methods = c(ain), gamma.0 = gamma.0)
    }

    for(method in methods){
      aout <- paste0(method, combination_pattern, ain)
      func <- norm_functions[[method]]
      if(method == "RobNorm"){
        se <- func(se, ain = ain, aout = aout, gamma.0 = gamma.0)
      } else {
        se <- func(se, ain = ain, aout = aout)
      }
      message(paste0(method, " normalization performed on ", ain, "-normalized data completed."))
    }
  }
  return(se)
}


#' Normalize SummarizedExperiment object using single normalization methods or specified combinations of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_all_normalization_methods())
#' @param combination_pattern String specifying how normalization methods are combined. For instance, methods = c("IRS", "Median_on_IRS"), combination_pattern = "_on_".
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
normalize_se <- function(se, methods, combination_pattern = "_on_", gamma.0 = 0.5){
  # extract combination of methods
  if(!is.null(combination_pattern)){
    comb_methods <- methods[stringr::str_detect(methods, combination_pattern)] #  combined methods
    sing_methods <- methods[!methods %in% comb_methods] # single methods
  } else {
    sing_methods <- methods
  }

  # single normalization
  if(length(sing_methods) > 0){
    se <- normalize_se_single(se, sing_methods, gamma.0 = gamma.0)
  }
  # combined normalization
  if(!is.null(combination_pattern)){
    if(length(comb_methods) > 0){
      for(m in comb_methods){
        method <- strsplit(m, combination_pattern)[[1]][1]
        ain <- strsplit(m, combination_pattern)[[1]][2]
        se <- normalize_se_combination(se, c(method), c(ain), combination_pattern, gamma.0 = gamma.0)
      }
    }
  }
  return(se)
}

