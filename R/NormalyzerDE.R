###### ------- NormalyzerDE Functions ------- ######

#' Calculate CV per replicate group and normalization technique
#'
#' Iterates through each normalization method and calculate average CV values
#' per replicate group.
#'
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return avgCVPerNormAndReplicates Matrix with group CVs as rows and
#'   normalization technique as columns
#' @keywords internal
calculateReplicateCV <- function(methodList, sampleReplicateGroups) {

  calculateFeatureCVs <- function(feature, groups) {
    featureCVs <- RcmdrMisc::numSummary(
      feature,
      statistics=c("cv"),
      groups=groups)

    featureCVs$table
  }

  calculateMethodReplicateCVs <- function(methodData, groups) {

    cvPerFeatureAndGroup <- apply(
      methodData,
      1,
      calculateFeatureCVs,
      groups=sampleReplicateGroups
    )

    if (length(unique(groups)) > 1) {
      # Transpose to get groups as column heads
      featureCondCVs <- t(cvPerFeatureAndGroup)
    }
    else if (length(unique(groups)) == 1) {
      # For one feature the CVs are dropped to a vector
      featureCondCVs <- data.frame(cvPerFeatureAndGroup)
      colnames(featureCondCVs) <- groups[1]
    }
    else {
      stop("Unknown state encountered for groups:", paste(groups, collapse=", "))
    }

    # Calculate mean CV for all features for each condition
    summedCondCVs <- colMeans(featureCondCVs, na.rm=TRUE)
    summedCondCVs * 100
  }

  avgCVPerNormAndReplicates <- vapply(
    methodList,
    calculateMethodReplicateCVs,
    rep(0, length(unique(sampleReplicateGroups))),
    groups=sampleReplicateGroups
  )

  # If only one group, this is needed to get data in correct shape and format
  if (length(unique(sampleReplicateGroups)) == 1) {
    avgCVPerNormAndReplicates <- t(as.matrix(avgCVPerNormAndReplicates))
  }

  avgCVPerNormAndReplicates
}

#' Calculate CV values for each feature. Iterates through each normalization
#' method and calculates a matrix of CV values where each column correspond to
#' a method and each row corresponds to a feature.
#'
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return methodFeatureCVMatrix Matrix with feature as rows and normalization
#'   method as columns
#' @keywords internal
calculateFeatureCV <- function(methodList) {

  methodCount <- length(methodList)
  numberFeatures <- nrow(methodList[[1]])

  calculateFeatureCVVector <- function(processedDataMatrix) {
    cv <- function(row) {
      stDev <- stats::sd(row, na.rm=TRUE)
      meanVal <- mean(unlist(row), na.rm=TRUE)
      stDev / mean(meanVal) * 100
    }
    featureCVs <- apply(processedDataMatrix, 1, cv)
    featureCVs
  }

  methodFeatureCVMatrix <- vapply(
    methodList,
    calculateFeatureCVVector,
    rep(0, numberFeatures)
  )

  colnames(methodFeatureCVMatrix) <- names(methodList)
  methodFeatureCVMatrix
}

#' Calculate average MAD (Median Absolute Deviation) for each feature in
#' each condition and then calculates the average for each replicate group
#'
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return condAvgMadMat Matrix with average MAD for each biological condition.
#' @keywords internal
calculateAvgMadMem <- function(methodList, sampleReplicateGroups) {

  groupIndexList <- getIndexList(sampleReplicateGroups)

  calculateAvgFeatureMadForGroup <- function(groupIndices, methodData) {
    groupData <- methodData[[groupIndices]] # .. added
    featureMAD <- matrixStats::rowMads(as.matrix(groupData), na.rm=TRUE) # as.matrix added
    featureMAD
  }

  calculateGroupMadForMethod <- function(methodData, groups, indexList) {

    # Extracts groups of replicates and calculate MAD for each feature
    featureMADMat <- vapply(
      indexList,
      calculateAvgFeatureMadForGroup,
      rep(0, nrow(methodData)),
      methodData=methodData
    )

    methodRepGroupMADMean <- colMeans(featureMADMat, na.rm=TRUE)
    methodRepGroupMADMean
  }

  condAvgMadMat <- vapply(
    methodList,
    calculateGroupMadForMethod,
    rep(0, length(unique(sampleReplicateGroups))),
    groups=sampleReplicateGroups,
    indexList=groupIndexList
  )

  condAvgMadMat
}

#' Calculate average variance for each feature in each condition and then
#' calculate the average for each replicate group
#'
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return avgVarianceMat Matrix with average variance for each biological
#' condition
#' @keywords internal
calculateAvgReplicateVariation <- function(methodList, sampleReplicateGroups) {

  groupIndexList <- getIndexList(sampleReplicateGroups)

  calculateReplicateGroupVariance <- function(groupIndices, methodData) {

    groupData <- methodData[[groupIndices]] # .. added
    groupData <- as.matrix(groupData) # line added
    rowNonNACount <- rowSums(!is.na(groupData)) - 1
    rowVariances <- rowNonNACount * matrixStats::rowVars(groupData, na.rm=TRUE)
    replicateGroupVariance <- sum(rowVariances, na.rm=TRUE) / sum(rowNonNACount, na.rm=TRUE)
    replicateGroupVariance
  }

  avgVarianceMat <- vapply(
    methodList,
    function(methodData) {
      replicateGroupVariance <- vapply(
        groupIndexList,
        calculateReplicateGroupVariance,
        0,
        methodData=methodData
      )
      replicateGroupVariance
    },
    rep(0, length(unique(sampleReplicateGroups)))
  )

  avgVarianceMat
}

#' General function for calculating percentage difference of average column
#' means in matrix
#'
#' @param targetMat Matrix for which column means should be compared
#' @return percDiffVector Vector with percentage difference, where first element
#'   always will be 100
#' @keywords internal
calculatePercentageAvgDiffInMat <- function(targetMat) {

  calculatePercDiff <- function (sampleIndex, mat) {
    mean(mat[, sampleIndex]) * 100 / mean(mat[,"log2"])
  }

  percDiffVector <- vapply(
    seq_len(ncol(targetMat)),
    calculatePercDiff,
    0,
    mat=targetMat)

  percDiffVector
}

#' Return list containing vector positions of values in string
#'
#' @param targetVector vector of sample groups
#' @return indexList List where key is condition level and values are indices
#'   for the condition
#' @keywords internal
getIndexList <- function(targetVector) {

  indexList <- list()
  uniqVals <- unique(targetVector)
  for (val in uniqVals) {
    indexList[[toString(val)]] <- which(targetVector == val)
  }
  indexList
}

#' Calculates correlation values between replicates for each condition matrix.
#' Finally returns a matrix containing the results for all dataset
#'
#' @param methodlist List containing normalized matrices for each normalization
#'   method
#' @param allReplicateGroups Vector with condition groups matching the columns
#'   found in the normalization methods
#' @param sampleGroupsWithReplicates Unique vector with condition groups
#'   present in two or more samples
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return avgCorSum Matrix with column corresponding to normalization
#' approaches and rows corresponding to replicate group
#' @keywords internal
calculateSummarizedCorrelationVector <- function(
    methodlist, allReplicateGroups, sampleGroupsWithReplicates, corrType) {

  validCorrTypes <- c("pearson", "spearman")
  if (!corrType %in% validCorrTypes) {
    stop("Unknown correlation type: ",
         corrType,
         " valid are: ",
         paste(validCorrTypes, collapse=", "))
  }

  corr_combination_count <- function(allReplicateGroups) {
    replicate_counts <- table(allReplicateGroups)
    sum(vapply(
      replicate_counts,
      function(count) { (count * (count-1)) / 2 },
      0))
  }

  avgCorSum <- vapply(
    methodlist,
    calculateCorrSum,
    rep(0, corr_combination_count(allReplicateGroups)),
    allReplicateGroups=allReplicateGroups,
    sampleGroupsWithReplicates=sampleGroupsWithReplicates,
    corrType=corrType
  )

  avgCorSum
}

#' Calculates internal correlations for each condition having at least two
#' samples and returns a vector with correlation values corresponding to each
#' condition
#'
#' @param methodData Expression data matrix
#' @param allReplicateGroups Full condition header corresponding to data tables
#'   columns
#' @param sampleGroupsWithReplicates Unique conditions where number of
#'   replicates exceeds one
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return corSums
#' @keywords internal
calculateCorrSum <- function(methodData, allReplicateGroups,
                             sampleGroupsWithReplicates, corrType) {
  methodData <- as.matrix(methodData) # line added
  corSums <- vector()
  for (groupNbr in seq_along(sampleGroupsWithReplicates)) {

    specificReplicateVals <- as.matrix(
      methodData[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
    class(specificReplicateVals) <- "numeric"
    corVals <- stats::cor(
      specificReplicateVals ,
      use="pairwise.complete.obs",
      method=corrType)

    for (index in seq_len(ncol(specificReplicateVals) - 1)) {
      corSums <- c(
        corSums,
        corVals[index, -(seq_len(index)), drop="FALSE"]
      )
    }
  }

  corSums
}
