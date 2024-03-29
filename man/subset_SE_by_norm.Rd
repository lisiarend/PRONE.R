% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SubsetSE.R
\name{subset_SE_by_norm}
\alias{subset_SE_by_norm}
\title{Subset SummarizedExperiment object by normalization assays}
\usage{
subset_SE_by_norm(se, ain)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{ain}{Character vector of assay names to keep in the SummarizedExperiment object}
}
\value{
SummarizedExperiment object with only the selected normalization assays
}
\description{
Subset SummarizedExperiment object by normalization assays
}
