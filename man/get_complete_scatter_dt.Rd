% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NormalizationPlots.R
\name{get_complete_scatter_dt}
\alias{get_complete_scatter_dt}
\title{Function to get a long data table of first and second sample values of all kind of normalization}
\usage{
get_complete_scatter_dt(se, ain = NULL)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomic dataset}

\item{ain}{Vector of strings which assay should be used as input (default NULL).
If NULL then all normalization of the se object are plotted next to each other.}
}
\value{
data table
}
\description{
Function to get a long data table of first and second sample values of all kind of normalization
}
