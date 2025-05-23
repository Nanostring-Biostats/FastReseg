% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_baselineCT}
\alias{example_baselineCT}
\title{Percentile profile for single cell distribution of transcript number and average tLLR score from example dataset for spatial transcriptional profiling of tissue}
\format{
An object of class \code{"list"} with 2 elements
\describe{
\item{span_score}{a matrix of average transcript tLLR score per cell for 22 distinct cell types in rows, percentile at (0\%, 25\%, 50\%, 75\%, 100\%) in columns}
\item{span_transNum}{a matrix of transcript number per cell for 22 distinct cell types in rows, percentile at (0\%, 25\%, 50\%, 75\%, 100\%) in columns}
}
}
\usage{
data(example_baselineCT)
}
\description{
The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. The single cell distribution was calculated for transcript number per cell and average transcript tLLR score per cell, and then grouped by the assigned cell types based on the original cell segmentation outcome to get the percentile profile for each cell type at 25\% increment step.
}
\keyword{datasets}
