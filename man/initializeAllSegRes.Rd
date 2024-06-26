% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastReseg_full_pipeline.R
\name{initializeAllSegRes}
\alias{initializeAllSegRes}
\title{initializeAllSegRes}
\usage{
initializeAllSegRes(
  prep_res,
  save_intermediates = TRUE,
  return_perCellData = TRUE
)
}
\arguments{
\item{prep_res}{nested list of elements for 'refProfiles', 'baselineData', 'cutoffs_list', outputs of \code{runPreprocess}}

\item{save_intermediates}{flag to save intermediate outputs into output folder, including data.frame for spatial modeling statistics of each cell,}

\item{return_perCellData}{flag to return and save to output folder for gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type}
}
\value{
a list
\describe{
\item{refProfiles}{a genes X clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
\item{baselineData}{a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; \code{span_score} is for the average per molecule transcript tLLR score of each cell, \code{span_transNum} is for the transcript number of each cell.}
\item{cutoffs_list}{a list of cutoffs used in resegmentation pipeline, including, \code{score_baseline}, \code{lowerCutoff_transNum}, \code{higherCutoff_transNum}, \code{cellular_distance_cutoff}, \code{molecular_distance_cutoff}}
\item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types}
\item{updated_perCellDT}{an empty list, return when \code{return_perCellData} = TRUE}
\item{updated_perCellExprs}{an empty list, return when \code{return_perCellData} = TRUE}
\item{reseg_actions}{a list of 4 empty elements, same output format as \code{decide_ReSegment_Operations}, return when \code{save_intermediates} = TRUE}
}
}
\description{
initialize result holder for \code{fastReseg_full_pipeline} based on preprocessed data
}
