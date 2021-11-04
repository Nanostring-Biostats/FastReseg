#' E step: estimate each cluster's mean profile (Ptolemy Estep function)
#' @title estimate_MeanProfile
#' @description  Given cell assignments (or posterior probabilities), estimate the mean profile of each cluster.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities of cells (rows) belonging to clusters (columns).
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @importFrom Matrix colMeans
#' @return A matrix of cluster profiles, genes * clusters
#' @details estimate each cluster's mean profile (Ptolemy Estep function)
estimate_MeanProfile <- function(counts, clust, s, bg) {
  
  # scale counts:
  counts = sweep(pmax(counts - bg, 0), 1, s, "/")
  
  # get cluster means:
  if (is.vector(clust)) {
    means = sapply(unique(clust), function(cl) {
      Matrix::colMeans(counts[clust == cl, , drop = FALSE])
    })
  }
  if (is.matrix(clust)) {
    means = apply(clust, 2, function(x) {
      wts = x / mean(x)
      Matrix::colMeans(sweep(counts, 1, wts, "*"))
    })
  }
  
  return(means)
}



# calculate the loglik of each gene in expression matrix
#' @title scoreGenesInRef
#' @description calculate log-likilhood score of each gene based on reference expression profiles
#' @param genes a vector of gene name to score
#' @param ref_profiles a gene x cell_type expression matrix for reference profiles
#' @return loglik, a gene x cell_type matrix of loglik score for each gene
scoreGenesInRef <- function(genes, ref_profiles){
  common_feats <- intersect(unique(genes), rownames(ref_profiles))
  if(length(common_feats) <1){
    stop("No common genes found in `genes` and `ref_profiles`. ref_profiles should be a gene x cell_type matrix.")
  }
  common_feats <- sort(common_feats)
  # filter, normalize and log transform the ref_profiles
  ref_profiles <- as.matrix(ref_profiles)
  ref_profiles <- ref_profiles[, order(colnames(ref_profiles))]
  ref_profiles <- ref_profiles[common_feats, ]
  libsize <- colSums(ref_profiles, na.rm = TRUE)
  norm_exprs <- Matrix::t(Matrix::t(ref_profiles)/ libsize)
  loglik <- log(norm_exprs)
  return(loglik)
}



# get new cell type assignment based on max score
#' @title getCellType_maxScore
#' @description get the cell type give maximum score
#' @param score_GeneMatrix a gene x cell-type score matrix
#' @param transcript_df the data.frame of transcript_ID and cell_ID
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param return_transMatrix logic flag whether to return the transcript level of score matrix (default = FALSE)
#' @return return a list
#' \enumerate{
#'    \item{cellType_DF, a data.frame of cell_ID and cell_type}
#'    \item{score_TransMatrix, a transcript x cell-type score matrix when return_transMatrix = TRUE}
#' }
getCellType_maxScore <- function(score_GeneMatrix, 
                                 transcript_df, 
                                 transID_coln = 'transcript_id',
                                 transGene_coln = "target",
                                 cellID_coln = "cell_ID", 
                                 return_transMatrix = FALSE){
  # check format of transcript_df
  if(any(!c(cellID_coln, transGene_coln, transID_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transGene_coln, transID_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transGene_coln, transID_coln)]
  
  
  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            unique(transcript_df[[transGene_coln]]))
  message(sprintf("Found %d common genes among transcript_df, and score_GeneMatrix. ", 
                  length(common_genes)))
  
  if(any(length(common_genes)<1)){
    stop("Too few common genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  transcript_df <- transcript_df[which(transcript_df[[transGene_coln]] %in% common_genes), ]
  
  
  # get score for each transcripts
  transcriptGeneScore <- score_GeneMatrix[transcript_df[[transGene_coln]], ]
  rownames(transcriptGeneScore) <- transcript_df[[transID_coln]]
  tmp_score <- as.data.frame(transcriptGeneScore)
  
  # sparse matrix gave bigger size, so not to use
  if(!return_transMatrix){
    # to save some memory 
    rm(transcriptGeneScore)
  }
  
  
  tmp_score[[cellID_coln]] <- transcript_df[[cellID_coln]]
  
  tmp_score <-data.table::setDT(tmp_score)[, lapply(.SD, sum), by = cellID_coln] 
  tmp_cellID <- tmp_score[[cellID_coln]]
  tmp_score[[cellID_coln]] <- NULL
  # assign cell type based on max values
  max_idx_1st <- max.col(tmp_score,ties.method="first")
  newCellTypes <- colnames(tmp_score)[max_idx_1st]
  
  cellType_DF <- data.frame(cell_ID = tmp_cellID, 
                            cell_type = newCellTypes)
  outputs <- list(cellType_DF = cellType_DF)
  
  if(return_transMatrix){
    outputs[[score_TransMatrix]] <- transcriptGeneScore 
  }
  
  return(outputs)
}




#' @title getScoreCellType_gene
#' @description get each transcript's score based on score matrix and chosen cell-type
#' @param score_GeneMatrix a gene x cell-type score matrix
#' @param transcript_df the data.frame of transcript_ID and cell_ID
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param celltype_coln the column name of cell type in transcript_df
#' @return score_df, a data.frame of "[transID_coln]" and "score_[celltype_coln]" column for chosen cell-type
getScoreCellType_gene <- function(score_GeneMatrix, transcript_df, 
                                  transID_coln = "transcript_id",
                                  transGene_coln = "target",
                                  celltype_coln = 'cell_type'){
  
  # check format of transcript_df
  if(any(!c(transID_coln, transGene_coln, celltype_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(transID_coln, transGene_coln, celltype_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(transID_coln, transGene_coln, celltype_coln)]
  
  common_genes <- intersect(transcript_df[[transGene_coln]], rownames(score_GeneMatrix))
  if(length(common_genes) ==0){
    stop("No common genes found between score_GeneMatrix and the transcripts of chosen_cells in transcript_df. Check if score_GeneMatrix in transcipt x cell-type format.")
  } 
  
  # get same transcript_id order 
  transcript_df <- transcript_df[which(transcript_df[[transGene_coln]] %in% common_genes), ]
  score_matrix <- score_GeneMatrix[transcript_df[[transGene_coln]], ]
  
  # loop over each cell type to get score
  all_celltypes <- unique(transcript_df[[celltype_coln]])
  score_df <- data.frame(transcript_id = transcript_df[[transID_coln]], 
                         celltype_NScore = rep(NA, nrow(transcript_df)))
  for (each_celltype in all_celltypes){
    rowidx <- which(transcript_df[[celltype_coln]] == each_celltype)
    score_df[['celltype_NScore']][rowidx] <- score_matrix[rowidx, each_celltype]
  }
  colnames(score_df) <- c(transID_coln, paste0('score_', celltype_coln))
  
  return(score_df)
}