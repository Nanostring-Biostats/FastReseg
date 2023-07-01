#' E step: estimate each cluster's mean profile (Ptolemy Estep function)
#' @title estimate_MeanProfile
#' @description  Given cell assignments (or posterior probabilities), estimate the mean profile of each cluster.
#' @param counts Counts matrix, cells X genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities of cells (rows) belonging to clusters (columns).
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @importFrom Matrix colMeans
#' @return A matrix of cluster profiles, genes X clusters
#' @details estimate each cluster's mean profile (Ptolemy Estep function)
#' @export
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
#' @description calculate log-likilhood score of each gene based on reference expression profiles and return the centered score matrix
#' @param genes a vector of gene name to score
#' @param ref_profiles a gene X cell_type expression matrix for reference profiles
#' @param flag_center flag to center the score matrix per gene before turn, default = TRUE 
#' @return loglik, a gene X cell_type matrix of centered loglik score for each gene
#' @export
scoreGenesInRef <- function(genes, ref_profiles, flag_center= TRUE){
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
  
  # tLLRv2 score, re-center on maximum per row/transcript
  if(flag_center){
    tmp_max <- apply(loglik, 1, max)
    loglik <- sweep(loglik, 1, tmp_max, '-')
  }
  
  
  return(loglik)
}



# get new cell type assignment based on max score
#' @title getCellType_maxScore
#' @description get the cell type give maximum score
#' @param score_GeneMatrix a gene x cell-type score matrix
#' @param transcript_df the data.frame of transcript_ID and cell_ID
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @return return a named vector with cell type in values and cell_ID in names
#' @export
getCellType_maxScore <- function(score_GeneMatrix, 
                                 transcript_df, 
                                 transGene_coln = "target",
                                 cellID_coln = "cell_ID"){
  # check format of transcript_df
  if(any(!c(cellID_coln, transGene_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transGene_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }

  # convert to cell x gene count matrix
  counts <- data.table::dcast(data.table::as.data.table(transcript_df), 
                              formula = as.formula(paste0(cellID_coln, " ~ ", transGene_coln)), 
                              drop = F, value.var = cellID_coln, fun.aggregate = length)
  counts <- Matrix::Matrix(as.matrix(counts[, -1]), 
                           dimnames = list(counts[[1]], colnames(counts)[-1]),
                           sparse = T)

  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            colnames(counts))
  message(sprintf("Found %d common genes among transcript_df and score_GeneMatrix. ", 
                  length(common_genes)))
  
  if(length(common_genes)<1){
    stop("Too few common genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  # get cell x cell-type score matrix 
  score_cellMatrix <- counts[, common_genes, drop = F] %*% score_GeneMatrix[common_genes, , drop = F]
  
  # assign cell type based on max values
  max_idx_1st <- max.col(score_cellMatrix,ties.method="first")
  newCellTypes <- colnames(score_cellMatrix)[max_idx_1st]
  names(newCellTypes) <- rownames(score_cellMatrix)
  
  return(newCellTypes)
}




#' @title getScoreCellType_gene
#' @description get each transcript's score based on score matrix and chosen cell-type
#' @param score_GeneMatrix a gene x cell-type score matrix
#' @param transcript_df the data.frame of transcript_ID and cell_ID
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param celltype_coln the column name of cell type in transcript_df
#' @return a named vector with score of given cell type in values and transcript_id in names
#' @export
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
  transcript_df <- as.data.frame(transcript_df)[, c(transID_coln, transGene_coln, celltype_coln)]
  
  common_genes <- intersect(transcript_df[[transGene_coln]], rownames(score_GeneMatrix))
  if(length(common_genes) ==0){
    stop("No common genes found between score_GeneMatrix and the transcripts of chosen_cells in transcript_df. Check if score_GeneMatrix in transcipt x cell-type format.")
  } 
  
  # get same transcript_id order 
  transcript_df <- transcript_df[which(transcript_df[[transGene_coln]] %in% common_genes), ]
  
  celltype_NScore <- score_GeneMatrix[as.matrix(transcript_df[, c(transGene_coln, celltype_coln)])]
  names(celltype_NScore) <- transcript_df[[transID_coln]]
  
  return(celltype_NScore)
}