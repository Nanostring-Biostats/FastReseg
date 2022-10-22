#' @title update_transDF_ResegActions
#' @description Update transcript data.frame based on resegmentation action, calculate the new cell type and mean per cell spatial coordinates
#' @param transcript_df the data.frame of transcript to be updated 
#' @param reseg_full_converter a named converter to update the cell ID in `transcript_df`, cell_ID in name would be converted to cell_ID in value; discard cell_ID with value = NA
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type, needed if using `LogLikeRatio` cell typing method (default  = NULL)
#' @param refProfiles A matrix of cluster profiles, genes X clusters, needed if using `NegBionomial` cell typing method (default  = NULL)
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param celltype_coln the column name of cell type in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in transcript_df 
#' @param return_perCellDF flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param celltype_method use either `LogLikeRatio` or `NegBinomial` method for quick cell typing and corresponding score_baseline calculation (default = LogLikeRatio)
#' @return a list 
#' \describe{
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{perCell_DT}{a per cell data.table with mean sptial coordinates and new cell type when return_perCellDF = TRUE}
#'    \item{perCell_expression}{a gene x cell count sparse matrix for updated transcript data.frame when return_perCellDF = TRUE}
#' }
#' @details Update transcript data.frame based on resegmentation action and get new cell type; when return_perCellDF = TRUE, return gene x cell count matrix and per cell data.frame with mean per cell spatial coordinates and new cell type.
#' @export
update_transDF_ResegActions <- function(transcript_df, 
                                        reseg_full_converter, 
                                        score_GeneMatrix = NULL,
                                        refProfiles = NULL, 
                                        transGene_coln = 'target',
                                        cellID_coln = 'cell_ID', 
                                        celltype_coln = 'cell_type',
                                        spatLocs_colns = c("x","y","z"),
                                        return_perCellDF = TRUE, 
                                        celltype_method = 'LogLikeRatio'){
  
  celltype_method <- match.arg(celltype_method, c('LogLikeRatio', 'NegBinomial')) 
  
  if(celltype_method == 'LogLikeRatio' & is.null(score_GeneMatrix)){
    stop("Must provided `score_GeneMatrix` when using log-likelihood ratio based cell typing method.")
  }
  
  if(celltype_method == 'NegBinomial' & is.null(refProfiles)){
    stop("Must provided `refProfiles` when using negative binomial cell typing method.")
  }
  
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transGene_coln, celltype_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transGene_coln, celltype_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  # get common cells
  common_cells <- intersect(unique(names(reseg_full_converter)), unique(transcript_df[[cellID_coln]]))
  
  # get common genes
  if(celltype_method == 'LogLikeRatio'){
    common_genes <- intersect(rownames(score_GeneMatrix), 
                              unique(transcript_df[[transGene_coln]]))
    score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  } else if (celltype_method =='NegBinomial'){
    common_genes <- intersect(rownames(refProfiles), 
                              unique(transcript_df[[transGene_coln]]))
    refProfiles <- refProfiles[common_genes, ]
  }
  
  
  message(sprintf("Found %d common cells and %d common genes among `names(reseg_full_converter)`, `transcript_df`, and `score_GeneMatrix` or `refProfiles. ", 
                  length(common_cells), length(common_genes)))
  
  if(any(length(common_cells) <1, length(common_genes)<1)){
    stop("Too few common cells or genes to proceed. Check if `score_GeneMatrix` or `refProfiles` is a gene x cell-type matrix.")
  }
  
  
  
  # split reseg_full_converter into different types of cells
  # get idx
  cells_to_discard <- which(is.na(reseg_full_converter))
  cells_to_keep <- which(reseg_full_converter == names(reseg_full_converter))
  cells_to_update <- setdiff(seq(1, length(reseg_full_converter)), c(cells_to_discard, cells_to_keep))
  # get values
  cells_to_update <- reseg_full_converter[cells_to_update]
  cells_to_keep <- names(reseg_full_converter)[cells_to_keep]
  cells_to_discard <- names(reseg_full_converter)[cells_to_discard]
  
  # discard cells
  transcript_df <- transcript_df[which(!(transcript_df[[cellID_coln]] %in% cells_to_discard)), ]
  
  # update cell ID
  transcript_df[['updated_cellID']] <- transcript_df[[cellID_coln]]
  tmp_idx <- which(transcript_df[[cellID_coln]] %in% names(cells_to_update))
  transcript_df[['updated_cellID']][tmp_idx] <- cells_to_update[transcript_df[[cellID_coln]][tmp_idx]]
  
  ## get new cell types for cells being updated ----
  subTransDF <- transcript_df[which(transcript_df[['updated_cellID']] %in% unique(cells_to_update) & transcript_df[[transGene_coln]] %in% common_genes), ]
  
  if(nrow(subTransDF)<1){
    newCellTypes <- NULL
  } else if(celltype_method == 'LogLikeRatio'){
    # get score for each transcripts
    transcriptGeneScore <- score_GeneMatrix[subTransDF[[transGene_coln]], ]
    
    tmp_score <- as.data.frame(transcriptGeneScore)
    tmp_score[['updated_cellID']] <- subTransDF[['updated_cellID']]
    
    tmp_score <- data.table::setDT(tmp_score)[, lapply(.SD, sum), by = 'updated_cellID'] 
    tmp_cellID <- tmp_score[['updated_cellID']]
    tmp_score[['updated_cellID']] <- NULL
    # assign cell type based on max values
    max_idx_1st <- max.col(tmp_score,ties.method="first")
    newCellTypes <- colnames(tmp_score)[max_idx_1st]
    names(newCellTypes) <- tmp_cellID
    
  } else if (celltype_method =='NegBinomial'){
    exprMat <- reshape2::acast(subTransDF, as.formula(paste('updated_cellID', '~', transGene_coln)), length)
    # fill missing genes that in refProfiles but not in current data as 0
    missingGenes <- setdiff(rownames(refProfiles), colnames(exprMat))
    exprMat <- cbind(exprMat, 
                     matrix(0, nrow = nrow(exprMat), ncol = length(missingGenes), 
                            dimnames = list(rownames(exprMat), missingGenes)))
    exprMat <- exprMat[, rownames(refProfiles), drop = FALSE]
    
    newCellTypes <- quick_celltype(exprMat, bg = 0, reference_profiles = refProfiles, align_genes = FALSE)[['clust']]
    rm(exprMat, missingGenes)
    
  }
  
  
  
  # update cell type
  transcript_df[['updated_celltype']] <- transcript_df[[celltype_coln]]
  tmp_idx <- which(transcript_df[['updated_cellID']] %in% unique(cells_to_update))
  transcript_df[['updated_celltype']][tmp_idx] <- newCellTypes[transcript_df[['updated_cellID']][tmp_idx]]
  
  outputs <- list(updated_transDF = transcript_df)
  
  # get per cell data frame
  if(return_perCellDF){
    transcript_df <- data.table::setDT(transcript_df)
    # get cell type
    perCell_DT <- transcript_df[, .SD[1], by = updated_cellID, .SDcols = 'updated_celltype']
    
    # get mean spatial coordinates
    if(!is.null(spatLocs_colns)){
      perCell_dt2 <- transcript_df[, lapply(.SD, mean), by = 'updated_cellID', .SDcols = spatLocs_colns] 
      perCell_DT <- merge(perCell_DT, perCell_dt2, by = 'updated_cellID')
    }
    
    # get per cell expression matrix
    exprs_tmp <- transcript_df[, .SD, .SDcols = c('updated_cellID', transGene_coln)]
    exprs_tmp <- reshape2::dcast(exprs_tmp, paste0("updated_cellID~", transGene_coln), 
                                 value.var = transGene_coln, 
                                 fun.aggregate = length, fill = 0)
    rownames(exprs_tmp) <- exprs_tmp[['updated_cellID']]
    exprs_tmp <- Matrix::Matrix(as.matrix(exprs_tmp[, -1]), sparse = TRUE)
    
    outputs[['perCell_DT']] <- perCell_DT
    outputs[['perCell_expression']] <- Matrix::t(exprs_tmp)
  }
  
  
  return(outputs)
}