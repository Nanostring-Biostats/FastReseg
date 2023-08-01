# modular wrapper to perform cell segmentation refinement
#' @title runSegRefinement
#' @description modular wrapper to evalute transcript groups in neighborhood, decide resegmentation operations and execute
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param chosen_cells the cell_ID of chosen cells need to be evaluate for re-segmentation
#' @param reseg_transcript_df the data.frame with transcript_id, target/geneName, x, y, cell_id for all transcript groups and the cell type of maximum transcript scores for each transcript group
#' @param reseg_cellID_coln the column name of cell_ID for all transcript groups in transcript_df
#' @param reseg_celltype_coln the column name of cell_type for all transcript groups in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param score_baseline a named vector of score baseline for all cell type listed in neighborhood_df such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param neighbor_distance_xy maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param distance_cutoff maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `neighbor_distance_xy`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param spatialMergeCheck_method use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space (default = "leidenCut")
#' @param cutoff_spatialMerge spatial constraint on a valid merging event between two source transcript groups, default = 0.5 for 50% cutoff, set to 0 to skip spatial constraint evaluation for merging.   
#' For `spatialMergeCheck_method = "leidenCut"`, this is the minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event.
#' For `spatialMergeCheck_method = "geometryDiff"`, this is the maximum percentage of white space change upon merging of query cell and neighbor cell for a valid merging event. 
#' @param leiden_config (leidenCut) a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.  
#' @param config_spatNW_transcript (leidenCut) configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config
#' @param return_intermediates flag to return intermediate outputs, including `neighborhoodDF_ToReseg` data.frame for neighborhood evaluation, `reseg_actions` list of resegmentation actions  
#' @param return_perCellData flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param includeAllRefGenes flag to include all genes in `score_GeneMatrix` in the returned `updated_perCellExprs` with missing genes of value 0 (default = FALSE)
#' @param seed_segRefine seed for transcript error correction step, default = NULL to skip the seed     
#' @return a list 
#' \describe{
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content` function, return when `return_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' @export
#' 
runSegRefinement <- function(score_GeneMatrix,  
                             chosen_cells = NULL, 
                             reseg_transcript_df, 
                             reseg_cellID_coln = "tmp_cellID", 
                             reseg_celltype_coln = "group_maxCellType", 
                             transID_coln = "transcript_id",
                             transGene_coln = "target", 
                             transSpatLocs_coln = c('x','y','z'),
                             score_baseline = NULL, 
                             lowerCutoff_transNum = NULL, 
                             higherCutoff_transNum= NULL, 
                             neighbor_distance_xy = NULL,
                             distance_cutoff = 2.7,
                             spatialMergeCheck_method = c("leidenCut", "geometryDiff"), 
                             cutoff_spatialMerge = 0.5, 
                             leiden_config = list(objective_function = c("CPM", "modularity"),
                                                  resolution_parameter = 1,
                                                  beta = 0.01,
                                                  n_iterations = 200),
                             config_spatNW_transcript = NULL,
                             return_intermediates = TRUE,
                             return_perCellData = TRUE,
                             includeAllRefGenes = FALSE, 
                             seed_segRefine = NULL){
  if(!is.null(seed_segRefine)){
    set.seed(seed_segRefine)
    }
  
  if(cutoff_spatialMerge <0 | cutoff_spatialMerge >1){
    stop(sprintf("The providied `cutoff_spatialMerge = %.3f`, must be within [0, 1].", cutoff_spatialMerge))
  } else if(cutoff_spatialMerge >0){
    # method of spatial constraint on merging
    spatialMergeCheck_method <- match.arg(spatialMergeCheck_method, c("leidenCut", "geometryDiff"))
    
    if(spatialMergeCheck_method == "leidenCut"){
      # check and set config for spatial network and leiden clustering 
      leiden_config <- check_config_leiden(leiden_config)
      config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript, 
                                                         spat_locs = transSpatLocs_coln)
       
    }
  }
  
  
  ### search within absolute distance, consider 25um in xy for cell level search and 15 pixel = 2.7um to be direct neighbor at transcript level.
  # using spatstat to locate neighbor cells and rank them by minimal molecular distance to query cell
  # `get_neighborhood_content` function returns a data.frame with each cell in row and its neighborhood information in columns
  neighborReSeg_df <- get_neighborhood_content(chosen_cells = chosen_cells,
                                                          score_GeneMatrix = score_GeneMatrix,
                                                          score_baseline = score_baseline,
                                                          neighbor_distance_xy = neighbor_distance_xy,
                                                          distance_cutoff = distance_cutoff,
                                                          transcript_df = reseg_transcript_df,
                                                          cellID_coln = reseg_cellID_coln,
                                                          celltype_coln = reseg_celltype_coln,
                                                          transID_coln = transID_coln,
                                                          transGene_coln = transGene_coln,
                                                          transSpatLocs_coln = transSpatLocs_coln)
  

  #### decide resegmentation operation: merge, new cell, or discard ----
  # # `decide_ReSegment_Operations` function returns a list containing the following 4 elements:
  # `cells_to_discard`: a vector of cell ID that should be discarded during resegmentation
  # `cells_to_update`: a named vector of cell ID whether the cell_ID in name would be replaced with cell_ID in value.
  # `cells_to_keep`: a vector of cell ID that should be kept as it is.
  # `reseg_full_converter`: a single named vector of cell ID to update the original cell ID, assign NA for cells_to_discard.
  reseg_actions <- decide_ReSegment_Operations(neighborhood_df = neighborReSeg_df, 
                                               score_baseline = score_baseline, 
                                               lowerCutoff_transNum = lowerCutoff_transNum, 
                                               higherCutoff_transNum = higherCutoff_transNum,
                                               transcript_df = reseg_transcript_df, 
                                               cellID_coln = reseg_cellID_coln,
                                               transID_coln = transID_coln,
                                               transSpatLocs_coln = transSpatLocs_coln, 
                                               spatialMergeCheck_method = spatialMergeCheck_method, 
                                               cutoff_spatialMerge = cutoff_spatialMerge, 
                                               leiden_config = leiden_config,
                                               config_spatNW_transcript = config_spatNW_transcript)
  
  #### update the transcript data.frame based on reseg_actions ----
  # update neighborReSeg_df_cleanSVM with resegmentaion actions
  neighborReSeg_df[['corrected_CellId']] <- reseg_actions[['reseg_full_converter']][neighborReSeg_df[['CellId']]]
  
 
  # # `update_transDF_ResegActions` function returns a list containing the following 3 elements:
  # `updated_transDF`, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter.
  # `perCell_DT`, a per cell data.table with mean spatial coordinates and new cell type when return_perCellDF = TRUE.
  # `perCell_expression`, a gene x cell count sparse matrix for updated transcript data.frame when return_perCellDF = TRUE.
  
  post_reseg_results <- update_transDF_ResegActions(transcript_df = reseg_transcript_df, 
                                                    reseg_full_converter = reseg_actions$reseg_full_converter, 
                                                    score_GeneMatrix = score_GeneMatrix, 
                                                    transGene_coln = transGene_coln, 
                                                    cellID_coln = reseg_cellID_coln, 
                                                    celltype_coln = reseg_celltype_coln, 
                                                    spatLocs_colns = transSpatLocs_coln, 
                                                    return_perCellDF = return_perCellData)
  
  
  # get tLLR score under updated cell type
  # `getScoreCellType_gene` function returns a named vector with score of given cell type in values and transcript_id in names
  score_transVector <- getScoreCellType_gene(score_GeneMatrix = score_GeneMatrix, 
                                             transcript_df = post_reseg_results$updated_transDF, 
                                             transID_coln = transID_coln,
                                             transGene_coln = transGene_coln,
                                             celltype_coln = 'updated_celltype')   
  post_reseg_results$updated_transDF[['score_updated_celltype']] <- score_transVector[post_reseg_results$updated_transDF[[transID_coln]]]
  post_reseg_results$updated_transDF <- post_reseg_results$updated_transDF[!is.na(post_reseg_results$updated_transDF[['score_updated_celltype']]), ]
  rm(score_transVector)
  
  # get perCell data
  if(return_perCellData){
    #### add in type of resegmentation action applied to each cells ----
    # cells altered by fastReseg pipeline
    altered_cells <- list(
      # original cell_ID with transcripts got discarded
      oriCells_trimmed = unique(sapply(strsplit(reseg_actions[['cells_to_discard']],"_g"),"[[",1)), 
      # original cell_ID with transcripts got split, the split groups were kept as new cells
      oriCells_split = unique(sapply(strsplit(reseg_actions[['cells_to_keep']],"_g"),"[[",1)), 
      # new cell_ID that received merge cells
      updatedCells_merged = unique(reseg_actions[['cells_to_update']]), 
      # new cell_ID that got split and kept as separate new cells
      updatedCells_kept = unique(reseg_actions[['cells_to_keep']]))
    
    # mark each cell with the type of resegmetaion actions applied to them
    post_reseg_results$perCell_DT[['reSeg_action']] <- 'none'
    
    tmp_idx <- which(post_reseg_results$perCell_DT[['updated_cellID']] %in% c(altered_cells[['oriCells_trimmed']], 
                                                                              altered_cells[['oriCells_split']]))
    post_reseg_results$perCell_DT[['reSeg_action']][tmp_idx] <- 'trim'
    
    tmp_idx <- which(post_reseg_results$perCell_DT[['updated_cellID']] %in% altered_cells[['updatedCells_kept']])
    post_reseg_results$perCell_DT[['reSeg_action']][tmp_idx] <- 'new'
    
    # some cells might have no change in cellID but still presence in this category, 
    # if they have been flagged for evaluation but does not satisfy the rules of merging 
    tmp_idx <- which(post_reseg_results$perCell_DT[['updated_cellID']] %in% altered_cells[['updatedCells_merged']])
    post_reseg_results$perCell_DT[['reSeg_action']][tmp_idx] <- 'merge_or_flagged'
    
    rm(tmp_idx)
   
    
    ## impute zero value for genes not in `post_reseg_results$perCell_expression` but in `score_GeneMatrix` ----
    if(includeAllRefGenes){
      missingGenes <- setdiff(rownames(score_GeneMatrix), rownames(post_reseg_results$perCell_expression))
      if(length(missingGenes)>0){
        message(sprintf("%d genes do not present in updated transcript data.frame, impute 0 for missing genes: `%s`.", 
                        length(missingGenes), paste0(missingGenes, collapse = '`, `')))
        mockExprs <- matrix(0, nrow = length(missingGenes), 
                            ncol = ncol(post_reseg_results$perCell_expression), 
                            dimnames = list(missingGenes, 
                                            colnames(post_reseg_results$perCell_expression)))
        mockExprs <- Matrix::Matrix(mockExprs, sparse = TRUE)
        post_reseg_results$perCell_expression <- rbind(post_reseg_results$perCell_expression, 
                                                       mockExprs)
      }
    }
  }
  
  
  #### return final results ----
  outs <- list(updated_transDF = as.data.frame(post_reseg_results$updated_transDF))
  
  if(return_intermediates){
    outs[['neighborhoodDF_ToReseg']] <- neighborReSeg_df
    outs[['reseg_actions']] <- reseg_actions
  }
  
  if(return_perCellData){
    outs[['updated_perCellDT']] <- post_reseg_results$perCell_DT
    outs[['updated_perCellExprs']] <- post_reseg_results$perCell_expression
  }
  
 return(outs)
}


#' @title update_transDF_ResegActions
#' @description Update transcript data.frame based on resegmentation action, calculate the new cell type and mean per cell spatial coordinates
#' @param transcript_df the data.frame of transcript to be updated 
#' @param reseg_full_converter a named converter to update the cell ID in `transcript_df`, cell_ID in name would be converted to cell_ID in value; discard cell_ID with value = NA
#' @param score_GeneMatrix a gene x cell-type score matrix
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param celltype_coln the column name of cell type in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in transcript_df 
#' @param return_perCellDF flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
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
                                        score_GeneMatrix, 
                                        transGene_coln = 'target',
                                        cellID_coln = 'cell_ID', 
                                        celltype_coln = 'cell_type',
                                        spatLocs_colns = c("x","y","z"),
                                        return_perCellDF = TRUE){
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transGene_coln, celltype_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transGene_coln, celltype_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  # get common cells
  common_cells <- intersect(unique(names(reseg_full_converter)), unique(transcript_df[[cellID_coln]]))
  
  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            unique(transcript_df[[transGene_coln]]))
  message(sprintf("Found %d common cells and %d common genes among `names(reseg_full_converter)`, `transcript_df`, and `score_GeneMatrix`. ", 
                  length(common_cells), length(common_genes)))
  
  if(any(length(common_cells) <1, length(common_genes)<1)){
    stop("Too few common cells or genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  
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
  
  # get score for each transcripts
  transcriptGeneScore <- score_GeneMatrix[subTransDF[[transGene_coln]], ]
  
  tmp_score <- as.data.frame(transcriptGeneScore)
  tmp_score[['updated_cellID']] <- subTransDF[['updated_cellID']]
  
  tmp_score <-data.table::setDT(tmp_score)[, lapply(.SD, sum), by = 'updated_cellID'] 
  tmp_cellID <- tmp_score[['updated_cellID']]
  tmp_score[['updated_cellID']] <- NULL
  # assign cell type based on max values
  max_idx_1st <- max.col(tmp_score,ties.method="first")
  newCellTypes <- colnames(tmp_score)[max_idx_1st]
  names(newCellTypes) <- tmp_cellID
  
  # update cell type
  transcript_df[['updated_celltype']] <- transcript_df[[celltype_coln]]
  tmp_idx <- which(transcript_df[['updated_cellID']] %in% unique(cells_to_update))
  transcript_df[['updated_celltype']][tmp_idx] <- newCellTypes[transcript_df[['updated_cellID']][tmp_idx]]
  
  outputs <- list()
  
  # get per cell data frame and expression matrix
  if(return_perCellDF){
    outputs <- transDF_to_perCell_data(transcript_df = transcript_df, 
                                       transGene_coln = transGene_coln,
                                       cellID_coln = "updated_cellID",
                                       spatLocs_colns = spatLocs_colns,
                                       celltype_coln = "updated_celltype",
                                       return_cellMeta = TRUE)

  }
  
  outputs[["updated_transDF"]] <- as.data.frame(transcript_df)
  
  return(outputs)
}



#' @title transDF_to_perCell_data
#' @description get per cell expression matrix and optional metadata like mean spatial coordinates and cell types from transcript data.frame
#' @param transcript_df transcript data.frame 
#' @param transGene_coln column name of target or gene name in `transcript_df`
#' @param cellID_coln column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param celltype_coln column name for cell type in `transcript_df`
#' @param return_cellMeta flag to return per cell data.table for metadata like mean spatial coordinates and cell types
#' @return a list 
#' \describe{
#'    \item{perCell_DT}{a per cell data.table with mean spatial coordinates, cell type, return when `return_cellMeta` = TRUE}
#'    \item{perCell_expression}{a gene x cell count sparse matrix derived from transcript data.frame after resegmentation}
#' }
transDF_to_perCell_data <- function(transcript_df, 
                                    transGene_coln = "target",
                                    cellID_coln = "updated_cellID",
                                    spatLocs_colns = c("x","y","z"),
                                    celltype_coln = "updated_celltype",
                                    return_cellMeta = FALSE){
  transcript_df <- data.table::setDT(transcript_df)
  
  # get per cell expression matrix
  exprs_tmp <- transcript_df[, .SD, .SDcols = c(cellID_coln, transGene_coln)]
  exprs_tmp <- reshape2::dcast(exprs_tmp, paste0(cellID_coln, "~", transGene_coln), 
                               value.var = transGene_coln, 
                               fun.aggregate = length, fill = 0)
  rownames(exprs_tmp) <- exprs_tmp[[cellID_coln]]
  exprs_tmp <- Matrix::Matrix(as.matrix(exprs_tmp[, -1, drop = F]), sparse = TRUE)
  
  outs <- list()
  
  if(return_cellMeta){
    # get cell type
    perCell_DT <- transcript_df[, .SD[1], by = cellID_coln, .SDcols = celltype_coln]
    
    # get mean spatial coordinates
    if(!is.null(spatLocs_colns)){
      perCell_dt2 <- transcript_df[, lapply(.SD, mean), by = cellID_coln, .SDcols = spatLocs_colns] 
      perCell_DT <- merge(perCell_DT, perCell_dt2, by = cellID_coln)
    }
    
    outs[['perCell_DT']] <- perCell_DT
  }
  
  outs[['perCell_expression']] <- Matrix::t(exprs_tmp)
  
  return(outs)
  
}