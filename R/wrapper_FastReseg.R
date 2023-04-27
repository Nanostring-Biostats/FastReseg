# core wrapper for resegmentation pipeline using external reference profiles and cutoffs
#' @title fastReseg_core_externalRef
#' @description core wrapper for resegmentation pipeline using external reference profiles and cutoffs. 
#' @param refProfiles A matrix of cluster profiles, genes X clusters
#' @param transcript_df the data.frame for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates.
#' @param transID_coln the column name of transcript_ID in `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param flagCell_lrtest_cutoff the cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param groupTranscripts_method use either 'delaunay' or 'dbscan' method to group transcripts in space (default = 'delaunay')
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param leiden_args a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.
#' @param flagMerge_sharedLeiden_cutoff minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event, default = 0.5 for 50% cutoff
#' @param return_intermediates flag to return intermediate outputs, including data.frame for spatial modeling statistics of each cell  
#' @param return_perCellData flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param includeAllRefGenes flag to include all genes in `refProfiles` in the returned `updated_perCellExprs` with missing genes of value 0 (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not present in `counts` or `refProfiles`; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @return a list 
#' \describe{
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, return when `return_intermediates` = TRUE}
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE}
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' @details The pipeline would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including triming to extracellular space, merging to neighbor cell or labeling as new cell. 
#' 
#' To account for genes missing in `refProfiles` but present in input transcript data.frame, genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. To avoid significant interference from those `ctrl_genes`, it's recommended to have total counts of those genes below 1% of total counts of all genes in each cell. 
#' 
#' @examples 
#' data(example_refProfiles)
#' data(mini_transcriptDF)
#' data(example_baselineCT)
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' score_baseline <- example_baselineCT[['span_score']][,"25%"]
#' lowerCutoff_transNum  <- example_baselineCT[['span_transNum']][,"25%"]
#' higherCutoff_transNum  <- example_baselineCT[['span_transNum']][,"50%"]
#' final_res <- fastReseg_core_externalRef(refProfiles = example_refProfiles, 
#'                                         transcript_df = mini_transcriptDF,
#'                                         extracellular_cellID = extracellular_cellID, 
#'                                         molecular_distance_cutoff = 2.7,
#'                                         cellular_distance_cutoff = 25,
#'                                         score_baseline = score_baseline, 
#'                                         lowerCutoff_transNum = lowerCutoff_transNum, 
#'                                         higherCutoff_transNum= higherCutoff_transNum)
#' @importFrom data.table as.data.table
#' @export
fastReseg_core_externalRef <- function(refProfiles, 
                                       transcript_df, 
                                       transID_coln = "transcript_id",
                                       transGene_coln = "target",
                                       cellID_coln = 'cell_ID', 
                                       spatLocs_colns = c('x','y','z'), 
                                       extracellular_cellID = NULL, 
                                       flagModel_TransNum_cutoff = 50, 
                                       flagCell_lrtest_cutoff = 5,
                                       svmClass_score_cutoff = -2, 
                                       svm_args = list(kernel = "radial", 
                                                       scale = FALSE, 
                                                       gamma = 0.4),
                                       groupTranscripts_method = c('delaunay', 'dbscan'),
                                       molecular_distance_cutoff = 2.7,
                                       cellular_distance_cutoff = NULL,
                                       score_baseline = NULL, 
                                       lowerCutoff_transNum = NULL, 
                                       higherCutoff_transNum= NULL, 
                                       leiden_args = list(objective_function = c("CPM", "modularity"),
                                                          resolution_parameter = 1,
                                                          beta = 0.01,
                                                          n_iterations = 200), 
                                       flagMerge_sharedLeiden_cutoff = 0.5,
                                       return_intermediates = TRUE,
                                       return_perCellData = TRUE, 
                                       includeAllRefGenes = FALSE, 
                                       ctrl_genes = NULL){
  
  groupTranscripts_method <- match.arg(groupTranscripts_method, c('delaunay', 'dbscan'))
    
  # final results
  final_res <- list()
  #### check inputs ----
  # check format of transcript_df
  if(any(!c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided `transcript_df`, missing columns include: `%s`.",
                 paste0(setdiff(c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  # extract the needed columns only
  transcript_df <- as.data.frame(transcript_df)
  transcript_df <- transcript_df[, c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln)]
  
  # check if cutoff for spatial modeling making sense
  if(flagModel_TransNum_cutoff < 2){
    stop(sprintf("`flagModel_TransNum_cutoff` must be no less than 2 to enable spatial modeling, current cutoff for per cell transcript number is %s", as.character(flagModel_TransNum_cutoff)))
  }
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("`spatLocs_colns` must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in `transcript_df`.")
  } else {
    message(sprintf("%d Dimension of spaital coordinates are provided.", d2_or_d3))
  }
  
  # numeric format or null
  if(!is.null(cellular_distance_cutoff)){
    if(!any(class(cellular_distance_cutoff) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for cell network, neighbor_distance_xy must be either NULL to use 2 times of average cell diameter or a numeric value to define the largest cell-to-cell distance.")
    } else if(cellular_distance_cutoff <0){
      stop("neighbor_distance_xy must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
    } 
  } 
  
  if(!is.null(molecular_distance_cutoff)){
    if(!any(class(molecular_distance_cutoff) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcript network, `molecular_distance_cutoff` must be either NULL to use 5 times of 90% quantile of minimal molecular distance within no more than 2500 randomly chosen cells or a numeric value to define the largest molecule-to-molecule distance.")
    } else if (molecular_distance_cutoff <= 0){
      stop("`molecular_distance_cutoff` must be either `NULL` or positive number")
    } else{
      message(sprintf('Use `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    }
  } 
  
  all_celltypes <- colnames(refProfiles)
  all_genes <- rownames(refProfiles)
  
  # check the format for all cutoff, make sure it's number and has names for all cell types used.
  for(cutoff_var in c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")){
    if(is.null(get(cutoff_var))){
      stop(sprintf("The `%s` must be provided as a named numeric vector for score cutoff under each cell type used in `refProfiles`.", cutoff_var))
    } else {
      if(!any(class(get(cutoff_var)) %in% c('numeric'))){
        stop(sprintf("The provided `%s` must be a named numeric vector. ", cutoff_var))
      }
      if(length(setdiff(all_celltypes, names(get(cutoff_var))))>0){
        stop(sprintf("The provided `%s` is missing for the following cell types used in `refProfiles`: `%s`.", 
                     cutoff_var, paste0(setdiff(all_celltypes, names(get(cutoff_var))), collapse="`, `")))
      }
    }
    
  }
  
  # make sure lowerCutoff is lower than higherCutoff such that the merge would overwrite the keep. 
  lowerCutoff_transNum <- lowerCutoff_transNum[match(all_celltypes, names(lowerCutoff_transNum))]
  higherCutoff_transNum <- higherCutoff_transNum[match(all_celltypes, names(higherCutoff_transNum))]
  if(any(lowerCutoff_transNum > higherCutoff_transNum)){
    stop(sprintf("`lowerCutoff_transNum` is larger than `higherCutoff_transNum` in cell type: `%s`.", 
                 paste0(names(lowerCutoff_transNum)[lowerCutoff_transNum > higherCutoff_transNum], collapse = "`, `")))
  }
  
  #### (0) prepare transcript_df for resegmentation ----
  ## (0.1) get transcript score matrix for each gene based on reference profile ----
  common_genes <- unique(transcript_df[[transGene_coln]])
  message(sprintf("%d unique genes are found in `transcript_df`.", length(common_genes)))
  common_genes <- intersect(common_genes, rownames(refProfiles))
  message(sprintf("%d unique genes are shared by the provided `refProfiles` that contains cluster profiles of %d genes for %d clusters.", 
                  length(common_genes), nrow(refProfiles), ncol(refProfiles)))
  if(length(common_genes)<2){
    stop("Too few common genes to proceed. Check if `refProfiles` is a gene x cell-cluster matrix.")
  }
  
  ## get tLL score matrix 
  meanCelltype_profiles <- pmax(refProfiles, 1e-5)
  # tLLR score, re-center on maximum per row/transcript
  tLLRv2_geneMatrix <- scoreGenesInRef(genes = common_genes, ref_profiles = meanCelltype_profiles)
  
  # set tLLR score for control genes, same as `svmClass_score_cutoff`
  if(!is.null(ctrl_genes)){
    message(sprintf("Include the following `ctrl_genes` in analysis: `%s`.\nIt's recommended to have total counts of those genes below 1%% of total counts of all genes in each cell.", 
                    paste0(ctrl_genes, collapse = "`, `")))
    
    if(any(ctrl_genes %in% rownames(refProfiles))){
      message(sprintf("Overwrite transcript score for %d `ctrl_genes` shared with `refProfiles`: `%s`.", 
                      sum(ctrl_genes %in% rownames(refProfiles)),
                      paste0(intersect(ctrl_genes, rownames(refProfiles)), collapse = "`, `")))
      
      tLLRv2_geneMatrix <- tLLRv2_geneMatrix[!(rownames(tLLRv2_geneMatrix) %in% ctrl_genes), ]
    }
    
    tLLRv2_geneMatrix <- rbind(tLLRv2_geneMatrix, 
                               matrix(svmClass_score_cutoff, 
                                      nrow = length(ctrl_genes), ncol = ncol(tLLRv2_geneMatrix),
                                      dimnames = list(ctrl_genes, colnames(tLLRv2_geneMatrix)))
    )
    
    # include ctrl_genes 
    all_genes <- unique(c(all_genes, ctrl_genes))
    
    
  }
  
  ## (0.2) remove extracellular transcript from transcript_df ----
  common_cells <- unique(transcript_df[[cellID_coln]])
  message(sprintf("Found %d transcript records and %d cells in input `transcript_df`.", nrow(transcript_df), length(common_cells)))
  
  if(!is.null(extracellular_cellID)){
    if(length(extracellular_cellID)>1){
      common_cells <- setdiff(common_cells, extracellular_cellID)
      transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells), ]
      
      message(sprintf("%d transcript records within %d cells remain after removal of extracellular transcripts within %d cell_ID in provided `extracellular_cellID` vector.", 
                      nrow(transcript_df), length(common_cells), length(extracellular_cellID)))
    }
    
  }
  
  
  ## (0.3) get distance cutoff for molecular-to-molecular and cell-to-cell ----
  ## get distance cutoff for neighborhood at transcript level and cell level
  # # get both distance cutoff with `choose_distance_cutoff` function
  # `cellular_distance_cutoff` is defined as maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact. 
  # The function calculates average 2D cell diameter from input data.frame and use 2 times of the mean cell diameter as `cellular_distance_cutoff`. 
  # `molecular_distance_cutoff` is defined as maximum molecule-to-molecule distance within connected transcript groups belonging to same source cells. 
  # When `run_molecularDist = TRUE`, the function would first randomly choose `sampleSize_cellNum` number of cells from `sampleSize_nROI`number of randomly picked ROIs
  # with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. 
  # The function would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. 
  # This calculation is slow and is not recommended for large transcript data.frame.
  
  if(is.null(molecular_distance_cutoff)){
    distCutoffs <- choose_distance_cutoff(transcript_df, 
                                          transID_coln = transID_coln,
                                          cellID_coln = cellID_coln, 
                                          spatLocs_colns = spatLocs_colns, 
                                          extracellular_cellID = NULL, 
                                          sampleSize_nROI = 10, 
                                          sampleSize_cellNum = 2500, 
                                          seed = 123, 
                                          run_molecularDist = TRUE)
    molecular_distance_cutoff <- distCutoffs[['molecular_distance_cutoff']]
    message(sprintf('Use `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    
    # get cellular_distance_cutoff from estimation outcomes if not NULL
    if(is.null(cellular_distance_cutoff)){
      cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
      message(sprintf("Use cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
    }
    
    rm(distCutoffs)
  }else {
    # get only cellular_distance_cutoff
    if(is.null(cellular_distance_cutoff)){
      distCutoffs <- choose_distance_cutoff(transcript_df, 
                                            transID_coln = transID_coln,
                                            cellID_coln = cellID_coln, 
                                            spatLocs_colns = spatLocs_colns, 
                                            extracellular_cellID = NULL, 
                                            sampleSize_nROI = 10, 
                                            sampleSize_cellNum = 2500, 
                                            seed = 123, 
                                            run_molecularDist = FALSE)
      cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
      message(sprintf("Use cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
      rm(distCutoffs)
    }
  }
  
  
  ## (0.4) for each cell, get new cell type based on maximum score ----
  # `getCellType_maxScore` function returns a list contains element `cellType_DF`, a data.frame with cell in row, cell_ID and cell_type in column.
  tmp_df <- getCellType_maxScore(score_GeneMatrix = tLLRv2_geneMatrix, 
                                 transcript_df = transcript_df, 
                                 transID_coln = transID_coln,
                                 transGene_coln = transGene_coln,
                                 cellID_coln = cellID_coln, 
                                 return_transMatrix = FALSE)
  
  select_cellmeta <- tmp_df[['cellType_DF']]
  colnames(select_cellmeta) <- c(cellID_coln,'tLLRv2_maxCellType')
  rm(tmp_df)
  
  transcript_df <- merge(transcript_df, select_cellmeta, by = cellID_coln)
  message(sprintf("Found %d cells and assigned cell type based on the provided 'refProfiles` cluster profiles.", nrow(select_cellmeta)))
  
  
  ## (0.5) for each transcript, calculate tLLR score based on the max cell type
  # `getScoreCellType_gene` function returns a data.frame with transcript in row and "[transID_coln]" and "score_[celltype_coln]" in column for chosen cell-type
  tmp_df <- getScoreCellType_gene(score_GeneMatrix = tLLRv2_geneMatrix, 
                                  transcript_df = transcript_df, 
                                  transID_coln = transID_coln,
                                  transGene_coln = transGene_coln,
                                  celltype_coln = 'tLLRv2_maxCellType')
  transcript_df <- merge(transcript_df, tmp_df, by = transID_coln)
  rm(tmp_df)
  
  ####re-segmentation based on tLLRv2 score (1) flag cells, identify transcript groups ----
  ## (1.1) spatial modeling of tLLR score profile within each cell to identify cells with strong spatail dependency 
  # `score_cell_segmentation_error` function returns a data.frame with cell in row and spatial modeling outcomes in columns
  tmp_df <- score_cell_segmentation_error(chosen_cells = common_cells, 
                                          transcript_df = transcript_df, 
                                          cellID_coln = cellID_coln, 
                                          transID_coln = transID_coln, 
                                          score_coln = 'score_tLLRv2_maxCellType',
                                          spatLocs_colns = spatLocs_colns, 
                                          model_cutoff = flagModel_TransNum_cutoff)
  
  # if no cells being evaluated due to too few counts per cell
  if(is.null(tmp_df)){
    if(return_intermediates){
      final_res[['modStats_ToFlagCells']] <- NULL
    }
    
    flagged_cells <- NULL
    
  } else{
    
    #-log10(P)
    tmp_df[['lrtest_-log10P']] <- -log10(tmp_df[['lrtest_Pr']])
    modStats_tLLRv2_3D <- merge(select_cellmeta, tmp_df, by.x = cellID_coln, by.y = 'cell_ID')
    rm(tmp_df)
    
    if(return_intermediates){
      final_res[['modStats_ToFlagCells']] <- modStats_tLLRv2_3D
    }
    
    
    ## (1.2) flag cells based on linear regression of tLLRv2, lrtest_-log10P
    flagged_cells <- modStats_tLLRv2_3D[[cellID_coln]][which(modStats_tLLRv2_3D[['lrtest_-log10P']] > flagCell_lrtest_cutoff )]
    message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > %.1f.", 
                    length(flagged_cells), length(flagged_cells)/nrow(modStats_tLLRv2_3D), flagCell_lrtest_cutoff))
    
    
  }
    
  ## (1.3) if no flagged cells for resegment, prepare the proper outputs and return the results right away 
  if(length(flagged_cells)<1){
    message("No cells being flagged for resegmentation, no further operation is performed on this dataset.")
    
    # return NULL for some intermediates
    if(return_intermediates){
      final_res[['groupDF_ToFlagTrans']] <- NULL

      final_res[['neighborhoodDF_ToReseg']] <- NULL
      final_res[['reseg_actions']] <- list(cells_to_discard = NULL, 
                                           cells_to_update = NULL, 
                                           cells_to_keep = NULL, 
                                           reseg_full_converter = NULL)
      
    }
    
    # updated transcript DF with missing columns
    # update the transcript_df with flagged transcript_group
    reseg_transcript_df <- transcript_df
    reseg_transcript_df[['connect_group']] <- 0
    reseg_transcript_df[['tmp_cellID']] <- reseg_transcript_df[[cellID_coln]]
    reseg_transcript_df[['group_maxCellType']] <- reseg_transcript_df[['tLLRv2_maxCellType']]
    
    # update transcript df with resegmentation outcomes
    reseg_transcript_df[['updated_cellID']] <- reseg_transcript_df[[cellID_coln]]
    reseg_transcript_df[['updated_celltype']] <- reseg_transcript_df[['tLLRv2_maxCellType']]
    reseg_transcript_df[['score_updated_celltype']] <- reseg_transcript_df[['score_tLLRv2_maxCellType']]
    
    final_res[['updated_transDF']] <- reseg_transcript_df
    
    # perCellDT
    # get cell type
    perCell_DT <- data.table::setDT(reseg_transcript_df)[, .SD[1], by = updated_cellID, .SDcols = 'updated_celltype']
    
    # get mean spatial coordinates
    if(!is.null(spatLocs_colns)){
      perCell_dt2 <- reseg_transcript_df[, lapply(.SD, mean), by = 'updated_cellID', .SDcols = spatLocs_colns] 
      perCell_DT <- merge(perCell_DT, perCell_dt2, by = 'updated_cellID')
    }
    
    # mark each cell with the type of resegmetaion actions applied to them
    perCell_DT[['reSeg_action']] <- 'none'
    
    final_res[['updated_perCellDT']] <- perCell_DT 
    
    
    if(return_perCellData){
      # get per cell expression matrix
      exprs_tmp <- data.table::setDT(reseg_transcript_df)[, .SD, .SDcols = c('updated_cellID', transGene_coln)]
      exprs_tmp <- reshape2::dcast(exprs_tmp, paste0("updated_cellID~", transGene_coln), 
                                   value.var = transGene_coln, 
                                   fun.aggregate = length, fill = 0)
      rownames(exprs_tmp) <- exprs_tmp[['updated_cellID']]
      exprs_tmp <- Matrix::Matrix(as.matrix(exprs_tmp[, -1]), sparse = TRUE)
      
      perCell_expression <- Matrix::t(exprs_tmp)
      
      
      # impute zero value for missing genes not in `perCell_expression` but in `all_genes`
      if(includeAllRefGenes){
        missingGenes <- setdiff(all_genes, rownames(perCell_expression))
        
        if(length(missingGenes)>0){
          message(sprintf("%d genes do not present in updated transcript data.frame, impute 0 for missing genes: `%s`.", 
                          length(missingGenes), paste0(missingGenes, collapse = '`, `')))
          mockExprs <- matrix(0, nrow = length(missingGenes), 
                              ncol = ncol(perCell_expression), 
                              dimnames = list(missingGenes, 
                                              colnames(perCell_expression)))
          mockExprs <- Matrix::Matrix(mockExprs, sparse = TRUE)
          perCell_expression <- rbind(perCell_expression, 
                                      mockExprs)
        }
      }
      final_res[['updated_perCellExprs']] <- perCell_expression
    }
    
    return(final_res)
    
  }
  
  ## (2) use SVM~hyperplane to identify the connected transcripts group based on tLLRv2 score ----
  # SVM can separate continuous low score transcript from the rest.
  # but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
  flagged_transDF3d <- transcript_df[which(transcript_df[[cellID_coln]] %in% flagged_cells),]
  
  # `flagTranscripts_SVM` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
  tmp_df <- flagTranscripts_SVM(chosen_cells = flagged_cells,
                                score_GeneMatrix = tLLRv2_geneMatrix,
                                transcript_df = flagged_transDF3d, 
                                cellID_coln = cellID_coln, 
                                transID_coln = transID_coln, 
                                score_coln = 'score_tLLRv2_maxCellType',
                                spatLocs_colns = spatLocs_colns, 
                                model_cutoff = flagModel_TransNum_cutoff, 
                                score_cutoff = svmClass_score_cutoff, 
                                svm_args = svm_args)
  
  # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
  flagged_transDF_SVM3 <- merge(flagged_transDF3d, 
                                as.data.frame(tmp_df)[, c(transID_coln,'DecVal','SVM_class','SVM_cell_type')], 
                                by = transID_coln)
  
  message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                  length(flagged_cells) - length(unique(flagged_transDF_SVM3[[cellID_coln]])), svmClass_score_cutoff))
  rm(tmp_df)
  
  
  # flagged transcript ID, character vector
  flaggedSVM_transID3d <- flagged_transDF_SVM3[flagged_transDF_SVM3[['SVM_class']] ==0, transID_coln]
  
  
  ## (3) do network analysis on flagged transcript in vectorized operation to see if more than 1 connected group ----
  ## (3.1) perform delaunay or dbscan on SVM-flagged transcripts within flagged cells ----
  # # config for spatial network for transcripts
  # this config would be used by two functions, `groupTranscripts_Delaunay` and `decide_ReSegment_Operations_leidenCut`, 
  # below to separate transcripts that would likely be from two different source cells based on their spatial clustering.
  
  config_spatNW <- list(
    # name for spatial network (default = 'Delaunay_network' or 'kNN_network' for method = 'Delaunay' or 'kNN', respectively if NULL)
    name = 'transcript_delaunay_network',
    # which spatial dimensions to use (default = all)
    dimensions = "all",
    # which method to use to create a spatial network, choose from c("Delaunay", "kNN"). (default = Delaunay)
    method = 'Delaunay',
    # minimum number of nearest neighbors if maximum_distance != NULL; used by both Delaunay and kNN methods
    minimum_k = 0,
    
    #### config for creating Delaunay network
    # Delaunay method to use, choose from c("deldir", "delaunayn_geometry", "RTriangle")
    delaunay_method = "delaunayn_geometry",
    # distance cuttoff for nearest neighbors to consider for Delaunay network. 
    maximum_distance_delaunay = "auto",
    
    ## only for delaunay_method = "delaunayn_geometry"
    # (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (http://www.qhull.org/html/qdelaun.htm) for the available options. (default = 'Pp', do not report precision problems)
    options = "Pp",
    
    ## only for delaunay_method = "RTriangle"
    # (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
    Y = TRUE,
    # (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
    j = TRUE,
    # (RTriangle) Specifies the maximum number of added Steiner points.
    S = 0
    
  )
  
  # `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` function returns a data.frame of connected transcripts among chosen_transcripts, 
  # with each transcript in row, the group ID for the connected transcript groups and the original cell ID, spatial coordinates in column.
  
  if(groupTranscripts_method == 'delaunay'){
    flaggedSVM_transGroupDF3d <- groupTranscripts_Delaunay(chosen_transcripts = flaggedSVM_transID3d, 
                                                           config_spatNW_transcript = config_spatNW, 
                                                           distance_cutoff = molecular_distance_cutoff,
                                                           transcript_df = flagged_transDF3d, 
                                                           cellID_coln = cellID_coln, 
                                                           transID_coln = transID_coln,
                                                           transSpatLocs_coln = spatLocs_colns)
  } else if (groupTranscripts_method == 'dbscan'){
    flaggedSVM_transGroupDF3d <- groupTranscripts_dbscan(chosen_transcripts = flaggedSVM_transID3d, 
                                                         distance_cutoff = molecular_distance_cutoff,
                                                         transcript_df = flagged_transDF3d, 
                                                         cellID_coln = cellID_coln, 
                                                         transID_coln = transID_coln,
                                                         transSpatLocs_coln = spatLocs_colns)
  }
  
  
  message(sprintf("SVM spatial model further identified %d cells with transcript score all in same class, exclude from transcript group analysis.", 
                  length(unique(flagged_transDF_SVM3[[cellID_coln]])) - length(unique(flaggedSVM_transGroupDF3d[[cellID_coln]]))))
  
  
  ## (3.2) generate tmp_cellID include group information and get max cell type for each group ----
  # assign group ID name for all transcripts
  # transcripts with SVM class = 1 are high-score molecules and would get `connect_group` = 0 and thus keep the original cell ID as the tmp_cellID
  # transcripts with SVM class = 0 are low-score molecules and would get `connect_group` to be same value as the group ID identified based on spatial network analysis and tmp_cellID modified from original cell_ID based on the `connect_group` value
  flagged_transDF_SVM3[['connect_group']] <- 1- as.numeric(as.character(flagged_transDF_SVM3[['SVM_class']]))
  group_converter <- flaggedSVM_transGroupDF3d[['transcript_group']] 
  names(group_converter) <- flaggedSVM_transGroupDF3d[[transID_coln]]
  
  tmp_idx <- which(flagged_transDF_SVM3[[transID_coln]] %in% flaggedSVM_transGroupDF3d[[transID_coln]])
  flagged_transDF_SVM3[['connect_group']][tmp_idx] <- group_converter[flagged_transDF_SVM3[[transID_coln]][tmp_idx]]
  
  rm(tmp_idx, group_converter)
  
  # get new cell_id and cell type for each group
  flagged_transDF_SVM3 <- data.table::as.data.table(flagged_transDF_SVM3)
  flagged_transDF_SVM3[, tmp_cellID := ifelse(connect_group == 0, get(cellID_coln), paste0(get(cellID_coln),'_g', connect_group))]
  
  # get new cell type of each group based on maximum
  # `getCellType_maxScore` function returns a list contains element `cellType_DF`, a data.frame with cell in row, cell_ID and cell_type in column.
  tmp_df <- getCellType_maxScore(score_GeneMatrix = tLLRv2_geneMatrix, 
                                 transcript_df = flagged_transDF_SVM3, 
                                 transID_coln = transID_coln,
                                 transGene_coln = transGene_coln,
                                 cellID_coln = "tmp_cellID", 
                                 return_transMatrix = FALSE)
  colnames(tmp_df[['cellType_DF']]) <- c('tmp_cellID','group_maxCellType')
  flagged_transDF_SVM3 <- merge(flagged_transDF_SVM3, tmp_df[['cellType_DF']], by = 'tmp_cellID', all.x = TRUE)
  flagged_transDF_SVM3 <- as.data.frame(flagged_transDF_SVM3)
  rm(tmp_df)
  
  if(return_intermediates){
    final_res[['groupDF_ToFlagTrans']] <- flagged_transDF_SVM3
  }
  
  #### (4) re-segmentation in neighborhood ----
  ## (4.1) get cutoff of transcript score and transcript number from baseline
  # the baseline is provided externally, no calculation needed
  
  ## (4.2) prepare transcript_df for re-segmentation ----
  # update the transcript_df with flagged transcript_group
  reseg_transcript_df <- merge(transcript_df, 
                               flagged_transDF_SVM3[, c(transID_coln, 'connect_group','tmp_cellID','group_maxCellType')], 
                               by = transID_coln, all.x = TRUE)
  # fill in the missing values for unflagged cells
  tmp_idx <- which(is.na(reseg_transcript_df[['connect_group']]))
  reseg_transcript_df[['connect_group']][tmp_idx]<-rep(0, length(tmp_idx))
  reseg_transcript_df[['tmp_cellID']][tmp_idx] <- reseg_transcript_df[[cellID_coln]][tmp_idx]
  reseg_transcript_df[['group_maxCellType']][tmp_idx] <- reseg_transcript_df[['tLLRv2_maxCellType']][tmp_idx]
  rm(tmp_idx)
  
  ## (4.3) evaluate the neighborhood of each group for re-segmentation ----
  cells_to_use <- unique(flagged_transDF_SVM3[which(flagged_transDF_SVM3[['connect_group']]!=0),][['tmp_cellID']])
  
  
  # did not consider extra cellular transcripts for neighbor identification. 
  
  ### search within absolute distance, consider 25um in xy for cell level search and 15 pixel = 2.7um to be direct neighbor at transcript level.
  # using spatstat to locate neighbor cells and rank them by minimal molecular distance to query cell
  # `neighborhood_for_resegment_spatstat` function returns a data.frame with each cell in row and its neighborhood information in columns
  neighborReSeg_df <- neighborhood_for_resegment_spatstat(chosen_cells = cells_to_use,
                                                          score_GeneMatrix = tLLRv2_geneMatrix,
                                                          score_baseline = score_baseline,
                                                          neighbor_distance_xy = cellular_distance_cutoff,
                                                          distance_cutoff = molecular_distance_cutoff,
                                                          transcript_df = reseg_transcript_df,
                                                          cellID_coln = "tmp_cellID",
                                                          celltype_coln = "group_maxCellType",
                                                          transID_coln = transID_coln,
                                                          transGene_coln = transGene_coln,
                                                          transSpatLocs_coln = spatLocs_colns)
  
  #### (4.4) decide resegmentation operation: merge, new cell, or discard ----
  # # `decide_ReSegment_Operations_leidenCut` function returns a list containing the following 4 elements:
  # `cells_to_discard`: a vector of cell ID that should be discarded during resegmentation
  # `cells_to_update`: a named vector of cell ID whether the cell_ID in name would be replaced with cell_ID in value.
  # `cells_to_keep`: a vector of cell ID that should be kept as it is.
  # `reseg_full_converter`: a single named vector of cell ID to update the original cell ID, assign NA for cells_to_discard.
  reseg_actions <- decide_ReSegment_Operations_leidenCut(neighborhood_df = neighborReSeg_df, 
                                                         score_baseline = score_baseline, 
                                                         lowerCutoff_transNum = lowerCutoff_transNum, 
                                                         higherCutoff_transNum = higherCutoff_transNum,
                                                         config_spatNW_transcript = config_spatNW,
                                                         transcript_df = reseg_transcript_df, 
                                                         cellID_coln = "tmp_cellID",
                                                         transID_coln = transID_coln,
                                                         transSpatLocs_coln = spatLocs_colns, 
                                                         leiden_config = leiden_args, 
                                                         cutoff_sharedLeiden = flagMerge_sharedLeiden_cutoff)
  
  #### (4.5) update the transcript data.frame based on reseg_actions ----
  # update neighborReSeg_df_cleanSVM with resegmentaion actions
  neighborReSeg_df[['corrected_CellId']] <- reseg_actions[['reseg_full_converter']][neighborReSeg_df[['CellId']]]
  if(return_intermediates){
    final_res[['neighborhoodDF_ToReseg']] <- neighborReSeg_df
    final_res[['reseg_actions']] <- reseg_actions
  } 
  
  # # `update_transDF_ResegActions` function returns a list containing the following 3 elements:
  # `updated_transDF`, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter.
  # `perCell_DT`, a per cell data.table with mean spatial coordinates and new cell type when return_perCellDF = TRUE.
  # `perCell_expression`, a gene x cell count sparse matrix for updated transcript data.frame when return_perCellDF = TRUE.
  
  post_reseg_results <- update_transDF_ResegActions(transcript_df = reseg_transcript_df, 
                                                    reseg_full_converter = reseg_actions$reseg_full_converter, 
                                                    score_GeneMatrix = tLLRv2_geneMatrix, 
                                                    transGene_coln = transGene_coln, 
                                                    cellID_coln = 'tmp_cellID', 
                                                    celltype_coln = 'group_maxCellType', 
                                                    spatLocs_colns = spatLocs_colns, 
                                                    return_perCellDF = return_perCellData)
  
  
  # get tLLRv2 score under updated cell type
  # `getScoreCellType_gene` function returns a data.frame with transcript in row and "[transID_coln]" and "score_[celltype_coln]" in column for chosen cell-type
  tmp_df <- getScoreCellType_gene(score_GeneMatrix = tLLRv2_geneMatrix, 
                                  transcript_df = post_reseg_results$updated_transDF, 
                                  transID_coln = transID_coln,
                                  transGene_coln = transGene_coln,
                                  celltype_coln = 'updated_celltype')   
  
  post_reseg_results$updated_transDF <- merge(post_reseg_results$updated_transDF, tmp_df, by = transID_coln)
  rm(tmp_df)
  
  
  #### (4.6) add in type of resegmentation action applied to each cells ----
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
  
  
  
  #### (5) return final results ----
  final_res[['updated_transDF']] <- post_reseg_results$updated_transDF
  final_res[['updated_perCellDT']] <- post_reseg_results$perCell_DT
  if(return_perCellData){
    # impute zero value for missing genes not in `post_reseg_results$perCell_expression` but in `all_genes`
    if(includeAllRefGenes){
      missingGenes <- setdiff(all_genes, rownames(post_reseg_results$perCell_expression))
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
    final_res[['updated_perCellExprs']] <- post_reseg_results$perCell_expression
  }
  
  return(final_res)
  
}


# wrapper for resegmentation pipeline using internal reference profiles and cutoffs
#' @title fastReseg_internalRef
#' @description wrapper for resegmentation pipeline using internal reference profiles and cutoffs. This function first estimates proper refernece profiles and cutoffs from the provided data and then use `fastReseg_core_externalRef` function to process each transcript data.frame. 
#' @param counts Counts matrix for entire dataset, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, when NULL to automatically assign the cell cluster for each cell based on maximum transcript score of given the provided `refProfiles`
#' @param refProfiles A matrix of cluster profiles, genes X clusters, default = NULL to use external cluster assignments
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire dataset; when NULL, use the provided `transcript_df` directly
#' @param filepath_fov_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
#' @param prefix_colns the column names of annotation in `transDF_fileInfo`, to be added to the CellId as prefix when creating unique cell_ID for entire dataset; set to NULL if use the original `transID_coln` or `cellID_coln` 
#' @param fovOffset_colns the column name of coordinate offsets in 1st and 2nd dimension for each per FOV transcript data.frame in `transDF_fileInfo`, unit in micron
#' Notice that some assays like SMI has XY axes swapped between stage and each FOV such that `fovOffset_colns` should be c("stage_Y", "stage_X").
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of each `transcript_df`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of each `transcript_df`
#' @param transcript_df the data.frame of transcript level information with unique CellId, default = NULL to read from the `transDF_fileInfo`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in each `transcript_df`; when `prefix_colns` != NULL, unique transcript_id would be generated from `prefix_colns` and `transID_coln` in each `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_colns` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param flagCell_lrtest_cutoff the cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param groupTranscripts_method use either 'delaunay' or 'dbscan' method to group transcripts in space (default = 'delaunay')
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, unit in micron (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence; default = NULL to calculate from `counts` and `refProfiles` 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is; default = NULL to calculate from `counts` and `refProfiles` 
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type; default = NULL to calculate from `counts` and `refProfiles` 
#' @param imputeFlag_missingCTs flag to impute `score_baseline`, `lowerCutoff_transNum`,`higherCutoff_transNum` for cell types present in `refProfiles` but missing in the provided transcript data files or the provided baseline and cutoffs; when TRUE, the median values of existing cell types would be used as the values for missing cell types.
#' @param leiden_args a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.
#' @param flagMerge_sharedLeiden_cutoff minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event, default = 0.5 for 50% cutoff
#' @param path_to_output the file path to output folder where the resegmentation data is saved; directory would be created by function if not exists; transcript data.frame `updated_transDF` is saved as individual csv files for each FOV, while cell data of all FOVs, `updated_perCellDT` and `updated_perCellExprs`, are combined to save as .RData object.
#' @param save_intermediates flag to save intermediate outputs into output folder, including data.frame for spatial modeling statistics of each cell,  
#' @param return_perCellData flag to return and save to output folder for gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param combine_extra flag to combine original extracellular transcripts and trimmed transcripts back to the updated transcript data.frame, slow process if many transcript in each FOV file. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not present in `counts` or `refProfiles`; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @return a list 
#' \describe{
#'    \item{refProfiles}{a genes X clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{cutoffs_list}{a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, save when `save_intermediates` = TRUE}
#' }
#' @details The pipeline would first estimate mean profile for each cell cluster based on the provided cell x gene count matrix and cluster assignment for entire data set. 
#' And then, the pipeline would use the estimated cluster-specific profile as reference profiles and calculate suitable cutoff for distance search, transcript number and score in first provided per FOV transcript data frame when those cutoffs are not provided. 
#' When transcript data.frame is provided as multiple file paths in `transDF_fileInfo` data.frame, the pipeline would further perform resegmentation on individual transcript data.frame using the baseline and cutoff defined globally. 
#' For each transcript data.frame, the pipeline would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, 
#' identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, 
#' evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including triming to extracellular space, merging to neighbor cell or labeling as new cell.
#' 
#' To account for genes missing in `refProfiles` but present in input transcript data.frame, genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. To avoid significant interference from those `ctrl_genes`, it's recommended to have total counts of those genes below 1% of total counts of all genes in each cell.
#' 
#' The pipeline would save the each per FOV output as individual file in `path_to_output` directory; `updated_transDF` would be saved as csv file. 
#' When save_intermediates = TRUE, all intermediate files and resegmenation outputs of each FOV would be saved as single .RData object in 1 list object `each_segRes` containing the following elements: 
#' \describe{
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, save when `save_intermediates` = TRUE}
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, save when `save_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, save when `save_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, save when `save_intermediates` = TRUE}
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' The pipeline would also combine per cell data for all FOVs, save and return the combined data when `return_perCellData` = TRUE; `updated_perCellDT` and `updated_perCellExprs` would be save as single .RData object in `path_to_output` directory.
#' \describe{
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' @examples 
#' # get example based on example dataset
#' data("mini_transcriptDF")
#' data("ori_RawExprs")
#' data("example_refProfiles")
#' data("example_baselineCT")
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' 
#' # case #1: provide `transcript_df` directly, 
#' # do auto-calculation of distance cutoff from data while using the provided cutoffs for score and transcript numbers. 
#' res1 <- fastReseg_internalRef(counts = ori_RawExprs, 
#'                               clust = NULL, 
#'                               refProfiles = example_refProfiles,
#'                               pixel_size = 1, 
#'                               zstep_size = 1, 
#'                               transcript_df = mini_transcriptDF, 
#'                               transID_coln = "transcript_id", 
#'                               transGene_coln = "target", 
#'                               cellID_coln = "cell_ID",
#'                               spatLocs_colns = c("x","y","z"), 
#'                               extracellular_cellID = extracellular_cellID, 
#'                               molecular_distance_cutoff = NULL,
#'                               cellular_distance_cutoff = NULL,
#'                               score_baseline = example_baselineCT[['score_baseline']], 
#'                               lowerCutoff_transNum = example_baselineCT[['lowerCutoff_transNum']], 
#'                               higherCutoff_transNum= example_baselineCT[['higherCutoff_transNum']],
#'                               imputeFlag_missingCTs = TRUE,
#'                               path_to_output = "res1_directDF")
#' 
#' # case #2: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV, 
#' # do auto-calculation of score and transcript number cutoffs from gene expression matrix, `counts`, and cluster assignment of each cell, `clust`, 
#' # do auto-calculation of distance cutoff from the 1st per FOV transcript data. 
#' data("example_CellGeneExpr")
#' data("example_clust")
#' # the example individual transcript files are stored under `data` directory of this package
#' # update your path accordingly
#' # Notice that some assays like SMI has XY axes swapped between stage and each FOV;
#' # coordinates for each FOV should have units in micron
#' fileInfo_DF <- data.frame(file_path = c("data/Run4104_FOV001__complete_code_cell_target_call_coord.csv", 
#'                                         "data/Run4104_FOV002__complete_code_cell_target_call_coord.csv"), 
#'                           slide = c(1, 1), 
#'                           fov = c(1,2), 
#'                           stage_X = 1000*c(5.13, -2.701), 
#'                           stage_Y = 1000*c(-0.452, 0.081))
#' 
#' res2 <- fastReseg_internalRef(counts = example_CellGeneExpr, 
#'                               clust = example_clust, 
#'                               refProfiles = NULL,
#'                               transDF_fileInfo =fileInfo_DF, 
#'                               filepath_coln = 'file_path', 
#'                               prefix_colns = c('slide','fov'), 
#'                               fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
#'                               pixel_size = 0.18, # 0.18 micron per pixel in transcript data
#'                               zstep_size = 0.8, # 0.8 micron per z step in transcript data
#'                               transcript_df = NULL, 
#'                               transID_coln = NULL, # row index as transcript_id
#'                               transGene_coln = "target", 
#'                               cellID_coln = "CellId",
#'                               spatLocs_colns = c("x","y","z"), 
#'                               extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data 
#'                               molecular_distance_cutoff = NULL,
#'                               cellular_distance_cutoff = NULL,
#'                               score_baseline = NULL, 
#'                               lowerCutoff_transNum = NULL, 
#'                               higherCutoff_transNum= NULL,
#'                               imputeFlag_missingCTs = TRUE,
#'                               path_to_output = "res2_multiFiles")
#' 
#' # case #3: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV, 
#' # do auto-calculation of score and transcript number cutoffs from gene expression matrix, `counts`, and cluster-specific reference profiles, `refProfiles`, 
#' # use the provided distance cutoff for `molecular_distance_cutoff` but calculate the `cellular_distance_cutoff`
#' res3 <- fastReseg_internalRef(counts = example_CellGeneExpr, 
#'                               clust = NULL, 
#'                               refProfiles = example_refProfiles,
#'                               transDF_fileInfo =fileInfo_DF, 
#'                               filepath_coln = 'file_path', 
#'                               prefix_colns = c('slide','fov'), 
#'                               fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
#'                               pixel_size = 0.18, # 0.18 micron per pixel in transcript data
#'                               zstep_size = 0.8, # 0.8 micron per z step in transcript data
#'                               transcript_df = NULL, 
#'                               transID_coln = NULL, # row index as transcript_id
#'                               transGene_coln = "target", 
#'                               cellID_coln = "CellId",
#'                               spatLocs_colns = c("x","y","z"), 
#'                               extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data 
#'                               molecular_distance_cutoff = 2.7,
#'                               cellular_distance_cutoff = NULL,
#'                               score_baseline = NULL, 
#'                               lowerCutoff_transNum = NULL, 
#'                               higherCutoff_transNum= NULL,
#'                               imputeFlag_missingCTs = TRUE,
#'                               path_to_output = "res3_multiFiles")
#' @importFrom fs path
#' @importFrom Matrix rowSums
#' @export 
#' 
fastReseg_internalRef <- function(counts, 
                                  clust, 
                                  refProfiles = NULL,
                                  transDF_fileInfo = NULL, 
                                  filepath_coln = 'file_path', 
                                  prefix_colns = c('slide','fov'), 
                                  fovOffset_colns = c('stage_X','stage_Y'), 
                                  pixel_size = 0.18, 
                                  zstep_size = 0.8,
                                  transcript_df = NULL, 
                                  transID_coln = NULL,
                                  transGene_coln = "target",
                                  cellID_coln = 'CellId', 
                                  spatLocs_colns = c('x','y','z'), 
                                  extracellular_cellID = NULL, 
                                  flagModel_TransNum_cutoff = 50, 
                                  flagCell_lrtest_cutoff = 5,
                                  svmClass_score_cutoff = -2, 
                                  svm_args = list(kernel = "radial", 
                                                  scale = FALSE, 
                                                  gamma = 0.4),
                                  groupTranscripts_method = c('delaunay', 'dbscan'),
                                  molecular_distance_cutoff = 2.7,
                                  cellular_distance_cutoff = NULL,
                                  score_baseline = NULL, 
                                  lowerCutoff_transNum = NULL, 
                                  higherCutoff_transNum= NULL, 
                                  imputeFlag_missingCTs = FALSE,
                                  leiden_args = list(objective_function = c("CPM", "modularity"),
                                                     resolution_parameter = 1,
                                                     beta = 0.01,
                                                     n_iterations = 200), 
                                  flagMerge_sharedLeiden_cutoff = 0.5,
                                  path_to_output = "reSeg_res", 
                                  save_intermediates = TRUE,
                                  return_perCellData = TRUE, 
                                  combine_extra = FALSE, 
                                  ctrl_genes = NULL){
  
  groupTranscripts_method <- match.arg(groupTranscripts_method, c('delaunay', 'dbscan'))
  message(sprintf("Use %s for grouping low-score transcripts within each cell in space. ", groupTranscripts_method))
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("`spatLocs_colns` must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in `transcript_df`.")
  } else {
    message(sprintf("%d Dimension of spaital coordinates are provided.", d2_or_d3))
  }
  
  # create output directory 
  if(!file.exists(path_to_output)) dir.create(path_to_output)
  
  ## check the format of transcript data.frame provided ----
  # a data.frame by itself or a data.frame with file path to each per FOV transcript data
  if(is.null(transDF_fileInfo) & is.null(transcript_df)){
    stop("Must provdie either `transcript_df` or `transDF_fileInfo`.")
  } 
  
  # a data.frame with file path to each per FOV transcript data
  if (!is.null(transDF_fileInfo)){
    if(! "data.frame" %in% class(transDF_fileInfo)){
      stop("The provided `transDF_fileInfo` is not a data.frame.")
    }
    
    message(sprintf("%d individual per FOV files are provided in `transDF_fileInfo`.", nrow(transDF_fileInfo)))
    # check if all needed information is present
    need_colns <- c(filepath_coln, fovOffset_colns)
    
    if(is.null(prefix_colns)){
      message("`prefix_colns` = NULL, use the `transID_coln` and `cellID_coln` as they are in each per FOV transcript_df.")
    } else {
      need_colns <- c(need_colns, prefix_colns)
      message(sprintf("`transID_coln` and `cellID_coln` of each per FOV transcript_df would be re-named based on `prefix_colns` = `%s`.", 
                      paste0(prefix_colns, collapse = "`,`")))
    }
    
    if(length(fovOffset_colns)!=2){
      stop("Must provide only 2 elements for the column names of fov coorindates offset in micrometer.")
    }
    
    # check format of transcript_df
    if(any(!need_colns %in% colnames(transDF_fileInfo))){
      stop(sprintf("Not all necessary columns can be found in provided `transDF_fileInfo`, missing columns include `%s`.",
                   paste0(setdiff(need_colns, colnames(transDF_fileInfo)),
                          collapse = "`, `")))
    }
    
  } else {
    # a data.frame by itself 
    if(! "data.frame" %in% class(transcript_df)){
      stop("The provided `transcript_df` is not a data.frame.")
    }
    message(sprintf('A single `transcript_df` is provided with unique `cellID_coln` = %s and `transID_coln` = %s (use row idx if NULL).', 
                    cellID_coln, transID_coln))
    
    # create `transDF_fileInfo` for the provided `transcript_df`
    transDF_fileInfo <- data.frame(file_path = NA, 
                                   stage_X = 0, 
                                   stage_Y = 0)
    filepath_coln = 'file_path'
    fovOffset_colns = c('stage_X', 'stage_Y')
    prefix_colns = NULL
    
  }
  
  
  
  ## get baseline based on cell x gene matrix and cluster assignments ----
  if(is.null(refProfiles) & is.null(clust)){
    stop('Must provide either `refProfiles` or `clust`.')
  } 
  
  # if provided cluster assignments
  if(!is.null(clust)){
    if(!is.vector(clust)){
      stop("The provided `clust` is not a vector of cluster assignment.")
    }
    if(length(clust) != nrow(counts)){
      stop("`clust` has different length from the row number of `counts`.")
    } 
    
    # if no reference profiles but only cluster assignments
    # get mean cluster-specific profiles from cell x gene expression matrix
    if(is.null(refProfiles)){
      # ignore background, use total count per cell as scaling factor
      # `estimate_MeanProfile` function returns a matrix of cluster profiles, genes X clusters
      refProfiles <- estimate_MeanProfile( counts = as.matrix(counts), 
                                           clust = as.character(clust), 
                                           s = Matrix::rowSums(as.matrix(counts)), 
                                           bg = rep(0, nrow(counts)))
    }
    
    # # `get_baselineCT` function gets cluster-specific quantile distribution of transcript number and per cell per molecule transcript score in the provided cell x gene expression matrix based on the reference profiles and cell cluster assignment. 
    # # The function returns a list containing the following elements:
    # span_score, a matrix of average transcript tLLR score per molecule per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns
    # span_transNum, a matrix of transcript number per cell for each distinct cell types in row, percentile at (0%, 25%, 50%, 75%, 100%) in columns
    # score_baseline, a named vector of 25% quantile of cluster-specific per cell transcript score, to be used as score baseline such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence
    # lowerCutoff_transNum, a named vector of 25% quantile of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that higher than the cutoff is required to keep query cell as it is
    # higherCutoff_transNum, a named vector of median value of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
    # clust_used,  a named vector of cluster assignments for each cell used in baseline calculation, cell_ID in `counts` as name
    baselineData <- get_baselineCT(refProfiles = refProfiles, counts = counts, clust = as.character(clust))
    
  } else {
    # reference profiles exists, but no cluster assignment
    baselineData <- get_baselineCT(refProfiles = refProfiles, counts = counts, clust = NULL)
    clust = baselineData[['clust_used']]
  }
  
  ## update and check format of baseline and cutoffs ----
  # update cutoffs for baseline only when not provided in inputs
  # check if any cell types not presented in cutoffs and baselines, but presents in refProfiles
  # fill in missing values based on median across all cell types when imputeFlag_missingCTs = TRUE
  myFun_updateCutoff <- function(name_cutoff, imputeFlag){
    val_cutoff <- get(name_cutoff)
    if(is.null(val_cutoff)){
      val_cutoff <- baselineData[[name_cutoff]]
    }
    
    # cutoff exists, check if cover all cell types in `refProfiles`
    missingCTs <- setdiff(colnames(refProfiles), names(val_cutoff))
    if(length(missingCTs)>0){
      if(imputeFlag){
        old_names <- names(val_cutoff)
        if(is.null(old_names)){
          # no name for cutoff, assign fake name
          old_names <- paste0("fakeName", seq_along(val_cutoff))
          message(sprintf("The provided `%s` has no name, median value would be assigned to all cell types presented in `refProfiles`.", 
                          name_cutoff))
        }
        # impute missing values based on median across all cell types
        val_cutoff <- c(val_cutoff, rep(median(val_cutoff), length(missingCTs)))
        names(val_cutoff)<- c(old_names, missingCTs)
        
      } else {
        stop(sprintf("The provided `%s` is missing for the following cell types used in `refProfiles`: `%s`. Try `imputeFlag_missingCTs = TRUE` for automatic imputation of values for missing cell types.",
                     name_cutoff, paste0(missingCTs, collapse = "`, `")))
      }
    }
    return(val_cutoff)
  }
  
  score_baseline <- myFun_updateCutoff('score_baseline', imputeFlag_missingCTs)
  lowerCutoff_transNum <- myFun_updateCutoff('lowerCutoff_transNum', imputeFlag_missingCTs)
  higherCutoff_transNum <- myFun_updateCutoff('higherCutoff_transNum', imputeFlag_missingCTs)
  
  
  ## get distance cutoff based on 1st provided FOV data ----
  # load and prepare first FOV transcript data
  path_to_transDF <- transDF_fileInfo[[filepath_coln]][1]
  if(!is.na(path_to_transDF)){
    transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
  }
  
  # fovOffset_colns must have XY axes of stage matched to XY axes of images
  # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
  transcript_df <- myFun_fov_prep_dropOrig(each_transDF = transcript_df, 
                                           fov_centerLocs = unlist(transDF_fileInfo[1, fovOffset_colns]),
                                           prefix_vals = unlist(transDF_fileInfo[1, prefix_colns]), 
                                           pixel_size = pixel_size, 
                                           zstep_size = zstep_size,
                                           transID_coln = transID_coln,
                                           transGene_coln = transGene_coln,
                                           cellID_coln = cellID_coln, 
                                           spatLocs_colns = spatLocs_colns, 
                                           extracellular_cellID = extracellular_cellID, 
                                           drop_original = TRUE)
  
  
  # # get both distance cutoff with `choose_distance_cutoff` function
  # `cellular_distance_cutoff` is defined as maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact. 
  # The function calculates average 2D cell diameter from input data.frame and use 2 times of the mean cell diameter as `cellular_distance_cutoff`. 
  # `molecular_distance_cutoff` is defined as maximum molecule-to-molecule distance within connected transcript groups belonging to same source cells. 
  # When `run_molecularDist = TRUE`, the function would first randomly choose `sampleSize_cellNum` number of cells from `sampleSize_nROI`number of randomly picked ROIs
  # with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. 
  # The function would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. 
  # This calculation is slow and is not recommended for large transcript data.frame.
  
  # get molecular_distance_cutoff
  if(is.null(molecular_distance_cutoff)){
    distCutoffs <- choose_distance_cutoff(transcript_df[['intraC']], 
                                          transID_coln = 'UMI_transID',
                                          cellID_coln = 'UMI_cellID', 
                                          spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                          extracellular_cellID = NULL, 
                                          sampleSize_nROI = 10, 
                                          sampleSize_cellNum = 2500, 
                                          seed = 123, 
                                          run_molecularDist = TRUE)
    molecular_distance_cutoff <- distCutoffs[['molecular_distance_cutoff']]
    message(sprintf('Use `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    
    # get cellular_distance_cutoff from estimation outcomes if not NULL
    if(is.null(cellular_distance_cutoff)){
      cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
      message(sprintf("Use cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
    }
    
    rm(distCutoffs)
    
  } else {
    # get only cellular_distance_cutoff
    if(is.null(cellular_distance_cutoff)){
      distCutoffs <- choose_distance_cutoff(transcript_df[['intraC']], 
                                            transID_coln = 'UMI_transID',
                                            cellID_coln = 'UMI_cellID', 
                                            spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                            extracellular_cellID = NULL, 
                                            sampleSize_nROI = 10, 
                                            sampleSize_cellNum = 2500, 
                                            seed = 123, 
                                            run_molecularDist = FALSE)
      
      cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
      message(sprintf("Use cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
      rm(distCutoffs)
    }
    
  }
  
  
  ## initialize list to collect each FOV outputs ----
  all_segRes <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  all_segRes[['refProfiles']] <- refProfiles
  all_segRes[['baselineData']] <- list(span_score = baselineData[['span_score']], 
                                       span_transNum = baselineData[['span_transNum']])
  all_segRes[['cutoffs_list']] <- list(score_baseline = score_baseline, 
                                       lowerCutoff_transNum = lowerCutoff_transNum, 
                                       higherCutoff_transNum = higherCutoff_transNum, 
                                       cellular_distance_cutoff = cellular_distance_cutoff, 
                                       molecular_distance_cutoff = molecular_distance_cutoff)
  
  rm(baselineData)
  
  # holder for perCell data and intermediate files that would be returned
  if(return_perCellData){
    all_segRes[['updated_perCellDT']] <- list()
    all_segRes[['updated_perCellExprs']] <- list()
  }
  
  if(save_intermediates){
    # `reseg_actions` would be combined for all FOVs before return
    all_segRes[['reseg_actions']] <- list(cells_to_discard = list(), 
                                          cells_to_update = list(), 
                                          cells_to_keep = list(), 
                                          reseg_full_converter = list())
  }
  
  # warning for special treatment on `ctrl_genes`
  if(!is.null(ctrl_genes)){
    message(sprintf("Include the following `ctrl_genes` in analysis: `%s`.\nIt's recommended to have total counts of those genes below 1%% of total counts of all genes in each cell.", 
                    paste0(ctrl_genes, collapse = "`, `")))
    
    if(any(ctrl_genes %in% rownames(refProfiles))){
      message(sprintf("Overwrite transcript score for %d `ctrl_genes` shared with `refProfiles`: `%s`.", 
                      sum(ctrl_genes %in% rownames(refProfiles)),
                      paste0(intersect(ctrl_genes, rownames(refProfiles)), collapse = "`, `")))
    }
    
    all_segRes[['ctrl_genes']] <- ctrl_genes
  }
  
  ## apply the fixed cutoffs settings to individual FOVs for resegmentation ----
  
  # save `updated_transDF` into csv file and intermeidates file into .RData for each FOV
  # but combine perCell data from all FOVs to return and save as single .RData object
  
  # function for processing each FOV for resegmentation
  # return `each_segRes` when complete, save results to disk
  myFun_reseg_eachFOV <- function(idx){
    if (idx !=1){
      # load and prep each FOV data 
      path_to_transDF <- transDF_fileInfo[[filepath_coln]][idx]
      if(!is.na(path_to_transDF)){
        transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
      }
      # fovOffset_colns must have XY axes of stage matched to XY axes of images
      # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
      transcript_df <- myFun_fov_prep_dropOrig(each_transDF = transcript_df, 
                                               fov_centerLocs = unlist(transDF_fileInfo[idx, fovOffset_colns]),
                                               prefix_vals = unlist(transDF_fileInfo[idx, prefix_colns]), 
                                               pixel_size = pixel_size, 
                                               zstep_size = zstep_size,
                                               transID_coln = transID_coln,
                                               transGene_coln = transGene_coln,
                                               cellID_coln = cellID_coln, 
                                               spatLocs_colns = spatLocs_colns, 
                                               extracellular_cellID = extracellular_cellID, 
                                               drop_original = TRUE)
    }
    
    # resegment current FOV
    message(sprintf("\n##############\nProcessing file `%d`: %s\n\n\n",
                    idx, path_to_transDF))
    timestamp()
    if(!is.null(transcript_df[['extraC']])){
      message(sprintf("Exclude %d extracellular transcripts from downstream, %.4f of total molecules.\n\n", 
                      nrow(transcript_df[['extraC']]), 
                      nrow(transcript_df[['extraC']])/(nrow(transcript_df[['extraC']]) + nrow(transcript_df[['intraC']]))))
    }
    
    
    # # `fastReseg_core_externalRef` function is a wrapper for resegmentation pipeline using external reference profiles and cutoffs. 
    # # The function returns a list containing the following elements:
    # modStats_ToFlagCells, a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, return when `return_intermediates` = TRUE
    # groupDF_ToFlagTrans, data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE
    # neighborhoodDF_ToReseg, a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE
    # reseg_actions, a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE
    # updated_transDF, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter
    # updated_perCellDT, a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE
    # updated_perCellExprs, a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE
    
    # set `includeAllRefGenes` to TRUE, to ensure return of all genes in `refProfiles` to `each_segRes[['updated_perCellExprs']]`
    # missing genes would have imputed value of zero
    
    # To account for genes missing in `refProfiles` but present in input transcript data.frame, 
    # genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. 
    # It's recommended to have total counts of those genes below 1% of total counts of all genes in each cell to avoid significant interference from those `ctrl_genes`,  
    each_segRes <- fastReseg_core_externalRef(refProfiles = refProfiles,
                                              transcript_df = transcript_df[['intraC']],
                                              extracellular_cellID = NULL,
                                              molecular_distance_cutoff = molecular_distance_cutoff,
                                              cellular_distance_cutoff = cellular_distance_cutoff,
                                              score_baseline = score_baseline,
                                              lowerCutoff_transNum = lowerCutoff_transNum,
                                              higherCutoff_transNum= higherCutoff_transNum,  
                                              transID_coln = "UMI_transID",
                                              transGene_coln = "target",
                                              cellID_coln = 'UMI_cellID', 
                                              spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                              flagModel_TransNum_cutoff = flagModel_TransNum_cutoff, 
                                              flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
                                              svmClass_score_cutoff = svmClass_score_cutoff, 
                                              svm_args = svm_args,
                                              groupTranscripts_method = groupTranscripts_method, 
                                              leiden_args = leiden_args, 
                                              flagMerge_sharedLeiden_cutoff = flagMerge_sharedLeiden_cutoff,
                                              return_intermediates = save_intermediates,
                                              return_perCellData = return_perCellData, 
                                              includeAllRefGenes = TRUE,
                                              ctrl_genes = ctrl_genes)
    
    # intracellular in original and updated segmentation
    each_segRes[['updated_transDF']][['transComp']] <- 'intraC' 
    
    # merge original transcript_df to updated results, NA for trimmed or extracellular transcripts
    # slow process when large transcript data.frame for given FOV
    if(combine_extra){
      # flag transcripts that got trimmed during resegmentation
      trimmed_transIDs <- setdiff(transcript_df[['intraC']][['UMI_transID']], 
                                  each_segRes[['updated_transDF']][['UMI_transID']])
      if(length(trimmed_transIDs)>0){
        message(sprintf("\n\nTrim %d transcripts during resegmentation, %.4f of all intracellular molecules.\nCombine extracellular and trimmed transcripts to the updated transcript data.frame.\n", 
                        length(trimmed_transIDs), 
                        length(trimmed_transIDs)/nrow(transcript_df[['intraC']])))
        trimmed_transDF <- transcript_df[['intraC']][transcript_df[['intraC']][['UMI_transID']] %in% trimmed_transIDs, ]
        trimmed_transDF[['transComp']] <- 'trimmed'
      }
      
      
      # combined trimmed transcripts with original extracellular
      if(!is.null(transcript_df[['extraC']])){
        transcript_df[['extraC']][['transComp']] <- 'extraC'
        if(length(trimmed_transIDs)>0){
          trimmed_transDF <- rbind(trimmed_transDF, transcript_df[['extraC']])
        } else {
          trimmed_transDF <- transcript_df[['extraC']]
        }
        
      } 
      
      
      # combine all transcript, fill missing value as NA
      if(any(length(trimmed_transIDs)>0, !is.null(transcript_df[['extraC']]))){
        colns_segRes_only <- setdiff(colnames(each_segRes[['updated_transDF']]), 
                                     colnames(trimmed_transDF))
        tmp_data <- matrix(NA, nrow = nrow(trimmed_transDF), ncol = length(colns_segRes_only))
        colnames(tmp_data) <- colns_segRes_only
        tmp_data <- as.data.frame(tmp_data)
        trimmed_transDF <- cbind(trimmed_transDF, tmp_data)
        each_segRes[['updated_transDF']] <- rbind(each_segRes[['updated_transDF']], trimmed_transDF)
        
        rm(trimmed_transIDs, colns_segRes_only, tmp_data, trimmed_transDF)
      }
    } else if(nrow(each_segRes[['updated_transDF']]) != nrow(transcript_df[['intraC']])){
      message(sprintf("\n\nTrim %d transcripts during resegmentation, %.4f of all intracellular molecules.\nThe updated transcript data.frame contains NO extracellular or trimmed transcripts.\n", 
                      nrow(transcript_df[['intraC']]) - nrow(each_segRes[['updated_transDF']]), 
                      1- nrow(each_segRes[['updated_transDF']])/ nrow(transcript_df[['intraC']])))
    }
    
    rm(transcript_df)
    
    # save `updated_transDF` into csv file for each FOV 
    write.csv(each_segRes[['updated_transDF']], 
              file = fs::path(path_to_output, paste0(idx, "_updated_transDF.csv")), 
              row.names = FALSE)
    
    # return only idx, perCell data and reseg_actions as a list
    res_to_return <- list(idx = idx)
    
    if(return_perCellData){
      res_to_return[['updated_perCellDT']] <- each_segRes[['updated_perCellDT']]
      res_to_return[['updated_perCellExprs']] <- each_segRes[['updated_perCellExprs']]
    }
    
    if(save_intermediates){
      res_to_return[['reseg_actions']] <- each_segRes[['reseg_actions']]
      
      # save intermediate files and all other outputs for current FOVs as single .RData
      save(each_segRes, 
           file = fs::path(path_to_output, paste0(idx, "_each_segRes.RData")))
      
    }
    
    rm(each_segRes)
    
    return(res_to_return)
  }
  
  # lapply() to get a list with each element for a list of results from each FOV
  process_outputs <- lapply(seq_len(nrow(transDF_fileInfo)), myFun_reseg_eachFOV)
  
  # combine data for each FOV
  if(return_perCellData){
    # extract perCell data from nested list of process_outputs
    updated_perCellDT <- lapply(process_outputs, '[[', 'updated_perCellDT')
    updated_perCellDT <- do.call(rbind, updated_perCellDT)
    
    # per cell gene expression in gene x cell sparse matrix format, do cbind
    updated_perCellExprs <- lapply(process_outputs, '[[', 'updated_perCellExprs')
    updated_perCellExprs <- do.call(cbind, updated_perCellExprs)
    
    # save combined perCell data into .RData
    save(updated_perCellDT, 
         updated_perCellExprs, 
         file = fs::path(path_to_output, "combined_updated_perCellDT_perCellExprs.RData"))
    
    all_segRes[['updated_perCellDT']] <- updated_perCellDT
    all_segRes[['updated_perCellExprs']] <- updated_perCellExprs
    
    rm(updated_perCellDT, updated_perCellExprs)
  }
  
  if(save_intermediates){
    # `reseg_actions` would be combined for all FOVs before return
    cells_to_discard <- unlist(lapply(lapply(process_outputs, '[[', 'reseg_actions'), '[[', 'cells_to_discard'))
    cells_to_update <- unlist(lapply(lapply(process_outputs, '[[', 'reseg_actions'), '[[', 'cells_to_update'))
    cells_to_keep <- unlist(lapply(lapply(process_outputs, '[[', 'reseg_actions'), '[[', 'cells_to_keep'))
    reseg_full_converter <- unlist(lapply(lapply(process_outputs, '[[', 'reseg_actions'), '[[', 'reseg_full_converter'))
    
    
    all_segRes[['reseg_actions']] <- list(cells_to_discard = cells_to_discard,
                                          cells_to_update = cells_to_update,
                                          cells_to_keep = cells_to_keep,
                                          reseg_full_converter = reseg_full_converter)
    
    rm(cells_to_discard, cells_to_update, cells_to_keep, reseg_full_converter)
  }
  
  rm(process_outputs)
  
  return(all_segRes)
  
}





#' @title findSegmentError_allFiles
#' @description Reformat the individual transcript data.frame to have unique IDs and a global coordinate system and save into disk, then score each cell for segmentation error and flag transcripts that have low goodness-of-fit to current cells.
#' @param counts Counts matrix for entire data set, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, when NULL to automatically assign the cell cluster for each cell based on maximum transcript score of given the provided `refProfiles`
#' @param refProfiles A matrix of cluster profiles, genes X clusters, default = NULL to use external cluster assignments
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set; when NULL, use the provided `transcript_df` directly
#' @param filepath_fov_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
#' @param prefix_colns the column names of annotation in `transDF_fileInfo`, to be added to the CellId as prefix when creating unique cell_ID for entire data set; set to NULL if use the original `transID_coln` or `cellID_coln` 
#' @param fovOffset_colns the column name of coordinate offsets in 1st and 2nd dimension for each per FOV transcript data.frame in `transDF_fileInfo`, unit in micron
#' Notice that some assays like SMI has XY axes swapped between stage and each FOV such that `fovOffset_colns` should be c("stage_Y", "stage_X").
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of each `transcript_df`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of each `transcript_df`
#' @param transcript_df the data.frame of transcript level information with unique CellId, default = NULL to read from the `transDF_fileInfo`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in each `transcript_df`; when `prefix_colns` != NULL, unique transcript_id would be generated from `prefix_colns` and `transID_coln` in each `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_colns` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param flagCell_lrtest_cutoff the cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param path_to_output the file path to output folder; directory would be created by function if not exists; `flagged_transDF`, the reformatted transcript data.frame with transcripts of low goodness-of-fit flagged by` SVM_class = 0`, and `modStats_ToFlagCells`, the per cell evaluation output of segmentation error, and `classDF_ToFlagTrans`, the class assignment of transcripts within each flagged cells are saved as individual csv files for each FOV, respectively.
#' @param combine_extra flag to combine original extracellular transcripts back to the flagged transcript data.frame. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not present in `counts` or `refProfiles`; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @return a list 
#' \describe{
#'    \item{refProfiles}{a genes * clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrice in cluster * percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.}
#'    \item{combined_modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell for all cells in the data set, output of `score_cell_segmentation_error` function}
#'    \item{combined_flaggedCells}{a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV}
#' }
#' @details The function would first estimate mean profile for each cell cluster based on the provided cell x gene count matrix and cluster assignment for entire data set. 
#' And then, the function would use the estimated cluster-specific profile as reference profiles when not provided. 
#' For each transcript data.frame, the function would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, and identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile. 
#' The function would save the each per FOV output as individual file in `path_to_output` directory; `flagged_transDF`, `modStats_ToFlagCells` and `classDF_ToFlagTrans` would be saved as csv file, respectively. 
#' \describe{
#'    \item{flagged_transDF}{a transcript data.frame for each FOV, with columns for unique IDs for transcripts `UMI_transID` and cells `UMI_cellID`, for global coordiante system `x`, `y`, `z`, and for the goodness-of-fit in original cell segment `SMI_class`; the original per FOV cell ID and pixel/index-based coordinates systems are saved under columns, `CellId`, `pixel_x`, `pixel_y`, `idx_z`}
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function}
#'    \item{classDF_ToFlagTrans}{data.frame for the class assignment of transcripts within putative wrongly segmented cells, output of `flagTranscripts_SVM` functions}
#' }
#' 
#' To account for genes missing in `refProfiles` but present in input transcript data.frame, genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. To avoid significant interference from those `ctrl_genes`, it's recommended to have total counts of those genes below 1% of total counts of all genes in each cell. 
#' 
#' @examples 
#' # get example based on example dataset
#' data("mini_transcriptDF")
#' data("ori_RawExprs")
#' data("example_refProfiles")
#' data("example_baselineCT")
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] #' cell_ID for extracellualr transcripts
#' 
#' # case #1: provide `transcript_df` directly,
#' # do auto cluster assignment of each cell based on gene expression matrix, `counts`, and cluster-specific reference profiles, `refProfiles`
#' res1 <- findSegmentError_allFiles(counts = ori_RawExprs,
#'                                   clust = NULL,
#'                                   refProfiles = example_refProfiles,
#'                                   pixel_size = 1,
#'                                   zstep_size = 1,
#'                                   transcript_df = mini_transcriptDF,
#'                                   transID_coln = "transcript_id",
#'                                   transGene_coln = "target",
#'                                   cellID_coln = "cell_ID",
#'                                   spatLocs_colns = c("x","y","z"),
#'                                   extracellular_cellID = extracellular_cellID,
#'                                   path_to_output = "res1f_directDF")
#' 
#' # case #2: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV,
#' # do auto-calculation of cluster-specific reference profiles from gene expression matrix, `counts`, and cluster assignment of each cell, `clust`.
#' data("example_CellGeneExpr")
#' data("example_clust")
#' # the example individual transcript files are stored under `data` directory of this package
#' # update your path accordingly
#' # Notice that some assays like SMI has XY axes swapped between stage and each FOV;
#' # coordinates for each FOV should have units in micron
#' fileInfo_DF <- data.frame(file_path = c("data/Run4104_FOV001__complete_code_cell_target_call_coord.csv",
#'                                         "data/Run4104_FOV002__complete_code_cell_target_call_coord.csv"),
#'                           slide = c(1, 1),
#'                           fov = c(1,2),
#'                           stage_X = 1000*c(5.13, -2.701),
#'                           stage_Y = 1000*c(-0.452, 0.081))
#' 
#' res2 <- findSegmentError_allFiles(counts = example_CellGeneExpr,
#'                                   clust = example_clust,
#'                                   refProfiles = NULL,
#'                                   transDF_fileInfo =fileInfo_DF,
#'                                   filepath_coln = 'file_path',
#'                                   prefix_colns = c('slide','fov'),
#'                                   fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
#'                                   pixel_size = 0.18, # 0.18 micron per pixel in transcript data
#'                                   zstep_size = 0.8, # 0.8 micron per z step in transcript data
#'                                   transcript_df = NULL,
#'                                   transID_coln = NULL, # row index as transcript_id
#'                                   transGene_coln = "target",
#'                                   cellID_coln = "CellId",
#'                                   spatLocs_colns = c("x","y","z"),
#'                                   extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data
#'                                   path_to_output = "res2f_multiFiles")
#' 
#' @importFrom fs path
#' @importFrom Matrix rowSums
#' @export 
#' 
findSegmentError_allFiles <- function(counts, 
                                      clust, 
                                      refProfiles = NULL,
                                      transDF_fileInfo = NULL, 
                                      filepath_coln = 'file_path', 
                                      prefix_colns = c('slide','fov'), 
                                      fovOffset_colns = c('stage_X','stage_Y'), 
                                      pixel_size = 0.18, 
                                      zstep_size = 0.8,
                                      transcript_df = NULL, 
                                      transID_coln = NULL,
                                      transGene_coln = "target",
                                      cellID_coln = 'CellId', 
                                      spatLocs_colns = c('x','y','z'), 
                                      extracellular_cellID = NULL, 
                                      flagModel_TransNum_cutoff = 50, 
                                      flagCell_lrtest_cutoff = 5,
                                      svmClass_score_cutoff = -2, 
                                      svm_args = list(kernel = "radial", 
                                                      scale = FALSE, 
                                                      gamma = 0.4),
                                      path_to_output = "reSeg_res", 
                                      combine_extra = FALSE, 
                                      ctrl_genes = NULL){
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("`spatLocs_colns` must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in `transcript_df`.")
  } else {
    message(sprintf("%d Dimension of spaital coordinates are provided.", d2_or_d3))
  }
  
  # create output directory 
  if(!file.exists(path_to_output)) dir.create(path_to_output)
  
  ## check the format of transcript data.frame provided ----
  # a data.frame by itself or a data.frame with file path to each per FOV transcript data
  if(is.null(transDF_fileInfo) & is.null(transcript_df)){
    stop("Must provdie either `transcript_df` or `transDF_fileInfo`.")
  } 
  
  # a data.frame with file path to each per FOV transcript data
  if (!is.null(transDF_fileInfo)){
    if(! "data.frame" %in% class(transDF_fileInfo)){
      stop("The provided `transDF_fileInfo` is not a data.frame.")
    }
    
    message(sprintf("%d individual per FOV files are provided in `transDF_fileInfo`.", nrow(transDF_fileInfo)))
    # check if all needed information is present
    need_colns <- c(filepath_coln, fovOffset_colns)
    
    if(is.null(prefix_colns)){
      message("`prefix_colns` = NULL, use the `transID_coln` and `cellID_coln` as they are in each per FOV transcript_df.")
    } else {
      need_colns <- c(need_colns, prefix_colns)
      message(sprintf("`transID_coln` and `cellID_coln` of each per FOV transcript_df would be re-named based on `prefix_colns` = `%s`.", 
                      paste0(prefix_colns, collapse = "`,`")))
    }
    
    if(length(fovOffset_colns)!=2){
      stop("Must provide only 2 elements for the column names of fov coorindates offset in micrometer.")
    }
    
    # check format of transcript_df
    if(any(!need_colns %in% colnames(transDF_fileInfo))){
      stop(sprintf("Not all necessary columns can be found in provided `transDF_fileInfo`, missing columns include `%s`.",
                   paste0(setdiff(need_colns, colnames(transDF_fileInfo)),
                          collapse = "`, `")))
    }
    
  } else {
    # a data.frame by itself 
    if(! "data.frame" %in% class(transcript_df)){
      stop("The provided `transcript_df` is not a data.frame.")
    }
    message(sprintf('A single `transcript_df` is provided with unique `cellID_coln` = %s and `transID_coln` = %s (use row idx if NULL).', 
                    cellID_coln, transID_coln))
    
    # create `transDF_fileInfo` for the provided `transcript_df`
    transDF_fileInfo <- data.frame(file_path = NA, 
                                   stage_X = 0, 
                                   stage_Y = 0)
    filepath_coln = 'file_path'
    fovOffset_colns = c('stage_X', 'stage_Y')
    prefix_colns = NULL
    
  }
  
  
  
  ## get baseline based on cell x gene matrix and cluster assignments ----
  if(is.null(refProfiles) & is.null(clust)){
    stop('Must provide either `refProfiles` or `clust`.')
  } 
  
  # if provided cluster assignments
  if(!is.null(clust)){
    if(!is.vector(clust)){
      stop("The provided `clust` is not a vector of cluster assignment.")
    }
    if(length(clust) != nrow(counts)){
      stop("`clust` has different length from the row number of `counts`.")
    } 
    
    # if no reference profiles but only cluster assignments
    # get mean cluster-specific profiles from cell x gene expression matrix
    if(is.null(refProfiles)){
      # ignore background, use total count per cell as scaling factor
      # `estimate_MeanProfile` function returns a matrix of cluster profiles, genes * clusters
      refProfiles <- estimate_MeanProfile( counts = as.matrix(counts), 
                                           clust = as.character(clust), 
                                           s = Matrix::rowSums(as.matrix(counts)), 
                                           bg = rep(0, nrow(counts)))
    }
    
    # # `get_baselineCT` function gets cluster-specific quantile distribution of transcript number and per cell per molecule transcript score in the provided cell x gene expression matrix based on the reference profiles and cell cluster assignment. 
    # # The function returns a list containing the following elements:
    # span_score, a matrix of average transcript tLLR score per molecule per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns
    # span_transNum, a matrix of transcript number per cell for each distinct cell types in row, percentile at (0%, 25%, 50%, 75%, 100%) in columns
    # score_baseline, a named vector of 25% quantile of cluster-specific per cell transcript score, to be used as score baseline such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence
    # lowerCutoff_transNum, a named vector of 25% quantile of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that higher than the cutoff is required to keep query cell as it is
    # higherCutoff_transNum, a named vector of median value of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
    # clust_used,  a named vector of cluster assignments for each cell used in baseline calculation, cell_ID in `counts` as name
    baselineData <- get_baselineCT(refProfiles = refProfiles, counts = counts, clust = as.character(clust))
    
  } else {
    # reference profiles exists, but no cluster assignment
    baselineData <- get_baselineCT(refProfiles = refProfiles, counts = counts, clust = NULL)
    clust = baselineData[['clust_used']]
  }
  
  ## common genes between reference profiles and count matrix
  common_genes <- intersect(rownames(refProfiles), colnames(counts))
  if(length(common_genes) <2){
    stop("Too few common genes between the `refProfiles` (genes X clusters) and `counts` (cells X genes), check if correct format. ")
  }
  
  ## initialize list to collect each FOV outputs ----
  all_segRes <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  all_segRes[['refProfiles']] <- refProfiles
  all_segRes[['baselineData']] <- list(span_score = baselineData[['span_score']], 
                                       span_transNum = baselineData[['span_transNum']])
  
  rm(baselineData)
  
  
  ## process individual FOVs for preprocessing, scoring segmentation error on cell basis and flagging transcripts with low goodness-of-fit ----
  # save `flagged_transDF`, `modStats_ToFlagCells`, and `classDF_ToFlagTrans` into csv file for each FOV
  # but also combine perCell data from all FOVs to return 
  
  ## (0) get transcript score matrix for each gene based on reference profile 
  tLLRv2_geneMatrix <- scoreGenesInRef(genes = common_genes, ref_profiles = pmax(refProfiles, 1e-5))
  
  # set tLLR score for control genes, same as `svmClass_score_cutoff`
  if(!is.null(ctrl_genes)){
    message(sprintf("Include the following `ctrl_genes` in analysis: `%s`.\nIt's recommended to have total counts of those genes below 1%% of total counts of all genes in each cell.", 
                    paste0(ctrl_genes, collapse = "`, `")))
    
    all_segRes[['ctrl_genes']] <- ctrl_genes
    
    if(any(ctrl_genes %in% rownames(refProfiles))){
      message(sprintf("Overwrite transcript score for %d `ctrl_genes` shared with `refProfiles`: `%s`.", 
                      sum(ctrl_genes %in% rownames(refProfiles)),
                      paste0(intersect(ctrl_genes, rownames(refProfiles)), collapse = "`, `")))
      
      tLLRv2_geneMatrix <- tLLRv2_geneMatrix[!(rownames(tLLRv2_geneMatrix) %in% ctrl_genes), ]
    }
    
    tLLRv2_geneMatrix <- rbind(tLLRv2_geneMatrix, 
                               matrix(svmClass_score_cutoff, 
                                      nrow = length(ctrl_genes), ncol = ncol(tLLRv2_geneMatrix),
                                      dimnames = list(ctrl_genes, colnames(tLLRv2_geneMatrix)))
                               )
    
    
  }
  
  
  # function for processing each FOV 
  # return `modStats_ToFlagCells` when complete, save results to disk
  myFun_flag_eachFOV <- function(idx){
    
    ## (1) load and prep each FOV data 
    path_to_transDF <- transDF_fileInfo[[filepath_coln]][idx]
    if(!is.na(path_to_transDF)){
      transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
    }
    # fovOffset_colns must have XY axes of stage matched to XY axes of images
    # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
    transcript_df <- myFun_fov_prep_dropOrig(each_transDF = transcript_df, 
                                             fov_centerLocs = unlist(transDF_fileInfo[idx, fovOffset_colns]),
                                             prefix_vals = unlist(transDF_fileInfo[idx, prefix_colns]), 
                                             pixel_size = pixel_size, 
                                             zstep_size = zstep_size,
                                             transID_coln = transID_coln,
                                             transGene_coln = transGene_coln,
                                             cellID_coln = cellID_coln, 
                                             spatLocs_colns = spatLocs_colns, 
                                             extracellular_cellID = extracellular_cellID, 
                                             drop_original = FALSE)
    
    # processing current FOV
    message(sprintf("\n##############\nProcessing file `%d`: %s\n\n\n",
                    idx, path_to_transDF))
    timestamp()
    if(!is.null(transcript_df[['extraC']])){
      message(sprintf("Exclude %d extracellular transcripts from downstream, %.4f of total molecules.\n\n", 
                      nrow(transcript_df[['extraC']]), 
                      nrow(transcript_df[['extraC']])/(nrow(transcript_df[['extraC']]) + nrow(transcript_df[['intraC']]))))
    }
    
    
    ## (2) for each cell, get new cell type based on maximum score ----
    # `getCellType_maxScore` function returns a list contains element `cellType_DF`, a data.frame with cell in row, cell_ID and cell_type in column.
    tmp_df <- getCellType_maxScore(score_GeneMatrix = tLLRv2_geneMatrix, 
                                   transcript_df = transcript_df[['intraC']], 
                                   transID_coln = 'UMI_transID',
                                   transGene_coln = 'target',
                                   cellID_coln = 'UMI_cellID', 
                                   return_transMatrix = FALSE)
    
    select_cellmeta <- tmp_df[['cellType_DF']]
    colnames(select_cellmeta) <- c('UMI_cellID','tLLRv2_maxCellType')
    rm(tmp_df)
    
    all_cells <- select_cellmeta[['UMI_cellID']]
    
    transcript_df[['intraC']] <- merge(transcript_df[['intraC']], select_cellmeta, by = 'UMI_cellID')
    message(sprintf("Found %d cells and assigned cell type based on the provided `refProfiles` cluster profiles.", nrow(select_cellmeta)))
    
    
    ## (3) for each transcript, calculate tLLR score based on the max cell type
    # `getScoreCellType_gene` function returns a data.frame with transcript in row and "[transID_coln]" and "score_[celltype_coln]" in column for chosen cell-type
    tmp_df <- getScoreCellType_gene(score_GeneMatrix = tLLRv2_geneMatrix, 
                                    transcript_df = transcript_df[['intraC']], 
                                    transID_coln = 'UMI_transID',
                                    transGene_coln = 'target',
                                    celltype_coln = 'tLLRv2_maxCellType')
    transcript_df[['intraC']] <- merge(transcript_df[['intraC']], tmp_df, by = 'UMI_transID')
    rm(tmp_df)
    
    
    
    ## (4.1) spatial modeling of tLLR score profile within each cell to identify cells with strong spatail dependency 
    # `score_cell_segmentation_error` function returns a data.frame with cell in row and spatial modeling outcomes in columns
    tmp_df <- score_cell_segmentation_error(chosen_cells = all_cells, 
                                            transcript_df = transcript_df[['intraC']], 
                                            cellID_coln = 'UMI_cellID', 
                                            transID_coln = 'UMI_transID', 
                                            score_coln = 'score_tLLRv2_maxCellType',
                                            spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                            model_cutoff = flagModel_TransNum_cutoff)
    
    if(is.null(tmp_df)){
      # if no cells with enough transcript per cell for model evaluation
      modStats_ToFlagCells <- NULL
      flagged_cells <- NULL
      
    } else {
      #-log10(P)
      tmp_df[['lrtest_-log10P']] <- -log10(tmp_df[['lrtest_Pr']])
      modStats_ToFlagCells <- merge(select_cellmeta, tmp_df, by.x = 'UMI_cellID', by.y = 'cell_ID')
      rm(tmp_df)
      
      
      ## (4.2) flag cells based on linear regression of tLLRv2, lrtest_-log10P
      modStats_ToFlagCells[['flagged']] <- (modStats_ToFlagCells[['lrtest_-log10P']] > flagCell_lrtest_cutoff )
      flagged_cells <- modStats_ToFlagCells[['UMI_cellID']][modStats_ToFlagCells[['flagged']]]
      message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > %.1f.", 
                      length(flagged_cells), length(flagged_cells)/nrow(modStats_ToFlagCells), flagCell_lrtest_cutoff))
      
      # write into disk
      # add idx as file idx
      modStats_ToFlagCells[['file_idx']] <- idx
      write.csv(modStats_ToFlagCells, file = fs::path(path_to_output, paste0(idx, '_modStats_ToFlagCells.csv')), row.names = FALSE)
      
      
    }
    
    
    
    ## (5) use SVM~hyperplane to identify the connected transcripts group based on tLLRv2 score ----
    # SVM can separate continuous low score transcript from the rest.
    # but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
    
    if(length(flagged_cells)>0){
      classDF_ToFlagTrans <- transcript_df[['intraC']][which(transcript_df[['intraC']][['UMI_cellID']] %in% flagged_cells),]
      
      # `flagTranscripts_SVM` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
      tmp_df <- flagTranscripts_SVM(chosen_cells = flagged_cells,
                                    score_GeneMatrix = tLLRv2_geneMatrix,
                                    transcript_df = classDF_ToFlagTrans, 
                                    cellID_coln = 'UMI_cellID', 
                                    transID_coln = 'UMI_transID', 
                                    score_coln = 'score_tLLRv2_maxCellType',
                                    spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                    model_cutoff = flagModel_TransNum_cutoff, 
                                    score_cutoff = svmClass_score_cutoff, 
                                    svm_args = svm_args)
      
      # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
      classDF_ToFlagTrans <- merge(classDF_ToFlagTrans, 
                                   as.data.frame(tmp_df)[, c('UMI_transID','DecVal','SVM_class','SVM_cell_type')], 
                                   by = 'UMI_transID')
      
      message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                      length(flagged_cells) - length(unique(classDF_ToFlagTrans[['UMI_cellID']])), svmClass_score_cutoff))
      rm(tmp_df)
      
      # write into disk
      write.csv(classDF_ToFlagTrans, file = fs::path(path_to_output, paste0(idx, '_classDF_ToFlagTrans.csv')), row.names = FALSE)
      
      
      # flagged transcript ID, character vector
      flaggedSVM_transID3d <- classDF_ToFlagTrans[classDF_ToFlagTrans[['SVM_class']] ==0, 'UMI_transID']
      # assign SVM_class =0 for transcripts with low goodness-of-fit
      transcript_df[['intraC']][['SVM_class']] <- 1- as.numeric(transcript_df[['intraC']][['UMI_transID']] %in% flaggedSVM_transID3d)
    } else {
      # no cells flaggged for resegmentation
      message("No cells being flagged for resegmentation, no SVM is performed on this dataset.")
      transcript_df[['intraC']][['SVM_class']] <- 1
    }

    
    # intracellular vs extracelluar compartment 
    transcript_df[['intraC']][['transComp']] <- 'intraC' 
    
    # combine extra cellular transcript back to complete transcript data.frame
    if(all(combine_extra, 
           !is.null(transcript_df[['extraC']]))){
      transcript_df[['extraC']][['transComp']] <- 'extraC'
      
      # add dummy values for extracellular transcripts
      colns_intraC_only <- setdiff(colnames(transcript_df[['intraC']]), 
                                   colnames(transcript_df[['extraC']]))
      tmp_data <- matrix(NA, nrow = nrow(transcript_df[['extraC']]), ncol = length(colns_intraC_only))
      colnames(tmp_data) <- colns_intraC_only
      
      transcript_df[['extraC']] <- cbind(transcript_df[['extraC']], as.data.frame(tmp_data))
      
      # combine intra and extra togehter
      transcript_df <- do.call(rbind, transcript_df)
      
      rm(tmp_data, colns_intraC_only)
      
      
    } else {
      
      # keep only intracellular transcripts
      transcript_df <- transcript_df[['intraC']]
    }
    
    # save `flagged_transDF` into csv file for each FOV 
    write.csv(transcript_df, 
              file = fs::path(path_to_output, paste0(idx, "_flagged_transDF.csv")), 
              row.names = FALSE)
    
    # return only idx, perCell data and flagged_cells as a list
    res_to_return <- list(file_idx = idx, 
                          flagged_cells = flagged_cells, 
                          modStats_ToFlagCells = modStats_ToFlagCells)
    
    return(res_to_return)
  }
  
  # lapply() to get a list with each element to be a list of results from each FOV
  process_outputs <- lapply(seq_len(nrow(transDF_fileInfo)), myFun_flag_eachFOV)
  
  # combine `modStats_ToFlagCells` data for each FOV
  combined_modStats_ToFlagCells <- lapply(process_outputs, '[[', 'modStats_ToFlagCells')
  combined_modStats_ToFlagCells <- do.call(rbind, 
                                           combined_modStats_ToFlagCells[!sapply(combined_modStats_ToFlagCells, is.null)])
  
  all_segRes[['combined_modStats_ToFlagCells']] <- combined_modStats_ToFlagCells
  all_segRes[['combined_flaggedCells']] <- lapply(process_outputs, '[[', 'flagged_cells')
  
  rm(process_outputs)
  
  return(all_segRes)
  
}





## supporting functions to prepare single FOV data


#' @title myFun_fov_load
#' @description supporting function for \code{fastReseg_internalRef} and \code{findSegmentError_allFiles}, to load transcript data.frame of each FOV from file path
#' @param path_to_fov
# function to load each FOV's transcript data.frame
myFun_fov_load <- function(path_to_fov){
  if(grepl(".RData$", path_to_fov, ignore.case = TRUE)){
    each_transDF <- get(load(path_to_fov))
  } else if (grepl(".csv$", path_to_fov)){
    each_transDF <- read.csv(path_to_fov, sep = ',', header = TRUE)
    
  } else if (grepl(".txt$", path_to_fov)){
    each_transDF <- read.csv(path_to_fov, sep = '\t', header = TRUE)
    
  } else {
    stop(sprintf('The per FOV transcript data.frame must be RData, csv, or txt file. Current `path_to_fov` = %s', 
                 path_to_fov))
  }
  
  return(each_transDF)
}



#' @title myFun_fov_prep_dropOrig
#' @description supporting function for \code{fastReseg_internalRef} and \code{findSegmentError_allFiles} to get unique IDs for cells and transcripts, and convert pixel coordinates to um; when `drop_original = FALSE, the function will also return original per FOV based cell ID and coordinates under columns `CellId`, `pixel_x`, `pixel_y`, `idx_z`.
#' @param each_transDF data.frame for raw transcript
#' @param fov_centerLocs a named vector of fov 2D coordinates
#' @param prefix_vals a named vector of values to be used as prefix in `UMI_transID` and `UMI_cellID`; when `prefix_vals` != NULL, unique transcript_id would be generated from `prefix_vals` and `transID_coln` in `each_transDF` 
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of `each_transDF`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of `each_transDF`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in `each_transDF`; when `prefix_vals` != NULL, unique transcript_id would be generated from `prefix_vals` and `transID_coln` in `each_transDF`
#' @param transGene_coln the column name of target or gene name in `each_transDF`
#' @param cellID_coln the column name of cell_ID in `each_transDF`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_vals` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `each_transDF` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param drop_original flag to drop original per FOV based cell ID and coordinates under columns `CellId`, `pixel_x`, `pixel_y`, `idx_z` (default = FALSE)
#' @return a list contains transcript_df for downstream process and extracellular transcript data.frame
#' ' \describe{
#'    \item{intraC}{a data.frame for intracellular transcript, `UMI_transID` and `UMI_cellID` as column names for unique transcript_id and cell_id, `target` as column name for target gene name}
#'    \item{extraC}{a data.frame for extracellular transcript, same structure as the `intraC` data.frame in returned list}
#' }
myFun_fov_prep_dropOrig <- function(each_transDF, 
                                    fov_centerLocs, 
                                    prefix_vals = NULL, 
                                    pixel_size = 0.18, 
                                    zstep_size = 0.8,
                                    transID_coln = NULL,
                                    transGene_coln = "target",
                                    cellID_coln = 'CellId', 
                                    spatLocs_colns = c('x','y','z'), 
                                    extracellular_cellID = NULL, 
                                    drop_original = FALSE){
  
  # check format of transcript_df
  if(any(!c(transID_coln, spatLocs_colns, transGene_coln, cellID_coln) %in% colnames(each_transDF))){
    stop(sprintf("Not all necessary columns can be found in provided `transcript_df`, missing columns include `%s`.",
                 paste0(setdiff(c(transID_coln, spatLocs_colns, transGene_coln, cellID_coln), 
                                colnames(each_transDF)), collapse = "`, `")))
  }
  each_transDF <- as.data.frame(each_transDF)[, c(transID_coln, spatLocs_colns, transGene_coln, cellID_coln)]
  d2_or_d3 <- length(spatLocs_colns)
  
  # add values for prefix_colns, prefix_vals is a named vector
  if(!is.null(prefix_vals)){
    tmp_df <- matrix(rep(prefix_vals, nrow(each_transDF)), byrow = TRUE, ncol = length(prefix_vals))
    colnames(tmp_df) <- names(prefix_vals)
    each_transDF <- cbind(each_transDF, as.data.frame(tmp_df))
    rm(tmp_df)
  }
  
  
  # use row idx as transcript id if transID_coln = NULL
  if(is.null(transID_coln)){
    tmp_transID_coln = "transcript_id"
    each_transDF[[tmp_transID_coln]] <- seq_len(nrow(each_transDF))
  }else {
    tmp_transID_coln = transID_coln
  }
  
  # generate unique IDs for whole data set based on prefix_vals
  if(!is.null(prefix_vals)){
    # get unique transcript_id
    each_transDF[['UMI_transID']] <- apply(each_transDF[, c(names(prefix_vals), tmp_transID_coln)], 
                                           MARGIN = 1, 
                                           function(x) paste0(c('t', x), collapse = '_'))
    # get unique cell_ID
    each_transDF[['UMI_cellID']] <- apply(each_transDF[, c(names(prefix_vals), cellID_coln)], 
                                          MARGIN = 1, 
                                          function(x) paste0(c('c', x), collapse = '_'))
    
  } else {
    # rename the existing transcript ID columns with UMI
    colnames(each_transDF)[which(colnames(each_transDF) == tmp_transID_coln)] <- 'UMI_transID'
    each_transDF[['UMI_transID']] <- as.character(each_transDF[['UMI_transID']])
    
    # keep the original copy of cellID_coln for extracelllar transcript filtering downstream
    each_transDF[['UMI_cellID']] <- as.character(each_transDF[[cellID_coln]])
  }
  
  
  # cleanup transcript data.frame
  each_transDF <- each_transDF[, c("UMI_cellID","UMI_transID", transGene_coln, cellID_coln, spatLocs_colns)]
  if(d2_or_d3 ==2){
    orig_spatLocs_colns <- c('pixel_x', 'pixel_y')
  } else {
    orig_spatLocs_colns <- c('pixel_x', 'pixel_y', 'idx_z')
  }
  
  colnames(each_transDF) <- c("UMI_cellID","UMI_transID", "target", 'CellId', orig_spatLocs_colns)
  
  # convert coordinates to um, include the center location of each FOV values
  raw_locs <- each_transDF[, orig_spatLocs_colns]
  
  # flip y coordinates (2nd) to have images shown from top to bottom
  raw_locs[[orig_spatLocs_colns[2]]] <- 0-raw_locs[[orig_spatLocs_colns[2]]]
  # place target coordinates in reference to whole slide 
  raw_locs[, 1:2] <- sweep(raw_locs[, 1:2] * pixel_size, 2, fov_centerLocs,"+")
  
  if(d2_or_d3 ==2){
    colnames(raw_locs) <-c('x','y')
  } else {
    raw_locs[, 3] <- raw_locs[, 3]*zstep_size
    colnames(raw_locs) <-c('x','y','z')
  }
  
  # add location in global coordinate
  each_transDF <- cbind(each_transDF, raw_locs)
  
  
  # initialize results with all transcripts as intracellular
  res <- list(intraC = each_transDF, 
              extraC = NULL)
  
  # remove extracellular transcript from each_transDF
  if(!is.null(extracellular_cellID)){
    if(length(extracellular_cellID)>0){
      extraC_idx <- which(each_transDF[[cellID_coln]] %in% extracellular_cellID)
      intraC_idx <- setdiff(seq_len(nrow(each_transDF)), extraC_idx)
      
      if(length(extraC_idx)>0){
        res <- list(intraC = each_transDF[intraC_idx, ],  
                    extraC = each_transDF[extraC_idx, ])
        
        # whether to drop original coordinates and CellId
        if(drop_original){
          res[['extraC']] <- res[['extraC']][, c("UMI_cellID","UMI_transID", "target", colnames(raw_locs))]
        }
        
      }
    }
  }
  
  # whether to drop original coordinates and CellId
  if(drop_original){
    res[['intraC']] <- res[['intraC']][, c("UMI_cellID","UMI_transID", "target", colnames(raw_locs))]
  }
  
  
  return(res)
  
}



