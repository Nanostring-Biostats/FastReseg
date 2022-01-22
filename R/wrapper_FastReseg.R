# wrapper for resegmentation pipeline using external reference profiles and cutoffs
#' @title fastReseg_externalRef
#' @description wrapper for resegmentation pipeline using external reference profiles and cutoffs. The pipeline would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including triming to extracellular space, merging to neighbor cell or labeling as new cell. 
#' @param refProfiles A matrix of cluster profiles, genes * clusters
#' @param transcript_df the data.frame for each transcript
#' @param transID_coln the column name of transcript_ID in `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param flagCell_lrtest_cutoff the cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 10 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param leiden_args a list of arguments to pass to reticulate and Giotto:::python_leiden function, including python path, resolution, partition_type, n_iterations, set_seed, seed_number. 
#' @param flagMerge_sharedLeiden_cutoff minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event, default = 0.5 for 50% cutoff
#' @param return_intermediates flag to return intermediate outputs, including data.frame for spatial modeling statistics of each cell,  
#' @param return_perCellData flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @return a list 
#' \enumerate{
#'    \item{modStats_ToFlagCells, a data.frame for spatial modeling statistics of each cell, output of `spatialModelScoreCell_hyperplane` function, return when return_intermediates = TRUE}
#'    \item{groupDF_ToFlagTrans, data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delanuay` functions, return when return_intermediates = TRUE}
#'    \item{neighborhoodDF_ToReseg, a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when return_intermediates = TRUE}
#'    \item{reseg_actions, a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when return_intermediates = TRUE}
#'    \item{updated_transDF, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{updated_perCellDT, a per cell data.table with mean sptial coordinates and new cell type after resegmentation, return when return_perCellData = TRUE}
#'    \item{updated_perCellExprs, a gene x cell count table for updated transcript data.frame after resegmentation, return when return_perCellData = TRUE}
#' }
#' @examples 
#' data(refProfiles)
#' data(transcriptDF)
#' data(baselineCT)
#' extracellular_cellID <- transcriptDF[which(transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' score_baseline <- baselineCT[['span_score']][,"25%"]
#' lowerCutoff_transNum  <- baselineCT[['span_transNum']][,"25%"]
#' higherCutoff_transNum  <- baselineCT[['span_transNum']][,"50%"]
#' final_res <- fastReseg_externalRef(refProfiles = refProfiles, 
#'                                    transcript_df = transcriptDF,
#'                                    extracellular_cellID = extracellular_cellID, 
#'                                    molecular_distance_cutoff = 2.7,
#'                                    cellular_distance_cutoff = 25,
#'                                    score_baseline = score_baseline, 
#'                                    lowerCutoff_transNum = lowerCutoff_transNum, 
#'                                    higherCutoff_transNum= higherCutoff_transNum)
#' @importFrom data.table as.data.table
#' @export
fastReseg_externalRef <- function(refProfiles, 
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
                                  molecular_distance_cutoff = 2.7,
                                  cellular_distance_cutoff = NULL,
                                  score_baseline = NULL, 
                                  lowerCutoff_transNum = NULL, 
                                  higherCutoff_transNum= NULL, 
                                  leiden_args = list(python_path = "/usr/bin/python3", 
                                                     partition_type = c("RBConfigurationVertexPartition", "ModularityVertexPartition"),
                                                     resolution =1,
                                                     n_iterations = 1000,
                                                     set_seed = T,
                                                     seed_number = 1234), 
                                  flagMerge_sharedLeiden_cutoff = 0.5,
                                  return_intermediates = TRUE,
                                  return_perCellData = TRUE){
  # final results
  final_res <- list()
  #### check inputs ----
  # check format of transcript_df
  if(any(!c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided `transcript_df`, missing columns include `%s`.",
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
      stop("To define the neighborhood to consider for transcript network, `molecular_distance_cutoff` must be either NULL to use 10 times of 90% quantile of minimal molecular distance within no more than 2500 randomly chosen cells or a numeric value to define the largest molecule-to-molecule distance.")
    } else if (molecular_distance_cutoff <= 0){
      stop("`molecular_distance_cutoff` must be either `NULL` or positive number")
    } else{
      message(sprintf('Use `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    }
  } 
  
  all_celltypes <- colnames(refProfiles)
  
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
  transcript_loglik <- scoreGenesInRef(genes = common_genes, ref_profiles = meanCelltype_profiles)
  
  # tLLRv2 score, re-center on maximum per row/transcript
  tmp_max <- apply(transcript_loglik, 1, max)
  tLLRv2_geneMatrix <- sweep(transcript_loglik, 1, tmp_max, '-')
  rm(tmp_max)
  
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
  if(is.null(cellular_distance_cutoff)){
    # # get both distance cutoff with `choose_distance_cutoff` function
    # `cellular_distance_cutoff` is defined as maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact. 
    # The function calculates average 2D cell diameter from input data.frame and use 2 times of the mean cell diameter as `cellular_distance_cutoff`. 
    # `molecular_distance_cutoff` is defined as maximum molecule-to-molecule distance within connected transcript groups belonging to same source cells. 
    # When `run_molecularDist = TRUE`, the function would first randomly choose `sampleSize_cellNum` number of cells from `sampleSize_nROI`number of randomly picked ROIs
    # with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. 
    # The function would further use the 10 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. 
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
      
    } else {
      # get only cellular distance cutoff
      distCutoffs <- choose_distance_cutoff(transcript_df, 
                                            transID_coln = transID_coln,
                                            cellID_coln = cellID_coln, 
                                            spatLocs_colns = spatLocs_colns, 
                                            extracellular_cellID = NULL, 
                                            sampleSize_nROI = 10, 
                                            sampleSize_cellNum = 2500, 
                                            seed = 123, 
                                            run_molecularDist = FALSE)
    }
    
    cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
    
    rm(distCutoffs)
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
  # `spatialModelScoreCell_hyperplane` function returns a data.frame with cell in row and spatial modeling outcomes in columns
  tmp_df <- spatialModelScoreCell_hyperplane(chosen_cells = common_cells, 
                                             transcript_df = transcript_df, 
                                             cellID_coln = cellID_coln, 
                                             transID_coln = transID_coln, 
                                             score_coln = 'score_tLLRv2_maxCellType',
                                             spatLocs_colns = spatLocs_colns, 
                                             model_cutoff = flagModel_TransNum_cutoff)
  
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
  
  ## (2) use SVM~hyperplane to identify the connected transcripts group based on tLLRv2 score ----
  # SVM can separate continuous low score transcript from the rest.
  # but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
  flagged_transDF3d <- transcript_df[which(transcript_df[[cellID_coln]] %in% flagged_cells),]
  
  # `flagTranscripts_SVM` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
  tmp_df <- flagTranscripts_SVM(chosen_cells = flagged_cells,
                                score_GeneMatrix = transcript_loglik,
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
  ## (3.1) perform delaunay on SVM-flagged transcripts within flagged cells ----
  # # config for spatial network for transcripts
  # this config would be used by two functions, `groupTranscripts_Delanuay` and `decide_ReSegment_Operations_leidenCut`, 
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
    
    #### Approach 1: create Delanuay network for HMRF module
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
    S = 0,
    
    #### Approach 2: create kNN network for leiden clustering
    
    # method to create kNN network
    knn_method = "dbscan",
    # number of nearest neighbors based on physical distance
    k = 6,
    # distance cuttoff for nearest neighbors to consider for kNN network
    maximum_distance_knn = NULL
  )
  
  # `groupTranscripts_Delanuay` function returns a data.frame of connected transcripts among chosen_transcripts, 
  # with each transcript in row, the group ID for the connected transcript groups and the original cell ID, spatial coordinates in column.
  flaggedSVM_transGroupDF3d <- groupTranscripts_Delanuay(chosen_transcripts = flaggedSVM_transID3d, 
                                                         config_spatNW_transcript = config_spatNW, 
                                                         distance_cutoff = molecular_distance_cutoff,
                                                         transcript_df = flagged_transDF3d, 
                                                         cellID_coln = cellID_coln, 
                                                         transID_coln = transID_coln,
                                                         transSpatLocs_coln = spatLocs_colns)
  
  message(sprintf("SVM spatial model further identified %d cells with transcript score all in same class, exclude from transcript group analysis.", 
                  length(unique(flagged_transDF_SVM3[[cellID_coln]])) - length(unique(flaggedSVM_transGroupDF3d[[cellID_coln]]))))
  
  
  ## (3.2) generate tmp_cellID include group information and get max cell type for each group ----
  # assign group ID name for all transcripts
  # transcripts with SVM class = 1 are high-score molecules and would get `connect_group` = 0 and thus keep the original cell ID as the tmp_cellID
  # transcripts with SVM class = 0 are low-score molecules and would get `connect_group` to be same value as the group ID identified based on spatial network analysis and tmp_cellID modified from original cell_ID based on the `connect_group` value
  flagged_transDF_SVM3[['connect_group']] <- 1- as.numeric(as.character(flagged_transDF_SVM3[['SVM_class']]))
  group_converter <- flaggedSVM_transGroupDF3d[['transcript_group']] 
  names(group_converter) <- flaggedSVM_transGroupDF3d[[transID_coln]]
  
  tmp_idx <- which(flagged_transDF_SVM3[['transcript_id']] %in% flaggedSVM_transGroupDF3d[[transID_coln]])
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
  # `perCell_DT`, a per cell data.table with mean sptial coordinates and new cell type when return_perCellDF = TRUE.
  # `perCell_expression`, a gene x cell count table for updated transcript data.frame when return_perCellDF = TRUE.
  
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
  
  final_res[['updated_transDF']] <- post_reseg_results$updated_transDF
  final_res[['updated_perCellDT']] <- post_reseg_results$perCell_DT
  if(return_perCellData){
    final_res[['updated_perCellExprs']] <- post_reseg_results$perCell_expression
  }
  
  return(final_res)
  
}
