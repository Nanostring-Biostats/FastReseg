# core wrapper for resegmentation pipeline on each FOV using transcript score matrix derived from external reference profiles and preset cutoffs
#' @title fastReseg_perFOV_full_process
#' @description core wrapper for resegmentation pipeline using transcript score matrix derived from external reference profiles and preset cutoffs
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
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
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `score_GeneMatrix` such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param groupTranscripts_method use either "dbscan" or "delaunay method" to group transcripts in space (default = "dbscan")
#' @param spatialMergeCheck_method use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space (default = "leidenCut")
#' @param cutoff_spatialMerge spatial constraint on a valid merging event between two source transcript groups, default = 0.5 for 50% cutoff, set to 0 to skip spatial constraint evaluation for merging.   
#' For `spatialMergeCheck_method = "leidenCut"`, this is the minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event.
#' For `spatialMergeCheck_method = "geometryDiff"`, this is the maximum percentage of white space change upon merging of query cell and neighbor cell for a valid merging event. 
#' @param leiden_config (leidenCut) a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.  
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config (default = NULL)
#' @param return_intermediates flag to return intermediate outputs, including data.frame for spatial modeling statistics of each cell  
#' @param return_perCellData flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param includeAllRefGenes flag to include all genes in `score_GeneMatrix` in the returned `updated_perCellExprs` with missing genes of value 0 (default = FALSE)
#' @param seed_process seed for per FOV processing, used in transcript error detection and correction steps, default = NULL to skip the seed  
#' @return a list 
#' \describe{
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, return when `return_intermediates` = TRUE}
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content` function, return when `return_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE}
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' @details The pipeline would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including triming to extracellular space, merging to neighbor cell or labeling as new cell. 
#' 
#' To account for genes missing in `score_GeneMatrix` but present in input transcript data.frame, genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. To avoid significant interference from those `ctrl_genes`, it's recommended to have total counts of those genes below 1% of total counts of all genes in each cell. 
#' 
#' @examples 
#' data(example_refProfiles)
#' data(mini_transcriptDF)
#' data(example_baselineCT)
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] #' cell_ID for extracellualr transcripts
#' score_baseline <- example_baselineCT[['span_score']][,"25%"]
#' lowerCutoff_transNum  <- example_baselineCT[['span_transNum']][,"25%"]
#' higherCutoff_transNum  <- example_baselineCT[['span_transNum']][,"50%"]
#' #' calculate log-likelihood of each gene under each cell type and center the score matrix on per gene basis
#' score_GeneMatrix <- scoreGenesInRef(genes = intersect(unique(mini_transcriptDF[["target"]]), rownames(example_refProfiles)),
#'                                     ref_profiles = pmax(example_refProfiles, 1e-5))
#' # case 1: run with default methods: "dbscan" for transcript grouping, "leidenCut" for merging check 
#' final_res1 <- fastReseg_perFOV_full_process(score_GeneMatrix= score_GeneMatrix,
#'                                            transcript_df = mini_transcriptDF,
#'                                            extracellular_cellID = extracellular_cellID,
#'                                            molecular_distance_cutoff = 2.7,
#'                                            cellular_distance_cutoff = 25,
#'                                            score_baseline = score_baseline,
#'                                            lowerCutoff_transNum = lowerCutoff_transNum,
#'                                            higherCutoff_transNum= higherCutoff_transNum)
#' 
#' 
#' # case 2: run with alternative methods: "delaunay" for transcript grouping, "geometryDiff" for merging check 
#' final_res2 <- fastReseg_perFOV_full_process(score_GeneMatrix= score_GeneMatrix,
#'                                            transcript_df = mini_transcriptDF,
#'                                            extracellular_cellID = extracellular_cellID,
#'                                            molecular_distance_cutoff = 2.7,
#'                                            cellular_distance_cutoff = 25,
#'                                            score_baseline = score_baseline,
#'                                            lowerCutoff_transNum = lowerCutoff_transNum,
#'                                            higherCutoff_transNum= higherCutoff_transNum,
#'                                            groupTranscripts_method = "delaunay",
#'                                            spatialMergeCheck_method = "geometryDiff")
#' 
#' 
#' # case 3: run with default "dbscan" for transcript grouping, but apply no spatial constraint on merging
#' final_res3 <- fastReseg_perFOV_full_process(score_GeneMatrix= score_GeneMatrix,
#'                                             transcript_df = mini_transcriptDF,
#'                                             extracellular_cellID = extracellular_cellID,
#'                                             molecular_distance_cutoff = 2.7,
#'                                             cellular_distance_cutoff = 25,
#'                                             score_baseline = score_baseline,
#'                                             lowerCutoff_transNum = lowerCutoff_transNum,
#'                                             higherCutoff_transNum= higherCutoff_transNum,
#'                                             cutoff_spatialMerge = 0)
#' @importFrom data.table as.data.table
#' @export
fastReseg_perFOV_full_process <- function(score_GeneMatrix, 
                                          transcript_df, 
                                          transID_coln = 'UMI_transID',
                                          transGene_coln = "target",
                                          cellID_coln = 'UMI_cellID', 
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
                                          groupTranscripts_method = c("dbscan", "delaunay"),
                                          spatialMergeCheck_method = c("leidenCut", "geometryDiff"), 
                                          cutoff_spatialMerge = 0.5, 
                                          leiden_config = list(objective_function = "CPM",
                                                               resolution_parameter = 1,
                                                               beta = 0.01,
                                                               n_iterations = 200), 
                                          config_spatNW_transcript = NULL, 
                                          return_intermediates = TRUE,
                                          return_perCellData = TRUE, 
                                          includeAllRefGenes = FALSE,
                                          seed_process = NULL){
  
  all_celltypes <- colnames(score_GeneMatrix)
  all_genes <- rownames(score_GeneMatrix)
  
  # check if any duplicate
  if(sum(duplicated(all_celltypes), duplicated(all_genes))>0){
    stop("`score_GeneMatrix` contains duplicated values in its column or row names!")
  }
  
  # check if all the inputs are valid and calculate distance cutoff from the current transcript_df if not provided
  # also check and set values for `config_spatNW_transcript`, which is required for merging evaluation regardless `groupTranscripts_method` used
  outs <- checkAndPrepInputs_perFOV(all_celltypes = all_celltypes,
                                    all_genes = all_genes, 
                                    transcript_df = transcript_df, 
                                    transID_coln = transID_coln,
                                    transGene_coln = transGene_coln,
                                    cellID_coln = cellID_coln, 
                                    spatLocs_colns = spatLocs_colns, 
                                    extracellular_cellID = extracellular_cellID, 
                                    flagModel_TransNum_cutoff = flagModel_TransNum_cutoff, 
                                    molecular_distance_cutoff = molecular_distance_cutoff,
                                    cellular_distance_cutoff = cellular_distance_cutoff,
                                    score_baseline = score_baseline, 
                                    lowerCutoff_transNum = lowerCutoff_transNum, 
                                    higherCutoff_transNum = higherCutoff_transNum,
                                    groupTranscripts_method = groupTranscripts_method,
                                    spatialMergeCheck_method = spatialMergeCheck_method, 
                                    cutoff_spatialMerge = cutoff_spatialMerge, 
                                    leiden_config = leiden_config, 
                                    config_spatNW_transcript = config_spatNW_transcript)
  transcript_df <- outs[["transcript_df"]]
  cellular_distance_cutoff <- outs[["cellular_distance_cutoff"]]
  molecular_distance_cutoff <- outs[["molecular_distance_cutoff"]]
  config_spatNW_transcript <- outs[["config_spatNW_transcript"]]
  leiden_config <- outs[["leiden_config"]]
  rm(outs)
  
  lowerCutoff_transNum <- lowerCutoff_transNum[match(all_celltypes, names(lowerCutoff_transNum))]
  higherCutoff_transNum <- higherCutoff_transNum[match(all_celltypes, names(higherCutoff_transNum))]
  
  # final results
  final_res <- list()
  
  ## (1) Flag cells with putative segmentation errors ----
  outs <- runSegErrorEvaluation(score_GeneMatrix= score_GeneMatrix, 
                                transcript_df = transcript_df, 
                                cellID_coln = cellID_coln, 
                                transID_coln = transID_coln,
                                transGene_coln = transGene_coln,
                                spatLocs_colns = spatLocs_colns,
                                # cutoff of transcript number to do spatial modeling
                                flagModel_TransNum_cutoff = flagModel_TransNum_cutoff) 
  
  modStats_ToFlagCells <- outs [['modStats_ToFlagCells']]
  transcript_df <- outs[['transcript_df']]
  rm(outs)
  
  # if no cells being evaluated due to too few counts per cell
  if(is.null(modStats_ToFlagCells)){
    if(return_intermediates){
      final_res[['modStats_ToFlagCells']] <- NULL
    }
    
    flagged_cells <- NULL
    
  } else{
    ## (1.2) flag cells based on linear regression of tLLR, lrtest_-log10P
    modStats_ToFlagCells[['flagged']] <- (modStats_ToFlagCells[['lrtest_-log10P']] > flagCell_lrtest_cutoff )
    flagged_cells <- modStats_ToFlagCells[['UMI_cellID']][modStats_ToFlagCells[['flagged']]]
    
    message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > %.1f.", 
                    length(flagged_cells), length(flagged_cells)/nrow(modStats_ToFlagCells), flagCell_lrtest_cutoff))
    
    if(return_intermediates){
      final_res[['modStats_ToFlagCells']] <- modStats_ToFlagCells
    }
    
  }
  
  # if no flagged cells for resegment, prepare the proper outputs and return the results right away 
  if(length(flagged_cells)<1){
    message("No cells being flagged for resegmentation, no further operation is performed on this dataset.")
    
    # create dummy outputs for all except 'modStats_ToFlagCells'
    outs <- makeDummyOuts_perFOV(all_genes = all_genes, 
                                 transcript_df = transcript_df, 
                                 transID_coln = transID_coln,
                                 transGene_coln = transGene_coln,
                                 cellID_coln = cellID_coln, 
                                 spatLocs_colns = spatLocs_colns, 
                                 return_intermediates = return_intermediates,
                                 return_perCellData = return_perCellData, 
                                 includeAllRefGenes = includeAllRefGenes)
    final_res <- c(final_res, outs)
    rm(outs)
    
    return(final_res)
    
  }
  
  ## (2) Identify wrongly segmented transcript groups using dbscan ----
  groupDF_ToFlagTrans <- runTranscriptErrorDetection(chosen_cells = flagged_cells,
                                                     score_GeneMatrix = score_GeneMatrix, 
                                                     transcript_df = transcript_df, 
                                                     cellID_coln = cellID_coln, 
                                                     transID_coln = transID_coln, 
                                                     # column for transcript score in current cell segment
                                                     score_coln = 'score_tLLR_maxCellType',
                                                     spatLocs_colns = spatLocs_colns,
                                                     model_cutoff = flagModel_TransNum_cutoff, 
                                                     score_cutoff = svmClass_score_cutoff, 
                                                     svm_args = svm_args,
                                                     groupTranscripts_method = groupTranscripts_method, 
                                                     # use the external provided `molecular_distance_cutoff` instead of `auto` which is recalculated based on `transcript_df`
                                                     distance_cutoff = molecular_distance_cutoff,
                                                     config_spatNW_transcript = config_spatNW_transcript,
                                                     seed_transError = seed_process)
  if(is.null(groupDF_ToFlagTrans)){
    message("No putative contaminating transcript groups detected, no further operation is performed on this dataset.")
    
    # create dummy outputs for all except 'modStats_ToFlagCells'
    outs <- makeDummyOuts_perFOV(all_genes = all_genes, 
                                 transcript_df = transcript_df, 
                                 transID_coln = transID_coln,
                                 transGene_coln = transGene_coln,
                                 cellID_coln = cellID_coln, 
                                 spatLocs_colns = spatLocs_colns, 
                                 return_intermediates = return_intermediates,
                                 return_perCellData = return_perCellData, 
                                 includeAllRefGenes = includeAllRefGenes)
    final_res <- c(final_res, outs)
    rm(outs)
    
    return(final_res)
    
  }
  
  if(return_intermediates){
    final_res[['groupDF_ToFlagTrans']] <- groupDF_ToFlagTrans
  }
  
  #### (3) get ready for resegmentation ----
  reSeg_ready_res <- prepResegDF(transcript_df = transcript_df, 
                                 groupDF_ToFlagTrans = groupDF_ToFlagTrans,
                                 cellID_coln = cellID_coln,
                                 transID_coln = transID_coln)
  
  #### (4) re-segmentation in neighborhood ----
  # did not consider extra cellular transcripts for neighbor identification. 
  
  outs <- runSegRefinement(score_GeneMatrix = score_GeneMatrix,  
                           chosen_cells =  reSeg_ready_res[["groups_to_reseg"]], 
                           reseg_transcript_df =  reSeg_ready_res[["reseg_transcript_df"]], 
                           reseg_cellID_coln = "tmp_cellID", 
                           reseg_celltype_coln = "group_maxCellType", 
                           transID_coln = transID_coln,
                           transGene_coln = transGene_coln, 
                           transSpatLocs_coln = spatLocs_colns,
                           score_baseline = score_baseline, 
                           lowerCutoff_transNum = lowerCutoff_transNum, 
                           higherCutoff_transNum= higherCutoff_transNum, 
                           neighbor_distance_xy = cellular_distance_cutoff,
                           distance_cutoff = molecular_distance_cutoff,
                           spatialMergeCheck_method = spatialMergeCheck_method, 
                           cutoff_spatialMerge = cutoff_spatialMerge, 
                           leiden_config = leiden_config, 
                           config_spatNW_transcript = config_spatNW_transcript,
                           return_intermediates = return_intermediates, 
                           return_perCellData = return_perCellData,  
                           includeAllRefGenes = includeAllRefGenes,
                           seed_segRefine = seed_process)
  final_res <- c(final_res, outs)
  rm(outs)
  
  return(final_res)
  
}


#' @title checkAndPrepInputs_perFOV
#' @description supporting function for \code{fastReseg_perFOV_full_process}, checks and preps inputs for full resegmentation pipeline on the provided `transcript_df` and calculates distance cutoffs and set values for `config_spatNW_transcript` and `leiden_config` if not provided. 
#' @param all_celltypes vector of all cell types consider in the analysis as listed in columns of `score_GeneMatrix` in parent function
#' @param all_genes vector of all genes consider in the analysis as listed in rows of `score_GeneMatrix` in parent function
#' @param transcript_df the data.frame for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates.
#' @param transID_coln the column name of transcript_ID in `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `score_GeneMatrix` such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param groupTranscripts_method use either "dbscan" or "delaunay method" to group transcripts in space (default = "dbscan")
#' @param spatialMergeCheck_method use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space (default = "leidenCut")
#' @param cutoff_spatialMerge spatial constraint on a valid merging event between two source transcript groups, default = 0.5 for 50% cutoff, set to 0 to skip spatial constraint evaluation for merging.   
#' For `spatialMergeCheck_method = "leidenCut"`, this is the minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event.
#' For `spatialMergeCheck_method = "geometryDiff"`, this is the maximum percentage of white space change upon merging of query cell and neighbor cell for a valid merging event. 
#' @param leiden_config (leidenCut) a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.  
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config (default = NULL)
#' @return a list 
#' \describe{
#'    \item{transcript_df}{transcript data.frame ready for downstream full resgmentation pipeline}
#'    \item{cellular_distance_cutoff}{maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate.}
#'    \item{molecular_distance_cutoff}{maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate}
#'    \item{config_spatNW_transcript}{configuration list to create spatial network at transcript level}
#'    \item{leiden_config}{configuration list to do leiden clustering on spatial network at transcript level for merge event evaluation}
#' }
#' 
checkAndPrepInputs_perFOV <- function(all_celltypes,
                                      all_genes, 
                                      transcript_df, 
                                      transID_coln = 'UMI_transID',
                                      transGene_coln = "target",
                                      cellID_coln = 'UMI_cellID', 
                                      spatLocs_colns = c('x','y','z'), 
                                      extracellular_cellID = NULL, 
                                      flagModel_TransNum_cutoff = 50, 
                                      molecular_distance_cutoff = 2.7,
                                      cellular_distance_cutoff = NULL,
                                      score_baseline = NULL, 
                                      lowerCutoff_transNum = NULL, 
                                      higherCutoff_transNum= NULL,
                                      groupTranscripts_method = c("dbscan", "delaunay"),
                                      spatialMergeCheck_method = c("leidenCut", "geometryDiff"), 
                                      cutoff_spatialMerge = 0.5, 
                                      leiden_config = NULL, 
                                      config_spatNW_transcript = NULL){
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("`spatLocs_colns` must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in `transcript_df`.")
  } else {
    message(sprintf("%d Dimension of spaital coordinates are provided.", d2_or_d3))
  }
  
  # config related to transcript level error detection and correction 
  groupTranscripts_method <- match.arg(groupTranscripts_method, c("dbscan", "delaunay"))
  if(groupTranscripts_method == 'delaunay'){
    config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript, 
                                                       spat_locs = spatLocs_colns)
  }
  
  if(cutoff_spatialMerge <0 | cutoff_spatialMerge >1){
    stop(sprintf("The providied `cutoff_spatialMerge = %.3f`, must be within [0, 1].", cutoff_spatialMerge))
  } else if(cutoff_spatialMerge > 0){
    spatialMergeCheck_method <- match.arg(spatialMergeCheck_method, c("leidenCut", "geometryDiff"))
    
    if(spatialMergeCheck_method == "leidenCut"){
      leiden_config <- check_config_leiden(leiden_config)
      if(groupTranscripts_method != 'delaunay'){
        config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript, 
                                                           spat_locs = spatLocs_colns)
      }
    }
  }
  
  
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
  
  
  # check the format for all cutoff, make sure it's number and has names for all cell types used.
  for(cutoff_var in c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")){
    if(is.null(get(cutoff_var))){
      stop(sprintf("The `%s` must be provided as a named numeric vector for score cutoff under each cell type.", cutoff_var))
    } else {
      if(!any(class(get(cutoff_var)) %in% c('numeric'))){
        stop(sprintf("The provided `%s` must be a named numeric vector. ", cutoff_var))
      }
      if(length(setdiff(all_celltypes, names(get(cutoff_var))))>0){
        stop(sprintf("The provided `%s` is missing for the following cell types: `%s`.", 
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
  common_genes <- unique(transcript_df[[transGene_coln]])
  message(sprintf("%d unique genes are found in `transcript_df`.", length(common_genes)))
  common_genes <- intersect(common_genes, all_genes)
  message(sprintf("%d unique genes are shared by the provided `score_GeneMatrix` that contains cluster profiles of %d genes for %d clusters.", 
                  length(common_genes), length(all_genes), length(all_celltypes)))
  if(length(common_genes)<2){
    stop("Too few common genes to proceed. Check if `score_GeneMatrix` is a gene x cell-cluster matrix.")
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
  if(is.null(molecular_distance_cutoff) | is.null(cellular_distance_cutoff)){
    if(is.null(molecular_distance_cutoff)){
      distCutoffs <- choose_distance_cutoff(transcript_df, 
                                            transID_coln = transID_coln,
                                            cellID_coln = cellID_coln, 
                                            spatLocs_colns = spatLocs_colns, 
                                            extracellular_cellID = NULL, # already removed from transcript_df
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
  } else {
    if(!any(c("numeric", "integer") %in% class(c(molecular_distance_cutoff, cellular_distance_cutoff)))){
      stop("Must provide numeric values for `molecular_distance_cutoff` and cellular_distance_cutoff` in microns.")
    }
    message(sprintf('Use the providied `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    message(sprintf("Use the providied `cellular_distance_cutoff` = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
    
  }
  
  outs <- list(transcript_df = transcript_df, 
               cellular_distance_cutoff = cellular_distance_cutoff, 
               molecular_distance_cutoff = molecular_distance_cutoff,
               config_spatNW_transcript = config_spatNW_transcript,
               leiden_config = leiden_config)
  
  return(outs)
}



#' @title makeDummyOuts_perFOV
#' @description supporting function for \code{fastReseg_perFOV_full_process}, makes dummy outputs based on input `transcript_df` in same format as the outputs of `fastReseg_perFOV_full_process()` but returns no `modStats_ToFlagCells` data.frame, used in case of no cells being flagged for cell segmentation error.
#' @param all_genes vector of all genes consider in the analysis as listed in rows of `score_GeneMatrix` in parent function
#' @param transcript_df the data.frame for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates, as well as `tLLR_maxCellType` and `score_tLLR_maxCellType` columns from earlier step in parent function
#' @param transID_coln the column name of transcript_ID in `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param return_intermediates flag to return intermediate outputs, including data.frame for spatial modeling statistics of each cell  
#' @param return_perCellData flag to return gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param includeAllRefGenes flag to include all genes in `score_GeneMatrix` in the returned `updated_perCellExprs` with missing genes of value 0 (default = FALSE)
#' @return a list 
#' \describe{
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content` function, return when `return_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE}
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
makeDummyOuts_perFOV <- function(all_genes, 
                                 transcript_df, 
                                 transID_coln = 'UMI_transID',
                                 transGene_coln = "target",
                                 cellID_coln = 'UMI_cellID', 
                                 spatLocs_colns = c('x','y','z'), 
                                 return_intermediates = TRUE,
                                 return_perCellData = TRUE, 
                                 includeAllRefGenes = FALSE){
  if(any(!c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln, 'tLLR_maxCellType', 'score_tLLR_maxCellType') %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided `transcript_df`, missing columns include: `%s`.",
                 paste0(setdiff(c(transID_coln, transGene_coln, spatLocs_colns, cellID_coln, 'tLLR_maxCellType', 'score_tLLR_maxCellType'), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  final_res <- list()
  
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
  reseg_transcript_df[['group_maxCellType']] <- reseg_transcript_df[['tLLR_maxCellType']]
  
  # update transcript df with resegmentation outcomes
  reseg_transcript_df[['updated_cellID']] <- reseg_transcript_df[[cellID_coln]]
  reseg_transcript_df[['updated_celltype']] <- reseg_transcript_df[['tLLR_maxCellType']]
  reseg_transcript_df[['score_updated_celltype']] <- reseg_transcript_df[['score_tLLR_maxCellType']]
  
  final_res[['updated_transDF']] <- as.data.frame(reseg_transcript_df)
  
  if(return_perCellData){
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

# prep transcript data.frame for resegmentation 
#' @title prepResegDF
#' @description supporting function for \code{fastReseg_perFOV_full_process}, combine `runTranscriptErrorDetection()` output with transcript data.frame to prep for `runSegRefinement()`
#' @param transcript_df the data.frame for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates and `tLLR_maxCellType`.
#' @param groupDF_ToFlagTrans data frame for transcripts in `chosen_cells` only, with columns for `connect_group`,`tmp_cellID`,`group_maxCellType`, output of `runTranscriptErrorDetection()`
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @return a list of two elements
#' \describe{
#'    \item{reseg_transcript_df}{data.frame with transcript_id, target/geneName, x, y, cell_id for all transcript groups `tmp_cellID` and the cell type of maximum transcript scores for each transcript group `group_maxCellType`}
#'    \item{groups_to_reseg}{vector of chosen transcript groups need to be evaluate for re-segmentation}
#' }
#' @export
prepResegDF <- function(transcript_df, 
                        groupDF_ToFlagTrans,
                        cellID_coln = "UMI_cellID",
                        transID_coln = "UMI_transID"){
  
  if(!all(c(cellID_coln, transID_coln, "tLLR_maxCellType") %in% colnames(transcript_df))){
    stop(sprintf("Missing necessary columns in `transcript_df`: `%s`.", 
                 paste0(setdiff(c(cellID_coln, transID_coln, "tLLR_maxCellType"), colnames(transcript_df)), 
                        collapse = "`, `")))
  }
  if(!all(c(transID_coln, "connect_group", "tmp_cellID", "group_maxCellType") %in% colnames(groupDF_ToFlagTrans))){
    stop(sprintf("Missing necessary columns in `groupDF_ToFlagTrans`: `%s`.", 
                 paste0(setdiff(c(transID_coln, "connect_group", "tmp_cellID", "group_maxCellType"), colnames(groupDF_ToFlagTrans)), 
                        collapse = "`, `")))
  }
  
  # update the transcript_df with flagged transcript_group
  reseg_transcript_df <- merge(transcript_df, 
                               groupDF_ToFlagTrans[, c(transID_coln, 'connect_group','tmp_cellID','group_maxCellType')], 
                               by = transID_coln, all.x = TRUE)
  # fill in the missing values for unflagged cells
  tmp_idx <- which(is.na(reseg_transcript_df[['connect_group']]))
  reseg_transcript_df[['connect_group']][tmp_idx]<-rep(0, length(tmp_idx))
  reseg_transcript_df[['tmp_cellID']][tmp_idx] <- reseg_transcript_df[[cellID_coln]][tmp_idx]
  reseg_transcript_df[['group_maxCellType']][tmp_idx] <- reseg_transcript_df[['tLLR_maxCellType']][tmp_idx]
  rm(tmp_idx)
  
  ## transcript groups need neighborhood evaluation 
  groups_to_reseg <- unique(groupDF_ToFlagTrans[which(groupDF_ToFlagTrans[['connect_group']]!=0),][['tmp_cellID']])
  
  outs <- list(reseg_transcript_df = reseg_transcript_df, 
               groups_to_reseg = groups_to_reseg)
  
  return(outs)
}

