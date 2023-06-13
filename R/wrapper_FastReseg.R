# wrapper to flag segmentation error in all files
#' @title fastReseg_flag_all_errors
#' @description Wrapper to process multiple files of one dataset for segmentation error detection in transcript level. The function reformats the individual transcript data.frame to have unique IDs and a global coordinate system and save into disk, then scores each cell for segmentation error and flags transcripts that have low goodness-of-fit to current cells.
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
#' @param flagCell_lrtest_cutoff the cutoff of `lrtest_-log10P` to identify putative wrongly segmented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param path_to_output the file path to output folder; directory would be created by function if not exists; `flagged_transDF`, the reformatted transcript data.frame with transcripts of low goodness-of-fit flagged by` SVM_class = 0`, and `modStats_ToFlagCells`, the per cell evaluation output of segmentation error, and `classDF_ToFlagTrans`, the class assignment of transcripts within each flagged cells are saved as individual csv files for each FOV, respectively.
#' @param combine_extra flag to combine original extracellular transcripts back to the flagged transcript data.frame. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not present in `counts` or `refProfiles`; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @param seed_transError seed for transcript error detection step, default = NULL to skip the seed   
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
#'    \item{flagged_transDF}{a transcript data.frame for each FOV, with columns for unique IDs of transcripts `UMI_transID` and cells `UMI_cellID`, for global coordiante system `x`, `y`, `z`, and for the goodness-of-fit in original cell segment `SMI_class`; the original per FOV cell ID and pixel/index-based coordinates systems are saved under columns, `CellId`, `pixel_x`, `pixel_y`, `idx_z`}
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function}
#'    \item{classDF_ToFlagTrans}{data.frame for the class assignment of transcripts within putative wrongly segmented cells, output of `flagTranscripts_SVM` functions}
#' }
#' 
#' To account for genes missing in `refProfiles` but present in input transcript data.frame, genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. To avoid significant interference from those `ctrl_genes`, it's recommended to have total counts of those genes below 1% of total counts of all genes in each cell. 
#' 
#' @examples 
#' data("mini_transcriptDF")
#' data("ori_RawExprs")
#' data("example_refProfiles")
#' data("example_baselineCT")
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' 
#' # case #'1: provide `transcript_df` directly,
#' # do auto cluster assignment of each cell based on gene expression matrix, `counts`, and cluster-specific reference profiles, `refProfiles`
#' res1 <- fastReseg_flag_all_errors(counts = ori_RawExprs,
#'                                   clust = NULL,
#'                                   refProfiles = example_refProfiles,
#'                                   pixel_size = 1,
#'                                   zstep_size = 1,
#'                                   transcript_df = mini_transcriptDF,
#'                                   transID_coln = "UMI_transID",
#'                                   transGene_coln = "target",
#'                                   cellID_coln = "UMI_cellID",
#'                                   spatLocs_colns = c("x","y","z"),
#'                                   extracellular_cellID = extracellular_cellID,
#'                                   path_to_output = "res1f_directDF")
#' 
#' # case #'2: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV,
#' # do auto-calculation of cluster-specific reference profiles from gene expression matrix, `counts`, and cluster assignment of each cell, `clust`.
#' data("example_CellGeneExpr")
#' data("example_clust")
#' 
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
#' res2 <- fastReseg_flag_all_errors(counts = example_CellGeneExpr,
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
#' @importFrom fs path
#' @importFrom Matrix rowSums
#' @export 
#' 
fastReseg_flag_all_errors <- function(counts, 
                                      clust = NULL, 
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
                                      ctrl_genes = NULL,
                                      seed_transError = NULL){
  
  # create output directory 
  if(!file.exists(path_to_output)) dir.create(path_to_output)
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  
  
  ## check inputs and then get baseline and cutoffs for counts
  # no need to get distance cutoff for error flagging, set values to skip calculation 
  prep_res <- runPreprocess(counts = counts, 
                            clust = clust, 
                            refProfiles = refProfiles,
                            score_baseline = NULL, 
                            lowerCutoff_transNum = NULL, 
                            higherCutoff_transNum= NULL, 
                            imputeFlag_missingCTs = TRUE,
                            ctrl_genes = ctrl_genes,
                            svmClass_score_cutoff = svmClass_score_cutoff,
                            molecular_distance_cutoff = 1,
                            cellular_distance_cutoff = 10,
                            transcript_df = NULL, 
                            transDF_fileInfo = NULL)
  
  
  ## check the format of transcript data.frame provided and load 1st fov 
  # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
  transcript_df <- checkTransFileInputsAndLoadFirst(transcript_df = transcript_df, 
                                                    transDF_fileInfo = transDF_fileInfo, 
                                                    filepath_coln = filepath_coln, 
                                                    prefix_colns = prefix_colns, 
                                                    fovOffset_colns = fovOffset_colns, 
                                                    pixel_size = pixel_size, 
                                                    zstep_size = zstep_size,
                                                    transID_coln = transID_coln,
                                                    transGene_coln = transGene_coln,
                                                    cellID_coln = cellID_coln, 
                                                    spatLocs_colns = spatLocs_colns, 
                                                    extracellular_cellID = extracellular_cellID)
  # when not a file list but direct transcript_df
  if(is.null(transDF_fileInfo)){
    # create `transDF_fileInfo` for the provided `transcript_df`
    transDF_fileInfo <- data.frame(file_path = NA, 
                                   stage_X = 0, 
                                   stage_Y = 0)
    filepath_coln = 'file_path'
    fovOffset_colns = c('stage_X', 'stage_Y')
    prefix_colns = NULL
  }
  
  ## initialize list to collect each FOV outputs ----
  all_segRes <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  all_segRes[['refProfiles']] <- prep_res[['refProfiles']] 
  all_segRes[['baselineData']] <- prep_res[['baselineData']]
  
  all_segRes[['ctrl_genes']] <- prep_res[['ctrl_genes']]
  
  ## transcript score matrix for each gene based on reference profile 
  ## where tLLR score for control genes is same as `svmClass_score_cutoff`
  tLLR_geneMatrix <- prep_res[['score_GeneMatrix']]
  
  
  
  # function for processing each FOV 
  # return `modStats_ToFlagCells` when complete, save results to disk
  # reset seed before transcript error detection for every FOV
  myFun_flag_eachFOV <- function(idx){
    
    ## (1) load and prep each FOV data 
    path_to_transDF <- transDF_fileInfo[[filepath_coln]][idx]
    
    # skip loading for the already loaded first FOV
    if(idx >1){
      if(!is.na(path_to_transDF)){
        transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
      }
      # fovOffset_colns must have XY axes of stage matched to XY axes of images
      # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
      transcript_df <- prepare_perFOV_transDF(each_transDF = transcript_df, 
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
    }
    
    
    
    # processing current FOV
    message(sprintf("\n##############\nProcessing file `%d`: %s\n\n\n",
                    idx, path_to_transDF))
    timestamp()
    if(!is.null(transcript_df[['extraC']])){
      message(sprintf("Exclude %d extracellular transcripts from downstream, %.4f of total molecules.\n\n", 
                      nrow(transcript_df[['extraC']]), 
                      nrow(transcript_df[['extraC']])/(nrow(transcript_df[['extraC']]) + nrow(transcript_df[['intraC']]))))
    }
    
    
    # evaluate segmentation errors
    outs <- runSegErrorEvaluation(score_GeneMatrix= tLLR_geneMatrix, 
                                  transcript_df = transcript_df[['intraC']], 
                                  cellID_coln = 'UMI_cellID', 
                                  transID_coln = 'UMI_transID',
                                  transGene_coln = 'target',
                                  spatLocs_colns = c('x','y','z')[1:d2_or_d3],
                                  flagModel_TransNum_cutoff = flagModel_TransNum_cutoff)
    modStats_ToFlagCells <- outs [['modStats_ToFlagCells']]
    transcript_df[['intraC']] <- outs[['transcript_df']]
    
    rm(outs)
    
    
    if(is.null(modStats_ToFlagCells)){
      # if no cells with enough transcript per cell for model evaluation
      flagged_cells <- NULL
      
    } else {
      ##  flag cells based on linear regression of tLLR, lrtest_-log10P
      modStats_ToFlagCells[['flagged']] <- (modStats_ToFlagCells[['lrtest_-log10P']] > flagCell_lrtest_cutoff )
      
      flagged_cells <- modStats_ToFlagCells[['UMI_cellID']][modStats_ToFlagCells[['flagged']]]
      message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > %.1f.", 
                      length(flagged_cells), length(flagged_cells)/nrow(modStats_ToFlagCells), flagCell_lrtest_cutoff))
      
      modStats_ToFlagCells[['file_idx']] <- idx
      write.csv(modStats_ToFlagCells, file = fs::path(path_to_output, paste0(idx, '_modStats_ToFlagCells.csv')), row.names = FALSE)
      
    }
    
    
    ## (5) use SVM~hyperplane to identify the connected transcripts group based on tLLR score ----
    # SVM can separate continuous low score transcript from the rest.
    # but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
    
    if(length(flagged_cells)>0){
      classDF_ToFlagTrans <- transcript_df[['intraC']][which(transcript_df[['intraC']][['UMI_cellID']] %in% flagged_cells),]
      
      if(!is.null(seed_transError)){
        set.seed(seed_transError)
      }
      # `flagTranscripts_SVM` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
      tmp_df <- flagTranscripts_SVM(chosen_cells = flagged_cells,
                                    score_GeneMatrix = tLLR_geneMatrix,
                                    transcript_df = classDF_ToFlagTrans, 
                                    cellID_coln = 'UMI_cellID', 
                                    transID_coln = 'UMI_transID', 
                                    score_coln = 'score_tLLR_maxCellType',
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
      flaggedSVM_transID <- classDF_ToFlagTrans[classDF_ToFlagTrans[['SVM_class']] ==0, 'UMI_transID']
      # assign SVM_class =0 for transcripts with low goodness-of-fit
      transcript_df[['intraC']][['SVM_class']] <- 1- as.numeric(transcript_df[['intraC']][['UMI_transID']] %in% flaggedSVM_transID)
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
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE}
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
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE}
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

  if(return_intermediates){
    final_res[['groupDF_ToFlagTrans']] <- groupDF_ToFlagTrans
  }
  
  #### (3) get ready for resegmentation ----
  ## (3.1) get cutoff of transcript score and transcript number from baseline
  # the baseline is provided externally, no calculation needed
  
  ## (3.2) prepare transcript_df for re-segmentation ----
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
  
  ## (3.3) evaluate the neighborhood of each group for re-segmentation ----
  groups_to_reseg <- unique(groupDF_ToFlagTrans[which(groupDF_ToFlagTrans[['connect_group']]!=0),][['tmp_cellID']])
  
  
  #### (4) re-segmentation in neighborhood ----
  # did not consider extra cellular transcripts for neighbor identification. 

  outs <- runSegRefinement(score_GeneMatrix = score_GeneMatrix,  
                           chosen_cells = groups_to_reseg, 
                           reseg_transcript_df = reseg_transcript_df, 
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


# wrapper for resegmentation pipeline using internal reference profiles and cutoffs
#' @title fastReseg_full_pipeline
#' @description wrapper for full resegmentation pipeline using internal reference profiles and cutoffs. This function first estimates proper reference profiles and cutoffs from the provided data and then use `fastReseg_perFOV_full_process` function to process each transcript data.frame. 
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
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, unit in micron (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param score_baseline a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence; default = NULL to calculate from `counts` and `refProfiles` 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is; default = NULL to calculate from `counts` and `refProfiles` 
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type; default = NULL to calculate from `counts` and `refProfiles` 
#' @param imputeFlag_missingCTs flag to impute `score_baseline`, `lowerCutoff_transNum`,`higherCutoff_transNum` for cell types present in `refProfiles` but missing in the provided transcript data files or the provided baseline and cutoffs; when TRUE, the median values of existing cell types would be used as the values for missing cell types.
#' @param groupTranscripts_method use either "dbscan" or "delaunay method" to group transcripts in space (default = "dbscan")
#' @param spatialMergeCheck_method use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space (default = "leidenCut")
#' @param cutoff_spatialMerge spatial constraint on a valid merging event between two source transcript groups, default = 0.5 for 50% cutoff, set to 0 to skip spatial constraint evaluation for merging.   
#' For `spatialMergeCheck_method = "leidenCut"`, this is the minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event.
#' For `spatialMergeCheck_method = "geometryDiff"`, this is the maximum percentage of white space change upon merging of query cell and neighbor cell for a valid merging event. 
#' @param leiden_config (leidenCut) a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.  
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config (default = NULL)
#' @param path_to_output the file path to output folder where the resegmentation data is saved; directory would be created by function if not exists; transcript data.frame `updated_transDF` is saved as individual csv files for each FOV, while cell data of all FOVs, `updated_perCellDT` and `updated_perCellExprs`, are combined to save as .RData object.
#' @param save_intermediates flag to save intermediate outputs into output folder, including data.frame for spatial modeling statistics of each cell,  
#' @param return_perCellData flag to return and save to output folder for gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @param combine_extra flag to combine original extracellular transcripts and trimmed transcripts back to the updated transcript data.frame, slow process if many transcript in each FOV file. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not in `refProfiles` and expect no cell type dependency, e.g. negative control probes; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @param seed_process seed for per FOV processing, used in transcript error detection and correction steps, default = NULL to skip the seed  
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
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations` function, save when `save_intermediates` = TRUE}
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
#' # case #'1: provide `transcript_df` directly,
#' # do auto-calculation of distance cutoff from data while using the provided cutoffs for score and transcript numbers.
#' res1 <- fastReseg_full_pipeline(counts = ori_RawExprs,
#'                                 clust = NULL,
#'                                 refProfiles = example_refProfiles,
#'                                 pixel_size = 1,
#'                                 zstep_size = 1,
#'                                 transcript_df = mini_transcriptDF,
#'                                 transID_coln = "UMI_transID",
#'                                 transGene_coln = "target",
#'                                 cellID_coln = "UMI_cellID",
#'                                 spatLocs_colns = c("x","y","z"),
#'                                 extracellular_cellID = extracellular_cellID,
#'                                 molecular_distance_cutoff = NULL,
#'                                 cellular_distance_cutoff = NULL,
#'                                 score_baseline = example_baselineCT[['score_baseline']],
#'                                 lowerCutoff_transNum = example_baselineCT[['lowerCutoff_transNum']],
#'                                 higherCutoff_transNum= example_baselineCT[['higherCutoff_transNum']],
#'                                 imputeFlag_missingCTs = TRUE,
#'                                 path_to_output = "res1_directDF")
#' 
#' # case #'2: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV,
#' # do auto-calculation of score and transcript number cutoffs from gene expression matrix, `counts`, and cluster assignment of each cell, `clust`,
#' # do auto-calculation of distance cutoff from the 1st per FOV transcript data.
#' data("example_CellGeneExpr")
#' data("example_clust")
#' 
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
#' res2 <- fastReseg_full_pipeline(counts = example_CellGeneExpr,
#'                                 clust = example_clust,
#'                                 refProfiles = NULL,
#'                                 transDF_fileInfo =fileInfo_DF,
#'                                 filepath_coln = 'file_path',
#'                                 prefix_colns = c('slide','fov'),
#'                                 fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
#'                                 pixel_size = 0.18, # 0.18 micron per pixel in transcript data
#'                                 zstep_size = 0.8, # 0.8 micron per z step in transcript data
#'                                 transcript_df = NULL,
#'                                 transID_coln = NULL, # row index as transcript_id
#'                                 transGene_coln = "target",
#'                                 cellID_coln = "CellId",
#'                                 spatLocs_colns = c("x","y","z"),
#'                                 extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data
#'                                 molecular_distance_cutoff = NULL,
#'                                 cellular_distance_cutoff = NULL,
#'                                 score_baseline = NULL,
#'                                 lowerCutoff_transNum = NULL,
#'                                 higherCutoff_transNum= NULL,
#'                                 imputeFlag_missingCTs = TRUE,
#'                                 path_to_output = "res2_multiFiles")
#' 
#' # case #'3: provide file paths to per FOV transcript data files and specify the spatial offset for each FOV,
#' # do auto-calculation of score and transcript number cutoffs from gene expression matrix, `counts`, and cluster-specific reference profiles, `refProfiles`,
#' # use the provided distance cutoff for `molecular_distance_cutoff` but calculate the `cellular_distance_cutoff`
#' res3 <- fastReseg_full_pipeline(counts = example_CellGeneExpr,
#'                                 clust = NULL,
#'                                 refProfiles = example_refProfiles,
#'                                 transDF_fileInfo =fileInfo_DF,
#'                                 filepath_coln = 'file_path',
#'                                 prefix_colns = c('slide','fov'),
#'                                 fovOffset_colns = c('stage_Y','stage_X'), 
#'                                 pixel_size = 0.18, 
#'                                 zstep_size = 0.8, 
#'                                 transcript_df = NULL,
#'                                 transID_coln = NULL, 
#'                                 transGene_coln = "target",
#'                                 cellID_coln = "CellId",
#'                                 spatLocs_colns = c("x","y","z"),
#'                                 extracellular_cellID = c(0), 
#'                                 molecular_distance_cutoff = 2.7,
#'                                 cellular_distance_cutoff = NULL,
#'                                 score_baseline = NULL,
#'                                 lowerCutoff_transNum = NULL,
#'                                 higherCutoff_transNum= NULL,
#'                                 imputeFlag_missingCTs = TRUE,
#'                                 path_to_output = "res3_multiFiles")
#' @importFrom fs path
#' @importFrom Matrix rowSums
#' @export 
#' 
fastReseg_full_pipeline <- function(counts, 
                                    clust = NULL, 
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
                                    molecular_distance_cutoff = 2.7,
                                    cellular_distance_cutoff = NULL,
                                    score_baseline = NULL, 
                                    lowerCutoff_transNum = NULL, 
                                    higherCutoff_transNum= NULL, 
                                    imputeFlag_missingCTs = TRUE,
                                    groupTranscripts_method = c("dbscan", "delaunay"),
                                    spatialMergeCheck_method = c("leidenCut", "geometryDiff"), 
                                    cutoff_spatialMerge = 0.5, 
                                    leiden_config = list(objective_function = "CPM",
                                                         resolution_parameter = 1,
                                                         beta = 0.01,
                                                         n_iterations = 200),
                                    config_spatNW_transcript = NULL,
                                    path_to_output = "reSeg_res", 
                                    save_intermediates = TRUE,
                                    return_perCellData = TRUE, 
                                    combine_extra = FALSE, 
                                    ctrl_genes = NULL,
                                    seed_process = NULL){
  
  # create output directory 
  if(!file.exists(path_to_output)) dir.create(path_to_output)
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  
  
  ## check inputs and then get baseline and cutoffs for counts
  ## if either distance cutoff is not provided, the function also checks the format of transcript data.frame provided and load 1st fov 
  prep_res <- runPreprocess(counts = counts, 
                            clust = clust, 
                            refProfiles = refProfiles,
                            score_baseline = score_baseline, 
                            lowerCutoff_transNum = lowerCutoff_transNum, 
                            higherCutoff_transNum= higherCutoff_transNum, 
                            imputeFlag_missingCTs = imputeFlag_missingCTs,
                            ctrl_genes = ctrl_genes,
                            svmClass_score_cutoff = svmClass_score_cutoff,
                            molecular_distance_cutoff = molecular_distance_cutoff,
                            cellular_distance_cutoff = cellular_distance_cutoff,
                            transcript_df = transcript_df, 
                            transDF_fileInfo = transDF_fileInfo, 
                            filepath_coln = filepath_coln, 
                            prefix_colns = prefix_colns, 
                            fovOffset_colns = fovOffset_colns, 
                            pixel_size = pixel_size, 
                            zstep_size = zstep_size,
                            transID_coln = transID_coln,
                            transGene_coln = transGene_coln,
                            cellID_coln = cellID_coln, 
                            spatLocs_colns = spatLocs_colns, 
                            extracellular_cellID = extracellular_cellID)
  
  # if both distance cutoff are provided, do checking and 1st fov loading in separate function
  if(is.null(prep_res[['processed_1st_transDF']])){
    ## check the format of transcript data.frame provided and load 1st fov 
    # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
    transcript_df <- checkTransFileInputsAndLoadFirst(transcript_df = transcript_df, 
                                                      transDF_fileInfo = transDF_fileInfo, 
                                                      filepath_coln = filepath_coln, 
                                                      prefix_colns = prefix_colns, 
                                                      fovOffset_colns = fovOffset_colns, 
                                                      pixel_size = pixel_size, 
                                                      zstep_size = zstep_size,
                                                      transID_coln = transID_coln,
                                                      transGene_coln = transGene_coln,
                                                      cellID_coln = cellID_coln, 
                                                      spatLocs_colns = spatLocs_colns, 
                                                      extracellular_cellID = extracellular_cellID)
  } else{
    transcript_df <- prep_res[['processed_1st_transDF']]
  }
  
  # when not a file list but direct transcript_df
  if(is.null(transDF_fileInfo)){
    # create `transDF_fileInfo` for the provided `transcript_df`
    transDF_fileInfo <- data.frame(file_path = NA, 
                                   stage_X = 0, 
                                   stage_Y = 0)
    filepath_coln = 'file_path'
    fovOffset_colns = c('stage_X', 'stage_Y')
    prefix_colns = NULL
  }
  
  
  
  ## initialize list to collect each FOV outputs ----
  all_segRes <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  all_segRes[['refProfiles']] <- prep_res[['refProfiles']] 
  all_segRes[['baselineData']] <- prep_res[['baselineData']]
  
  all_segRes[['cutoffs_list']] <- prep_res[['cutoffs_list']]
  
  
  score_baseline <- prep_res[['cutoffs_list']][['score_baseline']]
  lowerCutoff_transNum <- prep_res[['cutoffs_list']][['lowerCutoff_transNum']]
  higherCutoff_transNum <- prep_res[['cutoffs_list']][['higherCutoff_transNum']]
  molecular_distance_cutoff <- prep_res[['cutoffs_list']][['molecular_distance_cutoff']]
  cellular_distance_cutoff <- prep_res[['cutoffs_list']][['cellular_distance_cutoff']]
  
  
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
  
  all_segRes[['ctrl_genes']] <- prep_res[['ctrl_genes']]
  ## transcript score matrix for each gene based on reference profile 
  ## where tLLR score for control genes is same as `svmClass_score_cutoff`
  tLLR_geneMatrix <- prep_res[['score_GeneMatrix']]
  
  rm(prep_res)
  
  ## apply the fixed cutoffs settings to individual FOVs for resegmentation ----
  
  # save `updated_transDF` into csv file and intermeidates file into .RData for each FOV
  # but combine perCell data from all FOVs to return and save as single .RData object
  
  # function for processing each FOV for resegmentation
  # return `each_segRes` when complete, save results to disk
  myFun_reseg_eachFOV <- function(idx){
    # load and prep each FOV data 
    path_to_transDF <- transDF_fileInfo[[filepath_coln]][idx]
    
    if (idx !=1){
      if(!is.na(path_to_transDF)){
        transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
      }
      # fovOffset_colns must have XY axes of stage matched to XY axes of images
      # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
      transcript_df <- prepare_perFOV_transDF(each_transDF = transcript_df, 
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
    
    
    # # `fastReseg_perFOV_full_process` function is a wrapper for resegmentation pipeline using external reference profiles and cutoffs. 
    # # The function returns a list containing the following elements:
    # modStats_ToFlagCells, a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, return when `return_intermediates` = TRUE
    # groupDF_ToFlagTrans, data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flagTranscripts_SVM` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE
    # neighborhoodDF_ToReseg, a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE
    # reseg_actions, a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates` = TRUE
    # updated_transDF, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter
    # updated_perCellDT, a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE
    # updated_perCellExprs, a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE
    
    # set `includeAllRefGenes` to TRUE, to ensure return of all genes in `score_GeneMatrix` to `each_segRes[['updated_perCellExprs']]`
    # missing genes would have imputed value of zero
    
    # To account for genes missing in `refProfiles` but present in input transcript data.frame, 
    # genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. 
    # It's recommended to have total counts of those genes below 1% of total counts of all genes in each cell to avoid significant interference from those `ctrl_genes`,  
    each_segRes <- fastReseg_perFOV_full_process(score_GeneMatrix = tLLR_geneMatrix,
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
                                                 spatialMergeCheck_method = spatialMergeCheck_method, 
                                                 cutoff_spatialMerge = cutoff_spatialMerge, 
                                                 leiden_config = leiden_config, 
                                                 config_spatNW_transcript = config_spatNW_transcript,
                                                 return_intermediates = save_intermediates,
                                                 return_perCellData = return_perCellData, 
                                                 includeAllRefGenes = TRUE,
                                                 seed_process = seed_process)
    
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









