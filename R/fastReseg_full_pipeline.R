# wrapper for resegmentation pipeline using internal reference profiles and cutoffs
#' @title fastReseg_full_pipeline
#' @description wrapper for full resegmentation pipeline using internal reference profiles and cutoffs. This function first estimates proper reference profiles and cutoffs from the provided data and then use `fastReseg_perFOV_full_process` function to process each transcript data.frame. 
#' @param counts Counts matrix for entire dataset, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, when NULL to automatically assign the cell cluster for each cell based on maximum transcript score of given the provided `refProfiles`
#' @param refProfiles A matrix of cluster profiles, genes X clusters, default = NULL to use external cluster assignments
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire dataset; when NULL, use the provided `transcript_df` directly
#' @param filepath_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
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
#' @param invert_y flag to invert y axis of local coordinates during stitching (default = TRUE)
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @param flagCell_lrtest_cutoff the cutoff of lrtest_nlog10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
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
#' @param path_to_output the file path to output folder where the resegmentation data is saved; directory would be created by function if not exists; transcript data.frame `updated_transDF` is saved as individual csv files for each FOV, while cell data of all FOVs, `updated_perCellDT` and `updated_perCellExprs`, are combined to save as `.rds` object.
#' @param transDF_export_option option on how to export updated transcript_df, 0 for no export, 1 for write to `path_to_output` in disk as csv for each FOV, 2 for return to function as list (default = 1)
#' @param save_intermediates flag to save intermediate outputs into output folder, including data.frame for spatial modeling statistics of each cell,  
#' @param return_perCellData flag to return for gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type, also export to output folder in disk when `transDF_export_option = 1`.
#' @param combine_extra flag to combine original extracellular transcripts and trimmed transcripts back to the updated transcript data.frame, slow process if many transcript in each FOV file. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not in `refProfiles` and expect no cell type dependency, e.g. negative control probes; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @param seed_process seed for per FOV processing, used in transcript error detection and correction steps, default = NULL to skip the seed  
#' @param percentCores percent of cores to use for parallel processing (0-1] (default = 0.75)
#' @return a list 
#' \describe{
#'    \item{refProfiles}{a genes X clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{cutoffs_list}{a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations` function, return when `save_intermediates` = TRUE}
#'    \item{updated_transDF_list}{a list of per-FOV transcript data.frame with updated cell segmenation in `updated_cellID` and `updated_celltype` columns, return when `transDF_export_option = 2`}
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
#' When `save_intermediates` = TRUE, all intermediate files and resegmenation outputs of each FOV would be saved as single `.rds` object which is a list containing the following elements:  
#' \describe{
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function, save when `save_intermediates` = TRUE}
#'    \item{groupDF_ToFlagTrans}{data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, save when `save_intermediates` = TRUE}
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood environment of low-score transcript groups, output of `get_neighborhood_content` function, save when `save_intermediates` = TRUE}
#'    \item{reseg_actions}{a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations` function, save when `save_intermediates` = TRUE}
#'    \item{updated_transDF}{the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter, write to disk when `transDF_export_option =1`}
#'    \item{updated_perCellDT}{a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE}
#' }
#' The pipeline would also combine per cell data for all FOVs and return the combined data when `return_perCellData` = TRUE; `updated_perCellDT` and `updated_perCellExprs` would also be saved in a list as single `.rds` object in `path_to_output` directory when  `transDF_export_option = 1`.
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
#' # cell_ID for extracellualr transcripts
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] 
#' 
#' # case #'1: provide `transcript_df` directly,
#' # do auto-calculation of distance cutoff from data while using the provided 
#' # cutoffs for score and transcript numbers.
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
#' # case #'2: provide file paths to per FOV transcript data files and specify 
#' # the spatial offset for each FOV,
#' # do auto-calculation of score and transcript number cutoffs from gene 
#' # expression matrix, `counts`, and cluster assignment of each cell, `clust`,
#' # do auto-calculation of distance cutoff from the 1st per FOV transcript data.
#' data("example_CellGeneExpr")
#' data("example_clust")
#' 
#' # the example individual transcript files are stored under `data` directory of this package
#' # update your path accordingly
#' # Notice that some assays like SMI has XY axes swapped between stage and each FOV;
#' # coordinates for each FOV should have units in micron
#' dataDir <- system.file("extdata", package = "FastReseg")
#' fileInfo_DF <- data.frame(
#'   file_path = fs::path(dataDir,  
#'                        c("Run4104_FOV001__complete_code_cell_target_call_coord.csv",
#'                          "Run4104_FOV002__complete_code_cell_target_call_coord.csv")),
#'   slide = c(1, 1),
#'   fov = c(1,2),
#'   stage_X = 1000*c(5.13, -2.701),
#'   stage_Y = 1000*c(-0.452, 0.081))
#' 
#' res2 <- fastReseg_full_pipeline(counts = example_CellGeneExpr,
#'                                 clust = example_clust,
#'                                 refProfiles = NULL,
#'                                 transDF_fileInfo =fileInfo_DF,
#'                                 filepath_coln = 'file_path',
#'                                 prefix_colns = c('slide','fov'),
#'                                 
#'                                 # match XY axes between stage and each FOV
#'                                 fovOffset_colns = c('stage_Y','stage_X'), 
#'                                 # 0.18 micron per pixel in transcript data
#'                                 pixel_size = 0.18, 
#'                                 # 0.8 micron per z step in transcript data
#'                                 zstep_size = 0.8, 
#'                                 
#'                                 transcript_df = NULL,
#'                                 
#'                                 # row index as transcript_id
#'                                 transID_coln = NULL, 
#'                                 
#'                                 transGene_coln = "target",
#'                                 cellID_coln = "CellId",
#'                                 spatLocs_colns = c("x","y","z"),
#'                                 
#'                                 # CellId = 0 means extracelluar transcripts in raw data
#'                                 extracellular_cellID = c(0), 
#'                                 
#'                                 molecular_distance_cutoff = NULL,
#'                                 cellular_distance_cutoff = NULL,
#'                                 score_baseline = NULL,
#'                                 lowerCutoff_transNum = NULL,
#'                                 higherCutoff_transNum= NULL,
#'                                 imputeFlag_missingCTs = TRUE,
#'                                 path_to_output = "res2_multiFiles")
#' 
#' # case #'3: provide file paths to per FOV transcript data files and specify 
#' # the spatial offset for each FOV,
#' # do auto-calculation of score and transcript number cutoffs from gene 
#' # expression matrix, `counts`, and cluster-specific reference profiles, `refProfiles`,
#' # use the provided distance cutoff for `molecular_distance_cutoff` but 
#' # calculate the `cellular_distance_cutoff`
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
#' @importFrom utils timestamp write.csv
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
                                    invert_y = TRUE,
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
                                    transDF_export_option = c(1, 2, 0),
                                    save_intermediates = TRUE,
                                    return_perCellData = TRUE, 
                                    combine_extra = FALSE, 
                                    ctrl_genes = NULL,
                                    seed_process = NULL,
                                    percentCores = 0.75){
  transDF_export_option <- match.arg(as.character(transDF_export_option)[1], choices = c(1, 2, 0))
  if(transDF_export_option ==0){
    message("No transcript data.frame or per cell data would be exported with `transDF_export_option = 0`.")
  } else if (transDF_export_option ==1){
    message(sprintf("Per-FOV transcript data.frame with updated cell segmentation would be exported to disk at `path_to_output = '%s'`.", 
                    path_to_output))
  } else if (transDF_export_option ==2){
    message("Per-FOV transcript data.frame with updated cell segmentation would be returned to function in list `updated_transDF_list`.")
  } else {
    stop("Must specify how to export transcript data.frame in `transDF_export_option`.")
  }
  
  if(transDF_export_option ==1 || save_intermediates){
    if(!file.exists(path_to_output)) dir.create(path_to_output)
    if(save_intermediates || return_perCellData){
      message(sprintf("The %s would be exported to output directory at `path_to_output`.", 
                      paste0(c("intermediate results", "per cell data")[c(save_intermediates, return_perCellData)], 
                             collapse = " and ")))
    }
  }

  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  
  if(percentCores > 1 & percentCores <= 0){
    stop("percentCores is not a valid number, must be between 0-1")
  }
  
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
                            invert_y = invert_y,
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
                                                      invert_y = invert_y,
                                                      extracellular_cellID = extracellular_cellID)
  } else{
    transcript_df <- prep_res[['processed_1st_transDF']]
    prep_res[['processed_1st_transDF']] <- NULL
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
  }else{
    transDF_fileInfo <- as.data.frame(transDF_fileInfo)
  }
  

  ## initialize list to collect each FOV outputs
  all_segRes <- initializeAllSegRes(prep_res = prep_res, 
                                    save_intermediates = save_intermediates,
                                    return_perCellData = return_perCellData)
 
  ## apply the fixed cutoffs settings to individual FOVs for resegmentation ----
  
  # save `updated_transDF` into csv file and intermediates file into `.rds` for each FOV
  # but combine perCell data from all FOVs to return and save in a list as single `.rds` object
  
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
                                              invert_y = invert_y,
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
    # groupDF_ToFlagTrans, data.frame for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts` and `groupTranscripts_Delaunay` or `groupTranscripts_dbscan` functions, return when `return_intermediates` = TRUE
    # neighborhoodDF_ToReseg, a data.frame for neighborhood environment of low-score transcript groups, output of `get_neighborhood_content` function, return when `return_intermediates` = TRUE
    # reseg_actions, a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations` function, return when `return_intermediates` = TRUE
    # updated_transDF, the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter
    # updated_perCellDT, a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData` = TRUE
    # updated_perCellExprs, a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when `return_perCellData` = TRUE
    
    # set `includeAllRefGenes` to TRUE, to ensure return of all genes in `score_GeneMatrix` to `each_segRes[['updated_perCellExprs']]`
    # missing genes would have imputed value of zero
    
    # To account for genes missing in `refProfiles` but present in input transcript data.frame, 
    # genes in `ctrl_genes` would be assigned with goodness-of-fit score equal to `svmClass_score_cutoff` for all cell types to minimize the impact of those genes on the identification of low-score transcript groups via SVM. 
    # It's recommended to have total counts of those genes below 1% of total counts of all genes in each cell to avoid significant interference from those `ctrl_genes`,  
    each_segRes <- fastReseg_perFOV_full_process(score_GeneMatrix = prep_res[['score_GeneMatrix']],
                                                 transcript_df = transcript_df[['intraC']],
                                                 extracellular_cellID = NULL,
                                                 molecular_distance_cutoff = prep_res[['cutoffs_list']][['molecular_distance_cutoff']],
                                                 cellular_distance_cutoff = prep_res[['cutoffs_list']][['cellular_distance_cutoff']],
                                                 score_baseline = prep_res[['cutoffs_list']][['score_baseline']],
                                                 lowerCutoff_transNum = prep_res[['cutoffs_list']][['lowerCutoff_transNum']],
                                                 higherCutoff_transNum= prep_res[['cutoffs_list']][['higherCutoff_transNum']],  
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
    
    # add in transcript compartment label and 
    # optional to merge original transcript_df to updated results, NA for trimmed or extracellular transcripts
    each_segRes[['updated_transDF']] <- compartment_and_add_extra(updated_transDF = each_segRes[['updated_transDF']], 
                                                                  transcript_df = transcript_df, 
                                                                  updatedCellID_coln = "updated_cellID",
                                                                  compartment_coln = "transComp",
                                                                  combine_extra = combine_extra) 
    
    rm(transcript_df)
    
    # save `updated_transDF` into csv file for each FOV 
    if(transDF_export_option ==1){
      write.csv(each_segRes[['updated_transDF']], 
                file = fs::path(path_to_output, paste0(idx, "_updated_transDF.csv")), 
                row.names = FALSE)
    }
    
    
    # return only idx, perCell data and reseg_actions as a list
    res_to_return <- list(idx = idx)
    
    if(return_perCellData){
      res_to_return[['updated_perCellDT']] <- each_segRes[['updated_perCellDT']]
      res_to_return[['updated_perCellExprs']] <- each_segRes[['updated_perCellExprs']]
    }
    
    if(save_intermediates){
      res_to_return[['reseg_actions']] <- each_segRes[['reseg_actions']]
      
      # save intermediate files and all other outputs for current FOVs as single `.rds`
      saveRDS(each_segRes, 
              file = fs::path(path_to_output, paste0(idx, "_each_segRes.rds")))
    }
    
    if(transDF_export_option ==2){
      res_to_return[['updated_transDF']] <- each_segRes[['updated_transDF']]
    }
    
    rm(each_segRes)
    gc()
    
    return(res_to_return)
  }
  
  # processing each FOV in parallel
  process_outputs <- parallel::mclapply(X = seq_len(nrow(transDF_fileInfo)), 
                                        mc.allow.recursive = TRUE,
                                        mc.cores = numCores(percentCores = percentCores),
                                        FUN = myFun_reseg_eachFOV)
  
  ## combine data for each FOV ----
  if(return_perCellData){
    # extract perCell data from nested list of process_outputs
    updated_perCellDT <- lapply(process_outputs, '[[', 'updated_perCellDT')
    updated_perCellDT <- do.call(rbind, updated_perCellDT)
    
    # per cell gene expression in gene x cell sparse matrix format, do cbind
    updated_perCellExprs <- lapply(process_outputs, '[[', 'updated_perCellExprs')
    updated_perCellExprs <- combine_matrices_fast(matrix_list = updated_perCellExprs, 
                                                  bind = "cbind", fill = 0)
    
    # save combined perCell data into `.rds` object as a list
    if(transDF_export_option ==1){
      saveRDS(list(updated_perCellDT, 
                   updated_perCellExprs), 
              file = fs::path(path_to_output, "combined_updated_perCellDT_perCellExprs.rds"))
    }
    
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
  
  if(transDF_export_option ==2){
    all_segRes[['updated_transDF_list']] <- lapply(process_outputs, '[[', 'updated_transDF')
  }
  
  rm(process_outputs)
  gc()
  
  return(all_segRes)
  
}


#' @title initializeAllSegRes
#' @description initialize result holder for \code{fastReseg_full_pipeline} based on preprocessed data
#' @param prep_res nested list of elements for 'refProfiles', 'baselineData', 'cutoffs_list', outputs of \code{runPreprocess}
#' @param save_intermediates flag to save intermediate outputs into output folder, including data.frame for spatial modeling statistics of each cell,  
#' @param return_perCellData flag to return and save to output folder for gene x cell count matrix and per cell DF with updated mean spatial coordinates and new cell type
#' @return a list 
#' \describe{
#'    \item{refProfiles}{a genes X clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{cutoffs_list}{a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types}
#'    \item{updated_perCellDT}{an empty list, return when `return_perCellData` = TRUE}
#'    \item{updated_perCellExprs}{an empty list, return when `return_perCellData` = TRUE}
#'    \item{reseg_actions}{a list of 4 empty elements, same output format as \code{decide_ReSegment_Operations}, return when `save_intermediates` = TRUE}
#' }
initializeAllSegRes <- function(prep_res, 
                                save_intermediates = TRUE,
                                return_perCellData = TRUE){
  ## initialize list to collect each FOV outputs 
  all_segRes <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  all_segRes[['refProfiles']] <- prep_res[['refProfiles']] 
  all_segRes[['baselineData']] <- prep_res[['baselineData']]
  all_segRes[['cutoffs_list']] <- prep_res[['cutoffs_list']]
  
  all_segRes[['ctrl_genes']] <- prep_res[['ctrl_genes']]
  
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
  
  
  return(all_segRes)
}

#' @title compartment_and_add_extra
#' @description add compartment label to the updated transcript data.frame returned by \code{fastReseg_perFOV_full_process} and optional to combine it with transcripts got trimmed and originally extracellular transcripts
#' @param updated_transDF transcript data.frame with transcript id and updated cell segmentation, one of the outputs  of \code{fastReseg_perFOV_full_process}
#' @param transcript_df list of 2 elements, `intraC` and `extraC` for intracellular and extracellular transcript data.frame, respectively; columns contain transcript id and original cell segmentation, output of \code{prepare_perFOV_transDF}  
#' @param updatedCellID_coln column name of updated cell ID in `updated_transDF`, default = "updated_cellID"
#' @param compartment_coln column name for transcript compartment label, default = "transComp"
#' @param combine_extra flag to combine original extracellular transcripts and trimmed transcripts back to the updated transcript data.frame, slow process if many transcript in each FOV file. (default = FALSE)
#' @return a transcript data.frame containing all original transcripts with updated cell segmentation and transcript compartments for `intraC`, `extraC` and `trimmed`
compartment_and_add_extra <- function(updated_transDF, 
                                     transcript_df, 
                                     updatedCellID_coln = "updated_cellID",
                                     compartment_coln = "transComp",
                                     combine_extra = FALSE){
  if(!is.list(transcript_df)){
    stop("The provided `transcript_df` is not a list of 2 elements, `intraC` and `extraC`. The raw transcript data.frame needs to be processed by `prepare_perFOV_transDF()` function.")
  }

  # intracellular in original and updated segmentation
  updated_transDF[[compartment_coln]] <- 'intraC' 
  
  if(nrow(updated_transDF) != nrow(transcript_df[['intraC']])){
    message(sprintf("\n\nTrim %d transcripts during resegmentation, %.4f of all intracellular molecules.", 
                    nrow(transcript_df[['intraC']]) - nrow(updated_transDF), 
                    1- nrow(updated_transDF)/ nrow(transcript_df[['intraC']])))
    
    # merge original intracellular transcript_df to updated results, NA for trimmed transcripts
    # slow process when large transcript data.frame for given FOV
    if(combine_extra){
      shared_colns <- intersect(colnames(updated_transDF), colnames(transcript_df[['intraC']]))
      shared_colns <- setdiff(shared_colns, compartment_coln)
      updated_transDF <- merge(updated_transDF, 
                               transcript_df[['intraC']][, shared_colns], all = TRUE)
      
      updated_transDF[[compartment_coln]][is.na(updated_transDF[[ updatedCellID_coln]])] <- 'trimmed'
    }
  }
  
  # add original extracellular transcript_df to updated results, NA for extracellular transcripts
  if(all(combine_extra, 
         !is.null(transcript_df[['extraC']]))){
    message(sprintf("Combine newly trimmed and %d orignal extracellular transcripts to the updated transcript data.frame.\n", nrow(transcript_df[['extraC']]))) 
    
    # add dummy values for extracellular transcripts
    colns_to_add <- setdiff(colnames(updated_transDF), colnames(transcript_df[['extraC']]))
    dummy_data <- matrix(NA, nrow = nrow(transcript_df[['extraC']]), ncol = length(colns_to_add))
    colnames(dummy_data) <- colns_to_add
    dummy_data <- as.data.frame(dummy_data)
    dummy_data[[compartment_coln]]<- 'extraC'
    
    shared_colns <- intersect(colnames(updated_transDF), colnames(transcript_df[['extraC']]))
    transcript_df[['extraC']] <- cbind(transcript_df[['extraC']][, shared_colns], 
                                       dummy_data)
    
    updated_transDF <- rbind(updated_transDF, transcript_df[['extraC']])

    } 

  
  return(updated_transDF)
  
}


#' @title combine_matrices_fast
#' @description Efficiently combines a list of matrices by aligning along the necessary axis and then binding.
#' @param matrix_list A list of matrices to combine. All matrices must have non-NULL dimension names
#'   and be of same type, either sparse (`sparseMatrix` class from the Matrix package) or dense (`matrix`).
#' @param bind Character string, `"cbind"` or `"rbind"`. Controls the bind direction:
#'   - `"cbind"`: align rows (union of rownames) then column-bind.
#'   - `"rbind"`: align columns (union of colnames) then row-bind.
#' @param fill Value to fill missing entries introduced by alignment. Default `0`.
#'   - If `output_type = "sparse"` or all inputs are sparse and `fill = NA`, the function sets `fill = 0`.
#'   - If `fill != 0` and sparse output is requested, the function converts to dense and warns.
#' @param output_type One of `"auto"`, `"dense"`, `"sparse"`. Controls the output class.
#'   - `"auto"`: if all inputs are sparse → sparse; if all are dense → dense; if mixed → sparse **only when** `fill == 0`,
#'     otherwise dense.
#'   - `"dense"`: convert any sparse inputs to dense and return a dense matrix.
#'   - `"sparse"`: convert any dense inputs to sparse and return a sparse matrix. If `fill != 0`, converts to dense with warning.
#' @return A combined matrix:
#'   * If `output_type` resolves to sparse, returns a `Matrix::sparseMatrix` (typically `dgCMatrix`).
#'   * If `output_type` resolves to dense, returns a base R `matrix`.
#' @details
#' - Missing rows or columns are filled with the specified `fill` value.
#' - For sparse outputs, missing entries are naturally zeros. If `fill = NA`, it is replaced with `0` and a warning is issued.
#' - If `fill != 0` and sparse output requested, the function converts to dense and warns.
#' @export
combine_matrices_fast <- function(matrix_list, 
                                  bind = c("cbind", "rbind"),
                                  fill = 0,
                                  output_type = c("auto", "dense", "sparse")) {
  stopifnot(is.list(matrix_list))
  if(length(matrix_list) <2){
    return(matrix_list[[1]])
  }
  
  bind <- match.arg(bind, choices = c("cbind", "rbind"))
  output_type <- match.arg(output_type, choices =c("auto", "dense", "sparse"))
  
  # Class detection
  is_sparse_vec <- vapply(matrix_list, function(M) inherits(M, "sparseMatrix"), logical(1))
  is_dense_vec  <- vapply(matrix_list, is.matrix, logical(1))
  if (!all(is_sparse_vec | is_dense_vec)) {
    stop("All matrices must be either Matrix::sparseMatrix or base R dense matrices.")
  }
  
  # Check names presence
  if (!all(vapply(matrix_list, function(M) !is.null(rownames(M)) && !is.null(colnames(M)), logical(1)))) {
    stop("All matrices must have non-NULL rownames and colnames.")
  }
  
  # Decide target output type
  all_sparse <- all(is_sparse_vec)
  all_dense  <- all(is_dense_vec)
  target_sparse <- switch(output_type,
                          auto   = if (all_sparse) TRUE else if (all_dense) FALSE else fill == 0,
                          dense  = FALSE,
                          sparse = TRUE
  )
  
  # Handle fill for sparse output
  if (target_sparse) {
    if (is.na(fill)) {
      warning("Sparse output requested; NA fill replaced with 0 for efficiency.")
      fill <- 0
    }
    if (fill != 0) {
      warning("Sparse output requested but fill is non-zero. Converting to dense output.")
      target_sparse <- FALSE
    }
  }
  
  # Coerce inputs to target type if mixed or forced
  if (target_sparse) {
    matrix_list <- lapply(matrix_list, function(M) {
      if (inherits(M, "sparseMatrix")) M else Matrix::Matrix(M, sparse = TRUE)
    })
  } else {
    matrix_list <- lapply(matrix_list, function(M) {
      if (is.matrix(M)) M else as.matrix(M)
    })
  }
  
  # Axis universes
  all_rows <- if (bind == "cbind") Reduce(union, lapply(matrix_list, rownames)) else NULL
  all_cols <- if (bind == "rbind") Reduce(union, lapply(matrix_list, colnames)) else NULL
  
  # Alignment helpers
  align_dense_rows <- function(M, all_rows, fill) {
    out <- matrix(fill, nrow = length(all_rows), ncol = ncol(M),
                  dimnames = list(all_rows, colnames(M)))
    rn <- rownames(M); pos <- match(rn, all_rows)
    keep <- !is.na(pos)
    if (any(keep)) out[pos[keep], ] <- M[keep, , drop = FALSE]
    out
  }
  align_dense_cols <- function(M, all_cols, fill) {
    out <- matrix(fill, nrow = nrow(M), ncol = length(all_cols),
                  dimnames = list(rownames(M), all_cols))
    cn <- colnames(M); pos <- match(cn, all_cols)
    keep <- !is.na(pos)
    if (any(keep)) out[, pos[keep]] <- M[, keep, drop = FALSE]
    out
  }
  align_sparse_rows <- function(M, all_rows) {
    rn <- rownames(M)
    row_map <- setNames(seq_along(all_rows), all_rows)
    sm <- Matrix::summary(M)
    new_i <- row_map[rn[sm$i]]
    Matrix::sparseMatrix(i = new_i, j = sm$j, x = sm$x,
                         dims = c(length(all_rows), ncol(M)),
                         dimnames = list(all_rows, colnames(M)))
  }
  align_sparse_cols <- function(M, all_cols) {
    cn <- colnames(M)
    col_map <- setNames(seq_along(all_cols), all_cols)
    sm <- Matrix::summary(M)
    new_j <- col_map[cn[sm$j]]
    Matrix::sparseMatrix(i = sm$i, j = new_j, x = sm$x,
                         dims = c(nrow(M), length(all_cols)),
                         dimnames = list(rownames(M), all_cols))
  }
  
  if (bind == "cbind") {
    aligned <- if (target_sparse) {
      lapply(matrix_list, align_sparse_rows, all_rows = all_rows)
    } else {
      lapply(matrix_list, align_dense_rows, all_rows = all_rows, fill = fill)
    }
  } else { # bind == "rbind"
    aligned <- if (target_sparse) {
      lapply(matrix_list, align_sparse_cols, all_cols = all_cols)
    } else {
      lapply(matrix_list, align_dense_cols, all_cols = all_cols, fill = fill)
    }
  }
  
  
  # Bind
  out <- if (target_sparse) {
    res <- if (bind == "cbind") do.call(Matrix::cbind2, aligned) else do.call(Matrix::rbind2, aligned)
    Matrix::drop0(res)
  } else {
    if (bind == "cbind") do.call(cbind, aligned) else do.call(rbind, aligned)
  }
  
  return(out)
}

