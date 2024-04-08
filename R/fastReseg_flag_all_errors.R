# wrapper to flag segmentation error in all files
#' @title fastReseg_flag_all_errors
#' @description Wrapper to process multiple files of one dataset for segmentation error detection in transcript level. The function reformats the individual transcript data.frame to have unique IDs and a global coordinate system and save into disk, then scores each cell for segmentation error and flags transcripts that have low goodness-of-fit to current cells.
#' @param counts Counts matrix for entire data set, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, when NULL to automatically assign the cell cluster for each cell based on maximum transcript score of given the provided `refProfiles`
#' @param refProfiles A matrix of cluster profiles, genes X clusters, default = NULL to use external cluster assignments
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set; when NULL, use the provided `transcript_df` directly
#' @param filepath_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
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
#' @param flagCell_lrtest_cutoff the cutoff of `lrtest_nlog10P` to identify putative wrongly segmented cells with strong spatial dependency in transcript score profile
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
#' @param svm_args a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
#' @param path_to_output the file path to output folder; directory would be created by function if not exists; `flagged_transDF`, the reformatted transcript data.frame with transcripts of low goodness-of-fit flagged by` SVM_class = 0`, and `modStats_ToFlagCells`, the per cell evaluation output of segmentation error, and `classDF_ToFlagTrans`, the class assignment of transcripts within each flagged cells are saved as individual csv files for each FOV, respectively.
#' @param transDF_export_option option on how to export updated transcript_df, 0 for no export, 1 for write to `path_to_output` in disk as csv for each FOV, 2 for return to function as list (default = 1)
#' @param return_trimmed_perCell flag to return a gene x cell count sparse matrix where all putative contaminating transcripts are trimmed (default = FALSE)
#' @param combine_extra flag to combine original extracellular transcripts back to the flagged transcript data.frame. (default = FALSE)
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not present in `counts` or `refProfiles`; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @param percentCores percent of cores to use for parallel processing (0-1] (default = 0.75)
#' @param seed_transError seed for transcript error detection step, default = NULL to skip the seed   
#' @return a list 
#' \describe{
#'    \item{refProfiles}{a genes * clusters matrix of cluster-specific reference profiles used in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrices in cluster * percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.}
#'    \item{combined_modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell for all cells in the data set, output of `score_cell_segmentation_error` function}
#'    \item{combined_flaggedCells}{a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV}
#'    \item{trimmed_perCellExprs}{a gene x cell count sparse matrix where all putative contaminating transcripts are trimmed, return when `return_trimmed_perCell` = TRUE}
#'    \item{flagged_transDF_list}{a list of per-FOV transcript data.frame with flagging information in `SVM_class` column, return when `transDF_export_option = 2`}
#' }
#' @details The function would first estimate mean profile for each cell cluster based on the provided cell x gene count matrix and cluster assignment for entire data set. 
#' And then, the function would use the estimated cluster-specific profile as reference profiles when not provided. 
#' For each transcript data.frame, the function would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, and identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile. 
#' When `transDF_export_option =1`, the function would save the each per FOV output as individual file in `path_to_output` directory; `flagged_transDF`, `modStats_ToFlagCells` and `classDF_ToFlagTrans` would be saved as csv file, respectively. 
#' \describe{
#'    \item{flagged_transDF}{a transcript data.frame for each FOV, with columns for unique IDs of transcripts `UMI_transID` and cells `UMI_cellID`, for global coordinate system `x`, `y`, `z`, and for the goodness-of-fit in original cell segment `SMI_class`; the original per FOV cell ID and pixel/index-based coordinates systems are saved under columns, `CellId`, `pixel_x`, `pixel_y`, `idx_z`}
#'    \item{modStats_ToFlagCells}{a data.frame for spatial modeling statistics of each cell, output of `score_cell_segmentation_error` function}
#'    \item{classDF_ToFlagTrans}{data.frame for the class assignment of transcripts within putative wrongly segmented cells, output of `flag_bad_transcripts` functions}
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
#' dataDir <- system.file("extdata", package = "FastReseg")
#' fileInfo_DF <- data.frame(file_path = fs::path(dataDir,
#'                                                c("Run4104_FOV001__complete_code_cell_target_call_coord.csv",
#'                                                  "Run4104_FOV002__complete_code_cell_target_call_coord.csv")),
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
                                      transDF_export_option = c(1, 2, 0),
                                      return_trimmed_perCell = FALSE, 
                                      combine_extra = FALSE, 
                                      ctrl_genes = NULL,
                                      seed_transError = NULL,
                                      percentCores = 0.75){
  transDF_export_option <- match.arg(as.character(transDF_export_option)[1], choices = c(1, 2, 0))
  if(transDF_export_option ==0){
    message("No transcript data.frame or other per FOV outputs would be exported with `transDF_export_option = 0`.")
  } else if (transDF_export_option ==1){
    message(sprintf("Per-FOV outputs including transcript data.frame with flagging information would be exported to disk at `path_to_output = '%s'`.", 
                    path_to_output))
    # create output directory, for per-FOV outputs 
    if(!file.exists(path_to_output)) dir.create(path_to_output)
  } else if (transDF_export_option ==2){
    message("Per-FOV transcript data.frame with flagging information would be returned to function in list `flagged_transDF_list`.")
  } else {
    stop("Must specify how to export transcript data.frame in `transDF_export_option`.")
  }
  
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  
  if(percentCores > 1 & percentCores <= 0){
    stop("percentCores is not a valid number, must be between 0-1")
  }
  
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
    outs <- runSegErrorEvaluation(score_GeneMatrix= prep_res[['score_GeneMatrix']], 
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
      ##  flag cells based on linear regression of tLLR, lrtest_nlog10P
      modStats_ToFlagCells[['flagged']] <- (modStats_ToFlagCells[['lrtest_nlog10P']] > flagCell_lrtest_cutoff )
      
      flagged_cells <- modStats_ToFlagCells[['UMI_cellID']][modStats_ToFlagCells[['flagged']]]
      message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_nlog10P > %.1f.", 
                      length(flagged_cells), length(flagged_cells)/nrow(modStats_ToFlagCells), flagCell_lrtest_cutoff))
      
      modStats_ToFlagCells[['file_idx']] <- idx
      
      if(transDF_export_option ==1){
        write.csv(modStats_ToFlagCells, file = fs::path(path_to_output, paste0(idx, '_modStats_ToFlagCells.csv')), row.names = FALSE)
      }

    }
    
    
    ## (5) use SVM~hyperplane to identify the connected transcripts group based on tLLR score ----
    # SVM can separate continuous low score transcript from the rest.
    # but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
    
    transcript_df[['intraC']][['SVM_class']] <- 1
    
    if(length(flagged_cells)<1){
      message("No cells being flagged for resegmentation, no SVM is performed on this dataset.")
    }else {
      classDF_ToFlagTrans <- transcript_df[['intraC']][which(transcript_df[['intraC']][['UMI_cellID']] %in% flagged_cells),]
      
      if(!is.null(seed_transError)){
        set.seed(seed_transError)
      }
      # `flag_bad_transcripts` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
      tmp_df <- flag_bad_transcripts(chosen_cells = flagged_cells,
                                     score_GeneMatrix = prep_res[['score_GeneMatrix']],
                                     transcript_df = classDF_ToFlagTrans, 
                                     cellID_coln = 'UMI_cellID', 
                                     transID_coln = 'UMI_transID', 
                                     score_coln = 'score_tLLR_maxCellType',
                                     spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                     model_cutoff = flagModel_TransNum_cutoff, 
                                     score_cutoff = svmClass_score_cutoff, 
                                     svm_args = svm_args)
      if(!is.null(tmp_df)){
        # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
        classDF_ToFlagTrans <- merge(classDF_ToFlagTrans, 
                                     as.data.frame(tmp_df)[, c('UMI_transID','DecVal','SVM_class','SVM_cell_type')], 
                                     by = 'UMI_transID')
        
        message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                        length(flagged_cells) - length(unique(classDF_ToFlagTrans[['UMI_cellID']])), svmClass_score_cutoff))
        rm(tmp_df)
        
        # write into disk
        if(transDF_export_option ==1){
          write.csv(classDF_ToFlagTrans, file = fs::path(path_to_output, paste0(idx, '_classDF_ToFlagTrans.csv')), row.names = FALSE)
        }

        
        # flagged transcript ID, character vector
        flaggedSVM_transID <- classDF_ToFlagTrans[classDF_ToFlagTrans[['SVM_class']] ==0, 'UMI_transID']
        # assign SVM_class =0 for transcripts with low goodness-of-fit
        transcript_df[['intraC']][['SVM_class']] <- 1- as.numeric(transcript_df[['intraC']][['UMI_transID']] %in% flaggedSVM_transID)
      } 
      
    } 
    
    
    # intracellular vs extracelluar compartment 
    transcript_df[['intraC']][['transComp']] <- 'intraC' 
    
    # holder for per fov results to return 
    res_to_return <- list()

    # get gene x cell expression matrix after trimming all flagged transcripts
    if(return_trimmed_perCell){
      res_to_return <- transDF_to_perCell_data(
        transcript_df = data.table::as.data.table(transcript_df[['intraC']])[SVM_class == 1, ], 
        transGene_coln = 'target',
        cellID_coln = "UMI_cellID",
        spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
        celltype_coln = "SVM_cell_type",
        return_cellMeta = FALSE)
      
      ## impute zero value for genes not in `res_to_return$perCell_expression` but in `score_GeneMatrix` 
      missingGenes <- setdiff(rownames(prep_res[['score_GeneMatrix']]), 
                              rownames(res_to_return$perCell_expression))
      if(length(missingGenes)>0){
        mockExprs <- matrix(0, nrow = length(missingGenes), 
                            ncol = ncol(res_to_return$perCell_expression), 
                            dimnames = list(missingGenes, 
                                            colnames(res_to_return$perCell_expression)))
        mockExprs <- Matrix::Matrix(mockExprs, sparse = TRUE)
        res_to_return$perCell_expression <- rbind(res_to_return$perCell_expression, 
                                                  mockExprs)
      }
    }
    

    
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
    if(transDF_export_option ==1){
      write.csv(transcript_df, 
                file = fs::path(path_to_output, paste0(idx, "_flagged_transDF.csv")), 
                row.names = FALSE)
    }
    
    
    
    # return only idx, perCell data, flagged_cells and optional per cell expression after trimming as a list
    res_to_return <- c(list(file_idx = idx, 
                            flagged_cells = flagged_cells, 
                            modStats_ToFlagCells = modStats_ToFlagCells),
                       res_to_return)
    
    if(transDF_export_option ==2){
      res_to_return[["flagged_transDF"]] <- transcript_df
    }


    return(res_to_return)
  }
  
  # processing each FOV in parallel
  process_outputs <- parallel::mclapply(X = seq_len(nrow(transDF_fileInfo)), 
                                        mc.allow.recursive = TRUE,
                                        mc.cores = numCores(percentCores = percentCores),
                                        FUN = myFun_flag_eachFOV)
  
  # combine `modStats_ToFlagCells` data for each FOV
  combined_modStats_ToFlagCells <- lapply(process_outputs, '[[', 'modStats_ToFlagCells')
  combined_modStats_ToFlagCells <- do.call(rbind, 
                                           combined_modStats_ToFlagCells[!sapply(combined_modStats_ToFlagCells, is.null)])
  
  all_segRes[['combined_modStats_ToFlagCells']] <- combined_modStats_ToFlagCells
  all_segRes[['combined_flaggedCells']] <- lapply(process_outputs, '[[', 'flagged_cells')
  
  # per cell gene expression in gene x cell sparse matrix format, do cbind
  trimmed_perCellExprs <- lapply(process_outputs, '[[', 'perCell_expression')
  trimmed_perCellExprs <- do.call(cbind, trimmed_perCellExprs)
  all_segRes[['trimmed_perCellExprs']] <- trimmed_perCellExprs
  
  if(transDF_export_option ==2){
    all_segRes[['flagged_transDF_list']] <- lapply(process_outputs, '[[', 'flagged_transDF')
  }
  
  rm(process_outputs)
  
  return(all_segRes)
  
}

