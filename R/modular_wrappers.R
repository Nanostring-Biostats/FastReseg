# modular wrapper to preprocess entire dataset to get baseline and cutoffs
#' @title runPreprocess
#' @description modular wrapper to get baseline data and cutoffs from entire dataset
#' @param counts Counts matrix for entire dataset, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, when NULL to automatically assign the cell cluster for each cell based on maximum transcript score of given the provided `refProfiles`
#' @param refProfiles A matrix of cluster profiles, genes X clusters, default = NULL to use external cluster assignments.
#' Of note, when `refProfiles !=  NULL`, genes unique to `counts` but missing in `refProfiles` would be omitted from downstream analysis. 
#' @param score_baseline a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence; default = NULL to calculate from `counts` and `refProfiles` 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is; default = NULL to calculate from `counts` and `refProfiles` 
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type; default = NULL to calculate from `counts` and `refProfiles` 
#' @param imputeFlag_missingCTs flag to impute `score_baseline`, `lowerCutoff_transNum`,`higherCutoff_transNum` for cell types present in `refProfiles` but missing in the provided transcript data files or the provided baseline and cutoffs; when TRUE, the median values of existing cell types would be used as the values for missing cell types.
#' @param ctrl_genes a vector of control genes that are present in input transcript data.frame but not in `refProfiles` and expect no cell type dependency, e.g. negative control probes; the `ctrl_genes` would be included in FastReseg analysis. (default = NULL)
#' @param svmClass_score_cutoff the cutoff of transcript score to separate between high and low score transcripts in SVM, used as the score values for `ctrl_genes` (default = -2)
#' @param cellular_distance_cutoff maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. Default = NULL to use the 2 times of average 2D cell diameter.
#' @param molecular_distance_cutoff maximum molecule-to-molecule distance within connected transcript group, unit in micron (default = 2.7 micron). 
#' If set to NULL, the pipeline would first randomly choose no more than 2500 cells from up to 10 random picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The pipeline would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @param transcript_df the data.frame of transcript level information with unique CellId, default = NULL to read from the `transDF_fileInfo`
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire dataset; when NULL, use the provided `transcript_df` directly
#' @param filepath_fov_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
#' @param prefix_colns the column names of annotation in `transDF_fileInfo`, to be added to the CellId as prefix when creating unique cell_ID for entire dataset; set to NULL if use the original `transID_coln` or `cellID_coln` 
#' @param fovOffset_colns the column name of coordinate offsets in 1st and 2nd dimension for each per FOV transcript data.frame in `transDF_fileInfo`, unit in micron
#' Notice that some assays like SMI has XY axes swapped between stage and each FOV such that `fovOffset_colns` should be c("stage_Y", "stage_X").
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of each `transcript_df`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of each `transcript_df`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in each `transcript_df`; when `prefix_colns` != NULL, unique transcript_id would be generated from `prefix_colns` and `transID_coln` in each `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_colns` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @return a nested list 
#' \describe{
#'    \item{clust}{vector of cluster assignments for each cell in `counts`, used in caculating `baselineData`}
#'    \item{refProfiles}{a genes X clusters matrix of cluster-specific reference profiles to use in resegmenation pipeline}
#'    \item{baselineData}{a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.}
#'    \item{cutoffs_list}{a list of cutoffs to use in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`}
#'    \item{ctrl_genes}{a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.}
#'    \item{score_GeneMatrix}{a gene x cell-type score matrix to use in resegmenation pipeline, the scores for `ctrl_genes` are set to be the same as `svmClass_score_cutoff`}
#'    \item{processed_1st_transDF}{a list of 2 elements for the intracellular and extracellular transcript data.frame of the processed outcomes of 1st transcrip file}
#' }
#' The `cutoffs_list` is a list containing
#' \describe{ 
#'    \item{score_baseline}{a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence.}
#'    \item{lowerCutoff_transNum}{a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is.}
#'    \item{higherCutoff_transNum}{a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.}
#'    \item{cellular_distance_cutoff}{maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. }
#'    \item{molecular_distance_cutoff}{maximum molecule-to-molecule distance within connected transcript group, unit in micron.}
#'    }
#' @examples  
#' data("mini_transcriptDF")
#' data("example_CellGeneExpr")
#' data("example_clust")
#' data("example_refProfiles")
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' 
#' # case 1: use `clust` and `transcript_df` directly, with known distance cutoffs
#' prep_res1 <- runPreprocess(
#'   counts = example_CellGeneExpr,
#'   clust = example_clust,
#'   refProfiles = NULL,
#'   score_baseline = NULL,
#'   lowerCutoff_transNum = NULL,
#'   higherCutoff_transNum= NULL,
#'   imputeFlag_missingCTs = FALSE,
#'   ctrl_genes = NULL,
#'   svmClass_score_cutoff = -2,
#'   molecular_distance_cutoff = 2.7,
#'   cellular_distance_cutoff = 20,
#'   transcript_df = mini_transcriptDF, 
#'   transDF_fileInfo = NULL, 
#'   pixel_size = 0.18,
#'   zstep_size = 0.8, 
#'   transID_coln = NULL,
#'   transGene_coln = "target",
#'   cellID_coln = 'CellId',
#'   spatLocs_colns = c('x','y','z'),
#'   extracellular_cellID = 0 
#' )
#' # case 2: use `refProfiles` to get `clust`, use `transcript_df` directly, unknown distance cutoffs
#' prep_res2 <- runPreprocess(
#'   counts = example_CellGeneExpr,
#'   clust = NULL,
#'   refProfiles = example_refProfiles,
#'   score_baseline = NULL,
#'   lowerCutoff_transNum = NULL,
#'   higherCutoff_transNum= NULL,
#'   imputeFlag_missingCTs = TRUE, # impute for cell types missing in provided 'transcript_df' 
#'   ctrl_genes = NULL,
#'   svmClass_score_cutoff = -2,
#'   molecular_distance_cutoff = NULL,
#'   cellular_distance_cutoff = NULL,
#'   transcript_df = mini_transcriptDF, 
#'   transDF_fileInfo = NULL, 
#'   pixel_size = 0.18,
#'   zstep_size = 0.8, 
#'   transID_coln = NULL,
#'   transGene_coln = "target",
#'   cellID_coln = 'CellId',
#'   spatLocs_colns = c('x','y','z'),
#'   extracellular_cellID = 0 
#' )
#' # case 3: provide both `refProfiles` and `clust`, use transDF_fileInfo for multi-files, no known molecular distance cutoffs
#' transDF_fileInfo <- data.frame(file_path = c("data/Run4104_FOV001__complete_code_cell_target_call_coord.csv",
#'                                              "data/Run4104_FOV002__complete_code_cell_target_call_coord.csv"),
#'                                slide = c(1, 1),
#'                                fov = c(1,2),
#'                                stage_X = 1000*c(5.13, -2.701),
#'                                stage_Y = 1000*c(-0.452, 0.081))
#' prep_res3 <- runPreprocess(
#'   counts = example_CellGeneExpr,
#'   clust = example_clust,
#'   refProfiles = example_refProfiles,
#'   score_baseline = NULL,
#'   lowerCutoff_transNum = NULL,
#'   higherCutoff_transNum= NULL,
#'   imputeFlag_missingCTs = TRUE,
#'   ctrl_genes = NULL,
#'   svmClass_score_cutoff = -2,
#'   molecular_distance_cutoff = NULL,
#'   cellular_distance_cutoff = 20,
#'   transcript_df = NULL, 
#'   transDF_fileInfo = transDF_fileInfo, 
#'   filepath_coln = 'file_path', 
#'   prefix_colns = c('slide','fov'), 
#'   fovOffset_colns = c('stage_X','stage_Y'), 
#'   pixel_size = 0.18,
#'   zstep_size = 0.8, 
#'   transID_coln = NULL,
#'   transGene_coln = "target",
#'   cellID_coln = 'CellId',
#'   spatLocs_colns = c('x','y','z'),
#'   extracellular_cellID = 0 
#' )
#' @export
#' 
runPreprocess <- function(counts, 
                          clust = NULL, 
                          refProfiles = NULL,
                          score_baseline = NULL, 
                          lowerCutoff_transNum = NULL, 
                          higherCutoff_transNum= NULL, 
                          imputeFlag_missingCTs = TRUE,
                          ctrl_genes = NULL,
                          svmClass_score_cutoff = -2,
                          molecular_distance_cutoff = 2.7,
                          cellular_distance_cutoff = NULL,
                          transcript_df = NULL, 
                          transDF_fileInfo = NULL, 
                          filepath_coln = 'file_path', 
                          prefix_colns = c('slide','fov'), 
                          fovOffset_colns = c('stage_X','stage_Y'), 
                          pixel_size = 0.18, 
                          zstep_size = 0.8,
                          transID_coln = NULL,
                          transGene_coln = "target",
                          cellID_coln = 'CellId', 
                          spatLocs_colns = c('x','y','z'), 
                          extracellular_cellID = NULL){
  
  
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
    clust <- baselineData[['clust_used']]
  }
  
  ## common genes between reference profiles and count matrix
  common_genes <- intersect(rownames(refProfiles), colnames(counts))
  if(length(common_genes) <2){
    stop("Too few common genes between the `refProfiles` (genes X clusters) and `counts` (cells X genes), check if correct format. ")
  }
  ## genes in `counts` but missing in `refProfiles`, would be dropped in `score_GeneMatrix` and exclude from downstream
  missing_genes <- setdiff(colnames(counts), common_genes)
  if(length(missing_genes)>0){
    message(sprintf("%d genes in `counts` but not present in the provided `refProfiles` would be omitted from downstream analysis.\nTo include those genes, please provide `clust` while set `refProfiles` to NULL.\nThe genes omited are: `%s`.", 
                    length(missing_genes), paste0(missing_genes, collapse = "`, `")))
  }
  
  all_celltypes <- unique(colnames(refProfiles))
  
  ## update and check format of baseline and cutoffs ----
  # update cutoffs for baseline only when not provided in inputs
  # check if any cell types not presented in cutoffs and baselines, but presents in refProfiles
  # fill in missing values based on median across all cell types when imputeFlag_missingCTs = TRUE
  myFun_updateCutoff <- function(name_cutoff, imputeFlag){
    val_cutoff <- get(name_cutoff)
    if(is.null(val_cutoff)){
      val_cutoff <- baselineData[[name_cutoff]]
    }else{
      # check validity of cutoffs
      if(!('numeric' %in% class(val_cutoff))){
        stop(sprintf("The provided `%s` must be a named numeric vector. ", name_cutoff))
      }
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
    
    val_cutoff <- val_cutoff[match(all_celltypes, names(val_cutoff))]
    
    return(val_cutoff)
  }
  
  score_baseline <- myFun_updateCutoff('score_baseline', imputeFlag_missingCTs)
  lowerCutoff_transNum <- myFun_updateCutoff('lowerCutoff_transNum', imputeFlag_missingCTs)
  higherCutoff_transNum <- myFun_updateCutoff('higherCutoff_transNum', imputeFlag_missingCTs)
  
  # make sure lowerCutoff is lower than higherCutoff such that the merge would overwrite the keep. 
  if(any(lowerCutoff_transNum > higherCutoff_transNum)){
    stop(sprintf("`lowerCutoff_transNum` is larger than `higherCutoff_transNum` in cell type: `%s`.", 
                 paste0(names(lowerCutoff_transNum)[lowerCutoff_transNum > higherCutoff_transNum], collapse = "`, `")))
  }
  
  ## get distance cutoff based on 1st provided FOV data ----
  if(is.null(molecular_distance_cutoff) | is.null(cellular_distance_cutoff)){
    message('Extract distance cutoff from first input transcript data.')
    
    
    # spatial dimension
    d2_or_d3 <- length(spatLocs_colns)
    
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
                                            extracellular_cellID = NULL, # already remove in `transcript_df[['intraC']]`
                                            run_molecularDist = TRUE,
                                            sampleSize_nROI = 10, 
                                            sampleSize_cellNum = 2500, 
                                            seed = 123)
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
                                              run_molecularDist = FALSE, 
                                              sampleSize_nROI = 10, 
                                              sampleSize_cellNum = 2500, 
                                              seed = 123)
        
        cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
        message(sprintf("Use `cellular_distance_cutoff` = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
        rm(distCutoffs)
      }
      
    }
  } else {
    if(!any(c("numeric", "integer") %in% class(c(molecular_distance_cutoff, cellular_distance_cutoff)))){
      stop("Must provide numeric values for `molecular_distance_cutoff` and cellular_distance_cutoff` in microns.")
    }
    message(sprintf('Use the providied `molecular_distance_cutoff` = %.4f for defining direct neighbor cells based on molecule-to-molecule distance.', molecular_distance_cutoff))
    message(sprintf("Use the providied `cellular_distance_cutoff` = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
    
    # no loading of 1st fovs, set to NULL
    transcript_df <- NULL
  }
  
  

  
  ## initialize list to collect outputs ----
  outs <- list()
  
  # record the final cutoffs and reference profiles in use for all FOVs
  outs[['clust']] <- clust
  outs[['refProfiles']] <- refProfiles
  outs[['baselineData']] <- list(span_score = baselineData[['span_score']], 
                                 span_transNum = baselineData[['span_transNum']])
  outs[['cutoffs_list']] <- list(score_baseline = score_baseline, 
                                 lowerCutoff_transNum = lowerCutoff_transNum, 
                                 higherCutoff_transNum = higherCutoff_transNum, 
                                 cellular_distance_cutoff = cellular_distance_cutoff, 
                                 molecular_distance_cutoff = molecular_distance_cutoff)
  
  rm(baselineData)
  

  
  ## get transcript score matrix for each gene based on reference profile 
  tLLR_geneMatrix <- scoreGenesInRef(genes = common_genes, ref_profiles = pmax(refProfiles, 1e-5))
  
  # set tLLR score for control genes, same as `svmClass_score_cutoff`
  if(!is.null(ctrl_genes)){
    message(sprintf("Include the following `ctrl_genes` in analysis: `%s`.\nIt's recommended to have total counts of those genes below 1%% of total counts of all genes in each cell.", 
                    paste0(ctrl_genes, collapse = "`, `")))
    
    outs[['ctrl_genes']] <- ctrl_genes
    
    
    if(any(ctrl_genes %in% rownames(refProfiles))){
      message(sprintf("Overwrite transcript score for %d `ctrl_genes` shared with `refProfiles`: `%s`.", 
                      sum(ctrl_genes %in% rownames(refProfiles)),
                      paste0(intersect(ctrl_genes, rownames(refProfiles)), collapse = "`, `")))
      
      tLLR_geneMatrix <- tLLR_geneMatrix[!(rownames(tLLR_geneMatrix) %in% ctrl_genes), ]
    }
    
    tLLR_geneMatrix <- rbind(tLLR_geneMatrix, 
                             matrix(svmClass_score_cutoff, 
                                    nrow = length(ctrl_genes), ncol = ncol(tLLR_geneMatrix),
                                    dimnames = list(ctrl_genes, colnames(tLLR_geneMatrix)))
    )
  }
  outs[['score_GeneMatrix']] <- tLLR_geneMatrix
  outs[['processed_1st_transDF']] <- transcript_df
  
  return(outs)
  
}


# modular wrapper to flag cell segmentation error
#' @title runSegErrorEvaluation
#' @description modular wrapper to flag cell segmentation error 
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param transcript_df the data.frame of transcript_ID, cell_ID, score, spatial coordinates
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param score_coln the column name of score in transcript_df
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in transcript_df 
#' @param flagModel_TransNum_cutoff the cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
#' @return a list of two elements 
#' #' \enumerate{
#'    \item{modStats_ToFlagCells, a data.frame contains evaluation model statistics in columns for each cell's potential to have segmentation error}
#'    \item{transcript_df, transcript data.frame with 2 additional columns: `tLLR_maxCellType` for cell types of maxmium transcript score under current segments and `score_tLLR_maxCellType` for the corresponding transcript score for each transcript}
#' } 
#' @examples 
#' data("mini_transcriptDF")
#' data("example_CellGeneExpr")
#' data("example_refProfiles")
#' score_GeneMatrix <- scoreGenesInRef(
#'   genes = intersect(colnames(example_CellGeneExpr), rownames(example_refProfiles)), 
#'   ref_profiles = pmax(example_refProfiles, 1e-5))
#' 
#' res <- runSegErrorEvaluation(
#'   score_GeneMatrix= score_GeneMatrix, 
#'   transcript_df = mini_transcriptDF, 
#'   cellID_coln = 'UMI_cellID', 
#'   transID_coln = 'UMI_transID',
#'   transGene_coln = 'target',
#'   spatLocs_colns = c('x','y','z'),
#'   #' cutoff of transcript number to do spatial modeling
#'   flagModel_TransNum_cutoff = 50) 
#' @export
runSegErrorEvaluation <- function(score_GeneMatrix, 
                                  transcript_df,
                                  cellID_coln = 'UMI_cellID', 
                                  transID_coln = 'UMI_transID',
                                  transGene_coln = 'target',
                                  spatLocs_colns = c('x','y','z'),
                                  flagModel_TransNum_cutoff = 50){
  
  ## for each cell, get new cell type based on maximum score ----
  # `getCellType_maxScore` function returns a list contains element `cellType_DF`, a data.frame with cell in row, cell_ID and cell_type in column.
  select_cellmeta <- getCellType_maxScore(score_GeneMatrix = score_GeneMatrix, 
                                          transcript_df = transcript_df, 
                                          transID_coln = transID_coln,
                                          transGene_coln = transGene_coln,
                                          cellID_coln = cellID_coln, 
                                          return_transMatrix = FALSE)
  
  select_cellmeta <- select_cellmeta[['cellType_DF']]
  colnames(select_cellmeta) <- c(cellID_coln,'tLLR_maxCellType')
  
  
  all_cells <- select_cellmeta[[cellID_coln]]
  
  transcript_df <- merge(transcript_df, select_cellmeta, by = cellID_coln)
  message(sprintf("Found %d cells and assigned cell type based on the provided `refProfiles` cluster profiles.", nrow(select_cellmeta)))
  
  
  ##  for each transcript, calculate tLLR score based on the max cell type
  # `getScoreCellType_gene` function returns a data.frame with transcript in row and "[transID_coln]" and "score_[celltype_coln]" in column for chosen cell-type
  tmp_df <- getScoreCellType_gene(score_GeneMatrix = score_GeneMatrix, 
                                  transcript_df = transcript_df, 
                                  transID_coln = transID_coln,
                                  transGene_coln = transGene_coln,
                                  celltype_coln = 'tLLR_maxCellType')
  transcript_df <- merge(transcript_df, tmp_df, by = transID_coln)
  rm(tmp_df)
  
  
  
  ## spatial modeling of tLLR score profile within each cell to identify cells with strong spatial dependency 
  # `score_cell_segmentation_error` function returns a data.frame with cell in row and spatial modeling outcomes in columns
  modStats_ToFlagCells <- score_cell_segmentation_error(
    chosen_cells = all_cells, 
    transcript_df = transcript_df, 
    cellID_coln = cellID_coln, 
    transID_coln = transID_coln, 
    score_coln = 'score_tLLR_maxCellType',
    spatLocs_colns = spatLocs_colns, 
    model_cutoff = flagModel_TransNum_cutoff)
  
  if(!is.null(modStats_ToFlagCells)){
    #-log10(P)
    modStats_ToFlagCells[['lrtest_-log10P']] <- -log10(modStats_ToFlagCells[['lrtest_Pr']])
    modStats_ToFlagCells <- merge(select_cellmeta, modStats_ToFlagCells, by.x = cellID_coln, by.y = 'cell_ID')
    
  }
  
  outs <- list(modStats_ToFlagCells = modStats_ToFlagCells, 
               transcript_df = transcript_df)
  
  return(outs)
}



# modular wrapper to identify transcript groups of poor fit to current cell segments in space
#' @title runTranscriptErrorDetection
#' @description modular wrapper to identify transcript groups of poor fit to current cell segments in space
#' @param chosen_cells the cell_ID of chosen cells
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param transcript_df the data.frame of transcript_ID, cell_ID, score, spatial coordinates
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param score_coln the column name of score in transcript_df
#' @param spatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param model_cutoff the cutoff of transcript number to do spatial modeling (default = 50)
#' @param score_cutoff the cutoff of score to separate between high and low score transcripts (default = -2)
#' @param svm_args a list of arguments to pass to svm function, typically involve kernel, gamma, scale
#' @param groupTranscripts_method use either "dbscan" or "delaunay method" to group transcripts in space (default = "dbscan")
#' @param distance_cutoff maximum molecule-to-molecule distance within same transcript group (default = "auto")
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config
#' @param seed_transError seed for transcript error detection step, default = NULL to skip the seed   
#' @return data frame for transcripts in `chosen_cells` only, containing information for transcript score classifications and spatial group assignments as well as new cell/group ID for downstream resegmentation.
#' @export
#' 
runTranscriptErrorDetection <- function(chosen_cells, 
                                        score_GeneMatrix, 
                                        transcript_df, 
                                        cellID_coln = "CellId", 
                                        transID_coln = "transcript_id",
                                        transGene_coln = "target",
                                        score_coln = "score", 
                                        spatLocs_colns = c("x","y","z"),
                                        model_cutoff = 50, 
                                        score_cutoff = -2, 
                                        svm_args = list(kernel = "radial", 
                                                        scale = FALSE, 
                                                        gamma = 0.4),
                                        groupTranscripts_method = c("dbscan", "delaunay"),
                                        distance_cutoff = 'auto',
                                        config_spatNW_transcript = NULL, 
                                        seed_transError = NULL){
  
  groupTranscripts_method <- match.arg(groupTranscripts_method, c("dbscan", "delaunay"))
  if(groupTranscripts_method == 'delaunay'){
    # check `config_spatNW_transcript` and return default config if not provided 
    config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript, 
                                                       spat_locs = spatLocs_colns)
  }
  
  # transcript data.frame for flagged cells only, including column for transcript scores
  classDF_ToFlagTrans <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells),]
  
  if(nrow(classDF_ToFlagTrans)<1) {
    stop("Error: No transcripts within `chosen_cells`.")
  }
  
  if(!is.null(seed_transError)){
    set.seed(seed_transError)
  }
  
  # `flagTranscripts_SVM` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
  flagged_transDF_SVM <- flagTranscripts_SVM(chosen_cells = chosen_cells,
                                             score_GeneMatrix = score_GeneMatrix,
                                             transcript_df = classDF_ToFlagTrans, 
                                             cellID_coln = cellID_coln, 
                                             transID_coln = transID_coln, 
                                             transGene_coln = transGene_coln,
                                             score_coln =  score_coln,
                                             spatLocs_colns = c("x", "y", "z"), 
                                             model_cutoff = model_cutoff, 
                                             score_cutoff = score_cutoff, 
                                             svm_args = svm_args)
  
  
  # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
  flagged_transDF_SVM <- merge(classDF_ToFlagTrans, 
                               as.data.frame(flagged_transDF_SVM)[, c(transID_coln,'DecVal','SVM_class','SVM_cell_type')], 
                               by = transID_coln)
  
  message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                  length(chosen_cells) - length(unique(classDF_ToFlagTrans[[cellID_coln]])), score_cutoff))
  
  # flagged transcript ID, character vector
  flaggedSVM_transID <- flagged_transDF_SVM[flagged_transDF_SVM[['SVM_class']] ==0, transID_coln]
  
  # group the low-score transcripts in space using either `groupTranscripts_dbscan` or `groupTranscripts_Delaunay` function 
  # returns a data.frame of connected transcripts among chosen_transcripts, with each transcript in row, 
  # the group ID for the connected transcript groups and the original cell ID, spatial coordinates in column.
  
  if (groupTranscripts_method == 'dbscan'){
    flaggedSVM_transGroupDF <- groupTranscripts_dbscan(chosen_transcripts = flaggedSVM_transID, 
                                                       distance_cutoff = distance_cutoff, # for dbscan grouping and orphan transcript determination
                                                       transcript_df = classDF_ToFlagTrans, 
                                                       cellID_coln = cellID_coln, 
                                                       transID_coln = transID_coln,
                                                       transSpatLocs_coln = spatLocs_colns)
  } else if(groupTranscripts_method == 'delaunay'){
    flaggedSVM_transGroupDF <- groupTranscripts_Delaunay(chosen_transcripts = flaggedSVM_transID, 
                                                         config_spatNW_transcript = config_spatNW_transcript, # define spatial network including grouping in delaunay network 
                                                         distance_cutoff = distance_cutoff, # for orphan transcript determination only
                                                         transcript_df = classDF_ToFlagTrans,
                                                         cellID_coln = cellID_coln, 
                                                         transID_coln = transID_coln,
                                                         transSpatLocs_coln = spatLocs_colns)
  } 
  
  message(sprintf("SVM spatial model further identified %d cells with transcript score all in same class, exclude from transcript group analysis.", 
                  length(unique(flagged_transDF_SVM[[cellID_coln]])) - length(unique(flaggedSVM_transGroupDF[[cellID_coln]]))))
  
  ## (3.2) generate tmp_cellID include group information and get max cell type for each group ----
  # assign group ID name for all transcripts
  # transcripts with SVM class = 1 are high-score molecules and would get `connect_group` = 0 and thus keep the original cell ID as the tmp_cellID
  # transcripts with SVM class = 0 are low-score molecules and would get `connect_group` to be same value as the group ID identified based on spatial network analysis and tmp_cellID modified from original cell_ID based on the `connect_group` value
  flagged_transDF_SVM[['connect_group']] <- 1- as.numeric(as.character(flagged_transDF_SVM[['SVM_class']]))
  group_converter <- flaggedSVM_transGroupDF[['transcript_group']] 
  names(group_converter) <- flaggedSVM_transGroupDF[[transID_coln]]
  
  tmp_idx <- which(flagged_transDF_SVM[[transID_coln]] %in% flaggedSVM_transGroupDF[[transID_coln]])
  flagged_transDF_SVM[['connect_group']][tmp_idx] <- group_converter[flagged_transDF_SVM[[transID_coln]][tmp_idx]]
  
  rm(tmp_idx, group_converter)
  
  # get new cell_id and cell type for each group
  flagged_transDF_SVM <- data.table::as.data.table(flagged_transDF_SVM)
  flagged_transDF_SVM[, tmp_cellID := ifelse(connect_group == 0, get(cellID_coln), paste0(get(cellID_coln),'_g', connect_group))]
  
  # get new cell type of each group based on maximum
  # `getCellType_maxScore` function returns a list contains element `cellType_DF`, a data.frame with cell in row, cell_ID and cell_type in column.
  tmp_df <- getCellType_maxScore(score_GeneMatrix = score_GeneMatrix, 
                                 transcript_df = flagged_transDF_SVM, 
                                 transID_coln = transID_coln,
                                 transGene_coln = transGene_coln,
                                 cellID_coln = "tmp_cellID", 
                                 return_transMatrix = FALSE)
  colnames(tmp_df[['cellType_DF']]) <- c('tmp_cellID','group_maxCellType')
  flagged_transDF_SVM <- merge(flagged_transDF_SVM, tmp_df[['cellType_DF']], by = 'tmp_cellID', all.x = TRUE)
  flagged_transDF_SVM <- as.data.frame(flagged_transDF_SVM)
  rm(tmp_df)
  
  
  return(flagged_transDF_SVM)
}





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
#'    \item{neighborhoodDF_ToReseg}{a data.frame for neighborhood enviornment of low-score transcript groups, output of `neighborhood_for_resegment_spatstat` function, return when `return_intermediates` = TRUE}
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
  # `neighborhood_for_resegment_spatstat` function returns a data.frame with each cell in row and its neighborhood information in columns
  neighborReSeg_df <- neighborhood_for_resegment_spatstat(chosen_cells = chosen_cells,
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
  # `getScoreCellType_gene` function returns a data.frame with transcript in row and "[transID_coln]" and "score_[celltype_coln]" in column for chosen cell-type
  tmp_df <- getScoreCellType_gene(score_GeneMatrix = score_GeneMatrix, 
                                  transcript_df = post_reseg_results$updated_transDF, 
                                  transID_coln = transID_coln,
                                  transGene_coln = transGene_coln,
                                  celltype_coln = 'updated_celltype')   
  
  post_reseg_results$updated_transDF <- merge(post_reseg_results$updated_transDF, tmp_df, by = transID_coln)
  rm(tmp_df)
  
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
