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


