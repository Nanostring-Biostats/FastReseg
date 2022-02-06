# prepare resegmentation inputs from existing SMI object with multi-FOV multi-slide data
#' @title prepSMI_for_fastReseg
#' @description prepare resegmentation inputs from existing SMI object with multi-FOV multi-slide data 
#' @param path_to_SMIobject file path to \code{"Giotto"} object for existing SMI multi-FOV multi-slide data that was save in `.RData` format
#' @param config_loading a list holding arguments passed to `SMITAP::data_loading_from_config` when creating SMI object; this list should contain file path information for the raw data used for creating the corresponding multi-slide multi-FOV object.
#' @param refClus_coln the column name of cell cluster assignment saved in `cell_metadata` of the existing SMI object
#' @param cellClus_to_exclude a vector of cell cluster name to be excluded for calculation of reference profiles (default = NULL)
#' @param removeUnpaired flag to remove FOVs with unpaired target call files and fov position information; default = FALSE, to stop processing when missing target call files
#' @param blacklist_genes a vector of genes to be excluded from reference profile estimation (default = NULL)
#' @param pixel_size the micrometer size of image pixel listed in `Width` and `Height` dimension of each cell stored in `cell_metadata` of the existing SMI object (default = 0.18)
#' @return a list 
#' \describe{
#'    \item{counts}{a cells * genes count matrix for entire dataset, stored in SMI object.}
#'    \item{clust}{vector of cluster assignments for each cell in `counts`, stored in SMI object.}
#'    \item{refProfiles}{a genes * clusters matrix of cluster-specific reference profiles estimated from entire dataset.}
#'    \item{score_baseline}{a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence.}
#'    \item{lowerCutoff_transNum}{a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is.}
#'    \item{higherCutoff_transNum}{a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.}
#'    \item{cellular_distance_cutoff}{maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. Use the 2 times of average 2D cell diameter.}
#'    \item{transDF_fov_fileInfo}{a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire dataset, in colnames of `file_path`, `slide`,`slideName`,`fov`,`offset_x`, `offset_y`.}
#'    \item{sample_annot}{a data.frame of sample annotation with each row for each slide, columns for file path and metedata that are read from file path of `config_loading[["annotfile"]]`.}
#' }
#' @details The function requires the pre-built \code{"Giotto"} object that contains single cell typing results on original cell segmentation outcomes. This pre-built object may contain multi-slide multi-FOV data and should be built with raw data whose file path information is stored in `config_loading`.
#' `config_loading` should contain at least the following 6 elements:
#' \describe{
#'    \item{annotfile}{file path to the sample annotation file}
#'    \item{folderpathColumn}{Column name for full file path to run folder}
#'    \item{slidefoldersColumn}{Column name for base file path to slide folder in each run folder}
#'    \item{slidenameColumn}{Column name for slide_name}
#'    \item{votedfoldersColumn}{Column name for base file path to voted folders in each slide folder}
#'    \item{versionColumn}{Column name for the version of the target calling pipeline}
#' }
#' The function would first calculate reference profiles, `refProfiles`, and several cutoffs from the entire dataset, including the cutoffs for transcript number, `lowerCutoff_transNum` and `higherCutoff_transNum`, the cutoff for transcript tLLR score `score_baseline`, and the center-to-center distance cutoff for neighbor cells with direct contact, `cellular_distance_cutoff`.
#' Then, the function would check the existence of fov offset position and transcript data.frame for each FOV in each slide. 
#' Lastly, the function would prepare `transDF_fov_fileInfo` to be used with `fastReseg_internalRef` function for downstream resegmentation workflow which would use either NULL or a user-defined value for `molecular_distance_cutoff`, the molecule-to-molecule distance cutoff for neighbor molecules belonging to same cell. When `molecular_distance_cutoff` is set to NULL in downstream `fastReseg_internalRef` function, `molecular_distance_cutoff` would then be estimated from transcript data for the first FOV transcript data.frame. 
#' @importFrom Giotto pDataDT
#' @export 
prepSMI_for_fastReseg <- function(path_to_SMIobject, 
                                  config_loading, 
                                  refClus_coln = 'nb_clus',
                                  cellClus_to_exclude = NULL, 
                                  removeUnpaired = FALSE,
                                  blacklist_genes = NULL,
                                  pixel_size = 0.18){
  #### (1) prepare sample annotation file ----
  ## check config_loading
  colns_to_use <- c('folderpathColumn','slidefoldersColumn','slidenameColumn','votedfoldersColumn','versionColumn')
  if(!all(c(colns_to_use, 'annotfile') %in% names(config_loading))){
    stop(sprintf("The provided `config_loading` is missing elements for: `%s`.", 
                 paste0(setdiff(c(colns_to_use, 'annotfile'), names(config_loading)), 
                        collapse = '`,`')))
  }
  
  if(!file.exists(config_loading[['annotfile']])){
    stop(sprintf("File does not exist for `annotfile` in `config_loading`: %s", 
                 config_loading[['annotfile']]))
  }
  
  ## get sample annotaiton file
  ori_SampleAnnot <- read.csv(config_loading[['annotfile']])
  # keep a copy of original annotaiton file for results return
  
  ## clean up sample annotation file and rename colnames
  sample_annot <- ori_SampleAnnot 
  loading_colns <- unlist(config_loading)[names(config_loading) %in% colns_to_use]
  loading_colns <- loading_colns[match(colns_to_use, names(loading_colns))]
  sample_annot <- sample_annot[, loading_colns]
  colnames(sample_annot) <- colns_to_use
  
  sample_annot[['fullSlidePath']] <- fs::path(sample_annot[['folderpathColumn']], sample_annot[['slidefoldersColumn']])
  sample_annot[['fullVotedPath']] <- fs::path(sample_annot[['fullSlidePath']], sample_annot[['votedfoldersColumn']])
  sample_annot[['fovOffsetPath']] <- fs::path(sample_annot[['fullSlidePath']], "RunSummary","latest.fovs.csv")
  sample_annot[['versionColumn']] <- tolower(sample_annot[['versionColumn']])
  
  ## check if those files exist
  msg <- unlist(lapply(c('fullSlidePath', 'fullVotedPath', 'fovOffsetPath'), 
                       function(x) checkFileExists(sample_annot, x)))
  
  if ( length( msg ) > 0L ) {
    stop( "Annotation file issue:\n" , paste( msg , collapse = "\n" ) )
  }
  
  
  #### (2) get cell x gene expression matrix and cell typing info from SMI object ----
  ## load pre-built SMI giotto object
  if(grepl(".RData$", path_to_SMIobject, ignore.case = TRUE)){
    if(!file.exists(path_to_SMIobject)){
      stop(sprintf("`path_to_SMIobject` does NOT exist: %s", 
                   path_to_SMIobject))
    }
    
    gem <- get(load(path_to_SMIobject))
    
    if(!any(class(gem) %in% c('giotto'))){
      stop(sprintf("Data stored in `path_to_SMIobject` is not a giotto object, current path: %s",
                   path_to_SMIobject))
    }
    
  } else {
    stop(sprintf("`path_to_SMIobject` must be `.RData` file, current path: %s", 
                 path_to_SMIobject))
  }
  
  
  ## get cell typing and per cell info
  cell_metadata <- Giotto::pDataDT(gem)
  # cell_ID = c_[slide]_[FOV]_CellId
  celltype_metadata <- data.frame(cell_ID = cell_metadata$cell_ID, 
                                  FOV = cell_metadata$fov, 
                                  slide = cell_metadata$slide_ID_numeric, 
                                  slideName = cell_metadata[[config_loading$slidenameColumn]],
                                  cell_type = cell_metadata[[refClus_coln]])
  celltype_metadata[['CellId']] <- sapply(strsplit(celltype_metadata$cell_ID, '_'),'[[',4)
  
  ## get rna target only raw expression matrix, gene x cell
  exprs_tgrt <- Giotto:::select_expression_values(gem, feat_type = 'rna', values = 'raw')
  targets <- rownames(exprs_tgrt)
  
  # exclude blacklist genes from reference profile estimation
  if(!is.null(blacklist_genes)){
    targets <- setdiff(targets, blacklist_genes)
    exprs_tgrt <- exprs_tgrt[targets, ]
  }
  
  ## slide_ID_numeric in cell_metadata and cell_ID must align with sample_annot file
  sample_annot <- merge(sample_annot, 
                        unique(celltype_metadata[, c('slide','slideName')]), 
                        by.x = 'slidenameColumn', by.y = 'slideName', all.x = TRUE)
  
  if(any(is.na(sample_annot[['slide']]))){
    # some sample_annot does not present in giotto object
    message(sprintf("slide name = `%s` in annotation file is not included in `gem` object and thus excluded from `refProfiles` esimation.", 
                    paste0(sample_annot[['slidenameColumn']][is.na(sample_annot[['slide']])], collapse = "`, `")))
  }
  
  
  ## estimate `cellular_distance_cutoff`, maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. 
  ## Default to use the 2 times of average 2D cell diameter.
  if(all(c('Width','Height') %in% colnames(cell_metadata))){
    # giotto object `gem` has dimension of bounding box for each cell in `cell_metadata` already, unit in pixel
    cellsize_df <- cell_metadata[, c('Width', 'Height')]
    cellsize_df <- colMeans(cellsize_df)
    cellular_distance_cutoff <- 2*mean(cellsize_df)*pixel_size
    
    rm("cellsize_df")
  }else {
    ## if not exist, one could use the automatic estimation capacity of `fastReseg_externalRef` function
    cellular_distance_cutoff <- NULL
  }
  
  
  # remove additional information to free some space
  rm("gem","cell_metadata")
  
  
  #### (3) get reference profiles based on previous cell typing outcomes of entire dataset ----
  # remove cells with certain cell cluster assignment, e.g `NotDet`in case of nb_clus
  if(!is.null(cellClus_to_exclude)){
    tmp_idx <- which(!(grepl(paste0(unique(cellClus_to_exclude), collapse = '|'), celltype_metadata[["cell_type"]])))
    cells_to_keep <- celltype_metadata[['cell_ID']][tmp_idx]
    message(sprintf("Remove %d cells, %.2f%% of all, with cell type assignment of `%s` in `%s` from reference profile estimation. ", 
                    nrow(celltype_metadata) - length(cells_to_keep), 
                    100-length(cells_to_keep)/nrow(celltype_metadata)*100, 
                    paste0(unique(cellClus_to_exclude), collapse = '`, `'), 
                    refClus_coln))
    rm(tmp_idx)
  } else {
    cells_to_keep <- celltype_metadata[["cell_ID"]]
  }
  
  ## get same cell order for cell x gene counts matrix and cluster assignment
  # Counts matrix for entire dataset, cells * genes.
  counts <- Matrix::t(exprs_tgrt[, cells_to_keep])
  # Vector of cluster assignments for each cell in `counts`
  clust <- celltype_metadata[celltype_metadata[['cell_ID']] %in% cells_to_keep, c('cell_ID','cell_type')]
  clust <- clust [match(cells_to_keep, clust[['cell_ID']]),]
  clust <- clust[['cell_type']]
  
  # A matrix of cluster profiles, genes * clusters
  # ignore background, use total count per cell as scaling factor
  refProfiles <- estimate_MeanProfile( counts = as.matrix(counts), 
                                       clust = as.character(clust), 
                                       s = Matrix::rowSums(as.matrix(counts)), 
                                       bg = rep(0, nrow(counts)))
  # remove unused variables
  rm("exprs_tgrt","celltype_metadata")
  
  
  #### (4) get transcript number and tLLR transcript score cutoff based on entire dataset ----
  # # `get_baselineCT` function gets cluster-specific quantile distribution of transcript number and per cell per molecule transcript score in the provided cell x gene expression matrix based on the reference profiles and cell cluster assignment. 
  # # The function returns a list containing the following elements:
  # span_score, a matrix of average transcript tLLR score per molecule per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns
  # span_transNum, a matrix of transcript number per cell for each distinct cell types in row, percentile at (0%, 25%, 50%, 75%, 100%) in columns
  # score_baseline, a named vector of 25% quantile of cluster-specific per cell transcript score, to be used as score baseline such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence
  # lowerCutoff_transNum, a named vector of 25% quantile of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that higher than the cutoff is required to keep query cell as it is
  # higherCutoff_transNum, a named vector of median value of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
  # clust_used,  a named vector of cluster assignments for each cell used in baseline calculation, cell_ID in `counts` as name
  baselineData <- get_baselineCT(refProfiles = refProfiles, counts = counts, clust = as.character(clust))
  
  score_baseline <- baselineData[['score_baseline']]
  lowerCutoff_transNum <- baselineData[['lowerCutoff_transNum']]
  higherCutoff_transNum <- baselineData[['higherCutoff_transNum']]
  
  # remove unused variables
  rm("baselineData")
  
  
  
  #### (5) check existence of per FOV transcript data.frame and prepare data.frame for FOV file information ----
  transDF_fov_fileInfo <- list()
  
  for (idx in seq_len(nrow(sample_annot))){
    # locate transcript data.frame (target call files) for each FOV
    if ( sample_annot[['versionColumn']][idx] != "v1" ){
      exprs_files_list <- dir(path = sample_annot[['fullVotedPath']][idx], full.names = TRUE, recursive = TRUE, 
                              pattern = "^Run[0-9]+_FOV[0-9]+_?_complete_code_cell_target_call_coord.csv")
    } else {
      exprs_files_list <- dir(path = sample_annot[['fullVotedPath']][idx], full.names = TRUE, recursive = TRUE, 
                              pattern = "SUMMARY_Run[0-9]+_(.*)min[0-9]+_FOV[0-9]+_.xlsx")
    }
    
    if ( length(exprs_files_list[[idx]])==0 ){
      stop("No target call files are found.")
    }
    
    # extract FOV number from file name
    exprs_FOV <- sapply(unlist(stringr::str_extract_all(basename(exprs_files_list), "FOV[:digit:]+")), 
                        function(x){as.numeric(gsub("FOV", "", x))})
    
    fileInfo_to_use <- data.frame(file_path = exprs_files_list,
                                  fov = exprs_FOV)
    fileInfo_to_use[['slide']] <- sample_annot[['slide']][idx]
    fileInfo_to_use[['slideName']] <- sample_annot[['slidenameColumn']][idx]
    
    # fov offset position
    fov_offsets <- read.csv(sample_annot[['fovOffsetPath']][idx], header = FALSE, 
                            col.names = c("SlidePerRun", "X_mm", "Y_mm", "Z_mm", "ZOffset_mm", 
                                          "ROI", "fov"))
    # SMI assay has flipped axes between fov coordinates on stage and the pixel coordinates within individual FOV image
    # flip the xy axes to align stage coordinates to image coordinates
    fov_offsets[['offset_y']] <- fov_offsets[['X_mm']]*1000
    fov_offsets[['offset_x']] <- fov_offsets[['Y_mm']]*1000
    
    # add fov offset to file info data.frame
    fileInfo_to_use <- merge(fileInfo_to_use, fov_offsets[, c('fov','offset_x','offset_y')], by = 'fov')
    
    # warning for FOVs without expression files
    if(length(setdiff(fov_offsets[['fov']], exprs_FOV))>0){
      message(sprintf("Target Call Files for `%s` slide is missing for FOV = %s. ", 
                      sample_annot[['slidenameColumn']][idx], 
                      paste0(as.character(setdiff(fov_offsets[['fov']], exprs_FOV)), 
                             collapse = ", ")))
      if(removeUnpaired){
        message(sprintf("Proceed with remaining Target Call Files for FOV = %s. ", 
                        paste0(as.character(fileInfo_to_use[['fov']]), collapse = ", ")))
      } else {
        stop("Stop further processing, consider set `removeUnpaired` to TRUE to ignore those missing files.")
      }
    }
    
    transDF_fov_fileInfo[[idx]] <- fileInfo_to_use
    
  }
  
  transDF_fov_fileInfo <- do.call(rbind, transDF_fov_fileInfo)
  
  
  # # return results ----
  # counts: a cells * genes count matrix for entire dataset, stored in SMI object.
  # clust: vector of cluster assignments for each cell in `counts`, stored in SMI object.
  # refProfiles: a genes * clusters matrix of cluster-specific reference profiles estimated from entire dataset.
  # score_baseline: a named vector of score baseline under each cell type listed in `refProfiles` such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence.
  # lowerCutoff_transNum: a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is.
  # higherCutoff_transNum: a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
  # cellular_distance_cutoff: maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. Use the 2 times of average 2D cell diameter.
  # transDF_fov_fileInfo: a data.frame with each row for each individual file of per FOV transcript data.frame, columns for `file_path`, `slide`,`slideName`,`fov`,`offset_x`, `offset_y`.
  # sample_annot: a data.frame of sample annotation with each row for each slide, columns for file path and metedata that are read from file path of `config_loading[["annotfile"]]`.} 
  
  res <- list(counts = counts, # cells x genes count matrix stored in SMI object
              clust = clust, # cluster assignments stored in SMI object
              refProfiles = refProfiles, 
              score_baseline = score_baseline, 
              lowerCutoff_transNum = lowerCutoff_transNum,
              higherCutoff_transNum = higherCutoff_transNum, 
              cellular_distance_cutoff = cellular_distance_cutoff, 
              transDF_fov_fileInfo = transDF_fov_fileInfo)
  
  # clean up information stored in the original `annotfile`
  colns_to_remove <- c('folderpathColumn','slidefoldersColumn','votedfoldersColumn')
  loading_colns <- unlist(config_loading)[names(config_loading) %in% colns_to_remove]
  ori_SampleAnnot <- ori_SampleAnnot[, - which(colnames(ori_SampleAnnot) %in% loading_colns)]
  # add in slide_ID_numeric information stored in SMI giotto object
  slide_NameID <- sample_annot[, c('slide', 'slidenameColumn')]
  colnames(slide_NameID) <- c('slide_ID_numeric', config_loading[['slidenameColumn']])
  ori_SampleAnnot <- merge(slide_NameID, ori_SampleAnnot,
                        by = config_loading[['slidenameColumn']], all.x = TRUE)
  # append to result list fore return
  res[['sample_annot']] <- ori_SampleAnnot
  
  return(res)
  
}

#' @title checkFileExists
#' @description function to check if file is missing in the Annotation data.frame
#' @param Annotation the data.frame storing the sample annotation file 
#' @param name the column name in Annotation
#' @return a message about which file is missing
checkFileExists <- function( Annotation, name ){
  if ( all(file.exists(Annotation[[name]])) ){
    msg <- character()
  } else {
    msg <- sprintf( "\nFolder or file `%s` doesn't exist.", 
                    Annotation[[name]][ which(!file.exists(Annotation[[name]]))], 
                    name)
  }
  return( msg )
}