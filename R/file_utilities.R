# read in single FOV data
#' @title myFun_fov_load
#' @description supporting function for \code{runPreprocess}, \code{fastReseg_full_pipeline} and \code{fastReseg_flag_all_errors}, to load a single field-of-view (FOV) transcript data file into a data frame from file path.
#' @param path_to_fov file path to per fov transcript data.frame 
#' @return A data frame containing the transcript data.
#' @details This function automatically detects the file format based on its extension and supports common tabular and binary formats, including 
#' \itemize{
#'   \item `.csv`, `.txt`, `.tsv`: tabular text files (comma- or tab-delimited)
#'   \item `.csv.gz`, `.txt.gz`, `.tsv.gz`: gzip-compressed tabular files
#'   \item `.RData`: R workspace files containing a single data frame object
#'   \item `.rds`: serialized R object files
#' }
#' @examples
#' \dontrun{
#' df <- myFun_fov_load("data/fov1.csv")
#' df <- myFun_fov_load("data/fov2.tsv.gz")
#' df <- myFun_fov_load("data/fov3.RData")
#' }
#'
#' @importFrom utils read.csv
myFun_fov_load <- function(path_to_fov) {
  is_gz <- grepl("\\.gz$", path_to_fov, ignore.case = TRUE)
  base_path <- sub("\\.gz$", "", path_to_fov, ignore.case = TRUE)
  
  if (grepl("\\.RData$", base_path, ignore.case = TRUE)) {
    each_transDF <- get(load(path_to_fov))  # assume not gzipped
  } else if (grepl("\\.rds$", base_path, ignore.case = TRUE)) {
    each_transDF <- readRDS(path_to_fov)    # assume not gzipped
  } else if (grepl("\\.csv$", base_path, ignore.case = TRUE)) {
    con <- if (is_gz) gzfile(path_to_fov) else path_to_fov
    each_transDF <- read.csv(con, sep = ',', header = TRUE)
  } else if (grepl("\\.txt$|\\.tsv$", base_path, ignore.case = TRUE)) {
    con <- if (is_gz) gzfile(path_to_fov) else path_to_fov
    each_transDF <- read.csv(con, sep = '\t', header = TRUE)
  } else {
    stop(sprintf('Unsupported file type. Must be RData, rds, csv, txt, or tsv (optionally gzipped for tabular formats). Got: %s', path_to_fov))
  }
  
  return(each_transDF)
}


# function to load each FOV's transcript data.frame
#' @title prepare_perFOV_transDF
#' @description supporting function for \code{runPreprocess}, \code{fastReseg_internalRef} and \code{fastReseg_flag_all_errors} to get unique IDs for cells and transcripts, and convert pixel coordinates to um; when `drop_original = FALSE, the function will also return original per FOV based cell ID and coordinates under columns `CellId`, `pixel_x`, `pixel_y`, `idx_z`.
#' @param each_transDF data.frame for raw transcript
#' @param fov_centerLocs a named vector of fov 2D coordinates
#' @param prefix_vals a named vector of values to be used as prefix in `UMI_transID` and `UMI_cellID`; when `prefix_vals` != NULL, unique transcript_id would be generated from `prefix_vals` and `transID_coln` in `each_transDF` 
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of `each_transDF`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of `each_transDF`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in `each_transDF`; when `prefix_vals` != NULL, unique transcript_id would be generated from `prefix_vals` and `transID_coln` in `each_transDF`
#' @param transGene_coln the column name of target or gene name in `each_transDF`
#' @param cellID_coln the column name of cell_ID in `each_transDF`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_vals` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `each_transDF` 
#' @param invert_y flag to invert y axis of local coordinates during stitching (default = TRUE)
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param drop_original flag to drop original per FOV based cell ID and coordinates under columns `CellId`, `pixel_x`, `pixel_y`, `idx_z` (default = FALSE)
#' @return a list contains transcript_df for downstream process and extracellular transcript data.frame
#' ' \describe{
#'    \item{intraC}{a data.frame for intracellular transcript, `UMI_transID` and `UMI_cellID` as column names for unique transcript_id and cell_id, `target` as column name for target gene name}
#'    \item{extraC}{a data.frame for extracellular transcript, same structure as the `intraC` data.frame in returned list}
#' }
#' @export
prepare_perFOV_transDF <- function(each_transDF, 
                                   fov_centerLocs, 
                                   prefix_vals = NULL, 
                                   pixel_size = 0.18, 
                                   zstep_size = 0.8,
                                   transID_coln = NULL,
                                   transGene_coln = "target",
                                   cellID_coln = 'CellId', 
                                   spatLocs_colns = c('x','y','z'), 
                                   invert_y = TRUE,
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
  if(invert_y==TRUE){
    raw_locs[[orig_spatLocs_colns[2]]] <- 0-raw_locs[[orig_spatLocs_colns[2]]]
  }
  
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
      extraC_idx <- which(each_transDF[['CellId']] %in% extracellular_cellID)
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



#' @title checkTransFileInputsAndLoadFirst
#' @description check input formats for transcript data.frame file list and load 1st fov
#' @param transcript_df the data.frame of transcript level information with unique CellId, set to NULL if read from the `transDF_fileInfo`
#' @param transDF_fileInfo a data.frame with each row for each individual file of per FOV transcript data.frame within which the coordinates and CellId are unique, columns include the file path of per FOV transcript data.frame file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set; when NULL, use the provided `transcript_df` directly
#' @param filepath_coln the column name of each individual file of per FOV transcript data.frame in `transDF_fileInfo`
#' @param prefix_colns the column names of annotation in `transDF_fileInfo`, to be added to the CellId as prefix when creating unique cell_ID for entire data set; set to NULL if use the original `transID_coln` or `cellID_coln` 
#' @param fovOffset_colns the column name of coordinate offsets in 1st and 2nd dimension for each per FOV transcript data.frame in `transDF_fileInfo`, unit in micron
#' Notice that some assays like SMI has XY axes swapped between stage and each FOV such that `fovOffset_colns` should be c("stage_Y", "stage_X").
#' @param pixel_size the micrometer size of image pixel listed in 1st and 2nd dimension of `spatLocs_colns` of each `transcript_df`
#' @param zstep_size the micrometer size of z-step for the optional 3rd dimension of `spatLocs_colns` of each `transcript_df`
#' @param transID_coln the column name of transcript_ID in `transcript_df`, default = NULL to use row index of transcript in each `transcript_df`; when `prefix_colns` != NULL, unique transcript_id would be generated from `prefix_colns` and `transID_coln` in each `transcript_df`
#' @param transGene_coln the column name of target or gene name in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`; when `prefix_colns` != NULL, unique cell_ID would be generated from `prefix_colns` and `cellID_coln` in each `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param invert_y flag to invert y axis of local coordinates during stitching (default = TRUE)
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @return a list contains transcript_df for downstream process and extracellular transcript data.frame
#' ' \describe{
#'    \item{intraC}{a data.frame for intracellular transcript, `UMI_transID` and `UMI_cellID` as column names for unique transcript_id and cell_id, `target` as column name for target gene name}
#'    \item{extraC}{a data.frame for extracellular transcript, same structure as the `intraC` data.frame in returned list}
#' }
checkTransFileInputsAndLoadFirst <- function(transcript_df = NULL, 
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
                                             invert_y = TRUE,
                                             extracellular_cellID = NULL){
  # spatial dimension
  d2_or_d3 <- length(spatLocs_colns)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("`spatLocs_colns` must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in `transcript_df`.")
  } else {
    message(sprintf("%d Dimension of spaital coordinates are provided.", d2_or_d3))
  }
  
  
  ## check the format of transcript data.frame provided ----
  # a data.frame by itself or a data.frame with file path to each per FOV transcript data
  if(is.null(transDF_fileInfo) & is.null(transcript_df)){
    stop("Must provdie either `transcript_df` or `transDF_fileInfo` for distanct cutoff calculation.")
  } 
  
  # a data.frame with file path to each per FOV transcript data
  if (!is.null(transDF_fileInfo)){
    if(! "data.frame" %in% class(transDF_fileInfo)){
      stop("The provided `transDF_fileInfo` is not a data.frame.")
    }
    
    message(sprintf("%d individual per FOV files are provided in `transDF_fileInfo`, use the 1st file to calculate distance cutoffs", nrow(transDF_fileInfo)))
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
  
  
  
  # load and prepare first FOV transcript data
  path_to_transDF <- transDF_fileInfo[[filepath_coln]][1]
  if(!is.na(path_to_transDF)){
    transcript_df <- myFun_fov_load(path_to_fov = path_to_transDF)
  }
  
  # fovOffset_colns must have XY axes of stage matched to XY axes of images
  # return a list with two data.frame, `intraC` and `extraC` for intracelllular and extracellular transcripts, respectively
  transcript_df <- prepare_perFOV_transDF(each_transDF = transcript_df, 
                                          fov_centerLocs = unlist(transDF_fileInfo[1, fovOffset_colns]),
                                          prefix_vals = unlist(transDF_fileInfo[1, prefix_colns]), 
                                          pixel_size = pixel_size, 
                                          zstep_size = zstep_size,
                                          transID_coln = transID_coln,
                                          transGene_coln = transGene_coln,
                                          cellID_coln = cellID_coln, 
                                          spatLocs_colns = spatLocs_colns, 
                                          invert_y = invert_y,
                                          extracellular_cellID = extracellular_cellID, 
                                          drop_original = TRUE)
  
  return(transcript_df)
}


#' Get number of cores for parallelized operations
#'
#' @param percentCores percent of cores to use for parallelization, value range 0 to 1
#' @param minNotUsedCores minimum number of cores to leave for background processes
#' 
#' @return number of cores to use for mclapply
numCores <- function(percentCores = 0.9, minNotUsedCores = 2) {
  if(percentCores > 1 & percentCores <= 0){
    stop("percentCores is not a valid number, must be between 0-1")
  }
  
  num_cores <- 1
  if (.Platform$OS.type == "unix") {
    if (is.null(getOption("mc.cores"))) {
      num_cores <- parallel::detectCores()
      if(num_cores <= minNotUsedCores){
        stop("minNotUsedCores must be fewer than available cores")
      }
      num_cores <- min(floor(num_cores*percentCores), num_cores-minNotUsedCores)
    } else {
      num_cores <- getOption("mc.cores") 
    }
  }
  
  # limit the core number during R CMD check
  limitCores <- as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false"))
  ## Possible values: 'TRUE' 'false', 'warn', 'error'
  if(!is.na(limitCores)){
    if(limitCores){
      num_cores <- min(2, num_cores)
    }
  }
  
  return(num_cores)
}

