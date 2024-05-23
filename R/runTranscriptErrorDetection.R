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
#' @param spatLocs_colns the column names of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param model_cutoff the cutoff of transcript number to do spatial modeling (default = 50)
#' @param score_cutoff the cutoff of score to separate between high and low score transcripts (default = -2)
#' @param svm_args a list of arguments to pass to svm function, typically involve kernel, gamma, scale
#' @param groupTranscripts_method use either "dbscan" or "delaunay" method to group transcripts in space (default = "dbscan")
#' @param distance_cutoff maximum molecule-to-molecule distance within same transcript group (default = "auto")
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config
#' @param seed_transError seed for transcript error detection step, default = NULL to skip the seed   
#' @return data frame for transcripts in `chosen_cells` only, containing information for transcript score classifications and spatial group assignments as well as new cell/group ID for downstream resegmentation.
#' @importFrom data.table ':='
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
  classDF_ToFlagTrans <- as.data.frame(transcript_df)[which(transcript_df[[cellID_coln]] %in% chosen_cells),
                                                      setdiff(colnames(transcript_df), c('DecVal','SVM_class','SVM_cell_type'))]
  
  if(nrow(classDF_ToFlagTrans)<1) {
    stop("Error: No transcripts within `chosen_cells`.")
  }
  
  if(!is.null(seed_transError)){
    set.seed(seed_transError)
  }
  
  # `flag_bad_transcripts` function returns a data.frame with transcript in row, original cell_ID and SVM outcomes in column.
  flagged_transDF_SVM <- flag_bad_transcripts(chosen_cells = chosen_cells,
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
  if(is.null(flagged_transDF_SVM)){
    return(NULL)
  }
  
  # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
  flagged_transDF_SVM <- merge(classDF_ToFlagTrans, 
                               as.data.frame(flagged_transDF_SVM)[, c(transID_coln,'DecVal','SVM_class','SVM_cell_type')], 
                               by = transID_coln)
  
  message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                  length(chosen_cells) - length(unique(classDF_ToFlagTrans[[cellID_coln]])), score_cutoff))
  
  # flagged transcript ID, character vector
  flaggedSVM_transID <- flagged_transDF_SVM[flagged_transDF_SVM[['SVM_class']] ==0, transID_coln]
  
  if(length(flaggedSVM_transID)<1){
    return(NULL)
  }
  
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
  # `getCellType_maxScore` function returns a named vector with cell type in values and cell_ID in names
  celltype_cellVector <- getCellType_maxScore(score_GeneMatrix = score_GeneMatrix, 
                                              transcript_df = flagged_transDF_SVM, 
                                              transGene_coln = transGene_coln,
                                              cellID_coln = "tmp_cellID")
  flagged_transDF_SVM[['group_maxCellType']] <- celltype_cellVector[flagged_transDF_SVM[['tmp_cellID']]]
  flagged_transDF_SVM <- as.data.frame(flagged_transDF_SVM)
  rm(celltype_cellVector)
  
  return(flagged_transDF_SVM)
}

