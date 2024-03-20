#' @title flag_bad_transcripts
#' @description find out the spatially connected transcripts among chosen_transcripts based on SVM spatial model which scores each cell for how much their transcripts change their goodness-of-fit over space. 
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
#' @return a data.frame 
#' \enumerate{
#'    \item{cellID_coln, original cell id}
#'    \item{SVM_class, 0 for below cutoff, 1 for above cutoff}
#'    \item{transID_coln, original transcript_id}
#'    \item{transGene_coln, target gene of transcript}
#'    \item{score_coln, score in transcript_df}
#'    \item{spatLocs_colns for spatial coorindates of transcript}
#'    \item{DecVal, decision values of svm model output}
#'    \item{SVM_cell_type, new cell type for each transcript groups within each cells}
#' }
#' @details For score of transcripts within each cell, assign 0 or 1 label to each transcript based on whether the score is above score_cutoff; then run support vector machine on svm(above_cutoff ~ x + y) for 2D, svm(above_cutoff ~ x + y + z) for 3D, default to do radial kernal with scale = TRUE and gamma = 0.1.  
#' @export
flag_bad_transcripts <- function(chosen_cells, 
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
                                                 scale = TRUE, 
                                                 gamma = 0.1)){
  
  # check if model_cutoff making sense
  if(model_cutoff < 2){
    stop(sprintf("model_cutoff must be no less than 2 to enable spatial modeling, current model_cutoff = %s", as.character(model_cutoff)))
  }
  
  d2_or_d3 <- length(spatLocs_colns)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Run SVM in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns)]
  
  # get common cells
  common_cells <- intersect(unique(transcript_df[[cellID_coln]]), unique(chosen_cells))
  
  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            unique(transcript_df[[transGene_coln]]))
  message(sprintf("Found %d common cells and %d common genes among chosen_cells, transcript_df, and score_GeneMatrix. ", 
                  length(common_cells), length(common_genes)))
  
  if(any(length(common_cells) <1, length(common_genes)<1)){
    stop("Too few common cells or genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells & transcript_df[[transGene_coln]] %in% common_genes), ]
  
  
  ## operate in vector form to speed up the process
  # (1) filter cells based on transcript numbers
  count_df <- transcript_df[, .N,  by = .(get(cellID_coln))]
  chosen_cells2 <- count_df[N >= model_cutoff, get]
  
  warning(sprintf("Below model_cutoff = %s, skip %d cells with fewer transcripts. Move forward with remaining %d cells.", 
                  as.character(model_cutoff), length(chosen_cells) - length(chosen_cells2), length(chosen_cells2)))
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells2), ]
  
  if(nrow(transcript_df)<1){
    message("No single transcript left for the evaluation.")
    return(NULL)
  }
  
  # order by cell_ID before moving forward to keep same index
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  
  
  # (2) get new coordinate columns, class and cell_ID, scores, same order as transcript_df
  coord_df <- transcript_df[, .SD, .SDcols = c(transID_coln, cellID_coln, score_coln, spatLocs_colns)]
  
  
  # spatial SVM model at 2D and 3D for classification based on score cutoff
  if(d2_or_d3 ==2){
    colnames(coord_df) <- c('transcript_id', 'cell_ID','score','x','y')
    coord_colns <- c('x', 'y')
  } else {
    colnames(coord_df) <- c('transcript_id', 'cell_ID','score','x','y','z')
    coord_colns <- c('x', 'y', 'z')
  }
  
  # classify based on decision boundaries
  coord_df[['above_cutoff']] <- ifelse(coord_df[['score']]<score_cutoff, 0, 1)
  
  # skip the cell if no element in either group
  # count 0 or 1 for each cell
  count_df0 <- coord_df[above_cutoff ==0 , .N, by = .(cell_ID)]
  count_df1 <- coord_df[above_cutoff ==1 , .N, by = .(cell_ID)]
  count_df <- merge(count_df0, count_df1, by = 'cell_ID', all = TRUE)
  chosen_cells3 <- count_df[['cell_ID']][!as.logical(rowSums(is.na(count_df)))]
  
  warning(sprintf("Skip %d cells with all transcripts in same class given `score_cutoff = %s`. Move forward with remaining %d cells.", 
                  length(chosen_cells2) - length(chosen_cells3), as.character(score_cutoff), length(chosen_cells3)))
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells3), ]
  
  if(nrow(transcript_df)<1){
    message("No single transcript left for the evaluation.")
    return(NULL)
  }
  
  coord_df <- coord_df[which(coord_df[['cell_ID']] %in% chosen_cells3)]
  
  # (3) perform svm by group, not working for near-zero variance predictors
  my_fun <- function(data){
    current_args <- svm_args
    current_args[['x']] <- as.matrix(data[, .SD, .SDcols = c('transcript_id', coord_colns)], 
                                     rownames = 'transcript_id')
    current_args[['y']] <- as.factor(data[['above_cutoff']])
    mod_svm <- do.call(e1071::svm,current_args)
    
    # Make predictions
    outputs <- data.frame(transcript_id = names(mod_svm$fitted), 
                          SVM_class = mod_svm$fitted,
                          DecVal = mod_svm$decision.values[, 1])
    
    return(outputs)
  }
  
  model_stats <- by(coord_df, coord_df$cell_ID, my_fun)
  model_stats <- do.call(rbind, model_stats)
  
  # merge SVM_class and DecVal value to original transcript_df, by cell ID first, then transcript ID
  transcript_df <- merge(transcript_df, model_stats, by.x = transID_coln, by.y = 'transcript_id', all.x = TRUE)
  
  # order by cell_ID before moving forward to keep same index
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  
  # (4) assign cell type based on score matrix
  # get new cell type of each group based on maximum score
  
  # get score for each transcripts
  transcriptGeneScore <- score_GeneMatrix[transcript_df[[transGene_coln]], ]
  rownames(transcriptGeneScore) <- transcript_df[[transID_coln]]
  
  tmp_score <- as.data.frame(transcriptGeneScore)
  tmp_score[[cellID_coln]] <- transcript_df[[cellID_coln]]
  tmp_score[['SVM_class']] <- transcript_df[['SVM_class']]
  
  tmp_score <-data.table::setDT(tmp_score)[, lapply(.SD, sum), by = c(cellID_coln, 'SVM_class')] 
  cellType_DT <- tmp_score[, .SD, .SDcols = c(cellID_coln, 'SVM_class')]
  tmp_score[[cellID_coln]] <- NULL
  tmp_score[['SVM_class']] <- NULL
  # assign cell type based on max values
  max_idx_1st <- max.col(tmp_score,ties.method="first")
  cellType_DT[['SVM_cell_type']] <- colnames(tmp_score)[max_idx_1st]
  
  # merge cell type to original transcript_df
  transcript_df <- as.data.frame(merge(transcript_df, cellType_DT, by = c(cellID_coln, 'SVM_class')))
  
  
  return(transcript_df)
}



#' @title flagTranscripts_LDA_hyperplane
#' @description find out the spatially connected transcripts among chosen_transcripts based on LDA hyperplane spatial model which scores each cell for how much their transcripts change their goodness-of-fit over space. 
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
#' @return a data.frame 
#' \enumerate{
#'    \item{cellID_coln, original cell id}
#'    \item{LDA_class, 0 for below cutoff, 1 for above cutoff}
#'    \item{transID_coln, original transcript_id}
#'    \item{transGene_coln, target gene of transcript}
#'    \item{score_coln, score in transcript_df}
#'    \item{spatLocs_colns for spatial coorindates of transcript}
#'    \item{LD1, LD1 value of LDA model output}
#'    \item{LDA_cell_type, new cell type for each transcript groups within each cells}
#' }
#' @details For score of transcripts within each cell, assign 0 or 1 label to each transcript based on whether the score is above score_cutoff; then run linear discriminant analysis on lda(above_cutoff ~ x + y + xy + x^2 + y^2) for 2D, lda(above_cutoff ~ x + y + z+ xy +xz + yz + x^2 + y^2 + z^2) for 3D. Coordinate variables with variance less than 1e-8 would not be used. Z-step of 5E-4mm or similar often resulted in too small of variables for z2. 
#' @export
flagTranscripts_LDA_hyperplane <- function(chosen_cells, 
                                           score_GeneMatrix, 
                                           transcript_df, 
                                           cellID_coln = "CellId", 
                                           transID_coln = "transcript_id",
                                           transGene_coln = "target",
                                           score_coln = "score", 
                                           spatLocs_colns = c("x","y","z"),
                                           model_cutoff = 50, 
                                           score_cutoff = -2){
  
  # check if model_cutoff making sense
  if(model_cutoff < 2){
    stop(sprintf("model_cutoff must be no less than 2 to enable spatial modeling, current model_cutoff = %s", as.character(model_cutoff)))
  }
  
  d2_or_d3 <- length(spatLocs_colns)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Run LDA in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transGene_coln, score_coln, spatLocs_colns)]
  
  # get common cells
  common_cells <- intersect(unique(transcript_df[[cellID_coln]]), unique(chosen_cells))
  
  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            unique(transcript_df[[transGene_coln]]))
  message(sprintf("Found %d common cells and %d common genes among chosen_cells, transcript_df, and score_GeneMatrix. ", 
                  length(common_cells), length(common_genes)))
  
  if(any(length(common_cells) <1, length(common_genes)<1)){
    stop("Too few common cells or genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells & transcript_df[[transGene_coln]] %in% common_genes), ]
  
  
  
  ## operate in vector form to speed up the process
  # (1) filter cells based on transcript numbers
  count_df <- transcript_df[, .N,  by = .(get(cellID_coln))]
  chosen_cells2 <- count_df[N >= model_cutoff, get]
  
  warning(sprintf("Below model_cutoff = %s, skip %d cells with fewer transcripts. Move forward with remaining %d cells.", 
                  as.character(model_cutoff), length(chosen_cells) - length(chosen_cells2), length(chosen_cells2)))
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells2), ]
  if(nrow(transcript_df)<1){
    message("No single transcript left for the evaluation.")
    return(NULL)
  }
  # order by cell_ID before moving forward to keep same index
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  
  
  # (2) get new coordinate columns, class and cell_ID, scores, same order as transcript_df
  coord_df <- transcript_df[, .SD, .SDcols = c(transID_coln, cellID_coln, score_coln, spatLocs_colns)]
  
  
  # spatial LDA model at 2D and 3D for classification based on score cutoff
  if(d2_or_d3 ==2){
    colnames(coord_df) <- c('transcript_id', 'cell_ID','score','x','y')
    coord_colns <- c('x', 'y', 'x2','y2','xy')
  } else {
    colnames(coord_df) <- c('transcript_id', 'cell_ID','score','x','y','z')
    coord_colns <- c('x', 'y','z', 'x2','y2','z2','xy','xz','yz')
  }
  
  # classify based on decision boundaries
  coord_df[['above_cutoff']] <- ifelse(coord_df[['score']]<score_cutoff, 0, 1)
  
  # skip the cell if no element in either group
  # count 0 or 1 for each cell
  count_df0 <- coord_df[above_cutoff ==0 , .N, by = .(cell_ID)]
  count_df1 <- coord_df[above_cutoff ==1 , .N, by = .(cell_ID)]
  count_df <- merge(count_df0, count_df1, by = 'cell_ID', all = TRUE)
  chosen_cells3 <- count_df[['cell_ID']][!as.logical(rowSums(is.na(count_df)))]
  
  warning(sprintf("Skip %d cells with all transcripts in same class given `score_cutoff = %s`. Move forward with remaining %d cells.", 
                  length(chosen_cells2) - length(chosen_cells3), as.character(score_cutoff), length(chosen_cells3)))
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells3), ]
  if(nrow(transcript_df)<1){
    message("No single transcript left for the evaluation.")
    return(NULL)
  }
  coord_df <- coord_df[which(coord_df[['cell_ID']] %in% chosen_cells3)]
  
  # new spatial coordinate colns
  coord_df[['x2']] <- coord_df[['x']]^2
  coord_df[['y2']] <- coord_df[['y']]^2
  coord_df[['xy']] <- coord_df[['x']]*coord_df[['y']]
  
  if(d2_or_d3 ==3){
    coord_df[['z2']] <- coord_df[['z']]^2
    coord_df[['xz']] <- coord_df[['x']]*coord_df[['z']]
    coord_df[['yz']] <- coord_df[['y']]*coord_df[['z']]
  }
  
  
  # (3) perform lda by group, not working for near-zero variance predictors
  my_fun <- function(data){
    # remove colns with near zero variables before LDA, slow for large data set
    bad_colns <- names(which(data[, apply(data, MARGIN = 2, function(x) var(x, na.rm=TRUE))]<1e-8))
    mod_formula <- paste0('as.factor(above_cutoff) ~ ', paste0(setdiff(coord_colns, bad_colns), collapse = " + "))
    
    # linear discriminant analysis model, must have both groups in data
    # CV = TRUE to do leave-one-out cross-validation, default = FALSE to allow predict()
    mod_lda <- MASS::lda(as.formula(mod_formula), data = data)
    
    # Make predictions
    tmp_predictions <- predict(mod_lda, data = data)
    
    outputs <- data.frame(transcript_id = data[['transcript_id']], 
                          LDA_class = tmp_predictions$class,
                          LD1 = tmp_predictions$x)
    
    return(outputs)
  }
  
  model_stats <- by(coord_df, coord_df$cell_ID, my_fun)
  model_stats <- do.call(rbind, model_stats)
  
  # merge LDA_class and LD1 value to original transcript_df, by cell ID first, then transcript ID
  transcript_df <- merge(transcript_df, model_stats, by.x = transID_coln, by.y = 'transcript_id', all.x = TRUE)
  
  # order by cell_ID before moving forward to keep same index
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  
  # (4) assign cell type based on score matrix
  # get new cell type of each group based on maximum score
  
  # get score for each transcripts
  transcriptGeneScore <- score_GeneMatrix[transcript_df[[transGene_coln]], ]
  rownames(transcriptGeneScore) <- transcript_df[[transID_coln]]
  
  tmp_score <- as.data.frame(transcriptGeneScore)
  tmp_score[[cellID_coln]] <- transcript_df[[cellID_coln]]
  tmp_score[['LDA_class']] <- transcript_df[['LDA_class']]
  
  tmp_score <-data.table::setDT(tmp_score)[, lapply(.SD, sum), by = c(cellID_coln, 'LDA_class')] 
  cellType_DT <- tmp_score[, .SD, .SDcols = c(cellID_coln, 'LDA_class')]
  tmp_score[[cellID_coln]] <- NULL
  tmp_score[['LDA_class']] <- NULL
  # assign cell type based on max values
  max_idx_1st <- max.col(tmp_score,ties.method="first")
  cellType_DT[['LDA_cell_type']] <- colnames(tmp_score)[max_idx_1st]
  
  # merge cell type to original transcript_df
  transcript_df <- as.data.frame(merge(transcript_df, cellType_DT, by = c(cellID_coln, 'LDA_class')))
  
  return(transcript_df)
}

