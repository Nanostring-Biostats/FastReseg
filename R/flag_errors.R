#' @title score_cell_segmentation_error
#' @description Score each cell for how much their transcripts change their goodness-of-fit over space. 
#' @param chosen_cells the cell_ID of chosen cells
#' @param transcript_df the data.frame of transcript_ID, cell_ID, score, spatial coordinates
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param score_coln the column name of score in transcript_df
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in transcript_df 
#' @param model_cutoff the cutoff of transcript number to do spatial modeling (default = 50)
#' @return data.frame with columns for 
#' \enumerate{
#'    \item{cell_ID, cell id}
#'    \item{transcript_num, number of transcripts in given cell}
#'    \item{modAlt_rsq, summary(mod_alternative)$r.squared}
#'    \item{lrtest_ChiSq, lrtest chi-squared value}
#'    \item{lrtest_Pr, lrtest probability larger than chi-squared value, p-value}
#' }
#' @details For tLLRv2 score of transcripts within each cell,  run a quadratic model: mod_alternative = lm(tLLRv2 ~ x + y + x2 + y2 +xy) for 2D,  lm(tLLRv2 ~ x + y + z + x2 + y2 +z2 +xy + xz + yz) for 3D and a null model: mod_null = lm(tLLRv2 ~ 1); then run lmtest::lrtest(mod_alternative, mod_null). Return statistics for mod_alternative$fitted.values (standard deviation and minimal value), summary(mod_alternative)$r.squared and as well as lrtest chi-squared value.  
#' @export
score_cell_segmentation_error <- function(chosen_cells, transcript_df, 
                                          cellID_coln = "CellId", 
                                          transID_coln = "transcript_id",
                                          score_coln = "score", 
                                          spatLocs_colns = c("x","y","z"),
                                          model_cutoff = 50){
  
  # check if model_cutoff making sense
  if(model_cutoff < 2){
    stop(sprintf("model_cutoff must be no less than 2 to enable spatial modeling, current model_cutoff = %s", as.character(model_cutoff)))
  }
  
  
  d2_or_d3 <- length(spatLocs_colns)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Run linear regreassion in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, score_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, score_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, score_coln, spatLocs_colns)]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  ## operate in vector form to speed up the process
  # (1) filter cells based on transcript numbers
  count_df <- transcript_df[, .N,  by = .(get(cellID_coln))]
  chosen_cells2 <- count_df[N >= model_cutoff, get]
  
  warning(sprintf("Below model_cutoff = %s, skip %d cells with fewer transcripts. Move forward with remaining %d cells.", 
                  as.character(model_cutoff), length(chosen_cells) - length(chosen_cells2), length(chosen_cells2)))
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells2), ]
  
  # (2) get new coordinate columns and cell_ID, scores, same order as transcript_df
  coord_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, score_coln, spatLocs_colns)]
  
  # spatial quadratic model
  # lm(tLLRv2 ~ x + y + x2 + y2 + xy) for 2D,  lm(tLLRv2 ~ x + y + z + x2 + y2 +z2 +xy + xz + yz) for 3D
  if(d2_or_d3 ==2){
    colnames(coord_df) <- c('cell_ID','score','x','y')
    mod_formula <- 'score ~ x + y + x2 + y2 + xy'
  } else {
    colnames(coord_df) <- c('cell_ID','score','x','y','z')
    mod_formula <- 'score ~ x + y + z + x2 + y2 +z2 +xy + xz + yz'
  }
  
  coord_df[['x2']] <- coord_df[['x']]^2
  coord_df[['y2']] <- coord_df[['y']]^2
  coord_df[['xy']] <- coord_df[['x']]*coord_df[['y']]
  
  if(d2_or_d3 ==3){
    coord_df[['z2']] <- coord_df[['z']]^2
    coord_df[['xz']] <- coord_df[['x']]*coord_df[['z']]
    coord_df[['yz']] <- coord_df[['y']]*coord_df[['z']]
  }
  
  
  # (3) perform lm and lrtest by group
  my_fun <- function(data){
    # null linear model, lm(tLLRv2 ~ 1)
    mod_null <- lm(score~1, data = data)
    
    # spatial quadratic model
    # lm(tLLRv2 ~ x + y + x2 + y2 + xy) for 2D,  lm(tLLRv2 ~ x + y + z + x2 + y2 +z2 +xy + xz + yz) for 3D
    mod_alternative <- lm(as.formula(mod_formula), data = data)
    
    #likelihood ratio test of nested model
    lrtest_res <- lmtest::lrtest(mod_alternative, mod_null)
    
    outputs <- data.frame(transcript_num = nrow(data), 
                          modAlt_rsq = summary(mod_alternative)$r.squared,
                          lrtest_ChiSq = lrtest_res$Chisq[2],
                          lrtest_Pr = lrtest_res$`Pr(>Chisq)`[2])
    return(outputs)
  }
  
  model_stats <- by(coord_df, coord_df$cell_ID, my_fun)
  model_stats <- do.call(rbind, model_stats)
  model_stats[['cell_ID']] <- rownames(model_stats)
  
  return(model_stats)
  
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
#' @param spatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
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
  # model_stats[['cell_ID']] <- sapply(strsplit(rownames(model_stats), split = "[.]"),"[[", 1)
  
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
  transcript_df <- merge(transcript_df, cellType_DT, by = c(cellID_coln, 'LDA_class'))
  
  return(transcript_df)
}


#' @title flagTranscripts_SVM
#' @description find out the spatially connected transcripts among chosen_transcripts based on SVM spatial model which scores each cell for how much their transcripts change their goodness-of-fit over space. 
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
flagTranscripts_SVM <- function(chosen_cells, 
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
  # model_stats[['cell_ID']] <- sapply(strsplit(rownames(model_stats), split = "[.]"),"[[", 1)
  
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
  transcript_df <- merge(transcript_df, cellType_DT, by = c(cellID_coln, 'SVM_class'))
  
  return(transcript_df)
}

#' @title groupTranscripts_Delanuay
#' @description group the flagged transcript within each cell based on spatial connectivity of their transcript delaunay network
#' @param chosen_transcripts the transcript_id of chosen transcript
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level
#' @param distance_cutoff maximum distance within connected transcript group (default = "auto")
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @return data frame of connected transcripts among chosen_transcripts
#' #' \enumerate{
#'    \item{cellID_coln, orignal cell_ID}
#'    \item{transID_coln, connected transcripts among chosen_transcripts}
#'    \item{transSpatLocs_coln, spatial coordiantes of transcript}
#'    \item{transcript_group, group of chosen_transcripts}
#' } 
#' @details for query cell, build network on flagged transcripts only to identify groups. In case of no more than 3 transcripts, determine the grouping based on distance cutoff directly; when distance cutoff = 'auto', no additional edge filtering based on delaunay network output but use 20% average XY cell range as cutoff when no more than 3 transcript.  
#' @export
groupTranscripts_Delanuay <- function(chosen_transcripts = NULL, 
                                      config_spatNW_transcript, 
                                      distance_cutoff = "auto",
                                      transcript_df, 
                                      cellID_coln = "CellId", transID_coln = "transcript_id",
                                      transSpatLocs_coln = c('x','y','z'), 
                                      visualize_network = FALSE){
  
  
  if(distance_cutoff =='auto'){
    message("automatic cutoff for delauny network.")
  } else if (is.numeric(distance_cutoff)){
    if (distance_cutoff <= 0){
      stop("distance_cutoff must be either `auto` or positive number")
    }
  }else {
    stop("distance_cutoff must be either `auto` or positive number")
  }
  
  d2_or_d3 <- length(transSpatLocs_coln)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("transSpatLocs_coln must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Do spatial network analysis in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transSpatLocs_coln)]
  
  if(is.null(chosen_transcripts)){
    stop("chosen_transcripts = NULL. No resegmentation is done.")
  } else {
    chosen_transcripts <- intersect(chosen_transcripts, transcript_df[[transID_coln]])
    if(length(chosen_transcripts) <1){
      stop("Common cells contain no provided chosen_transcripts. No resegmentation is done.")
    } else {
      message(sprintf("%d chosen_transcripts are found in common cells.", length(chosen_transcripts)))
    }
  }
  
  # check config_spatNW_transcript
  config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript,
                                                     spat_locs = transcript_df)
  # which cells contained the chosen_transcripts
  # only chosen transcript
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  chosenTrans_df <- transcript_df[which(transcript_df[[transID_coln]] %in% chosen_transcripts), ]
  chosen_cells <- unique(chosenTrans_df[[cellID_coln]])
  # all transcripts of chosen_cells
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  ## operate in vector form to speed up the process
  # (1) get number of flagged transcripts and distance cutoff 
  count_df <- chosenTrans_df[, .N,  by = .(get(cellID_coln))]
  
  # get distance cutoff based on cell size if "auto"
  if(distance_cutoff == 'auto'){
    # size range of cells
    cellsize_df <- transcript_df[, lapply(.SD, function(x) (range(x)[2] - range(x)[1])),  by = .(get(cellID_coln)), .SDcols = transSpatLocs_coln]
    # cutoff for absolute distance between orphan transcript, defined as 20% of xy cell size range as a cube
    orphan_cutoff <- 0.2*rowMeans(cellsize_df[, 2:3])
    names(orphan_cutoff) <- cellsize_df[[1]]
    
  }else{
    orphan_cutoff <- rep(distance_cutoff, nrow(count_df))
    names(orphan_cutoff) <- count_df[[1]]
    
  }
  
  
  # (2) assign group for 1 flagged transcript cases
  transcript_df[['transcript_group']] <- as.numeric(transcript_df[[transID_coln]] %in% chosen_transcripts)
  
  # single flagged transcripts, as it is
  
  # (3) assign group for 2 flagged transcript cases
  # two flagged transcripts; evaluate distance on flagged transcript only
  chosenTrans_df2 <-  chosenTrans_df[which(chosenTrans_df[[cellID_coln]] %in% count_df[N ==2, get]), ]
  if(nrow(chosenTrans_df2) >0){
    data.table::setkeyv(chosenTrans_df2, c(cellID_coln, transID_coln))
    # distance between 2 points
    dist_df <- cbind(chosenTrans_df2[seq(1, nrow(chosenTrans_df2), by = 2), .SD, .SDcols = transSpatLocs_coln], 
                     chosenTrans_df2[seq(2, nrow(chosenTrans_df2), by = 2), .SD, .SDcols = transSpatLocs_coln])
    colnames(dist_df) <- c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))
    dist_df[['from']] <- chosenTrans_df2[[transID_coln]][seq(1, nrow(chosenTrans_df2), by = 2)]
    dist_df[['to']] <- chosenTrans_df2[[transID_coln]][seq(2, nrow(chosenTrans_df2), by = 2)]
    dist_df[['cell_ID']] <- chosenTrans_df2[[cellID_coln]][seq(1, nrow(chosenTrans_df2), by = 2)]
    dist_df <- data.table::setDT(dist_df)
    
    dist_df[, `:=`(distance, stats::dist(x = matrix(.SD, nrow = 2, byrow = T))), 
            by = 1:nrow(dist_df), 
            .SDcols = c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))]
    
    dist_df[, `:=`(distance, as.numeric(distance))]
    #dist_df[, `:=`(weight, 1/distance)]
    dist_df[, `:=` (separate, distance > orphan_cutoff[cell_ID])]
    tmp_transID <- dist_df[separate == TRUE, to]
    transcript_df[['transcript_group']][which(transcript_df[[transID_coln]] %in% tmp_transID)] <- 2
  }
  
  
  # (4) assign group for 3 flagged transcript cases
  # for each cell of 3 flagged transcripts in case of 3D
  my_fun_3points <- function(df_subset){
    # distance between 3 points, must be 3D
    dist_df <- cbind(df_subset[1:3, .SD, .SDcols = transSpatLocs_coln], df_subset[c(2,3,1), .SD, .SDcols = transSpatLocs_coln])
    colnames(dist_df) <- c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))
    dist_df[['from']] <- df_subset[[transID_coln]]
    dist_df[['to']] <- df_subset[[transID_coln]][c(2,3,1)]
    dist_df[['cell_ID']] <- df_subset[[cellID_coln]]
    dist_df <- data.table::setDT(dist_df)
    
    dist_df[, `:=`(distance, stats::dist(x = matrix(.SD, nrow = 2, byrow = T))), 
            by = 1:nrow(dist_df), 
            .SDcols = c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))]
    dist_df[, `:=`(distance, as.numeric(distance))]
    dist_df[, `:=`(weight, 1/distance)]
    dist_df[, `:=` (separate, distance > orphan_cutoff[cell_ID])]
    dist_df[['group']] <- rep(1, nrow(dist_df))
    if(dist_df[['separate']][1]){
      # separate 1 and 2
      dist_df[['group']][2] <- 2
    }
    if(dist_df[['separate']][2]){
      # separate 2 and 3
      dist_df[['group']][3] <- dist_df[['group']][2]+1
    }
    
    if(!dist_df[['separate']][3]){
      # keep 1 and 3
      dist_df[['group']][3] <- dist_df[['group']][1]
    }
    
    df_subset[['transcript_group']] <- dist_df[['group']]
    return(df_subset)
  }
  
  chosenTrans_df3 <-  chosenTrans_df[which(chosenTrans_df[[cellID_coln]] %in% count_df[N ==3, get]), ]
  if(d2_or_d3 ==3 & nrow(chosenTrans_df3) >0){
    tmp_group <- by(chosenTrans_df3, chosenTrans_df3[[cellID_coln]], my_fun_3points)
    tmp_group <- do.call(rbind, tmp_group)
    group_converter <- tmp_group[['transcript_group']] 
    names(group_converter) <- tmp_group[[transID_coln]]
    tmp_idx <- which(transcript_df[[transID_coln]] %in% tmp_group[[transID_coln]])
    transcript_df[['transcript_group']][tmp_idx] <- group_converter[transcript_df[[transID_coln]][tmp_idx]]
  }
  
  # (5) assign group for multiple flagged transcript cases using delanuay network analysis
  # function for each cell based on delaunay
  my_fun_delanuay <- function(df_subset){
    each_cell <- df_subset[[cellID_coln]][1]
    
    # some cells may have multiple transcripts in exactly same location
    # this would result error in delaunay network generation
    dfCoord_subset <- unique(df_subset[, .SD, .SDcols = transSpatLocs_coln])
    
    # 1 point
    if(nrow(dfCoord_subset) ==1){
      dfCoord_subset[['transcript_group']] = 1
    } else if(nrow(dfCoord_subset) ==2){
      # 2 points, group by direct cutoff
      message(sprintf("%s has %d unique coordinates, grouped by direct cutoff %.4f. ", 
                      each_cell, nrow(dfCoord_subset), orphan_cutoff[each_cell]))
      
      if(stats::dist(dfCoord_subset) > orphan_cutoff[each_cell]){
        dfCoord_subset[['transcript_group']] <- seq_len(nrow(dfCoord_subset))
      } else {
        dfCoord_subset[['transcript_group']] = 1
      }
      
      
    } else if (nrow(dfCoord_subset) ==3 & d2_or_d3 ==3){
      # 3 points in 3D
      message(sprintf("%s has %d unique coordinates, grouped by direct cutoff %.4f. ", 
                      each_cell, nrow(dfCoord_subset), orphan_cutoff[each_cell]))
      
      dfCoord_subset <- myFun_3point_singleCell(dfCoord_subset, 
                                                transSpatLocs_coln = transSpatLocs_coln, 
                                                distance_cutoff = orphan_cutoff[each_cell], 
                                                startGroup = 1)
      
    } else {
      # delaunay network for 3+ points in 2D and 4+ points in 3D
      dfCoord_subset[['tmp_transID']] <- seq_len(nrow(dfCoord_subset))
      
      if(d2_or_d3 ==2){
        # run 2D
        delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                spatLocs_df = dfCoord_subset, 
                                                                ID_column = 'tmp_transID',
                                                                spatLocs_column = transSpatLocs_coln)
      } else {
        # 3D specify, check if all selected transcript in same z plane
        z_num <- unique(dfCoord_subset[[transSpatLocs_coln[3]]])
        if(length(z_num) <2){
          # single z-plane, run as 2D
          message(sprintf("%s with all %d flagged transcripts of unique transcripts in same z plane, run 2D network analysis.", 
                          each_cell, nrow(dfCoord_subset)))
          delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                  spatLocs_df = dfCoord_subset, 
                                                                  ID_column = 'tmp_transID',
                                                                  spatLocs_column = transSpatLocs_coln[1:2])
        } else {
          # multiple z-plane, run as 3D
          delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                  spatLocs_df = dfCoord_subset, 
                                                                  ID_column = 'tmp_transID',
                                                                  spatLocs_column = transSpatLocs_coln)
        }
      }
      
      # if failed to give delaunay network, assign separate group to exact same coordinates
      if(!is.null(delaunayNW_Obj)){
        # make decision based on network
        transNetWorkDT <- delaunayNW_Obj$networkDT
        
        ## convert to igraph object
        all_index = unique(x = c(transNetWorkDT$from, transNetWorkDT$to))
        network_igraph = igraph::graph_from_data_frame(transNetWorkDT[,.(from, to, weight, distance)], directed = TRUE, vertices = all_index)
        
        # remove edges that are longer than cutoff, check the number of unconnected group
        if(distance_cutoff !='auto'){
          # better to use auto, if not, use distance_cutoff >10 pixel = 1.8um, cutoff = 15 pixel = 2.7um works best 
          network_igraph <- igraph::delete.edges(network_igraph, igraph::E(network_igraph)[distance > distance_cutoff])
        } 
        
        # identify groups
        group_vector <- igraph::components(network_igraph, mode = "weak")[['membership']]
        dfCoord_subset[['transcript_group']] <- 0
        dfCoord_subset[['transcript_group']] <- group_vector[dfCoord_subset[['tmp_transID']]]
        
        # some solo transcript would have membership/transcript_group = NA, assign new group based on distance
        solo_transcripts <- dfCoord_subset[is.na(dfCoord_subset[['transcript_group']]), ]
        if(nrow(solo_transcripts) >0){
          
          groupNum <- length(unique(dfCoord_subset[['transcript_group']])) # include na
          message(sprintf("%s has %d solo transcripts, flag based on distance cutoff = %.4f.", 
                          each_cell, nrow(solo_transcripts), orphan_cutoff[each_cell]))
          
          if(nrow(solo_transcripts) ==1){
            # 1 transcript
            solo_transcripts[['transcript_group']]<- groupNum
          } else if(nrow(solo_transcripts) ==2){
            # 2 transcripts
            if(stats::dist(solo_transcripts[, .SD, .SDcols = transSpatLocs_coln]) > orphan_cutoff[each_cell]){
              solo_transcripts[['transcript_group']] <- c(groupNum,groupNum+1)
            }else {
              solo_transcripts[['transcript_group']] <- groupNum
            }
          } else if(nrow(solo_transcripts)== 3){
            # distance between 3 points, must be 3D
            solo_transcripts <- myFun_3point_singleCell(solo_transcripts, 
                                                        transSpatLocs_coln = transSpatLocs_coln, 
                                                        distance_cutoff = orphan_cutoff[each_cell], 
                                                        startGroup = groupNum)
            
          } else {
            # too many solo transcripts, all listed as individual transcript group
            solo_transcripts[['transcript_group']]<- seq(groupNum, groupNum + nrow(solo_transcripts) -1)
          }
          
          # update dfCoord_subset with solo_transcripts
          solo_converter <- solo_transcripts[['transcript_group']]
          names(solo_converter) <- solo_transcripts[['tmp_transID']]
          dfCoord_subset[['transcript_group']][is.na(dfCoord_subset[['transcript_group']])] <- solo_converter[dfCoord_subset[['tmp_transID']][is.na(dfCoord_subset[['transcript_group']])]]
        }
        
      }else{
        # if still NULL network, get assign separate groupID
        dfCoord_subset[['transcript_group']] <- seq_len(nrow(dfCoord_subset))
        message(sprintf("%s return NULL delaunay network, has %d unique coordinates as separate group. ", 
                        each_cell, nrow(dfCoord_subset)))
        
      }
      # data. table
      dfCoord_subset <- dfCoord_subset[, .SD, .SDcols = c(transSpatLocs_coln, 'transcript_group')]
    }
    # assign directly if all unique
    if(nrow(df_subset) == nrow(dfCoord_subset)){
      df_subset[['transcript_group']] <- dfCoord_subset[['transcript_group']]
    } else {
      # merge by coordinate if duplicate coordinates exist
      df_subset <- merge(df_subset, dfCoord_subset, by = transSpatLocs_coln)
    }
    
    return(df_subset)
  }
  
  chosenTrans_df4 <-  chosenTrans_df[which(chosenTrans_df[[cellID_coln]] %in% count_df[N >3, get]), ]
  if(nrow(chosenTrans_df4) >0){
    tmp_group <- by(chosenTrans_df4, chosenTrans_df4[[cellID_coln]], my_fun_delanuay)
    tmp_group <- do.call(rbind, tmp_group)
    group_converter <- tmp_group[['transcript_group']] 
    names(group_converter) <- tmp_group[[transID_coln]]
    tmp_idx <- which(transcript_df[[transID_coln]] %in% tmp_group[[transID_coln]])
    transcript_df[['transcript_group']][tmp_idx] <- group_converter[transcript_df[[transID_coln]][tmp_idx]]
  }
  
  return(transcript_df)
}


#' @title myFun_3point_singleCell
#' @description supporting function for \code{groupTranscripts_Delanuay}, assign group ID for 3 transcripts in single cell in 3D based on distant cutoff
#' @param dfCoord_subset transcript data.table for single cell with only 3 transcripts in rows
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param distance_cutoff maximum distance within connected transcript group 
#' @return a data.table with `transcript_group` column added to original input data.table
#' @importFrom data.table as.data.table setDT
# function for 3 point case in a single cell only
myFun_3point_singleCell <- function(dfCoord_subset, 
                                    transSpatLocs_coln = c('x','y','z'),
                                    distance_cutoff = 2.7, 
                                    startGroup = 1){
  
  dfCoord_subset <- data.table::as.data.table(dfCoord_subset)
  dist_df <- cbind(dfCoord_subset[1:3, .SD, .SDcols = transSpatLocs_coln], dfCoord_subset[c(2,3,1), .SD, .SDcols = transSpatLocs_coln])
  colnames(dist_df) <- c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))
  dist_df[['from']] <- as.character(seq_len(nrow(dfCoord_subset)))
  dist_df[['to']] <- as.character(seq_len(nrow(dfCoord_subset)))[c(2,3,1)]
  dist_df <- data.table::setDT(dist_df)
  
  dist_df[, `:=`(distance, stats::dist(x = matrix(.SD, nrow = 2, byrow = T))), 
          by = 1:nrow(dist_df), 
          .SDcols = c(paste0(transSpatLocs_coln,'_begin'), paste0(transSpatLocs_coln,'_end'))]
  dist_df[, `:=`(distance, as.numeric(distance))]
  dist_df[, `:=`(weight, 1/distance)]
  dist_df[, `:=` (separate, distance >distance_cutoff)]
  dist_df[['group']] <- startGroup
  if(dist_df[['separate']][1]){
    # separate 1 and 2
    dist_df[['group']][2] <- startGroup+1
  }
  if(dist_df[['separate']][2]){
    # separate 2 and 3
    dist_df[['group']][3] <- dist_df[['group']][2]+1
  }
  
  if(!dist_df[['separate']][3]){
    # keep 1 and 3
    dist_df[['group']][3] <- dist_df[['group']][1]
  }
  
  dfCoord_subset[['transcript_group']] <- dist_df[['group']]
  
  return(dfCoord_subset)
}


#' @title groupTranscripts_dbscan
#' @description group the flagged transcript within each cell based on spatial clustering using dbscan
#' @param chosen_transcripts the transcript_id of chosen transcript
#' @param distance_cutoff maximum molecule-to-molecule distance within same transcript group (default = "auto")
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @return data frame of connected transcripts among chosen_transcripts
#' #' \enumerate{
#'    \item{cellID_coln, orignal cell_ID}
#'    \item{transID_coln, connected transcripts among chosen_transcripts}
#'    \item{transSpatLocs_coln, spatial coordiantes of transcript}
#'    \item{transcript_group, group of chosen_transcripts}
#' } 
#' @details For query cell, group flagged transcripts only based on their molecular distance to each other. When distance cutoff = 'auto', use 20% average XY cell range as cutoff. In case of no more than 3 flagged transcripts per cell, determine the grouping based on distance cutoff directly. In case of more transcripts per cell, use \code{dbscan} to group transcripts with distance_cutoff as `eps` and `minPts = 1`.  
#' @importFrom dbscan dbscan
#' @export
groupTranscripts_dbscan <- function(chosen_transcripts = NULL, 
                                    distance_cutoff = 'auto',
                                    transcript_df, 
                                    cellID_coln = "CellId", transID_coln = "transcript_id",
                                    transSpatLocs_coln = c('x','y','z')){
  
  
  if(distance_cutoff =='auto'){
    message("Automatic cutoff as 20% diameter of query cell.")
  } else if (is.numeric(distance_cutoff)){
    if (distance_cutoff <= 0){
      stop("distance_cutoff must be positive number")
    }
  }else {
    stop("distance_cutoff must be positive number")
  }
  
  d2_or_d3 <- length(transSpatLocs_coln)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("transSpatLocs_coln must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Do spatial network analysis in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transSpatLocs_coln)]
  
  if(is.null(chosen_transcripts)){
    stop("chosen_transcripts = NULL. No resegmentation is done.")
  } else {
    chosen_transcripts <- intersect(chosen_transcripts, transcript_df[[transID_coln]])
    if(length(chosen_transcripts) <1){
      stop("Common cells contain no provided chosen_transcripts. No resegmentation is done.")
    } else {
      message(sprintf("%d chosen_transcripts are found in common cells.", length(chosen_transcripts)))
    }
  }
  
  
  # which cells contained the chosen_transcripts
  # only cells with chosen transcript
  chosen_cells <- unique(transcript_df[which(transcript_df[[transID_coln]] %in% chosen_transcripts), get(cellID_coln)])
  
  # all transcripts of chosen_cells
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  ## operate in vector form to speed up the process
  # (1) get distance cutoff based on cell size if "auto"
  if(distance_cutoff == 'auto'){
    # size range of cells
    cellsize_df <- transcript_df[, lapply(.SD, function(x) (range(x)[2] - range(x)[1])),  by = .(get(cellID_coln)), .SDcols = transSpatLocs_coln]
    # cutoff for absolute distance between orphan transcript, defined as 20% of xy cell size range as a cube
    orphan_cutoff <- 0.2*rowMeans(cellsize_df[, 2:3])
    names(orphan_cutoff) <- cellsize_df[[1]]
    
    transcript_df[['distCutoff']] <- orphan_cutoff[transcript_df[[cellID_coln]]]
    
  }else{
    transcript_df[['distCutoff']] <- distance_cutoff
    
  }
  
  # (2) apply dbscan to the flagged transcripts of each cells
  idx_flagTrans <- which(transcript_df[[transID_coln]] %in% chosen_transcripts)
  groupDF <- by(transcript_df[idx_flagTrans, ], 
                transcript_df[[cellID_coln]][idx_flagTrans], 
                function(df_subset){
                  df_subset[['transcript_group']] <- dbscan::dbscan(as.data.frame(df_subset)[, transSpatLocs_coln], 
                                                                    eps = df_subset[['distCutoff']][1], minPts = 1)$cluster
                  return(df_subset)
                })
  groupDF <- do.call(rbind, groupDF)
  # combine back to original transcript DF, assign 0 for unflagged transcripts
  idx_otherTrans <- setdiff(seq_len(nrow(transcript_df)), idx_flagTrans)
  transcript_df <- transcript_df[idx_otherTrans, ]
  transcript_df[['transcript_group']] <- 0
  transcript_df <- rbind(transcript_df, groupDF)
  
  transcript_df[['distCutoff']] <- NULL
  
  return(transcript_df)
}

