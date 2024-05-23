# modular wrapper to flag cell segmentation error
#' @title runSegErrorEvaluation
#' @description modular wrapper to flag cell segmentation error 
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param transcript_df the data.frame of transcript_ID, cell_ID, score, spatial coordinates
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
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
  # `getCellType_maxScore` function returns a named vector with cell type in values and cell_ID in names
  celltype_cellVector <- getCellType_maxScore(score_GeneMatrix = score_GeneMatrix, 
                                              transcript_df = transcript_df, 
                                              transGene_coln = transGene_coln,
                                              cellID_coln = cellID_coln)
  
  transcript_df[['tLLR_maxCellType']] <- celltype_cellVector[transcript_df[[cellID_coln]]]
  transcript_df <- transcript_df[!is.na(transcript_df[['tLLR_maxCellType']]), ]
  message(sprintf("Found %d cells and assigned cell type based on the provided `refProfiles` cluster profiles.", length(celltype_cellVector)))
  
  
  ##  for each transcript, calculate tLLR score based on the max cell type
  # `getScoreCellType_gene` function returns a named vector with score of given cell type in values and transcript_id in names
  score_transVector <- getScoreCellType_gene(score_GeneMatrix = score_GeneMatrix, 
                                             transcript_df = transcript_df, 
                                             transID_coln = transID_coln,
                                             transGene_coln = transGene_coln,
                                             celltype_coln = 'tLLR_maxCellType')
  transcript_df[['score_tLLR_maxCellType']] <- score_transVector[transcript_df[[transID_coln]]]
  transcript_df <- transcript_df[!is.na(transcript_df[['score_tLLR_maxCellType']]), ]
  rm(score_transVector)
  
  ## spatial modeling of tLLR score profile within each cell to identify cells with strong spatial dependency 
  # `score_cell_segmentation_error` function returns a data.frame with cell in row and spatial modeling outcomes in columns
  modStats_ToFlagCells <- score_cell_segmentation_error(
    chosen_cells = names(celltype_cellVector), 
    transcript_df = transcript_df, 
    cellID_coln = cellID_coln, 
    transID_coln = transID_coln, 
    score_coln = 'score_tLLR_maxCellType',
    spatLocs_colns = spatLocs_colns, 
    model_cutoff = flagModel_TransNum_cutoff)
  
  if(!is.null(modStats_ToFlagCells)){
    #-log10(P)
    modStats_ToFlagCells[['lrtest_nlog10P']] <- (-log10(modStats_ToFlagCells[['lrtest_Pr']]))
    modStats_ToFlagCells[['tLLR_maxCellType']] <- celltype_cellVector[modStats_ToFlagCells[['cell_ID']]]
    if(cellID_coln != 'cell_ID'){
      colnames(modStats_ToFlagCells)[which(colnames(modStats_ToFlagCells) == 'cell_ID')] <- cellID_coln
    }
  }
  
  outs <- list(modStats_ToFlagCells = modStats_ToFlagCells, 
               transcript_df = transcript_df)
  
  return(outs)
}


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
#' @importFrom stats as.formula lm
#' @importFrom data.table .N .SD
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
  
  # check if any cell required evaluation 
  if(nrow(transcript_df)>0){
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
  } else {
    # no cells for evaluation, return NULL for model_stats
    message(sprintf("No single cell with transcript number above model_cutoff = %s, skip the evaluation.", 
                    as.character(model_cutoff)))
    model_stats <- NULL
    
  }
  
  
  return(model_stats)
  
}