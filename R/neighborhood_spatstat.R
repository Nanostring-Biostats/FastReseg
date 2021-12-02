
#' @title neighborhood_for_resegment_spatstat
#' @description find neighbor cells with transcripts that are direct neighbor of chosen_cell, check tLLRv2 score under neighbor cell type, return neighborhood information
#' @param chosen_cells the cell_ID of chosen cells need to be evaluate for re-segmentation
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param score_baseline a named vector of score baseline for all cell type listed in score_GeneMatrix
#' @param neighbor_distance_xy the distance in x, y from the center of each chosen cell to the center of neighborhood network building at transcript level. Default = NULL to use the 2 times of average cell diameter.
#' @param distance_cutoff maximum distance within connected transcript group. Default = NULL to use the 10 times of 90% quantile of minimal molecular distance within all chosen cells.
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param celltype_coln the column name of cell_type in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @importFrom spatstat.geom ppp subset.ppp nncross
#' @return a data.frame 
#' #' \enumerate{
#'    \item{CellId, original cell id of chosen cells}
#'    \item{cell_type, original cell type of chosen cells}
#'    \item{transcript_num, number of transcripts in chosen cells}
#'    \item{self_celltype, cell type give maximum score for query cell only}
#'    \item{score_under_self, score in query cell under its own maximum celltype}
#'    \item{neighbor_CellId, cell id of neighbor cell whose cell type gives maximum score in query cell among all neighbors, not including query cell itself}
#'    \item{neighbor_celltype, cell type that gives maximum score in query cell among all non-self neighbor cells}
#'    \item{score_under_neighbor, score in query cell under neighbor_celltype}
#' }
#' @details Locate neighbor cells of each query cell in 1st and 2nd dimension, firstly via cell-to-cell distance within neighbor_distance_xy, then via molecule-to-molecule distance within distance_cutoff. If no neighbor cells found for query cell, use the cell id and cell type of query cell to fill in the columns for neighbor cells in returned data.frame
#' @export
neighborhood_for_resegment_spatstat <- function(chosen_cells = NULL, 
                                                score_GeneMatrix,  
                                                score_baseline = NULL, 
                                                neighbor_distance_xy = NULL,
                                                distance_cutoff = NULL,
                                                transcript_df, 
                                                cellID_coln = "CellId", 
                                                celltype_coln = "cell_type", 
                                                transID_coln = "transcript_id",
                                                transGene_coln = "target", 
                                                transSpatLocs_coln = c('x','y','z')){
                                       
  if(is.null(chosen_cells)){
    stop("Must define chosen_cells to start resegmentation evaluation in neighborhood of each chosen cell.")
  }
  
  if(!is.null(score_baseline)){
    if(!any(class(score_baseline) %in% c('numeric'))){
      stop("The provided score_baseline must be either a named numeric vector or NUll which would disable comparision with score_baseline. ")
    }
    if(length(setdiff(colnames(score_GeneMatrix), names(score_baseline)))>0){
      stop(sprintf("The provided score_baseline is missing for the following cell types used in score_GeneMatrix: `%s`.", 
                   paste0(setdiff(colnames(score_GeneMatrix), names(score_baseline)), collapse ="`, `")))
    }
  }
  
  
  
  # numeric format or null
  if(!is.null(neighbor_distance_xy)){
    if(!any(class(neighbor_distance_xy) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for cell network, neighbor_distance_xy must be either NULL to use 2 times of average cell diameter or a numeric value to define the largest cell-to-cell distance.")
    } else if(neighbor_distance_xy <0){
      stop("neighbor_distance_xy must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
    } 
  } 
  
  if(!is.null(distance_cutoff)){
    if(!any(class(distance_cutoff) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcript network, distance_cutoff must be either NULL to use 10 times of 90% quantile of minimal molecular distance within all chosen cells or a numeric value to define the largest molecule-to-molecule distance.")
    } else if (distance_cutoff <= 0){
      stop("distance_cutoff must be either `auto` or positive number")
    } else{
      message(sprintf('Use distance_cutoff = %.4f for .', distance_cutoff))
    }
  } 
  
  
  d2_or_d3 <- length(transSpatLocs_coln)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Use the 1st and 2nd Dimension out of the provided %dD data for neighborhood search.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, celltype_coln, transID_coln, transGene_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, celltype_coln, transID_coln, transGene_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, celltype_coln,transID_coln, transGene_coln, transSpatLocs_coln)]
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
    
  # get common cells
  common_cells <- unique(transcript_df[[cellID_coln]])
  
  # get common genes
  common_genes <- intersect(rownames(score_GeneMatrix), 
                            unique(transcript_df[[transGene_coln]]))
  message(sprintf("Found %d common cells and %d common genes among transcript_df, cell_networkDT, and score_GeneMatrix. ", 
                  length(common_cells), length(common_genes)))
  
  if(any(length(common_cells) <2, length(common_genes)<1)){
    stop("Too few common cells or genes to proceed. Check if score_GeneMatrix is a gene x cell-type matrix.")
  }
  
  if(is.null(chosen_cells)){
    message("chosen_cells = NULL, analyze all common cells shared between transcript_df and cell_networkDT.")
    chosen_cells <- common_cells
  } else {
    chosen_cells <- intersect(chosen_cells, common_cells)
    if(length(chosen_cells) <1){
      stop("Common cells contain no provided chosen_cells. No calculation is done.")
    } else {
      message(sprintf("%d chosen_cells are found in common cells.", length(chosen_cells)))
    }
  }
  
  score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells & transcript_df[[transGene_coln]] %in% common_genes), ]
  
  # get per cell dataframe and search range if neighbor_distance_xy = NULL
  perCell_coordM <- transcript_df[, list(CenterX = mean(get(transSpatLocs_coln[1])), 
                                         CenterY = mean(get(transSpatLocs_coln[2])),
                                         Width = diff(range(get(transSpatLocs_coln[1]))),
                                         Height = diff(range(get(transSpatLocs_coln[2])))), 
                                  by = get(cellID_coln)]
  colnames(perCell_coordM)[1] <- cellID_coln
  if(is.null(neighbor_distance_xy)){
    neighbor_distance_xy <- max(colMeans(perCell_coordM[, c('Width','Height')]))*2
    message(sprintf("Use 2 times of average cell diameter as neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
  }
  
  # get neighbors cells for all chosen cells
  perCell_pp <- spatstat.geom::ppp(x = perCell_coordM[['CenterX']], 
                                   y = perCell_coordM[['CenterY']], 
                                   range(perCell_coordM[['CenterX']]), 
                                   range(perCell_coordM[['CenterY']]), 
                                   marks = factor(perCell_coordM[[cellID_coln]]), 
                                   unitname = c("um","um"))
  # subset to get ppp for all query cells
  query_pp <- spatstat.geom::subset.ppp(perCell_pp, marks %in% chosen_cells)
  
  # get closest 9 neighbors (exclude itself) data.frame, return what = which as index in original ppp  
  neighbor_matrix <- as.matrix(spatstat.geom::nncross(query_pp, perCell_pp, what = "which", k = 2:10))
  rownames(neighbor_matrix) <- perCell_pp$marks[perCell_pp$marks %in% chosen_cells]
  # filter based on distance
  distFlag <- which(spatstat.geom::nncross(query_pp, perCell_pp, what = "dist", k = 2:10) > neighbor_distance_xy)
  neighbor_matrix[distFlag] <- NA
  
  # data for chosen cells only
  chosen_transDF <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  # get distance cutoff for all chosen cells at transcript level
  if(is.null(distance_cutoff)){
    queryTrans_pp <- spatstat.geom::ppp(x = chosen_transDF[[transSpatLocs_coln[1]]], 
                                        y = chosen_transDF[[transSpatLocs_coln[2]]], 
                                        range(chosen_transDF[[transSpatLocs_coln[1]]]), 
                                        range(chosen_transDF[[transSpatLocs_coln[2]]]), 
                                        marks = factor(chosen_transDF[[cellID_coln]]), 
                                        unitname = c("um","um"))
    # get distribution of minimal molecule-to-molecule distance for each transcript in query cell
    dist_profle <- quantile(spatstat.geom::nndist(queryTrans_pp), seq(0,1,by=0.1))
    message(sprintf("Distribution of minimal molecular distance within all chosen cells is %s, at quantile = %s.", 
                    paste0(round(dist_profile, 2), collapse = ", "), 
                    paste0(names(dist_profile), collapse = ", ")))
    # define cutoff as 10 times of 90% quantile value
    distance_cutoff <- 10*dist_profile[['90%']]
    message(sprintf("Use 10 times of 90% quantile of minimal molecular distance within all query cells as distance_cutoff = %.4f for defining direct neighbor cells.", distance_cutoff))
    rm(queryTrans_pp)
    }
  
                                   
  
  # functions for each cell
  my_fun_eachCell <- function(query_df){
    each_cell <- query_df[[cellID_coln]][1]
    query_transID <- query_df[[transID_coln]] 
    
    # get neighbor cells
    neighborCells_NT <- perCell_pp$marks[neighbor_matrix[each_cell,]]
    neighborCells_NT <- as.character(neighborCells_NT[!is.na(neighborCells_NT)])
    
    # for neighbor cells only
    df_subsetNT <- transcript_df[get(cellID_coln) %in% neighborCells_NT, ]
    
    # for query cell only
    df_subsetQ <- transcript_df[get(cellID_coln) == each_cell, ]
    
    ### use spatstat to find direct neighbor cells based on transcript-to-transcript distance in xy
    neighborhood_pp <- spatstat.geom::ppp(x = df_subsetNT[[transSpatLocs_coln[1]]], 
                                          y = df_subsetNT[[transSpatLocs_coln[2]]], 
                                          range(df_subsetNT[[transSpatLocs_coln[1]]]), 
                                          range(df_subsetNT[[transSpatLocs_coln[2]]]), 
                                          marks = factor(df_subsetNT[[cellID_coln]]), 
                                          unitname = c("um","um"))
    neiQuery_pp <- spatstat.geom::ppp(x = df_subsetQ[[transSpatLocs_coln[1]]], 
                                      y = df_subsetQ[[transSpatLocs_coln[2]]], 
                                      range(df_subsetQ[[transSpatLocs_coln[1]]]), 
                                      range(df_subsetQ[[transSpatLocs_coln[2]]]), 
                                      marks = factor(df_subsetQ[[cellID_coln]]), 
                                      unitname = c("um","um"))
    
    # get minimal distance of each neighborhood cell to any transcript inside query cell
    neighbor_distDF <- aggregate(spatstat.geom::nncross(neighborhood_pp, neiQuery_pp, what = "dist", k = 1), 
                                 list(neighborhood_pp$marks), FUN = min)
    colnames(neighbor_distDF) <- c('neighbor_cellID', 'min_dist')
    # re-order based on minimal molecular distance of neighbor cells to query cell
    neighbor_distDF <- neighbor_distDF[order(neighbor_distDF$min_dist), ]
    # direct cell neighbors should be within the distance_cutoff
    directCell_neighbors <- neighbor_distDF$neighbor_cellID[neighbor_distDF$min_dist <= distance_cutoff]
    
    # get score matrix for each transcript in query cell
    cell_score <- score_GeneMatrix[query_df[[transGene_coln]], ]
    if(nrow(query_df) >1 ){
      cell_score <- colSums(cell_score)
    }
    cell_score <- matrix(cell_score, nrow = 1, dimnames = list(each_cell, colnames(score_GeneMatrix)))
    max_idx_1st <- max.col(cell_score, ties.method="first")
    
    
    queryPerCell_df <- data.frame(CellId = each_cell, 
                                  cell_type = query_df[[celltype_coln]][1], 
                                  transcript_num = nrow(query_df), 
                                  self_celltype = colnames(cell_score)[max_idx_1st], 
                                  score_under_self = cell_score[each_cell, max_idx_1st]/nrow(query_df))
    
    
    # get neighbor cell type and check the query transcript cell type against it
    if(!is.na(directCell_neighbors)){
      # found connected neighbor cells with close enough transcripts
      neighbors_df <- df_subsetNT[which(df_subsetNT[[cellID_coln]] %in% directCell_neighbors),]
      neighbors_df <- as.data.frame(neighbors_df)[, c(cellID_coln, celltype_coln)]
      neighbors_df <- unique(neighbors_df)
      neighbors_df[['score_under_neighbor']] <- cell_score[each_cell, neighbors_df[[celltype_coln]]]/nrow(query_df)
      
      # if score baseline is provided, compare the net score above baseline to choose the consistent neighbor cell types
      if(!is.null(score_baseline)){
        neighbors_df[['baseline']] <- score_baseline[neighbors_df[[celltype_coln]]]
        neighbors_df[['net_score']] <- neighbors_df[['score_under_neighbor']]  - neighbors_df[['baseline']] 
      }
      
      # if more than 1 neighbor
      if(length(directCell_neighbors) >1){
        # ideally, to rank the neighbor cells based on their shortest edge to transcript of query cells
        # such that for multiple neighbor cells of same consistent cell type, choose the neighbor with shortest edge
        # re-order based on directCell_neighbors
        neighbors_df <- neighbors_df[match(directCell_neighbors, neighbors_df[[cellID_coln]]), ]
        
        # if score baseline is provided, compare the net score above baseline to choose the consistent neighbor cell types
        if(is.null(score_baseline)){
          tmp_score_neighbor <- matrix(neighbors_df[['score_under_neighbor']], nrow = 1, dimnames = list(each_cell, neighbors_df[[cellID_coln]]))
        } else {
          tmp_score_neighbor <- matrix(neighbors_df[['net_score']], nrow = 1, dimnames = list(each_cell, neighbors_df[[cellID_coln]]))
        }
        
        max_idx_neighbor <- max.col(tmp_score_neighbor, ties.method="first")
        
        
      } else {
        max_idx_neighbor = 1
      }
      
      neighborInfo_df <- data.frame(neighbor_CellId = neighbors_df[[cellID_coln]][max_idx_neighbor], 
                                    neighbor_celltype = neighbors_df[[celltype_coln]][max_idx_neighbor], 
                                    score_under_neighbor = neighbors_df[['score_under_neighbor']][max_idx_neighbor])
      
    } else{
      # no neighbor cells, keep original cell type
      neighborInfo_df <- data.frame(neighbor_CellId = each_cell, 
                                    neighbor_celltype = queryPerCell_df[['self_celltype']], 
                                    score_under_neighbor = queryPerCell_df[['score_under_self']])
    }
    
    outputs <- cbind(queryPerCell_df, neighborInfo_df)
    return(outputs)
  }
  
  # for each chosen cell, find out its neighborhood
  resegmented_df <- by(chosen_transDF, chosen_transDF[[cellID_coln]], my_fun_eachCell)
  resegmented_df <- do.call(rbind, resegmented_df)
  return(resegmented_df)
}


#' get neighborhood transcript data frame for visualization
#' @title getNeighbors_transDF_spatstat
#' @description find neighbor cells of chosen_cells and return the relevant transcript data.frame for both query 
#' @param chosen_cells the cell_ID of chosen cells need to be evaluate for re-segmentation
#' @param neighbor_distance_xy the distance in x, y from the center of each chosen cell to the center of neighborhood network building at transcript level. Default = NULL to use the 2 times of average cell diameter.
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @importFrom spatstat.geom ppp subset.ppp nncross
#' @return a data.frame 
#' \enumerate{
#'    \item{transcript_id, transcript id of all transcript near chosen cell's neighborhood}
#'    \item{transSpatLocs_coln, coordinates}
#'    \item{in_query_cell, flag whether current transcript is in the query cell, false means in neighbor cells}
#'    \item{CellId, original cell id of each transcript}
#'    \item{query_CellId, original query cell id of transcript's neighborhood}
#' }
#' @details Locate neighbor cells of each query cell in 1st and 2nd dimension, firstly via cell-to-cell distance within neighbor_distance_xy. If no neighbor cells found for query cell, return query cell information only. Do not consider extracellular transcripts.
#' @importFrom dplyr between
#' @export
getNeighbors_transDF_spatstat <- function(chosen_cells = NULL, 
                                           neighbor_distance_xy = NULL,
                                           transcript_df, 
                                           cellID_coln = "CellId", 
                                           transID_coln = "transcript_id",
                                           transSpatLocs_coln = c('x','y','z')){
                                 
  if(is.null(chosen_cells)){
    stop("Must define chosen_cells to start resegmentation evaluation in neighborhood of each chosen cell.")
  }
  
  
  # numeric format or null
  if(!is.null(neighbor_distance_xy)){
    if(!any(class(neighbor_distance_xy) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcriptnetwork, neighbor_distance_xy must be either NULL to use 2 times of average cell diameter or a numeric value to define the coordiante range.")
    } else if(neighbor_distance_xy <0){
      stop("neighbor_distance_xy must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
    } 
  } 
  
  
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transSpatLocs_coln)]
  data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))

  common_cells <- unique(transcript_df[[cellID_coln]])
  message(sprintf("Found %d common cells among transcript_df, cell_networkDT. ", 
                  length(common_cells)))
  
  if(any(length(common_cells) <2)){
    stop("Too few common cells to proceed.")
  }
  
  if(is.null(chosen_cells)){
    message("chosen_cells = NULL, analyze all common cells shared between transcript_df and cell_networkDT.")
    chosen_cells <- common_cells
  } else {
    chosen_cells <- intersect(chosen_cells, common_cells)
    if(length(chosen_cells) <1){
      stop("Common cells contain no provided chosen_cells. No calculation is done.")
    } else {
      message(sprintf("%d chosen_cells are found in common cells.", length(chosen_cells)))
    }
  }
  
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells), ]
  # get per cell dataframe and search range if neighbor_distance_xy = NULL
  perCell_coordM <- transcript_df[, list(CenterX = mean(get(transSpatLocs_coln[1])), 
                                         CenterY = mean(get(transSpatLocs_coln[2])),
                                         Width = diff(range(get(transSpatLocs_coln[1]))),
                                         Height = diff(range(get(transSpatLocs_coln[2])))), 
                                  by = get(cellID_coln)]
  colnames(perCell_coordM)[1] <- cellID_coln
  if(is.null(neighbor_distance_xy)){
    neighbor_distance_xy <- max(colMeans(perCell_coordM[, c('Width','Height')]))*2
    message(sprintf("Use 2 times of average cell diameter as neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
  }
  
  # data for chosen cells only
  chosen_transDF <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  # get neighbors cells for all chosen cells
  perCell_pp <- spatstat.geom::ppp(x = perCell_coordM[['CenterX']], 
                                   y = perCell_coordM[['CenterY']], 
                                   range(perCell_coordM[['CenterX']]), 
                                   range(perCell_coordM[['CenterY']]), 
                                   marks = factor(perCell_coordM[[cellID_coln]]), 
                                   unitname = c("um","um"))
  # subset to get ppp for all query cells
  query_pp <- spatstat.geom::subset.ppp(perCell_pp, marks %in% chosen_cells)
  
  # get closest 9 neighbors (exclude itself) data.frame, return what = which as index in orignal ppp  
  neighbor_matrix <- as.matrix(spatstat.geom::nncross(query_pp, perCell_pp, what = "which", k = 2:10))
  rownames(neighbor_matrix) <- perCell_pp$marks[perCell_pp$marks %in% chosen_cells]
  # filter based on distance
  distFlag <- which(spatstat.geom::nncross(query_pp, perCell_pp, what = "dist", k = 2:10) > neighbor_distance_xy)
  neighbor_matrix[distFlag] <- NA
  
  # functions for each cell
  my_fun_eachCell <- function(query_df){
    each_cell <- query_df[[cellID_coln]][1]
    query_transID <- query_df[[transID_coln]] 
    
    neighborCells_NT <- perCell_pp$marks[neighbor_matrix[each_cell,]]
    neighborCells_NT <- as.character(neighborCells_NT[!is.na(neighborCells_NT)])
    # data frame for each cell and its neighbor cells
    df_subset <- transcript_df[get(cellID_coln) %in% unique(c(each_cell, neighborCells_NT)), ]
    
    df_subset[['query_CellId']] <- rep(each_cell, nrow(df_subset))
    df_subset[['in_query_cell']] <- (df_subset[[cellID_coln]] == df_subset[['query_CellId']])
    
    return(df_subset)
    
  }
  
  # for each chosen cell, find out its neighborhood
  neighborhood_df <- by(chosen_transDF, chosen_transDF[[cellID_coln]], my_fun_eachCell)
  neighborhood_df <- do.call(rbind, neighborhood_df)
  return(neighborhood_df)
}