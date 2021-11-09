
#' @title neighborhood_for_resegment
#' @description find neighbor cells with transcripts that are direct neighbor of chosen_cell, check tLLRv2 score under neighbor cell type, return neighborhood information
#' @param chosen_cells the cell_ID of chosen cells need to be evaluate for re-segmentation
#' @param score_GeneMatrix the gene x cell-type matrix of log-like score of gene in each cell type
#' @param score_baseline a named vector of score baseline for all cell type listed in score_GeneMatrix
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level
#' @param cell_networkDT the data.table component of cell-level delaunay network object, would be used for neighborhood search if neighbor_distance_xy is NULL.
#' @param network_cellID_split the deliminator used to convert cellID_coln in transcript_df to the cell ID in cell_networkDT
#' @param neighbor_distance_xy the distance in x, y from the outermost transcripts of each chosen cell to the furthest non-self transcripts for neighborhood network building at transcript level. Default = NULL to use the cell network for selection of neighborhood.
#' @param neighbor_distance_z the distance in z direction from the outermost transcripts of each chosen cell to the furthest non-self transcripts for neighborhood network building at transcript level. Default = NULL to use the cell network for selection of neighborhood.
#' @param distance_cutoff maximum distance within connected transcript group (default = "auto").
#' @param transNum_cutoff minimal number of transcripts in query cell and its neighborhood for delaunay network analysis (default = 4).
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param celltype_coln the column name of cell_type in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transGene_coln the column name of target or gene name in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param gridSpat_coln an optional vector of column names in transcript_df that could act as grid to separate different transcripts in space, such as columns for FOV, slideID when using direct distance threshold to define neighborhood search range.
#' @import dplyr
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
#' @details If no neighbor cells found for query cell, use the cell id and cell type of query cell to fill in the columns for neighbor cells in returned data.frame
#' @export
neighborhood_for_resegment <- function(chosen_cells = NULL, 
                                       score_GeneMatrix,  
                                       score_baseline = NULL, 
                                       cell_networkDT = NULL, 
                                       network_cellID_split = '_g',
                                       config_spatNW_transcript, 
                                       neighbor_distance_xy = NULL,
                                       neighbor_distance_z = NULL,
                                       distance_cutoff = 'auto',
                                       transNum_cutoff = 4,
                                       transcript_df, 
                                       cellID_coln = "CellId", 
                                       celltype_coln = "cell_type", 
                                       transID_coln = "transcript_id",
                                       transGene_coln = "target", 
                                       transSpatLocs_coln = c('x','y','z'), 
                                       gridSpat_coln = c('fov', 'slide')){
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
  
  if(transNum_cutoff <4){
    stop('transNum_cutoff must be no smaller than 4 to enable delanuay network analysis within neighhood of chosen_cells.')
  }
  
  if(is.null(neighbor_distance_xy) & is.null(cell_networkDT)){
    stop("Both neighbor_distance_xy and cell_networkDT are NULL, please provide either one to define the search area for neighborhood analysis.")
  }
  
  # numeric format or null
  if(!is.null(neighbor_distance_xy)){
    if(!any(class(neighbor_distance_xy) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcriptnetwork, neighbor_distance_xy must be either NULL to use cell_networkDT or a numeric value to define the coordiante range.")
    } else if(neighbor_distance_xy <0){
      stop("neighbor_distance_xy must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
      cell_networkDT <- NULL
      
      # check if the provided gridSpat_coln exists in transcript_df
      if(!is.null(gridSpat_coln)){
        gridSpat_coln <- intersect(gridSpat_coln, colnames(transcript_df))
        if(length(gridSpat_coln) ==0){
          message("The provided gridSpat_coln do not exist in transcript_df. Ignore spatial filtering with gridSpat_coln.")
          gridSpat_coln <- NULL
        } else {
          message(sprintf("Use gridSpat_coln = `%s` for spaital filtering when searching neighbor cells. ", 
                          paste0(gridSpat_coln, collapse = "`, `")))
        }
      }
    }
  } 
  
  
  if(!is.null(neighbor_distance_z)){
    if(!any(class(neighbor_distance_z) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcriptnetwork, neighbor_distance_z must be either NULL to use cell_networkDT or a numeric value to define the coordiante range.")
    }else if(neighbor_distance_z <0){
      stop("neighbor_distance_z must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_z = %.4f for searching of neighbor cells.", neighbor_distance_z))
    }
  }
  
  # check if right format of cell_networkDT
  if(!is.null(cell_networkDT)){
    if(!any(grepl('data.table', class(cell_networkDT)))){
      stop("The provided cell_networkDT must be a data.table.")
    }
    if(any(!c("from","to") %in% colnames(cell_networkDT))){
      stop("The provided cell_networkDT is not a spatial network object with both `from` and `to` columns.")
    }
    # not to use column based spaital filerting when cell_networkDT is provided. 
    gridSpat_coln <- NULL
  }
  
  
  if(distance_cutoff =='auto'){
    message("automatic cutoff for delauny network of transcript connectivity analysis.")
  } else if (is.numeric(distance_cutoff)){
    if (distance_cutoff <= 0){
      stop("distance_cutoff must be either `auto` or positive number")
    } else{
      message(sprintf('use distance_cutoff = %.4f for delauny network of transcript connectivity analysis.', distance_cutoff))
    }
  }else {
    stop("distance_cutoff must be either `auto` or positive number")
  }
  
  d2_or_d3 <- length(transSpatLocs_coln)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Run delanuay network in %d Dimension.", d2_or_d3))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, celltype_coln, transID_coln, transGene_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, celltype_coln, transID_coln, transGene_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  transcript_df <- data.table::as.data.table(transcript_df)
  # when using absolute distance for search, include gridSpat_coln
  if(!is.null(gridSpat_coln)){
    transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, celltype_coln,transID_coln, transGene_coln, transSpatLocs_coln, gridSpat_coln)]
    
    # add uid based on grid colns
    transcript_df[, gridAll_UID := apply(.SD, 1, paste0, collapse = '_'), .SDcols = gridSpat_coln]
    
    data.table::setkeyv(transcript_df, c('gridAll_UID', cellID_coln, transID_coln))
    
  } else {
    transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, celltype_coln,transID_coln, transGene_coln, transSpatLocs_coln)]
    data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  }
  
  
  # get common cells
  if(!is.null(cell_networkDT)){
    # get cellID in network
    transcript_df[['network_cellID']] <- sapply(strsplit(transcript_df[[cellID_coln]], split = network_cellID_split),"[[", 1)
    common_cells2 <- intersect(unique(transcript_df[['network_cellID']]), 
                               unique(c(cell_networkDT$from, cell_networkDT$to)))
    common_cells <- unique(transcript_df[[cellID_coln]][which(transcript_df[['network_cellID']] %in% common_cells2)])
  } else {
    common_cells <- unique(transcript_df[[cellID_coln]])
  } 
  
  
  
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
  
  if(!is.null(cell_networkDT)){
    cell_networkDT <- cell_networkDT[to %in% common_cells2 & from %in% common_cells2]
  }
  
  score_GeneMatrix <- score_GeneMatrix[common_genes, ]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells & transcript_df[[transGene_coln]] %in% common_genes), ]
  
  # check config_spatNW_transcript
  config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript,
                                                     spat_locs = transcript_df)
  # data for chosen cells only
  chosen_transDF <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  
  # function to search based on distance 
  my_function_searchDistance <- function(current_cell, current_transDT){
    # find range of coordinate for query cell
    df_subset <- as.data.frame(current_transDT)[which(current_transDT[[cellID_coln]] == current_cell), transSpatLocs_coln]
    cellsize_range <- apply(df_subset, 2, FUN = range)
    # add in the distance defined by neighbor_distance for search_range
    search_range <- cellsize_range + c(- neighbor_distance_xy, neighbor_distance_xy)
    if(!is.null(neighbor_distance_z)){
      # different range for 3rd/ z dimension
      if(length(transSpatLocs_coln) >2){
        # 3D coordinates
        search_range[,transSpatLocs_coln[3]] <- cellsize_range[,transSpatLocs_coln[3]] + c(- neighbor_distance_z, neighbor_distance_z)
      }
    }
    
    # find all transcripts within this search_range
    # operate with values
    # # need import entire dplyr to use %between%
    # chosen_transcripts <- lapply(seq_len(ncol(search_range)),
    #                              function(i) current_transDT[get(transSpatLocs_coln[i]) %between% search_range[,i], get(transID_coln)])
    
    # without import entire dplyr
    chosen_transcripts <- lapply(seq_len(ncol(search_range)),
                                 function(i) current_transDT[dplyr::between(get(transSpatLocs_coln[i]), 
                                                                            search_range[1, i], search_range[2, i]), 
                                                             get(transID_coln)])
    
    chosen_transcripts <- Reduce(intersect, chosen_transcripts)
    df_subset <- current_transDT[which(current_transDT[[transID_coln]] %in% chosen_transcripts), ]
    # # operate with TRUE and FALSE, slower
    # chosen_transcripts <- lapply(seq_len(ncol(search_range)),
    #                              function(i) current_transDT[[transSpatLocs_coln[i]]] %between% search_range[, i])
    # df_subset <- current_transDT[Reduce('&', chosen_transcripts), ]
    
    return(df_subset)
  }
  
  
  # functions for each cell
  my_fun_eachCell <- function(query_df){
    each_cell <- query_df[[cellID_coln]][1]
    query_transID <- query_df[[transID_coln]] 
    
    if(is.null(neighbor_distance_xy) & !is.null(cell_networkDT)){
      # define neighborhood using cell spatial network
      eachCell_NT <- query_df[['network_cellID']][1]
      neighborCells_NT <- c(cell_networkDT[from == eachCell_NT]$to, cell_networkDT[to == eachCell_NT]$from)
      neighborCells_NT <- unique(c(eachCell_NT, neighborCells_NT))
      df_subset <- transcript_df[which(transcript_df[['network_cellID']] %in% neighborCells_NT), ]
    } else {
      # define neighborhood based on absolute distance to query cell
      # subset data based on gridSpat_coln
      if(!is.null(gridSpat_coln)){
        # # current_grid as a vector
        # current_grid <- as.matrix(query_df[1, .SD, .SDcols = gridSpat_coln])[1,]
        # # operate with index
        # subTransIdx <- lapply(seq_len(length(gridSpat_coln)),
        #                       function(i) which(transcript_df[[gridSpat_coln[i]]] == current_grid[i]))
        # subTransIdx <- Reduce(intersect, subTransIdx)
        # subTransDT <- transcript_df[subTransIdx, ]
        # # operate with TRUE and FALSE, slower
        # subTransIdx <- lapply(seq_len(length(gridSpat_coln)),
        #                       function(i) transcript_df[[gridSpat_coln[i]]] == current_grid[i])
        # subTransDT <- transcript_df[Reduce('&', subTransIdx), ]
        
        # # current_grid to locate directly
        # subTransDT <- transcript_df[gridAll_UID == query_df[['gridAll_UID']][1], ]
        
        # use grid index to locate 
        current_grid <- query_df[['gridAll_UID']][1]
        subTransDT <- transcript_df[which(transcript_df[['gridAll_UID']] == current_grid), ]
        
        df_subset <- my_function_searchDistance(current_cell = each_cell, current_transDT = subTransDT)
        
      } else {
        # no subset, search directly
        df_subset <- my_function_searchDistance(current_cell = each_cell, current_transDT = transcript_df)
      }
      
    }
    
    
    # find connected neighborhood environment, be careful for the case with too few transcripts for delaunay 
    if(nrow(df_subset) < transNum_cutoff){
      # too few query and neighbor transcripts for network analysis and cell typing
      # treat as no valid direct neighbor
      connected_neighbor <- NULL
    }else{
      if(d2_or_d3 ==2){
        # run 2D
        delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                spatLocs_df = df_subset, 
                                                                ID_column = transID_coln,
                                                                spatLocs_column = transSpatLocs_coln)
      } else {
        # 3D specify, check if all selected transcript in same z plane
        z_num <- unique(df_subset[[transSpatLocs_coln[3]]])
        if(length(z_num) <2){
          # single z-plane, run as 2D
          message(sprintf("%s with all %d neighborhood transcripts in same z plane, run 2D network analysis.", each_cell, nrow(df_subset)))
          delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                  spatLocs_df = df_subset, 
                                                                  ID_column = transID_coln,
                                                                  spatLocs_column = transSpatLocs_coln[1:2])
        } else {
          # multiple z-plane, run as 3D
          delaunayNW_Obj <- createSpatialDelaunayNW_from_spatLocs(config_spatNW = config_spatNW_transcript, 
                                                                  spatLocs_df = df_subset, 
                                                                  ID_column = transID_coln,
                                                                  spatLocs_column = transSpatLocs_coln)
        }
        
        
      }
      
      transNetWorkDT <- delaunayNW_Obj$networkDT
      # remove edges that are longer than cutoff
      if(distance_cutoff !='auto'){
        # better to use auto, if not, use distance_cutoff >10 pixel = 1.8um, cutoff = 15 pixel = 2.7um works best 
        transNetWorkDT <- transNetWorkDT[distance <= distance_cutoff]
      }
      
      # does any neighbor belong to chosen_transcripts
      connect_flag <- tapply(transNetWorkDT$from, transNetWorkDT$to, function(x) any(x %in% query_transID))
      connected_all <- names(connect_flag)[connect_flag]
      connect_flag <- tapply(transNetWorkDT$to, transNetWorkDT$from, function(x) any(x %in% query_transID))
      connected_all <- unique(c(names(connect_flag)[connect_flag], connected_all))
      
      # transcripts in neighbor cells but connected to chosen_transcripts
      connected_neighbor <- setdiff(connected_all, query_transID)
      
    }
    
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
    if(length(connected_neighbor) >0){
      # found connected neighbor transcripts
      neighbors_df <- df_subset[which(df_subset[[transID_coln]] %in% connected_neighbor),]
      # ideally, to rank the neighbor cells based on their shortest edge to transcript of query cells
      # such that for multiple neighbor cells of same consistent cell type, choose the neighbor with shortest edge
      neighbors_df <- as.data.frame(neighbors_df)[, c(cellID_coln, celltype_coln)]
      neighbors_df <- unique(neighbors_df)
      neighbors_df[['score_under_neighbor']] <- cell_score[each_cell, neighbors_df[[celltype_coln]]]/nrow(query_df)
      
      # if score baseline is provided, compare the net score above baseline to choose the consistent neighbor cell types
      if(!is.null(score_baseline)){
        neighbors_df[['baseline']] <- score_baseline[neighbors_df[[celltype_coln]]]
        neighbors_df[['net_score']] <- neighbors_df[['score_under_neighbor']]  - neighbors_df[['baseline']] 
      }
      
      # if more than 1 neighbor
      if(nrow(neighbors_df) >1){
        # ideally, to rank the neighbor cells based on their shortest edge to transcript of query cells
        # such that for multiple neighbor cells of same consistent cell type, choose the neighbor with shortest edge
        idx1 <- which((transNetWorkDT$from %in% query_transID) & (transNetWorkDT$to %in% connected_neighbor))
        edge_df1 <- transNetWorkDT[idx1, ]
        edge_df1[['neighborCell']] <- df_subset[match(edge_df1$to, df_subset[[transID_coln]]),][[cellID_coln]]
        idx2 <- which((transNetWorkDT$to %in% query_transID) & (transNetWorkDT$from %in% connected_neighbor))  
        edge_df2 <- transNetWorkDT[idx2, ]
        edge_df2[['neighborCell']] <- df_subset[match(edge_df2$from, df_subset[[transID_coln]]),][[cellID_coln]]
        edge_df <- rbind(edge_df1, edge_df2)
        # get minimal distance
        edge_df <- aggregate(edge_df$distance, by = list(edge_df$neighborCell), FUN = "min")
        colnames(edge_df) <- c(cellID_coln, 'distance')
        
        neighbors_df <- merge(neighbors_df, edge_df, by = cellID_coln, all.x = TRUE)
        neighbors_df <- neighbors_df[order(neighbors_df$distance), ]
        
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
#' @title getNeighbors_transDF
#' @description find neighbor cells of chosen_cells and return the relevant transcript data.frame for both query 
#' @param chosen_cells the cell_ID of chosen cells need to be evaluate for re-segmentation
#' @param cell_networkDT the data.table component of cell-level delaunay network object, would be used for neighborhood search if neighbor_distance_xy is NULL.
#' @param network_cellID_split the deliminator used to convert cellID_coln in transcript_df to the cell ID in cell_networkDT
#' @param neighbor_distance_xy the distance in x, y from the outermost transcripts of each chosen cell to the furthest non-self transcripts for neighborhood network building at transcript level. Default = NULL to use the cell network for selection of neighborhood.
#' @param neighbor_distance_z the distance in z direction from the outermost transcripts of each chosen cell to the furthest non-self transcripts for neighborhood network building at transcript level. Default = NULL to use the cell network for selection of neighborhood.
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param gridSpat_coln an optional vector of column names in transcript_df that could act as grid to separate different transcripts in space, such as columns for FOV, slideID when using direct distance threshold to define neighborhood search range.
#' @return a data.frame 
#' \enumerate{
#'    \item{transcript_id, transcript id of all transcript near chosen cell's neighborhood}
#'    \item{transSpatLocs_coln, coordinates}
#'    \item{in_query_cell, flag whether current transcript is in the query cell, false means in neighbor cells}
#'    \item{CellId, original cell id of each transcript}
#'    \item{query_CellId, original query cell id of transcript's neighborhood}
#' }
#' @details If no neighbor cells found for query cell, return query cell information only. Do not consider extracellular transcripts.
#' @export
getNeighbors_transDF <- function(chosen_cells = NULL, 
                                 cell_networkDT = NULL, 
                                 network_cellID_split = '_g',
                                 neighbor_distance_xy = NULL,
                                 neighbor_distance_z = NULL,
                                 transcript_df, 
                                 cellID_coln = "CellId", 
                                 transID_coln = "transcript_id",
                                 transSpatLocs_coln = c('x','y','z'), 
                                 gridSpat_coln = c('fov', 'slide')){
  if(is.null(chosen_cells)){
    stop("Must define chosen_cells to start resegmentation evaluation in neighborhood of each chosen cell.")
  }
  
  
  if(is.null(neighbor_distance_xy) & is.null(cell_networkDT)){
    stop("Both neighbor_distance_xy and cell_networkDT are NULL, please provide either one to define the search area for neighborhood analysis.")
  }
  
  # numeric format or null
  if(!is.null(neighbor_distance_xy)){
    if(!any(class(neighbor_distance_xy) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcriptnetwork, neighbor_distance_xy must be either NULL to use cell_networkDT or a numeric value to define the coordiante range.")
    } else if(neighbor_distance_xy <0){
      stop("neighbor_distance_xy must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_xy = %.4f for searching of neighbor cells.", neighbor_distance_xy))
      cell_networkDT <- NULL
      
      # check if the provided gridSpat_coln exists in transcript_df
      if(!is.null(gridSpat_coln)){
        gridSpat_coln <- intersect(gridSpat_coln, colnames(transcript_df))
        if(length(gridSpat_coln) ==0){
          message("The provided gridSpat_coln do not exist in transcript_df. Ignore spatial filtering with gridSpat_coln.")
          gridSpat_coln <- NULL
        } else {
          message(sprintf("Use gridSpat_coln = `%s` for spaital filtering when searching neighbor cells. ", 
                          paste0(gridSpat_coln, collapse = "`, `")))
        }
      }
    }
  } 
  
  
  if(!is.null(neighbor_distance_z)){
    if(!any(class(neighbor_distance_z) %in% c('numeric','integer'))){
      stop("To define the neighborhood to consider for transcriptnetwork, neighbor_distance_z must be either NULL to use cell_networkDT or a numeric value to define the coordiante range.")
    }else if(neighbor_distance_z <0){
      stop("neighbor_distance_z must be positive number to define the neighborhood.")
    } else {
      message(sprintf("Use neighbor_distance_z = %.4f for searching of neighbor cells.", neighbor_distance_z))
    }
  }
  
  # check if right format of cell_networkDT
  if(!is.null(cell_networkDT)){
    if(!any(grepl('data.table', class(cell_networkDT)))){
      stop("The provided cell_networkDT must be a data.table.")
    }
    if(any(!c("from","to") %in% colnames(cell_networkDT))){
      stop("The provided cell_networkDT is not a spatial network object with both `from` and `to` columns.")
    }
    # not to use column based spaital filerting when cell_networkDT is provided. 
    gridSpat_coln <- NULL
  }
  
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  transcript_df <- data.table::as.data.table(transcript_df)
  # when using absolute distance for search, include gridSpat_coln
  if(!is.null(gridSpat_coln)){
    transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transSpatLocs_coln, gridSpat_coln)]
    
    # add uid based on grid colns
    transcript_df[, gridAll_UID := apply(.SD, 1, paste0, collapse = '_'), .SDcols = gridSpat_coln]
    
    data.table::setkeyv(transcript_df, c('gridAll_UID', cellID_coln, transID_coln))
    
  } else {
    transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln, transSpatLocs_coln)]
    data.table::setkeyv(transcript_df, c(cellID_coln, transID_coln))
  }
  
  
  # get common cells
  if(!is.null(cell_networkDT)){
    # get cellID in network
    transcript_df[['network_cellID']] <- sapply(strsplit(transcript_df[[cellID_coln]], split = network_cellID_split),"[[", 1)
    common_cells2 <- intersect(unique(transcript_df[['network_cellID']]), 
                               unique(c(cell_networkDT$from, cell_networkDT$to)))
    common_cells <- unique(transcript_df[[cellID_coln]][which(transcript_df[['network_cellID']] %in% common_cells2)])
  } else {
    common_cells <- unique(transcript_df[[cellID_coln]])
  } 
  
  
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
  
  if(!is.null(cell_networkDT)){
    cell_networkDT <- cell_networkDT[to %in% common_cells2 & from %in% common_cells2]
  }
  
  
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells), ]
  
  
  # data for chosen cells only
  chosen_transDF <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  
  # function to search based on distance 
  my_function_searchDistance <- function(current_cell, current_transDT){
    # find range of coordinate for query cell
    df_subset <- as.data.frame(current_transDT)[which(current_transDT[[cellID_coln]] == current_cell), transSpatLocs_coln]
    cellsize_range <- apply(df_subset, 2, FUN = range)
    # add in the distance defined by neighbor_distance for search_range
    search_range <- cellsize_range + c(- neighbor_distance_xy, neighbor_distance_xy)
    if(!is.null(neighbor_distance_z)){
      # different range for 3rd/ z dimension
      if(length(transSpatLocs_coln) >2){
        # 3D coordinates
        search_range[,transSpatLocs_coln[3]] <- cellsize_range[,transSpatLocs_coln[3]] + c(- neighbor_distance_z, neighbor_distance_z)
      }
    }
    
    # find all transcripts within this search_range
    # operate with values
    # # need import entire dplyr to use %between%
    # chosen_transcripts <- lapply(seq_len(ncol(search_range)),
    #                              function(i) current_transDT[get(transSpatLocs_coln[i]) %between% search_range[,i], get(transID_coln)])
    
    # without import entire dplyr
    chosen_transcripts <- lapply(seq_len(ncol(search_range)),
                                 function(i) current_transDT[dplyr::between(get(transSpatLocs_coln[i]), 
                                                                            search_range[1, i], search_range[2, i]), 
                                                             get(transID_coln)])
    
    chosen_transcripts <- Reduce(intersect, chosen_transcripts)
    df_subset <- current_transDT[which(current_transDT[[transID_coln]] %in% chosen_transcripts), ]
    # # operate with TRUE and FALSE, slower
    # chosen_transcripts <- lapply(seq_len(ncol(search_range)),
    #                              function(i) current_transDT[[transSpatLocs_coln[i]]] %between% search_range[, i])
    # df_subset <- current_transDT[Reduce('&', chosen_transcripts), ]
    
    return(df_subset)
  }
  
  
  # functions for each cell
  my_fun_eachCell <- function(query_df){
    each_cell <- query_df[[cellID_coln]][1]
    query_transID <- query_df[[transID_coln]] 
    
    if(is.null(neighbor_distance_xy) & !is.null(cell_networkDT)){
      # define neighborhood using cell spatial network
      eachCell_NT <- query_df[['network_cellID']][1]
      neighborCells_NT <- c(cell_networkDT[from == eachCell_NT]$to, cell_networkDT[to == eachCell_NT]$from)
      neighborCells_NT <- unique(c(eachCell_NT, neighborCells_NT))
      df_subset <- transcript_df[which(transcript_df[['network_cellID']] %in% neighborCells_NT), ]
    } else {
      # define neighborhood based on absolute distance to query cell
      # subset data based on gridSpat_coln
      if(!is.null(gridSpat_coln)){
        # # current_grid to locate directly
        # subTransDT <- transcript_df[gridAll_UID == query_df[['gridAll_UID']][1], ]
        
        # use grid index to locate 
        current_grid <- query_df[['gridAll_UID']][1]
        subTransDT <- transcript_df[which(transcript_df[['gridAll_UID']] == current_grid), ]
        
        df_subset <- my_function_searchDistance(current_cell = each_cell, current_transDT = subTransDT)
        
      } else {
        # no subset, search directly
        df_subset <- my_function_searchDistance(current_cell = each_cell, current_transDT = transcript_df)
      }
      
    }
    
    
    df_subset[['query_CellId']] <- rep(each_cell, nrow(df_subset))
    df_subset[['in_query_cell']] <- (df_subset[[cellID_coln]] == df_subset[['query_CellId']])
    
    return(df_subset)
    
  }
  
  # for each chosen cell, find out its neighborhood
  neighborhood_df <- by(chosen_transDF, chosen_transDF[[cellID_coln]], my_fun_eachCell)
  neighborhood_df <- do.call(rbind, neighborhood_df)
  return(neighborhood_df)
}