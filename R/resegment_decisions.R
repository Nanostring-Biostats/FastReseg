#' @title decide_ReSegment_Operations
#' @description Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations. 
#' @param neighborhood_df the data.frame containing neighborhood information for each query cells, expected to be output of neighborhood_for_resegment function.
#' @param selfcellID_coln the column name of cell_ID of query cell in neighborhood_df 
#' @param transNum_coln the column name of transcript number of query cell in neighborhood_df
#' @param selfCellType_coln the column name of cell_type under query cell in neighborhood_df 
#' @param selfScore_coln the column name of average transcript score under query cell in neighborhood_df 
#' @param neighborcellID_coln the column name of cell_ID of neighbor cell in neighborhood_df 
#' @param neighborCellType_coln the column name of cell_type under neighbor cell in neighborhood_df 
#' @param neighborScore_coln the column name of average transcript score under neighbor cell in neighborhood_df 
#' @param score_baseline a named vector of score baseline for all cell type listed in neighborhood_df
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @return a list 
#' \enumerate{
#'    \item{cells_to_discard, a vector of cell ID that should be discarded during resegmentation}
#'    \item{cells_to_update, a named vector of cell ID whether the cell_ID in name would be replaced with cell_ID in value.}
#'    \item{cells_to_keep, a vector of cell ID that should be kept as it is.}
#'    \item{reseg_full_converter, a single named vector of cell ID to update the original cell ID, assign NA for cells_to_discard.}
#' }
#' @details Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations.1) merge query to neighbor if consist cell type and fewer than average transcript number cutoff, higherCutoff_transNum; 2) keep query as new cell id if no consist neighbor cell type, but high self score and higher than minimal transcript number, lowerCutoff_transNum; 3) discard the rest of query cells that have no consistent neighbor cell type, fewer transcript number based on lowerCutoff_transNum, and/or low self score. Use network component analysis to resolve any conflic due to merging multiple query cells into one.  
#' @export
decide_ReSegment_Operations <- function(neighborhood_df,
                                        selfcellID_coln = 'CellId', 
                                        transNum_coln = 'transcript_num', 
                                        selfCellType_coln = 'self_celltype',
                                        selfScore_coln = 'score_under_self', 
                                        neighborcellID_coln = 'neighbor_CellId', 
                                        neighborCellType_coln = 'neighbor_celltype', 
                                        neighborScore_coln = 'score_under_neighbor',
                                        score_baseline = NULL, 
                                        lowerCutoff_transNum = NULL, 
                                        higherCutoff_transNum= NULL){
  
  # check format of transcript_df
  if(any(!c(selfcellID_coln, transNum_coln, selfCellType_coln, selfScore_coln, 
            neighborcellID_coln, neighborCellType_coln, neighborScore_coln) %in% colnames(neighborhood_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(selfcellID_coln, transNum_coln, selfCellType_coln, selfScore_coln, 
                                  neighborcellID_coln, neighborCellType_coln, neighborScore_coln), 
                                colnames(neighborhood_df)), collapse = "`, `")))
  }
  
  all_celltypes <- unique(c(neighborhood_df[[selfCellType_coln]], neighborhood_df[[neighborCellType_coln]]))
  
  # check the format for all cutoff, make sure it's number and has names for all cell types used.
  for(cutoff_var in c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")){
    if(is.null(get(cutoff_var))){
      stop(sprintf("The `%s` must be provided as a named numeric vector for score cutoff under each cell type used in neighborhood_df.", cutoff_var))
    } else {
      if(!any(class(get(cutoff_var)) %in% c('numeric'))){
        stop(sprintf("The provided `%s` must be a named numeric vector. ", cutoff_var))
      }
      if(length(setdiff(all_celltypes, names(get(cutoff_var))))>0){
        stop(sprintf("The provided `%s` is missing for the following cell types used in neighborhood_df: `%s`.", 
                     cutoff_var, paste0(setdiff(all_celltypes, names(get(cutoff_var))), collapse="`, `")))
      }
    }
    
  }
  
  # make sure lowerCutoff is lower than higherCutoff such that the merge would overwrite the keep. 
  lowerCutoff_transNum <- lowerCutoff_transNum[match(all_celltypes, names(lowerCutoff_transNum))]
  higherCutoff_transNum <- higherCutoff_transNum[match(all_celltypes, names(higherCutoff_transNum))]
  if(any(lowerCutoff_transNum > higherCutoff_transNum)){
    stop(sprintf("`lowerCutoff_transNum` is larger than `higherCutoff_transNum` in cell type: `%s`.", 
                 paste0(names(lowerCutoff_transNum)[lowerCutoff_transNum > higherCutoff_transNum], collapse = "`, `")))
  }
  
  
  # (1) whether consistent neighbor
  neighborhood_df[['consist_neighbor']] <- (neighborhood_df[[selfCellType_coln]] == neighborhood_df[[neighborCellType_coln]])
  
  # (2) whether high enough transcript number on self cell type for evaluate merging or keep as new cell
  neighborhood_df[['enough_to_keep']] <-( neighborhood_df[[transNum_coln]] > higherCutoff_transNum[neighborhood_df[[selfCellType_coln]]])
  neighborhood_df[['enough_to_merge']] <-( neighborhood_df[[transNum_coln]] > lowerCutoff_transNum[neighborhood_df[[selfCellType_coln]]])
  
  # (3) high enough score under self
  neighborhood_df[['high_selfscore']] <- (neighborhood_df[[selfScore_coln]] > score_baseline[neighborhood_df[[selfCellType_coln]]])
  
  # assign new cell id, NA to discard the cells with no consistent neighbor, few transcripts and low self score
  neighborhood_df[['corrected_cellID']] <- rep(NA, nrow(neighborhood_df))
  # merge if consist neighbors and fewer than average transcript number
  idx_to_merge <- which(neighborhood_df[['consist_neighbor']] & !(neighborhood_df[['enough_to_keep']]))
  # keep as new cell id if no consist neighbor, but high self score and higher than minimal transcript number
  idx_to_keep <- which(neighborhood_df[['high_selfscore']] & neighborhood_df[['enough_to_merge']])
  
  neighborhood_df[['corrected_cellID']][idx_to_keep] <- neighborhood_df[[selfcellID_coln]][idx_to_keep]
  neighborhood_df[['corrected_cellID']][idx_to_merge] <- neighborhood_df[[neighborcellID_coln]][idx_to_merge]
  
  #### clean up cells that need update and resolve circular referencing using network component analysis----
  cells_to_discard <- neighborhood_df[[selfcellID_coln]][is.na(neighborhood_df[['corrected_cellID']])]
  cells_to_update <- neighborhood_df[['corrected_cellID']][!(is.na(neighborhood_df[['corrected_cellID']]))]
  names(cells_to_update) <- neighborhood_df[[selfcellID_coln]][!(is.na(neighborhood_df[['corrected_cellID']]))]
  
  # remove the ones in both names and values
  cells_to_update <- cells_to_update[which(cells_to_update != names(cells_to_update))]
  # deal with pairs and groups, take advantage of the component network analysis to find membership between cells_to_update
  tmp_networkDF <- data.frame(from = names(cells_to_update), 
                              to = cells_to_update, 
                              distance = rep(1, length(cells_to_update)))
  
  all_index = unique(x = c(tmp_networkDF$from, tmp_networkDF$to))
  network_igraph = igraph::graph_from_data_frame(tmp_networkDF, directed = TRUE, vertices = all_index)
  group_vector <- igraph::components(network_igraph, mode = "weak")[['membership']]
  cells_to_update <- NULL
  for (groupID in unique(group_vector)){
    nodes <- names(group_vector)[which(group_vector == groupID)]
    # sort by iD
    nodes <- sort(nodes)
    current_vector <- rep(nodes[1], length(nodes)-1)
    names(current_vector) <- nodes[2: length(nodes)]
    cells_to_update <- c(cells_to_update, current_vector)
  }
  
  # update cells_to_discard, make sure not in the values or names of cells_to_udpate
  cells_to_discard <- cells_to_discard[which(!(cells_to_discard %in% cells_to_update))]
  cells_to_discard <- cells_to_discard[which(!(cells_to_discard %in% names(cells_to_update)))]
  
  # cells to keep as it is
  cells_to_keep <- setdiff(neighborhood_df[[selfcellID_coln]], c(cells_to_discard, names(cells_to_update)))
  
  # compile the single cellID converter for all operations, assign NA for cells_to_discard
  reseg_full_converter <- c(cells_to_update, cells_to_keep, rep(NA, length(cells_to_discard)))
  names(reseg_full_converter) <- c(names(cells_to_update), cells_to_keep, cells_to_discard)
  
  reseg_operations <- list(cells_to_discard = cells_to_discard, 
                           cells_to_update = cells_to_update, 
                           cells_to_keep = cells_to_keep, 
                           reseg_full_converter = reseg_full_converter)
  
  return(reseg_operations)
}

#' @title decide_ReSegment_Operations_leidenCut
#' @description Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations. 
#' @param neighborhood_df the data.frame containing neighborhood information for each query cells, expected to be output of neighborhood_for_resegment function.Use leiden clustering to determine whether a merge event is allowed.
#' @param selfcellID_coln the column name of cell_ID of query cell in neighborhood_df 
#' @param transNum_coln the column name of transcript number of query cell in neighborhood_df
#' @param selfCellType_coln the column name of cell_type under query cell in neighborhood_df 
#' @param selfScore_coln the column name of average transcript score under query cell in neighborhood_df 
#' @param neighborcellID_coln the column name of cell_ID of neighbor cell in neighborhood_df 
#' @param neighborCellType_coln the column name of cell_type under neighbor cell in neighborhood_df 
#' @param neighborScore_coln the column name of average transcript score under neighbor cell in neighborhood_df 
#' @param score_baseline a named vector of score baseline for all cell type listed in neighborhood_df
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param leiden_config a list of configuration to pass to reticulate and Giotto:::python_leiden function, including python path, resolution, partition_type, n_iterations, set_seed, seed_number. 
#' @param cutoff_sharedLeiden minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event, default = 0.5 for 50% cutoff
#' @return a list 
#' \enumerate{
#'    \item{cells_to_discard, a vector of cell ID that should be discarded during resegmentation}
#'    \item{cells_to_update, a named vector of cell ID whether the cell_ID in name would be replaced with cell_ID in value.}
#'    \item{cells_to_keep, a vector of cell ID that should be kept as it is.}
#'    \item{reseg_full_converter, a single named vector of cell ID to update the original cell ID, assign NA for cells_to_discard.}
#' }
#' @details Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations.1) merge query to neighbor if consist cell type and fewer than average transcript number cutoff, higherCutoff_transNum; 2) keep query as new cell id if no consist neighbor cell type, but high self score and higher than minimal transcript number, lowerCutoff_transNum; 3) discard the rest of query cells that have no consistent neighbor cell type, fewer transcript number based on lowerCutoff_transNum, and/or low self score. Use network component analysis to resolve any conflic due to merging multiple query cells into one. In case of merging into neighbor cell, leiden clustering on transcript level of spatial network is performed to decide whether the merge should be allowed. 
#' @export
decide_ReSegment_Operations_leidenCut <- function(neighborhood_df,
                                                  selfcellID_coln = 'CellId', 
                                                  transNum_coln = 'transcript_num', 
                                                  selfCellType_coln = 'self_celltype',
                                                  selfScore_coln = 'score_under_self', 
                                                  neighborcellID_coln = 'neighbor_CellId', 
                                                  neighborCellType_coln = 'neighbor_celltype', 
                                                  neighborScore_coln = 'score_under_neighbor',
                                                  score_baseline = NULL, 
                                                  lowerCutoff_transNum = NULL, 
                                                  higherCutoff_transNum= NULL, 
                                                  config_spatNW_transcript,
                                                  transcript_df, 
                                                  cellID_coln = "CellId",
                                                  transID_coln = "transcript_id",
                                                  transSpatLocs_coln = c('x','y','z'), 
                                                  leiden_config = list(python_path = "/usr/bin/python3", 
                                                                       partition_type = c("RBConfigurationVertexPartition", "ModularityVertexPartition"),
                                                                       resolution =1,
                                                                       n_iterations = 1000,
                                                                       set_seed = T,
                                                                       seed_number = 1234), 
                                                  cutoff_sharedLeiden = 0.5){
  
  
  # check format of neighborhood_df
  if(any(!c(selfcellID_coln, transNum_coln, selfCellType_coln, selfScore_coln, 
            neighborcellID_coln, neighborCellType_coln, neighborScore_coln) %in% colnames(neighborhood_df))){
    stop(sprintf("Not all necessary columns can be found in provided neighborhood_df, missing columns include `%s`.",
                 paste0(setdiff(c(selfcellID_coln, transNum_coln, selfCellType_coln, selfScore_coln, 
                                  neighborcellID_coln, neighborCellType_coln, neighborScore_coln), 
                                colnames(neighborhood_df)), collapse = "`, `")))
  }
  
  all_celltypes <- unique(c(neighborhood_df[[selfCellType_coln]], neighborhood_df[[neighborCellType_coln]]))
  
  # check the format for all cutoff, make sure it's number and has names for all cell types used.
  for(cutoff_var in c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")){
    if(is.null(get(cutoff_var))){
      stop(sprintf("The `%s` must be provided as a named numeric vector for score cutoff under each cell type used in neighborhood_df.", cutoff_var))
    } else {
      if(!any(class(get(cutoff_var)) %in% c('numeric'))){
        stop(sprintf("The provided `%s` must be a named numeric vector. ", cutoff_var))
      }
      if(length(setdiff(all_celltypes, names(get(cutoff_var))))>0){
        stop(sprintf("The provided `%s` is missing for the following cell types used in neighborhood_df: `%s`.", 
                     cutoff_var, paste0(setdiff(all_celltypes, names(get(cutoff_var))), collapse="`, `")))
      }
    }
    
  }
  
  # make sure lowerCutoff is lower than higherCutoff such that the merge would overwrite the keep. 
  lowerCutoff_transNum <- lowerCutoff_transNum[match(all_celltypes, names(lowerCutoff_transNum))]
  higherCutoff_transNum <- higherCutoff_transNum[match(all_celltypes, names(higherCutoff_transNum))]
  if(any(lowerCutoff_transNum > higherCutoff_transNum)){
    stop(sprintf("`lowerCutoff_transNum` is larger than `higherCutoff_transNum` in cell type: `%s`.", 
                 paste0(names(lowerCutoff_transNum)[lowerCutoff_transNum > higherCutoff_transNum], collapse = "`, `")))
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, transSpatLocs_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, transSpatLocs_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  
  d2_or_d3 <- length(transSpatLocs_coln)
  
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } else {
    message(sprintf("Run delanuay network in %d Dimension.", d2_or_d3))
  }
  
  # get transcript information for cells in neighborhood_df
  all_cellIDs <- unique(c(neighborhood_df[[selfcellID_coln]], neighborhood_df[[neighborcellID_coln]]))
  transcript_df <- as.data.frame(transcript_df)
  transcript_df <- transcript_df[transcript_df[[cellID_coln]] %in% all_cellIDs, ]
  
  # initialize leiden cluster configurations
  if(cutoff_sharedLeiden <=0 | cutoff_sharedLeiden >1){
    stop(sprintf("The providied cutoff_sharedLeiden = %.3f, must be within (0, 1].", cutoff_sharedLeiden))
  } else {
    message(sprintf("A valid merging event must have query cell with %.3f transcript shared same membership as neighbor cell of consistent cell type. ", cutoff_sharedLeiden))
  }
  
  if(is.null(leiden_config$python_path)){
    leiden_config$python_path <- "/usr/bin/python3"
  }
  if(!file.exists(leiden_config$python_path)){
    stop(sprintf("The provided python path does not exist, %s. ", leiden_config$python_path))
  } else {
    reticulate::use_python(required = T, python = leiden_config$python_path)
    python_leiden_function = system.file("python", "python_leiden.py", 
                                         package = "Giotto")
    reticulate::source_python(file = python_leiden_function)
  }
  if (leiden_config$set_seed == TRUE && !is.null(leiden_config$seed_number)) {
    leiden_config$seed_number = as.integer(leiden_config$seed_number)
  } else {
    leiden_config$seed_number = as.integer(sample(x = 1:10000, size = 1))
  }
  reticulate::py_set_seed(seed = leiden_config$seed_number, disable_hash_randomization = TRUE)
  
  if(is.null(leiden_config$partition_type)){
    leiden_config$partition_type <-  "RBConfigurationVertexPartition"
  }else {
    leiden_config$partition_type <- match.arg(leiden_config$partition_type, 
                                              choices = c("RBConfigurationVertexPartition", 
                                                          "ModularityVertexPartition"))
  }
  if(is.null(leiden_config$n_iterations)){
    leiden_config$n_iterations = 1000
  }
  if(is.null(leiden_config$resolution)){
    leiden_config$resolution = 1
  } else {
    if(leiden_config$resolution >1 | leiden_config$resolution <1e-10){
      stop(sprintf("The provided resolution for leiden_config = %.3f, which is outside of (0,1].", leiden_config$resolution))
    }
  }
  message(sprintf("Perform leiden clustering at resolution = %.3f.", leiden_config$resolution))
  
  
  # (1) whether consistent neighbor
  neighborhood_df[['consist_neighbor']] <- (neighborhood_df[[selfCellType_coln]] == neighborhood_df[[neighborCellType_coln]])
  
  # (2) whether high enough transcript number on self cell type for evaluate merging or keep as new cell
  neighborhood_df[['enough_to_keep']] <-( neighborhood_df[[transNum_coln]] > higherCutoff_transNum[neighborhood_df[[selfCellType_coln]]])
  neighborhood_df[['enough_to_merge']] <-( neighborhood_df[[transNum_coln]] > lowerCutoff_transNum[neighborhood_df[[selfCellType_coln]]])
  
  # (3) high enough score under self
  neighborhood_df[['high_selfscore']] <- (neighborhood_df[[selfScore_coln]] > score_baseline[neighborhood_df[[selfCellType_coln]]])
  
  # merge if consist neighbors and fewer than average transcript number
  idx_to_merge <- which(neighborhood_df[['consist_neighbor']] & !(neighborhood_df[['enough_to_keep']]))
  # keep as new cell id if no consist neighbor, but high self score and higher than minimal transcript number
  idx_to_keep <- which(neighborhood_df[['high_selfscore']] & neighborhood_df[['enough_to_merge']])
  
  # (4) for cells that might be merged, check spatial network to decide if merging is allowed by leiden clustering
  
  # for each pair of cells that would be merged
  myfun_NWclustering <- function(neighborDF_eachCell){
    query_cellID <- neighborDF_eachCell[[selfcellID_coln]][1]
    neighbor_cellID <- neighborDF_eachCell[[neighborcellID_coln]][1]
    outputs <- data.frame(query_cellID = query_cellID, 
                          neighbor_cellID = neighbor_cellID, 
                          origin_Idx = neighborDF_eachCell[['origin_Idx']][1])
    
    merged_pair <- unique(c(query_cellID, neighbor_cellID))
    if (length(merged_pair) !=2) {
      outputs[['merge']] <- FALSE
    } else {
      df_subset <- transcript_df[which(transcript_df[[cellID_coln]] %in% merged_pair), ]
      
      # find connected neighborhood environment, be careful for the case with too few transcripts for delaunay 
      if(nrow(df_subset) < 4){
        # too few query and neighbor transcripts for network analysis and cell typing
        # treat as no valid direct neighbor
        outputs[['merge']] <- FALSE
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
            message(sprintf("(`%s`) cell pair with all %d transcripts in same z plane, run 2D network analysis.", 
                            paste0(merged_pair, collapse = "``, `"), nrow(df_subset)))
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
        
        # if failed delanuay network, keep cell pair separated without merging
        if(is.null(delaunayNW_Obj)){
          outputs[['merge']] <- FALSE
        }else{
          # run leiden clustering on igraph
          network_edge_dt <- delaunayNW_Obj$networkDT[, c("from", "to", "weight"), with = F]
          # assign init membership based on source cell_ID, query cell = 2, neighbor cell =1
          init_membership <- as.integer(df_subset[[cellID_coln]] == query_cellID) +1
          names(init_membership) <- df_subset[[transID_coln]]
          
          # for unknown reason init_membership does not always have consistent length with leiden cluster
          # use NULL instead
          
          pyth_leid_result = python_leiden(df = network_edge_dt, partition_type = leiden_config$partition_type, 
                                           initial_membership = NULL, weights = "weight", 
                                           n_iterations = leiden_config$n_iterations, 
                                           seed = leiden_config$seed_number, 
                                           resolution_parameter = leiden_config$resolution)
          
          ident_clusters_DT = data.table::data.table(transcript_id = pyth_leid_result[[1]], 
                                                     clusterID = pyth_leid_result[[2]])
          # check if how many transcripts in query cell shared same clusterID as neighbor cell
          ident_clusters_DT[['init_cluster']] <- init_membership[ident_clusters_DT[['transcript_id']]]
          
          # use orignal input for total transcript since some transcripts would be missing in spatial network
          n_query <- sum(df_subset[[cellID_coln]] == query_cellID)
          lc_in_neighbor <- unique(ident_clusters_DT[['clusterID']][which(ident_clusters_DT[['init_cluster']] ==1)])
          n_sharedLC <- sum(ident_clusters_DT[['clusterID']][which(ident_clusters_DT[['init_cluster']] ==2)] %in% lc_in_neighbor)
          
          if(n_sharedLC/n_query < cutoff_sharedLeiden){
            outputs[['merge']] <- FALSE
          } else {
            outputs[['merge']] <- TRUE
          } 
        }
        
        
        
        
      }
      
    }
    return(outputs)
  }
  
  neighDF_for_mergeCheck <- neighborhood_df[idx_to_merge, ]
  neighDF_for_mergeCheck[['origin_Idx']] <- idx_to_merge
  
  message(sprintf("Perform ledien clustering on %d potential merging events. ", nrow(neighDF_for_mergeCheck)))
  
  mergeCheck_res <- by(neighDF_for_mergeCheck, neighDF_for_mergeCheck[[selfcellID_coln]], myfun_NWclustering)
  mergeCheck_res <- do.call(rbind, mergeCheck_res)
  
  # update the idx_to_merge after merge checking with leiden clustering
  idx_to_merge <- mergeCheck_res[['origin_Idx']][which(mergeCheck_res[['merge']] == TRUE)]
  
  # assign new cell id, NA to discard the cells with no consistent neighbor, few transcripts and low self score
  neighborhood_df[['corrected_cellID']] <- rep(NA, nrow(neighborhood_df))
  neighborhood_df[['corrected_cellID']][idx_to_keep] <- neighborhood_df[[selfcellID_coln]][idx_to_keep]
  neighborhood_df[['corrected_cellID']][idx_to_merge] <- neighborhood_df[[neighborcellID_coln]][idx_to_merge]
  
  #### clean up cells that need update and resolve circular referencing using network component analysis----
  cells_to_discard <- neighborhood_df[[selfcellID_coln]][is.na(neighborhood_df[['corrected_cellID']])]
  cells_to_update <- neighborhood_df[['corrected_cellID']][!(is.na(neighborhood_df[['corrected_cellID']]))]
  names(cells_to_update) <- neighborhood_df[[selfcellID_coln]][!(is.na(neighborhood_df[['corrected_cellID']]))]
  
  # remove the ones in both names and values
  cells_to_update <- cells_to_update[which(cells_to_update != names(cells_to_update))]
  # deal with pairs and groups, take advantage of the component network analysis to find membership between cells_to_update
  tmp_networkDF <- data.frame(from = names(cells_to_update), 
                              to = cells_to_update, 
                              distance = rep(1, length(cells_to_update)))
  
  all_index = unique(x = c(tmp_networkDF$from, tmp_networkDF$to))
  network_igraph = igraph::graph_from_data_frame(tmp_networkDF, directed = TRUE, vertices = all_index)
  group_vector <- igraph::components(network_igraph, mode = "weak")[['membership']]
  cells_to_update <- NULL
  for (groupID in unique(group_vector)){
    nodes <- names(group_vector)[which(group_vector == groupID)]
    # sort by iD
    nodes <- sort(nodes)
    current_vector <- rep(nodes[1], length(nodes)-1)
    names(current_vector) <- nodes[2: length(nodes)]
    cells_to_update <- c(cells_to_update, current_vector)
  }
  
  # update cells_to_discard, make sure not in the values or names of cells_to_udpate
  cells_to_discard <- cells_to_discard[which(!(cells_to_discard %in% cells_to_update))]
  cells_to_discard <- cells_to_discard[which(!(cells_to_discard %in% names(cells_to_update)))]
  
  # cells to keep as it is
  cells_to_keep <- setdiff(neighborhood_df[[selfcellID_coln]], c(cells_to_discard, names(cells_to_update)))
  
  # compile the single cellID converter for all operations, assign NA for cells_to_discard
  reseg_full_converter <- c(cells_to_update, cells_to_keep, rep(NA, length(cells_to_discard)))
  names(reseg_full_converter) <- c(names(cells_to_update), cells_to_keep, cells_to_discard)
  
  reseg_operations <- list(cells_to_discard = cells_to_discard, 
                           cells_to_update = cells_to_update, 
                           cells_to_keep = cells_to_keep, 
                           reseg_full_converter = reseg_full_converter)
  
  return(reseg_operations)
}
