#' @title decide_ReSegment_Operations
#' @description Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations. Use either leiden clustering or geometry statistics to determine whether a merge event is allowed.
#' @param neighborhood_df the data.frame containing neighborhood information for each query cells, expected to be output of get_neighborhood_content function.
#' @param selfcellID_coln the column name of cell_ID of query cell in neighborhood_df 
#' @param transNum_coln the column name of transcript number of query cell in neighborhood_df
#' @param selfCellType_coln the column name of cell_type under query cell in neighborhood_df 
#' @param selfScore_coln the column name of average transcript score under query cell in neighborhood_df 
#' @param neighborcellID_coln the column name of cell_ID of neighbor cell in neighborhood_df 
#' @param neighborCellType_coln the column name of cell_type under neighbor cell in neighborhood_df 
#' @param neighborScore_coln the column name of average transcript score under neighbor cell in neighborhood_df 
#' @param score_baseline a named vector of score baseline for all cell type listed in neighborhood_df such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence 
#' @param lowerCutoff_transNum a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is
#' @param higherCutoff_transNum a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param spatialMergeCheck_method use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space (default = "leidenCut")
#' @param cutoff_spatialMerge spatial constraint on a valid merging event between two source transcript groups, default = 0.5 for 50% cutoff, set to 0 to skip spatial constraint evaluation for merging.   
#' For `spatialMergeCheck_method = "leidenCut"`, this is the minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event.
#' For `spatialMergeCheck_method = "geometryDiff"`, this is the maximum percentage of white space change upon merging of query cell and neighbor cell for a valid merging event. 
#' @param leiden_config (leidenCut) a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.  
#' @param config_spatNW_transcript (leidenCut) configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details, set to NULL to use default config
#' @return a list 
#' \enumerate{
#'    \item{cells_to_discard, a vector of cell ID that should be discarded during resegmentation}
#'    \item{cells_to_update, a named vector of cell ID where the cell_ID in name would be replaced with cell_ID in value.}
#'    \item{cells_to_keep, a vector of cell ID that should be kept as it is.}
#'    \item{reseg_full_converter, a single named vector of cell ID to update the original cell ID, assign NA for cells_to_discard.}
#' }
#' @details Evaluate neighborhood information against score and transcript number cutoff to decide the resegmetation operations like the following: 
#'    * merge query to neighbor if consist cell type and fewer than average transcript number cutoff, higherCutoff_transNum; 
#'    * keep query as new cell id if no consist neighbor cell type, but high self score and higher than minimal transcript number, lowerCutoff_transNum; 
#'    * discard the rest of query cells that have no consistent neighbor cell type, fewer transcript number based on lowerCutoff_transNum, and/or low self score. 
#' The function uses network component analysis to resolve any conflict due to merging multiple query cells into one. 
#' When `cutoff_spatialMerge > 0`, the function applies additional spatial constraint on a valid merging event of query cell into neighbor cell. 
#'    * In case of `spatialMergeCheck_method = "leidenCut"`, the function builds spatial network at transcript level, does leiden clustering  on the spatial network, and then decides whether the merge should be allowed based on the observed shared leiden membership of the two source transcript groups for a putative merging event; the provided `cutoff_spatialMerge` gives the minimal values of shared leiden memberhsip for a valid merging event. 
#'    * In case of `spatialMergeCheck_method = "geometryDiff"`, the function would first calculate white space, i.e. the area difference between convex and concave hulls, respectively, for query cell, neighbor cell, and the corresponding merged cell; and then calculate the white space difference between the merged cell and two separate cells and normalize that value with respect to the concave area of query and neighbor cells, respectively; lastly, allow a valid merging when the normalized white space difference upon merging for both query and neighbor cells are smaller than the provided `cutoff_spatialMerge`.
#' @importFrom igraph cluster_leiden graph_from_data_frame membership components
#' @importFrom spatstat.geom area.owin convexhull.xy
#' @importFrom concaveman concaveman
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
                                        higherCutoff_transNum= NULL, 
                                        transcript_df, 
                                        cellID_coln = "CellId",
                                        transID_coln = "transcript_id",
                                        transSpatLocs_coln = c('x','y','z'), 
                                        spatialMergeCheck_method = c("leidenCut", "geometryDiff"), 
                                        cutoff_spatialMerge = 0.5, 
                                        leiden_config = list(objective_function = c("CPM", "modularity"),
                                                             resolution_parameter = 1,
                                                             beta = 0.01,
                                                             n_iterations = 200),
                                        config_spatNW_transcript = NULL){
  
  d2_or_d3 <- length(transSpatLocs_coln)
  if(!(d2_or_d3 %in% c(2,3))){
    stop("spatLocs_colns must be the column names for 1st, 2nd, optional 3rd dimension of spatial coordinates in transcript_df.")
  } 
  
  if(cutoff_spatialMerge <0 | cutoff_spatialMerge >1){
    stop(sprintf("The providied `cutoff_spatialMerge = %.3f`, must be within [0, 1].", cutoff_spatialMerge))
  } else if(cutoff_spatialMerge ==0){
    message(sprintf('The provided `cutoff_spatialMerge = 0`, no spatial evaluation would be done for potential cell merging events. '))
  }else {
    # method of spatial constraint on merging
    spatialMergeCheck_method <- match.arg(spatialMergeCheck_method, c("leidenCut", "geometryDiff"))
    message(sprintf('Use `%s` method to evaluate putative merging event in space. ', spatialMergeCheck_method))
    
    if(spatialMergeCheck_method == "leidenCut"){
      # check and set config for spatial network and leiden clustering 
      leiden_config <- check_config_leiden(leiden_config)
      config_spatNW_transcript <- check_config_spatialNW(config = config_spatNW_transcript, 
                                                         spat_locs = transSpatLocs_coln)
      
      message(sprintf("A valid merging event must have query cell with %.3f transcript shared same membership as neighbor cell of consistent cell type. ", cutoff_spatialMerge))
      message(sprintf("Run delanuay network in %d Dimension.", d2_or_d3))
      
    }else {
      message(sprintf("A valid merging event to neighbor cell of consistent cell type must have no more than %.3f area change in white space upon merging with respect to the concave area of either source cells. ", cutoff_spatialMerge))
      message(sprintf("Perform geometry analysis in 2D for potential merging events despite the provided %d Dimension data.", d2_or_d3))
    }
  }
  
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
  
  
  # get transcript information for cells in neighborhood_df
  all_cellIDs <- unique(c(neighborhood_df[[selfcellID_coln]], neighborhood_df[[neighborcellID_coln]]))
  transcript_df <- as.data.frame(transcript_df)
  transcript_df <- transcript_df[transcript_df[[cellID_coln]] %in% all_cellIDs, ]

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

  # check merging interface if any merging candidate and cutoff >0 
  if(length(idx_to_merge) >0 & cutoff_spatialMerge >0){
    neighDF_for_mergeCheck <- neighborhood_df[idx_to_merge, ]
    neighDF_for_mergeCheck[['origin_Idx']] <- idx_to_merge
    
    if(spatialMergeCheck_method == "leidenCut"){
    message(sprintf("Perform ledien clustering on %d potential merging events. ", length(idx_to_merge)))
    
    # for each pair of cells that would be merged, build spatial network and do leiden clustering
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
                              paste0(merged_pair, collapse = "`, `"), nrow(df_subset)))
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
            network_edge_dt <- delaunayNW_Obj@networkDT[, c("from", "to", "weight"), with = F]
            # assign init membership based on source cell_ID, query cell = 2, neighbor cell =1
            init_membership <- as.integer(df_subset[[cellID_coln]] == query_cellID) +1
            names(init_membership) <- df_subset[[transID_coln]]
            
            # undirected graph needed for igraph:::cluster_leiden, igraph >= 1.2.7
            all_index = unique(c(network_edge_dt$from, network_edge_dt$to))
            network_igraph = igraph::graph_from_data_frame(as.data.frame(network_edge_dt), directed = FALSE, vertices = all_index)
            # igraph::community object
            leid_result <- igraph:::cluster_leiden(network_igraph, 
                                                   objective_function = leiden_config$objective_function,
                                                   resolution_parameter = leiden_config$resolution_parameter,
                                                   beta = leiden_config$ beta,
                                                   initial_membership = NULL,
                                                   n_iterations = leiden_config$n_iterations)
            # igraph::membership object
            leid_result <- igraph::membership(leid_result)
            
            ident_clusters_DF <- data.frame(transcript_id = names(leid_result), 
                                            clusterID = as.vector(leid_result))
            
            # check if how many transcripts in query cell shared same clusterID as neighbor cell
            ident_clusters_DF[['init_cluster']] <- init_membership[ident_clusters_DF[['transcript_id']]]
            
            # use original input for total transcript since some transcripts would be missing in spatial network
            n_query <- sum(df_subset[[cellID_coln]] == query_cellID)
            lc_in_neighbor <- unique(ident_clusters_DF[['clusterID']][which(ident_clusters_DF[['init_cluster']] ==1)])
            n_sharedLC <- sum(ident_clusters_DF[['clusterID']][which(ident_clusters_DF[['init_cluster']] ==2)] %in% lc_in_neighbor)
            
            if(n_sharedLC/n_query < cutoff_spatialMerge){
              outputs[['merge']] <- FALSE
            } else {
              outputs[['merge']] <- TRUE
            } 
          }
          
          
        }
        
      }
      return(outputs)
    }
    
    
    } else if(spatialMergeCheck_method == "geometryDiff"){
      message(sprintf("Perform geometry analysis on %d potential merging events. ", length(idx_to_merge)))
      
      # for each pair of cells that would be merged, do geometry analysis 
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
          queryDF <- as.data.frame(transcript_df[which(transcript_df[[cellID_coln]] == query_cellID), ])
          neighDF <- as.data.frame(transcript_df[which(transcript_df[[cellID_coln]] == neighbor_cellID), ])
          
          # if either cell has only 1~2 transcripts of unique coordinates, merge as it is
          flag_geometry <- sapply(list(query = queryDF, 
                                       neighbor = neighDF), 
                                  function(eachDF){
                                    nrow(unique(eachDF[, transSpatLocs_coln[1:2]])) <3
                                  })
          
          if(any(flag_geometry)){
            outputs[['merge']] <- TRUE
          } else {
            # geometry analysis on 2D
            areaRes <- sapply(list(query = queryDF, 
                                   neighbor = neighDF, 
                                   merged = rbind(queryDF, neighDF)),
                              function(eachDF){
                                # convex hull area
                                convexArea <- spatstat.geom::area.owin(spatstat.geom::convexhull.xy(eachDF[, transSpatLocs_coln[1:2]]))
                                
                                # concave hull area
                                polyClockwise <- concaveman::concaveman(as.matrix(eachDF[, transSpatLocs_coln[1:2]]))
                                # need anti-clockwise polygonal as external boundary 
                                concaveArea <- spatstat.geom::area.owin(spatstat.geom::owin(poly=list(x=polyClockwise[nrow(polyClockwise):1, 1],y=polyClockwise[nrow(polyClockwise):1, 2])))
                                
                                areaStats <- c(convexArea, concaveArea, convexArea - concaveArea)
                                names(areaStats) <- c('convex', 'concave', 'whitespace')
                                
                                return(areaStats)
                              })
            
            # difference in white space due to merged, negative value means merged cell has high convexity 
            diffWB <- areaRes['whitespace', 'merged'] - areaRes['whitespace', 'query'] - areaRes['whitespace', 'neighbor']
            
            # normalized by the the concave area of either source cell, check if maximum normlaized value less than cutoff
            if(max(diffWB/areaRes['concave', c('query', 'neighbor')]) < cutoff_spatialMerge){
              outputs[['merge']] <- TRUE
            } else {
              outputs[['merge']] <- FALSE
            } 
            
          }
          
          
        }
        return(outputs)
      }
      
      
    }
    

    mergeCheck_res <- by(neighDF_for_mergeCheck, neighDF_for_mergeCheck[[selfcellID_coln]], myfun_NWclustering)
    mergeCheck_res <- do.call(rbind, mergeCheck_res)
    
    # update the idx_to_merge after merge checking with leiden clustering
    idx_to_merge <- mergeCheck_res[['origin_Idx']][which(mergeCheck_res[['merge']] == TRUE)]
    
  }  
  
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


#' @title check_config_leiden
#' @description check config used for leiden clustering, assign default values if missing arguments
#' @param config a list of config would be used to create spatial network using delaunay method only via GiottoClass::createSpatialNetwork
#' @return the corrected config list
check_config_leiden <- function(config){
  if(is.null(config$objective_function)){
    config$objective_function <-  "CPM"
  }else {
    config$objective_function <- match.arg(config$objective_function, 
                                           choices = c("CPM", "modularity"))
  }
  
  # positive integer
  if(is.null(config$n_iterations)){
    config$n_iterations = 200
  }else{
    msg <- checkTypeLengthValue(config, "n_iterations", 
                                expect_type = c("numeric","integer"), 
                                expect_len = 1, expect_range = "larger", 
                                expect_value = 0)
    if(length(msg) > 0L){
      stop( "Configuration Issues:\n" , paste( msg , collapse = "\n" ) )
    }
  }
  
  # Parameter affecting the randomness in the Leiden algorithm. This affects only the refinement step of the algorithm.
  if(is.null(config$beta)){
    config$beta = 0.01
  } else {
    if(config$beta >1 | config$beta <1e-10){
      stop(sprintf("The provided `beta` for leiden_config = %.3f, which is outside of (0,1].", config$beta))
    }
  }
  
  if(is.null(config$resolution_parameter)){
    config$resolution_parameter = 1
  } else {
    if(config$resolution_parameter >1 | config$resolution_parameter <1e-10){
      stop(sprintf("The provided `resolution_parameter` for leiden_config = %.3f, which is outside of (0,1].", config$resolution_parameter))
    }
  }
  message(sprintf("Perform leiden clustering at resolution_parameter = %.3f.", config$resolution_parameter))
  
  return(config)
}



