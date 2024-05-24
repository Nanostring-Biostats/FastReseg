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
#'    \item{transSpatLocs_coln, spatial coordinates of transcript}
#'    \item{transcript_group, group of chosen_transcripts}
#' } 
#' @details For query cell, group flagged transcripts only based on their molecular distance to each other. When distance cutoff = 'auto', use 20% average XY cell range as cutoff. In case of no more than 3 flagged transcripts per cell, determine the grouping based on distance cutoff directly. In case of more transcripts per cell, use \code{dbscan} to group transcripts with distance_cutoff as `eps` and `minPts = 1`.  
#' @importFrom dbscan dbscan
#' @importFrom data.table .SD
#' @export
groupTranscripts_dbscan <- function(chosen_transcripts = NULL, 
                                    distance_cutoff = 'auto',
                                    transcript_df, 
                                    cellID_coln = "CellId", transID_coln = "transcript_id",
                                    transSpatLocs_coln = c('x','y','z')){
  
  # `distance_cutoff` is for orphan transcript determination and overall grouping with `dbscan`
  if(distance_cutoff =='auto'){
    message("Automatic cutoff as 20% diameter of query cell.")
  } else if (is.numeric(distance_cutoff)){
    if (distance_cutoff <= 0){
      stop("`distance_cutoff` must be either `auto` or a positive number")
    }
  }else {
    stop("`distance_cutoff` must be either `auto` or a positive number")
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



#' @title groupTranscripts_Delaunay
#' @description group the flagged transcript within each cell based on spatial connectivity of their transcript delaunay network
#' @param chosen_transcripts the transcript_id of chosen transcript
#' @param config_spatNW_transcript configuration list to create spatial network at transcript level, see manual for \code{createSpatialDelaunayNW_from_spatLocs} for more details
#' @param distance_cutoff maximum distance within connected transcript group (default = "auto")
#' @param transcript_df the data.frame with transcript_id, target/geneName, x, y and cell_id
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @return data frame of connected transcripts among chosen_transcripts
#' #' \enumerate{
#'    \item{cellID_coln, orignal cell_ID}
#'    \item{transID_coln, connected transcripts among chosen_transcripts}
#'    \item{transSpatLocs_coln, spatial coordinates of transcript}
#'    \item{transcript_group, group of chosen_transcripts}
#' } 
#' @details for query cell, build network on flagged transcripts only to identify groups. In case of no more than 3 transcripts, determine the grouping based on distance cutoff directly; when distance cutoff = 'auto', no additional edge filtering based on delaunay network output but use 20% average XY cell range as cutoff when no more than 3 transcript.  
#' @importFrom data.table ':=' .N .SD
#' @export
groupTranscripts_Delaunay <- function(chosen_transcripts = NULL, 
                                      config_spatNW_transcript = list(name = 'transcript_delaunay_network',
                                                                      dimensions = "all",
                                                                      method = 'Delaunay',
                                                                      minimum_k = 0,
                                                                      delaunay_method = "delaunayn_geometry",
                                                                      maximum_distance_delaunay = "auto",
                                                                      options = "Pp",
                                                                      Y = TRUE,
                                                                      j = TRUE,
                                                                      S = 0),
                                      distance_cutoff = "auto",
                                      transcript_df, 
                                      cellID_coln = "CellId", 
                                      transID_coln = "transcript_id",
                                      transSpatLocs_coln = c('x','y','z')){
  
  # `distance_cutoff` is for orphan transcript determination, not for spatial network building which is defined in `config_spatNW_transcript` separately
  if(distance_cutoff =='auto'){
    message("automatic cutoff for delauny network.")
  } else if (is.numeric(distance_cutoff)){
    if (distance_cutoff <= 0){
      stop("`distance_cutoff` must be either `auto` or a positive number")
    }
  }else {
    stop("`distance_cutoff` must be either `auto` or a positive number")
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
  # (1) get number of flagged transcripts and distance cutoff which is used for orphan transcript determination but not in spatial network building
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
    transcript_df[['transcript_group']][tmp_idx] <- group_converter[as.character(transcript_df[[transID_coln]][tmp_idx])]
  }
  
  # (5) assign group for multiple flagged transcript cases using delaunay network analysis
  # function for each cell based on delaunay
  my_fun_delaunay <- function(df_subset){
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
        transNetWorkDT <- delaunayNW_Obj@networkDT
        
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
          names(solo_converter) <- as.character(solo_transcripts[['tmp_transID']])
          dfCoord_subset[['transcript_group']][is.na(dfCoord_subset[['transcript_group']])] <- solo_converter[as.character(dfCoord_subset[['tmp_transID']][is.na(dfCoord_subset[['transcript_group']])])]
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
    tmp_group <- by(chosenTrans_df4, chosenTrans_df4[[cellID_coln]], my_fun_delaunay)
    tmp_group <- do.call(rbind, tmp_group)
    group_converter <- tmp_group[['transcript_group']] 
    names(group_converter) <- tmp_group[[transID_coln]]
    tmp_idx <- which(transcript_df[[transID_coln]] %in% tmp_group[[transID_coln]])
    transcript_df[['transcript_group']][tmp_idx] <- group_converter[as.character(transcript_df[[transID_coln]][tmp_idx])]
  }
  
  return(transcript_df)
}


#' @title myFun_3point_singleCell
#' @description supporting function for \code{groupTranscripts_Delaunay}, assign group ID for 3 transcripts in single cell in 3D based on distant cutoff
#' @param dfCoord_subset transcript data.table for single cell with only 3 transcripts in rows
#' @param transSpatLocs_coln the column name of 1st, 2nd, optional 3rd spatial dimension of each transcript in transcript_df
#' @param distance_cutoff maximum distance within connected transcript group 
#' @param startGroup the index of starting group ID
#' @return a data.table with `transcript_group` column added to original input data.table
#' @importFrom data.table as.data.table setDT .SD ':='
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


