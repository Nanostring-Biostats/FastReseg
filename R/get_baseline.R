#' @title choose_distance_cutoff
#' @description Choose appropriate cellular distance cutoff and molecular distance cutoff based on input transcript data.frame for downstream resegmentation; cellular distance cutoff is defined as the search radius of direct neighbor cell, while molecular distance cutoff is defined as the maximum distance between two neighbor transcripts from same source cells.
#' @param transcript_df the data.frame for each transcript
#' @param transID_coln the column name of transcript_ID in `transcript_df`
#' @param cellID_coln the column name of cell_ID in `transcript_df`
#' @param spatLocs_colns column names for 1st, 2nd and optional 3rd dimension of spatial coordinates in `transcript_df` 
#' @param extracellular_cellID a vector of cell_ID for extracellular transcripts which would be removed from the resegmention pipeline (default = NULL)
#' @param sampleSize_nROI number of ROIs randomly picked from data for molecular distance cutoff estimation
#' @param sampleSize_cellNum maximum number of cells from the picked ROIs for molecular distance cutoff estimation
#' @param seed a random seed for sub-sampling cells from whole dataset for molecular distance cutoff estimation
#' @param run_molecularDist flag to run molecular distant cutoff estimation, default = TRUE 
#' @return a list
#' \describe{
#'    \item{cellular_distance_cutoff}{maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate. }
#'    \item{perCell_coordDT}{a data.table with cell in row, spatial XY coordiantes of centroid and dimensions of bounding box in column}
#'    \item{molecular_distance_cutoff}{maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate; return if run_molecularDist = TRUE}
#'    \item{distance_profile}{a named vector for the quantile profile of minimal molecular distance between transcripts belong to different cells at step size of 10% quantile; return if run_molecularDist = TRUE}
#' }
#' @details `cellular_distance_cutoff` is defined as maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact. The function calculates average 2D cell diameter from input data.frame and use 2 times of the mean cell diameter as `cellular_distance_cutoff`. `molecular_distance_cutoff` is defined as maximum molecule-to-molecule distance within connected transcript groups belonging to same source cells. The function would first randomly choose `sampleSize_cellNum` number of cells from `sampleSize_nROI` number of randomly picked ROIs with search radius to be 5 times of `cellular_distance_cutoff`, and then calculate the minimal molecular distance between picked cells. The function would further use the 5 times of 90% quantile of minimal molecular distance as `molecular_distance_cutoff`. This calculation is slow and is not recommended for large transcript data.frame.
#' @examples 
#' data(mini_transcriptDF)
#' extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID'] # cell_ID for extracellualr transcripts
#' distCutoffs <- choose_distance_cutoff(mini_transcriptDF,
#'                                       extracellular_cellID = extracellular_cellID)
#' @importFrom data.table as.data.table
#' @importFrom dplyr between
#' @importFrom spatstat.geom ppp pp3 box3 nndist
#' @export
choose_distance_cutoff <- function(transcript_df, 
                                   transID_coln = "transcript_id",
                                   cellID_coln = 'cell_ID', 
                                   spatLocs_colns = c('x','y','z'), 
                                   extracellular_cellID = NULL, 
                                   sampleSize_nROI = 10, 
                                   sampleSize_cellNum = 2500, 
                                   seed = 123, 
                                   run_molecularDist = TRUE){
  
  #### check inputs ----
  # check format of transcript_df
  if(any(!c(transID_coln, spatLocs_colns, cellID_coln) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided `transcript_df`, missing columns include `%s`.",
                 paste0(setdiff(c(transID_coln, spatLocs_colns, cellID_coln), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  # extract the needed columns only
  transcript_df <- as.data.frame(transcript_df)
  transcript_df <- transcript_df[, c(transID_coln, spatLocs_colns, cellID_coln)]
  
  ## remove extracellular transcript from transcript_df ----
  if(!is.null(extracellular_cellID)){
    if(length(extracellular_cellID)>1){
      common_cells <- unique(transcript_df[[cellID_coln]])
      common_cells <- setdiff(common_cells, extracellular_cellID)
      transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% common_cells), ]
      
      message(sprintf("%d transcript records within %d cells remain after removal of extracellular transcripts within %d cell_ID in provided `extracellular_cellID` vector.", 
                      nrow(transcript_df), length(common_cells), length(extracellular_cellID)))
    }
    
  }

 
  ## get cellular distance cutoff ---
  # get per cell data.table for centroid spatial coordinates and geometry information of bounding box
  perCell_coordM <- data.table::as.data.table(transcript_df)[, list(CenterX = mean(get(spatLocs_colns[1])), 
                                                                    CenterY = mean(get(spatLocs_colns[2])),
                                                                    Width = diff(range(get(spatLocs_colns[1]))),
                                                                    Height = diff(range(get(spatLocs_colns[2])))), 
                                                             by = get(cellID_coln)]
  colnames(perCell_coordM)[1] <- cellID_coln
  
  # Use 2 times of average 2D cell diameter as cellular_distance_cutoff for searching of neighbor cells.
  cellular_distance_cutoff <- max(colMeans(perCell_coordM[, c('Width','Height')]))*2
  message(sprintf("Use 2 times of average 2D cell diameter as cellular_distance_cutoff = %.4f for searching of neighbor cells.", cellular_distance_cutoff))
  
  # final results list
  res <- list(cellular_distance_cutoff = cellular_distance_cutoff, 
              perCell_coordDT = perCell_coordM)
  
  
  ## get molecular distance cutoff ---
  if(run_molecularDist){
    if(nrow(perCell_coordM)> sampleSize_cellNum){
      # random subset ROIs from whole transcript_df
      set.seed(seed)
      
      radius_ROIs = 5*cellular_distance_cutoff
      spatial_locs <- as.matrix(perCell_coordM[, -1])
      rownames(spatial_locs) <- perCell_coordM[[cellID_coln]]
      centroid_colns <- c("CenterX", "CenterY")
      
      # indices of cells to keep:
      cells_to_keep = c()
      for (i in seq_len(sampleSize_nROI)) {
        # cells that haven't been picked yet:
        available_cell_indices <- setdiff(seq_along(perCell_coordM[[cellID_coln]]), cells_to_keep)
        if(length(available_cell_indices) ==0) break
        
        # pick a center cell:
        randomStart <- sample(available_cell_indices, 1)
        center <- as.vector(t(spatial_locs[randomStart, centroid_colns]))
        search_range <- matrix(c(center - radius_ROIs, center + radius_ROIs), 
                               nrow = 2, byrow = TRUE, dimnames = list(c(), centroid_colns))
        # find cell within bounding box of radius_ROIs in length
        closest_cells <- lapply(seq_len(ncol(search_range)),
                                function(i) perCell_coordM[dplyr::between(get(centroid_colns[i]), 
                                                                          search_range[1, i], search_range[2, i]), 
                                                           which = TRUE])
        
        closest_cells <- Reduce(intersect, closest_cells)
        cells_to_keep <- unique(c(cells_to_keep, closest_cells))
        
        if(length(cells_to_keep)>sampleSize_cellNum) break
      }
      cells_to_keep <- perCell_coordM[[cellID_coln]][cells_to_keep]
      
    } else {
      cells_to_keep <- perCell_coordM[[cellID_coln]]
    }
    
    cutoff_transDF <- transcript_df[which(transcript_df[[cellID_coln]] %in% cells_to_keep), ]
    # drop the dimension without variance in coordinates
    spatLocs_to_use <- spatLocs_colns[apply(cutoff_transDF[, spatLocs_colns], 
                                            2, function(x) diff(range(x)))>0]
    message(sprintf("Identified %dD coordinates with variance. ", length(spatLocs_to_use)))
    
    if(length(spatLocs_to_use) ==2){
      queryTrans_pp <- spatstat.geom::ppp(x = cutoff_transDF[[spatLocs_to_use[1]]], 
                                          y = cutoff_transDF[[spatLocs_to_use[2]]], 
                                          range(cutoff_transDF[[spatLocs_to_use[1]]]), 
                                          range(cutoff_transDF[[spatLocs_to_use[2]]]), 
                                          marks = factor(cutoff_transDF[[cellID_coln]]), 
                                          unitname = c("um","um"))
    } else {
      queryTrans_pp <- spatstat.geom::pp3(x = cutoff_transDF[[spatLocs_colns[1]]], 
                                          y = cutoff_transDF[[spatLocs_colns[2]]], 
                                          z = cutoff_transDF[[spatLocs_colns[3]]],
                                          spatstat.geom::box3(range(cutoff_transDF[[spatLocs_colns[1]]]), 
                                                              range(cutoff_transDF[[spatLocs_colns[2]]]), 
                                                              range(cutoff_transDF[[spatLocs_colns[3]]]),
                                                              unitname = "um"),
                                          marks = factor(cutoff_transDF[[cellID_coln]]))
    }
    
    
    # get distribution of minimal molecule-to-molecule distance for each transcript in query cell
    dist_profile <- quantile(spatstat.geom::nndist(queryTrans_pp), seq(0,1,by=0.1))
    message(sprintf("Distribution of minimal molecular distance between %d cells: %s, at quantile = %s.", 
                    length(cells_to_keep), 
                    paste0(round(dist_profile, 2), collapse = ", "), 
                    paste0(names(dist_profile), collapse = ", ")))
    # define cutoff as 5 times of 90% quantile value
    molecular_distance_cutoff <- 5*dist_profile[['90%']]
    message(sprintf("Use 5 times of 90%% quantile of minimal %dD molecular distance between picked cells as `molecular_distance_cutoff` = %.4f for defining direct neighbor cells.", 
                    length(spatLocs_to_use), molecular_distance_cutoff))
    rm(queryTrans_pp, cutoff_transDF)
    # final results list
    res[['molecular_distance_cutoff']] <- molecular_distance_cutoff
    res[['distance_profile']] <- dist_profile
  }
  
  return(res)

}



#' @title get_baselineCT
#' @description get cluster-specific quantile distribution of transcript number and per cell per molecule transcript score in the provided cell x gene expression matrix based on the reference profiles and cell cluster assignment
#' @param refProfiles A matrix of cluster profiles, genes X clusters
#' @param counts Counts matrix, cells X genes.
#' @param clust Vector of cluster assignments for each cell in `counts`, default = NULL to automatically assign the cell cluster for each cell based on maximum transcript score  
#' @return a list
#' \enumerate{
#'    \item{span_score, a matrix of average transcript tLLR score per molecule per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns}
#'    \item{span_transNum, a matrix of transcript number per cell for each distinct cell types in row, percentile at (0%, 25%, 50%, 75%, 100%) in columns}
#'    \item{score_baseline, a named vector of 25% quantile of cluster-specific per cell transcript score, to be used as score baseline such that  per cell transcript score higher than the baseline is required to call a cell type of high enough confidence} 
#'    \item{lowerCutoff_transNum, a named vector of 25% quantile of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that higher than the cutoff is required to keep query cell as it is}
#'    \item{higherCutoff_transNum, a named vector of median value of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.}
#'    \item{clust_used,  a named vector of cluster assignments for each cell used in baseline calculation, cell_ID in `counts` as name}
#' }
#' @details Calculate average per molecule transcript score for each cell in `counts` expression matrix based on the provided cluster profiles `refProfiles` and cluster assignment for each cell `clust`; then get the quantile distribution of transcript number and per molecule per cell transcript score under each cluster. The function would also recommend the cutoff for transcript score and transcript number to be used in re-segmentation pipeline based on the calculated quantile distribution. 
#' @examples 
#' data(example_refProfiles)
#' data(ori_RawExprs)
#' baselineData <- get_baselineCT(refProfiles = example_refProfiles, counts = ori_RawExprs, clust = NULL)
#' @export
get_baselineCT <- function(refProfiles, 
                           counts, 
                           clust = NULL){
  # get common genes
  common_genes <- intersect(rownames(refProfiles), colnames(counts))

  if(length(common_genes)<1){
    stop("Too few common genes to proceed. Check if `refProfiles` is a gene x cell-type matrix and `counts` is a cell x gene matrix.")
  } else {
    message(sprintf("Found %d common genes among `refProfiles` and `counts`. ", 
                    length(common_genes)))
  }
  
  if(!is.null(clust)){
    if(!is.vector(clust)){
      stop("The provided `clust` is not a vector of cluster assignment.")
    }
    
    if(length(clust) != nrow(counts)){
      message("`clust` has different length from the row number of `counts`.")
      clust = NULL
    } 
  } 
  
  common_celltypes <- intersect(unique(clust), colnames(refProfiles))
  if(length(common_celltypes) ==0){
    message("No common cell types/clusters found between `clust` and `refProfiles`.")
    clust = NULL
    common_celltypes <- colnames(refProfiles)
  } else {
    # drop the cells with clusters not in refProfiles for baseline calculation
    # cells with common clusters
    tmp_idx <- which(clust %in% colnames(refProfiles))
    if(length(tmp_idx) != length(clust)){
      message(sprintf("%d cells with assigned `clust` not presented in `refProfiles`: `%s`; exclude those cells from baseline calcuation.", 
                      length(clust) - length(tmp_idx),
                      paste0(setdiff(unique(clust), colnames(refProfiles)), collapse = "`, `")))
      clust <- clust[tmp_idx]
      counts <- as.matrix(counts)[tmp_idx, ] 
    }

    rm(tmp_idx)

  }
  
  # filter and re-order data
  refProfiles <- as.matrix(refProfiles)[common_genes, common_celltypes]
  counts <- as.matrix(counts)[, common_genes]
  
  # get score matrix based on refProfiles for each gene and cell ----
  # replace zero in mean profiles with 1E-5
  refProfiles <- pmax(refProfiles, 1e-5)
  # tLLR score, re-center on maximum per row/transcript
  tLLRv2_geneMatrix <- scoreGenesInRef(genes = common_genes, ref_profiles = refProfiles)
  
  
  # get cell x cell-cluster score matrix = counts (cell x gene) %*% tLLR_score (gene x cell-cluster)
  tLLRv2_cellMatrix <- counts %*% tLLRv2_geneMatrix
  
  # assign cell type for each cell if not provided ----
  if(is.null(clust)){
    message('Perform cluster assignment based on maximum transcript score given the provided `refProfiles`.')

    # assign cell type based on max values
    max_idx_1st <- max.col(tLLRv2_cellMatrix, ties.method="first")
    clust <- colnames(tLLRv2_cellMatrix)[max_idx_1st]

    common_celltypes <- unique(clust)
    rm(max_idx_1st)
  }
  
  # get transcript number quantile profile ---
  all_transNum <- rowSums(counts)
  span_transNum_CellType <- tapply(all_transNum, 
                                   clust, 
                                   function(x) quantile(x, probs = seq(0, 1, 0.25)))
  span_transNum_CellType <- do.call(rbind, span_transNum_CellType)
  
  # get transcript score quantile profile ----
  # loop over each cell type to get score based on assigned clusters
  all_tLLRv2 <- rep(NA, length(clust))
  for (each_celltype in common_celltypes){
    rowidx <- which(clust == each_celltype)
    all_tLLRv2[rowidx] <- tLLRv2_cellMatrix[rowidx, each_celltype]
  }
  # normalized by transcript number to get per molecule transcript score for each cell
  all_tLLRv2 <- all_tLLRv2/all_transNum
  span_tLLRv2_CellType <- tapply(all_tLLRv2, 
                                 clust, 
                                 function(x) quantile(x, probs = seq(0, 1, 0.25)))
  span_tLLRv2_CellType <- do.call(rbind, span_tLLRv2_CellType)
  
  # return final results ---
  names(clust) <- rownames(counts)
  final_res <- list(span_score = span_tLLRv2_CellType, 
                    span_transNum = span_transNum_CellType, 
                    score_baseline = span_tLLRv2_CellType[, "25%"],
                    lowerCutoff_transNum = span_transNum_CellType[, "25%"], 
                    higherCutoff_transNum = span_transNum_CellType[, "50%"], 
                    clust_used = clust)
  
  return(final_res)

}
