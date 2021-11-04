#' @title check_config_spatialNW
#' @description check config used to create spatial network, assign default values if missing arguments
#' @param config a list of config would be used to create spatial network via Giotto::createSpatialNetwork
#' @param spat_locs a data.frame with spatial location info
#' @importFrom SMITAP checkTypeLengthValue
#' @return the corrected config list
check_config_spatialNW <- function(config, spat_locs){
  # check config, assign default value if NULL
  if(is.null(config$method)){
    config$method <- 'Delaunay'
    message("Create Delanay network when config$method is NULL.")
  }
  config$method <- match.arg(config$method, c("Delaunay", "kNN"))
  msg <- character()
  
  # check common config for both methods
  if(!is.null(config$name)){
    msg <- c(msg, SMITAP::checkTypeLengthValue(config, "name", expect_type = "character", 
                                               expect_len = 1))
  } else{
    # assign default name based on method
    config$name <- paste0(config$method,'_network')
    message(sprintf("Name the spatial network based on method as `%s` when config$name is NULL.", 
                    config$name))
  }
  
  # check dimensions in spat_locs
  if(!is.null(config$dimensions)){
    # if integer vector, get column name based on index
    if(is.numeric(config$dimensions) | is.integer(config$dimensions)){
      locs_colnames <- colnames(spat_locs)[config$dimenions]
      # must be within  c("sdimx", "sdimy","sdimz")
      if(any(!(locs_colnames %in% c("sdimx","sdimy","sdimz")))){
        msg <- c(msg, "dimensions in config have values other than c(\"sdimx\",\"sdimy\",\"sdimz\"). ")
      } 
    } else {
      # if character, check column name
      msg <- c(msg, SMITAP::checkTypeLengthValue(config, "dimensions", expect_type = "character",
                                                 expect_range = "within", 
                                                 expect_value = c("all","sdimx","sdimy","sdimz")))
    }
    
  } else{
    # use 'all'
    config$dimensions <- 'all'
  }
  # check minimum_k, non-negative integer
  if(!is.null(config$minimum_k)){
    msg <- c(msg, SMITAP::checkTypeLengthValue(config, "minimum_k", expect_type = c("numeric","integer"), 
                                               expect_len = 1, expect_range = "larger", 
                                               expect_value = (-1e-23)))
  }else {
    config$minimum_k <- 0
  }
  
  # check config for delanay network
  if(config$method == 'Delaunay'){
    if(!is.null(config$delaunay_method)){
      msg <- c(msg, SMITAP::checkTypeLengthValue(config, "delaunay_method", 
                                                 expect_type = "character", 
                                                 expect_len = 1, expect_range = "within",
                                                 expect_value = c("deldir", "delaunayn_geometry", "RTriangle")))
    } else {
      config$delaunay_method <- "delaunayn_geometry" # only method works for 3D
    }
    # distant cutoff = "auto" or non-negative number
    if(!is.null(config$maximum_distance_delaunay)){
      # check for character
      msg1 <- SMITAP::checkTypeLengthValue(config,"maximum_distance_delaunay", 
                                           expect_type = "character", expect_len = 1, 
                                           expect_range = "within", expect_value = "auto")
      if(length(msg1) > 0L){
        # not character, check for numeric values
        msg2 <- SMITAP::checkTypeLengthValue(config, "maximum_distance_delaunay",
                                             expect_type = "numeric", expect_len = 1, 
                                             expect_range = "larger", 
                                             expect_value = (-1e-23))
        # if both not satisfied, output both messages
        if(length(msg2) > 0L){
          msg <- c(msg, msg1, msg2)
        }
      }
      
    } else {
      config$maximum_distance_delaunay <- 'auto'
    }
    
    # options for "delaunayn_geometry", too many possible options to list here. 
    # see details in http://www.qhull.org/html/qdelaun.htm
    if(!is.null(config$options) & config$delaunay_method == "delaunayn_geometry"){
      msg <- c(msg, SMITAP::checkTypeLengthValue(config, "options", expect_type = "character"))
    } else {
      config$options <- 'Pp'
    }
    
    if(!is.null(config$Y) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,SMITAP::checkTypeLengthValue(config,"Y", expect_type = "logical", 
                                                expect_len = 1))
    } else {
      config$Y <- TRUE
    }
    if(!is.null(config$j) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,SMITAP::checkTypeLengthValue(config,"j", expect_type = "logical", 
                                                expect_len = 1))
    } else {
      config$j <- TRUE
    }
    if(!is.null(config$S) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,SMITAP::checkTypeLengthValue(config,"S", expect_type = c("numeric","integer"), 
                                                expect_len = 1, expect_range = "larger", 
                                                expect_value = (-1e-23)))
    } else {
      config$S <- 0
    }
    
  }
  
  # check config for kNN network
  if(config$method == 'kNN'){
    if(!is.null(config$knn_method)){
      msg <- c(msg, SMITAP::checkTypeLengthValue(config, "knn_method", expect_type = "character", 
                                                 expect_len = 1, expect_range = "within", 
                                                 expect_value = "dbscan"))
    } else {
      config$knn_method <- 'dbscan'
    }
    if(!is.null(config$k)){
      # positive integer
      msg <- c(msg, SMITAP::checkTypeLengthValue(config,"k", expect_type = c("numeric","integer"), 
                                                 expect_len = 1, expect_range = "larger", 
                                                 expect_value = 0))
    } else {
      config$k <- 4
    }
    if(!is.null(config$maximum_distance_knn)){
      # NULL or non-negative number
      msg <- c(msg, SMITAP::checkTypeLengthValue(config, "maximum_distance_knn",
                                                 expect_type = "numeric", expect_len = 1, 
                                                 expect_range = "larger", 
                                                 expect_value = (-1e-23)))
    } else {
      config$maximum_distance_knn <- NULL
    }
    
  }
  
  if(length(msg) > 0L){
    stop( "Configuration Issues:\n" , paste( msg , collapse = "\n" ) )
  }
  return(config)
}


#' @title createSpatialDelaunayNW_from_spatLocs
#' @description generate delaunay network based on provided config and spatial location using Giotto functions
#' @param config_spatNW configuration list 
#' @param spatLocs_df data.frame for spatial location of each entry for cell or transcript
#' @param ID_column column name for entry ID in spatLocs_df
#' @param spatLocs_column column name for 1st, 2nd, optional 3rd dimension of spatial coordinates in spatLocs_df 
#' @importFrom data.table as.data.table
#' @return delaunay_network_Obj, a spatial network object created by Giotto function
createSpatialDelaunayNW_from_spatLocs <- function(config_spatNW, 
                                                  spatLocs_df, 
                                                  ID_column = 'cell_ID',
                                                  spatLocs_column = c("sdimx","sdimy","sdimz")){
  if(!'data.frame' %in% class(spatLocs_df)){
    stop("spatLocs_df is not a data.frame with columns for transcript ID and coordinates.")
  }
  if(!ID_column %in% colnames(spatLocs_df)){
    stop(sprintf("ID_column = `%s` is not found in the column names of spatLocs_df.", 
                 ID_column))
  }
  if(any(!spatLocs_column %in% colnames(spatLocs_df))){
    stop(sprintf("spatLocs_column = `%s` is not found in the column names of spatLocs_df.", 
                 paste0(spatLocs_column, collapse ='`, `' )))
  }
  if(!length(spatLocs_column) %in% c(2,3)){
    stop("spatLocs_column must have 2 or 3 elements to define 1st, 2nd, 3rd spatial dimension in spatLocs_df.")
  }
  
  # reformat the column name to "cell_ID", "sdimx","sdimy","sdimz"
  spatLocs_df <- data.table::as.data.table(spatLocs_df)
  spatLocs_df <- spatLocs_df[,.SD, .SDcols = c(ID_column, spatLocs_column)]
  colname_converter <- colnames(spatLocs_df)
  names(colname_converter) <- colnames(spatLocs_df)
  colname_converter[[ID_column]] <- "cell_ID"
  colname_converter[[spatLocs_column[1]]] <- 'sdimx'
  colname_converter[[spatLocs_column[2]]] <- 'sdimy'
  if(length(spatLocs_column) ==3){
    # 3rd dimension
    colname_converter[[spatLocs_column[3]]] <- 'sdimz'
  } 
  
  colnames(spatLocs_df) <- unname(colname_converter)
  
  # check config_spatNW
  config_spatNW <- check_config_spatialNW(config = config_spatNW,
                                          spat_locs = spatLocs_df)
  if(config_spatNW$method != 'Delaunay'){
    stop(sprintf("This function only creates Delaunay network. But the provided config_spatNW$method = `%s`.", config_spatNW$method))
  }
  
  #### create delaunay spatial network: ---------------------------
  method <- match.arg(config_spatNW$delaunay_method, c("deldir", "delaunayn_geometry", 
                                                       "RTriangle"))
  spatial_locations <- data.table::as.data.table(spatLocs_df)
  # net matrix for spatial location
  spatLocs_matrix <- spatial_locations[, grepl("sdim", colnames(spatial_locations)), 
                                       with = F]
  spatLocs_matrix <- as.matrix(spatLocs_matrix)
  d2_or_d3 = dim(spatLocs_matrix)[2]
  
  if (d2_or_d3 == 2) {
    # 2D network
    first_dimension = colnames(spatLocs_matrix)[[1]]
    second_dimension = colnames(spatLocs_matrix)[[2]]
    spatLoc_in_use <- spatial_locations[, c("cell_ID", first_dimension, second_dimension), with = F]
    if (method == "RTriangle") {
      delaunay_output = Giotto:::create_delaunayNetwork_RTriangle(spatial_locations = spatLoc_in_use, 
                                                                  sdimx = first_dimension, sdimy = second_dimension, 
                                                                  Y = config_spatNW$Y, j = config_spatNW$j, S = config_spatNW$S)
      outputObj = delaunay_output$geometry_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(maximum_distance = config_spatNW$maximum_distance_delaunay, 
                        minimum_k = config_spatNW$minimum_k, 
                        Y = config_spatNW$Y, j = config_spatNW$j, S = config_spatNW$S)
      outputObj = outputObj
      
    }
    else if (method == "deldir") {
      delaunay_output = Giotto:::create_delaunayNetwork_deldir(spatial_locations = spatLoc_in_use, 
                                                               sdimx = first_dimension, sdimy = second_dimension)
      outputObj = delaunay_output$geometry_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(maximum_distance = config_spatNW$maximum_distance_delaunay, 
                        minimum_k = config_spatNW$minimum_k)
      outputObj = outputObj
    }
    else if (method == "delaunayn_geometry") {
      delaunay_output = Giotto:::create_delaunayNetwork_geometry(spatial_locations = spatLoc_in_use, 
                                                                 sdimx = first_dimension, sdimy = second_dimension, 
                                                                 options = config_spatNW$options)
      outputObj = delaunay_output$geometry_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(options = config_spatNW$options)
      outputObj = outputObj
    }
    # if return empty networkDT, return NULL directly
    if(nrow(delaunay_network_DT) == 0 ){
      return(NULL)
    }
    
    delaunay_network_DT = Giotto:::calculate_distance_and_weight(delaunay_network_DT, 
                                                                 sdimx = first_dimension, 
                                                                 sdimy = second_dimension, 
                                                                 d2_or_d3 = 2)
  }else if (d2_or_d3 == 3) {
    if (method != "delaunayn_geometry") {
      stop(method, " method only applies to 2D data, use delaunayn_geometry, see details \n")
    }
    else {
      # 3D network
      first_dimension = colnames(spatLocs_matrix)[[1]]
      second_dimension = colnames(spatLocs_matrix)[[2]]
      third_dimension = colnames(spatLocs_matrix)[[3]]
      spatLoc_in_use <- spatial_locations[, c("cell_ID", first_dimension, second_dimension,third_dimension), with = F]
      delaunay_output <- Giotto:::create_delaunayNetwork_geometry_3D(spatial_locations = spatLoc_in_use, 
                                                                     sdimx = first_dimension, sdimy = second_dimension, 
                                                                     sdimz = third_dimension, options = config_spatNW$options)
      outputObj = delaunay_output$geometry_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(options = config_spatNW$options)
      outputObj = outputObj
    }
    
    # if return empty networkDT, return NULL directly
    if(nrow(delaunay_network_DT) == 0 ){
      return(NULL)
    }
    
    delaunay_network_DT = Giotto:::calculate_distance_and_weight(delaunay_network_DT, 
                                                                 sdimx = first_dimension, sdimy = second_dimension,
                                                                 sdimz = third_dimension, d2_or_d3 = 3)
    
    
  }
  
  networkDT_before_filter = delaunay_network_DT
  delaunay_network_DT = Giotto:::filter_network(networkDT_before_filter, 
                                                maximum_distance = config_spatNW$maximum_distance_delaunay, 
                                                minimum_k = config_spatNW$minimum_k)
  meanCellDistance = Giotto:::get_distance(delaunay_network_DT, method = "mean")
  medianCellDistance = Giotto:::get_distance(delaunay_network_DT, method = "median")
  cellShapeObj = list(meanCellDistance = meanCellDistance, 
                      medianCellDistance = medianCellDistance)
  delaunay_network_Obj = Giotto:::create_spatialNetworkObject(name = config_spatNW$name, 
                                                              method = method, parameters = parameters, outputObj = outputObj, 
                                                              networkDT = delaunay_network_DT, networkDT_before_filter = networkDT_before_filter, 
                                                              cellShapeObj = cellShapeObj, misc = NULL)
  return(delaunay_network_Obj)
  
}
