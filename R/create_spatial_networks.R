#' @title createSpatialDelaunayNW_from_spatLocs
#' @description generate delaunay network based on provided config and spatial location using `GiottoClass` functions
#' @param config_spatNW configuration list 
#' @param spatLocs_df data.frame for spatial location of each entry for cell or transcript
#' @param ID_column column name for entry ID in `spatLocs_df`
#' @param spatLocs_column column name for 1st, 2nd, optional 3rd dimension of spatial coordinates in `spatLocs_df` 
#' @importFrom data.table as.data.table .SD
#' @importFrom methods new
#' @return delaunay_network_Obj, a spatial network object created by `GiottoClass` functions
#' @details This function leverages `GiottoClass` package to create spatial networks from spatial coordinates. An example `config_spatNW` list is shown below with possible options on controlling the spatial network generation. For more details, see the manual for `GiottoClass::createSpatialNetwork`. 
#'#' ' \describe{
#'    \item{name}{spatial network name; default = 'spatial_network'}
#'    \item{dimensions}{a vector for which spatial dimensions to use, default = 'all' to use all dimentions}
#'    \item{method}{method name for creating a spatial network, default = 'Delaunay'}
#'    \item{minimum_k}{minimum number of nearest neighbors if maximum_distance != NULL}
#'    \item{delaunay_method}{Delaunay method to use, choose from c("delaunayn_geometry", "deldir", "RTriangle"), default = "delaunayn_geometry"}
#'    \item{maximum_distance_delaunay}{distance cuttoff for nearest neighbors to consider for Delaunay network, default = "auto"}
#'    \item{options}{(geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options; default = `Pp`, do not report precision problems)}
#'    \item{Y}{(RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary}
#'    \item{j}{(RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.}
#'    \item{S}{(RTriangle) Specifies the maximum number of added Steiner points.}
#' }
#' @export
createSpatialDelaunayNW_from_spatLocs <- function(config_spatNW = list(name = 'spatial_network',
                                                                       dimensions = "all",
                                                                       method = 'Delaunay',
                                                                       minimum_k = 0,
                                                                       delaunay_method = "delaunayn_geometry",
                                                                       maximum_distance_delaunay = "auto",
                                                                       options = "Pp",
                                                                       Y = TRUE,
                                                                       j = TRUE,
                                                                       S = 0), 
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
      delaunay_output = GiottoClass:::.create_delaunaynetwork_RTriangle(
        spatial_locations = spatLoc_in_use, 
        sdimx = first_dimension, sdimy = second_dimension, 
        Y = config_spatNW$Y, j = config_spatNW$j, S = config_spatNW$S)
      
      outputObj = delaunay_output$RTriangle_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(maximum_distance = config_spatNW$maximum_distance_delaunay, 
                        minimum_k = config_spatNW$minimum_k, 
                        Y = config_spatNW$Y, j = config_spatNW$j, S = config_spatNW$S)
      outputObj = outputObj
      
    }
    else if (method == "deldir") {
      delaunay_output = GiottoClass:::.create_delaunaynetwork_deldir(
        spatial_locations = spatLoc_in_use, 
        sdimx = first_dimension, sdimy = second_dimension)
      
      outputObj = delaunay_output$deldir_obj
      delaunay_network_DT = delaunay_output$delaunay_network_DT
      parameters = list(maximum_distance = config_spatNW$maximum_distance_delaunay, 
                        minimum_k = config_spatNW$minimum_k)
      outputObj = outputObj
    }
    else if (method == "delaunayn_geometry") {
      delaunay_output = GiottoClass:::.create_delaunaynetwork_geometry(
        spatial_locations = spatLoc_in_use, 
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
    
    delaunay_network_DT = GiottoClass:::.calculate_distance_and_weight(
      delaunay_network_DT, 
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
      delaunay_output <- GiottoClass:::.create_delaunaynetwork_geometry_3d(
        spatial_locations = spatLoc_in_use, 
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
    
    delaunay_network_DT = GiottoClass:::.calculate_distance_and_weight(
      delaunay_network_DT, 
      sdimx = first_dimension, sdimy = second_dimension,
      sdimz = third_dimension, d2_or_d3 = 3)
    
    
  }
  
  networkDT_before_filter = delaunay_network_DT
  delaunay_network_DT = GiottoClass:::.filter_network(networkDT_before_filter, 
                                                      maximum_distance = config_spatNW$maximum_distance_delaunay, 
                                                      minimum_k = config_spatNW$minimum_k)
  meanCellDistance = GiottoClass:::get_distance(delaunay_network_DT, method = "mean")
  medianCellDistance = GiottoClass:::get_distance(delaunay_network_DT, method = "median")
  cellShapeObj = list(meanCellDistance = meanCellDistance, 
                      medianCellDistance = medianCellDistance)
  delaunay_network_Obj = new("spatialNetworkObj", 
                             name = config_spatNW$name, 
                             method = method, parameters = parameters, outputObj = outputObj, 
                             networkDT = delaunay_network_DT, networkDT_before_filter = networkDT_before_filter, 
                             cellShapeObj = cellShapeObj, crossSectionObjects = NULL, 
                             spat_unit = "cell", provenance = NULL, misc = NULL)
  return(delaunay_network_Obj)
  
}


#' @title checkTypeLengthValue
#' @description checking whether a single value in config have correct data type, length and value range
#' @param config the list storing config
#' @param name parameter name in config list
#' @param expect_type expected data type for vector, e.g. c("numeric","integer") for integer values, "character", "logical", etc
#' @param expect_len expected length, NULL for any length
#' @param expect_range expect value range; use "larger","smaller", "equal", "within" for numeric values, but use only "within" for character
#' @param expect_value a single numeric value or a character vector that would be used to check against with expect_range, NULL for any value
#' @return a message if config[[name]] does not satisfy the criteria
#' @examples 
#' config <- list(pos_Integer = 1, neg_value = -0.2, flag = TRUE,
#'                length2character = c("a","b"))
#' # check if positive integer of any length
#' FastReseg:::checkTypeLengthValue(config, "pos_Integer",
#'                                  expect_type = c("numeric","integer"),
#'                                  expect_range = "larger", expect_value = 0)
#' # check if negative value of any length
#' FastReseg:::checkTypeLengthValue(config, "neg_value",
#'                                  expect_type = "numeric",
#'                                  expect_range = "smaller", expect_value = 0)
#' # check if logical value
#' FastReseg:::checkTypeLengthValue(config, "flag", expect_type = "logical")
#' # check if character has 2 elements within c("a","b")
#' FastReseg:::checkTypeLengthValue(config, "length2character",
#'                                  expect_type = "character", expect_len = 2,
#'                                  expect_range = "within", expect_value = c("a","b"))
checkTypeLengthValue <- function(config, name, 
                                 expect_type, 
                                 expect_len = NULL,
                                 expect_range = c("equal","larger","smaller","within"),
                                 expect_value = NULL){
  
  expect_range = match.arg(expect_range, c("equal","larger","smaller","within"))
  
  msg <- character()
  
  # check class in general
  if(!(class(config[[name]]) %in% expect_type)){
    msg <- c(msg, sprintf("`%s` value in config is either not defined or not in class of `%s`. ", 
                          name, paste0(expect_type, collapse = "`, `")))
  }
  
  # if expect integer in numeric form, use expect_type = c("numeric","integer")
  # check for the integer case
  if(("integer" %in% expect_type) & (length(msg) == 0L)){
    if(any(round(config[[name]]) != config[[name]])){
      msg <- c(msg, sprintf( "`%s` value in config is not integer.", name))
    }
  }
  
  # check length 
  if(!is.null(expect_len)){
    if(length(expect_len) >1){
      stop("expect_len has more than 1 element.")
    }
    if(!is.numeric(expect_len)){
      stop("expect_len is not a number.")
    }
    # when expect_len is a single number, do length checking
    if(length(config[[name]]) != expect_len){
      msg <- c(msg, sprintf("`%s` value in config is not in expected length of %d.", name, round(expect_len)))
    }
  }
  
  
  # check value only when correct type and length
  if(!is.null(expect_value) & (length(msg) == 0L)){
    
    # if numeric values check if satisfy expect_range and values
    if("numeric" %in% expect_type){
      if(expect_range == "larger"){
        if(any(config[[name]] <= expect_value)){
          msg <- c(msg, sprintf( "`%s` value in config is no larger than %d.", 
                                 name, round(expect_value)))
        }
      } else if (expect_range == "smaller"){
        if(any(config[[name]] >= expect_value)){
          msg <- c(msg, sprintf( "`%s` value in config is no smaller than %d.", 
                                 name, round(expect_value)))
        }
      } else if (expect_range == "equal"){
        if(any(config[[name]] != expect_value)){
          msg <- c(msg, sprintf( "`%s` value in config is no equal to %d.", 
                                 name, round(expect_value)))
        } 
      } else if(expect_range == "within"){
        if(any(!(config[[name]] %in% expect_value))){
          msg <- c(msg, sprintf( "`%s` value in config contains value other than `%s`.", 
                                 name, paste0(expect_value, collapse = "`, `")))
        }
      } else {
        stop(sprintf("`%s` cannot be used for expect_type = `%s`.", 
                     expect_range, paste0(expect_type, collapse = "`, `")))
      }
      
    }
    
    # if character values check if within expect_values
    if("character" %in% expect_type){
      if(expect_range == "within"){
        if(any(!(config[[name]] %in% expect_value))){
          msg <- c(msg, sprintf( "`%s` value in config contains value other than `%s`.", 
                                 name, paste0(expect_value, collapse = "`, `")))
        }
      } else {
        stop(sprintf("`%s` cannot be used for expect_type = `%s`.", 
                     expect_range, paste0(expect_type, collapse = "`, `")))
      }
    }
    
  }
  
  return(msg)
}



#' @title check_config_spatialNW
#' @description check config used to create spatial network, assign default values if missing arguments
#' @param config a list of config would be used to create spatial network using delaunay method only via GiottoClass::createSpatialNetwork
#' @param spat_locs a data.frame with spatial location info
#' @importFrom GiottoClass createSpatialNetwork
#' @return the corrected config list
check_config_spatialNW <- function(config, spat_locs){
  # check config, assign default value if NULL
  if(is.null(config$method)){
    config$method <- 'Delaunay'
    message("Create Delanay network when config$method is NULL.")
  }
  # allow only 'Delaunay'
  config$method <- match.arg(config$method, c("Delaunay"))
  msg <- character()
  
  # check common config for both methods
  if(!is.null(config$name)){
    msg <- c(msg, checkTypeLengthValue(config, "name", expect_type = "character", 
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
      msg <- c(msg, checkTypeLengthValue(config, "dimensions", expect_type = "character",
                                                 expect_range = "within", 
                                                 expect_value = c("all","sdimx","sdimy","sdimz")))
    }
    
  } else{
    # use 'all'
    config$dimensions <- 'all'
  }
  # check minimum_k, non-negative integer
  if(!is.null(config$minimum_k)){
    msg <- c(msg, checkTypeLengthValue(config, "minimum_k", expect_type = c("numeric","integer"), 
                                               expect_len = 1, expect_range = "larger", 
                                               expect_value = (-1e-23)))
  }else {
    config$minimum_k <- 0
  }
  
  # check config for delanay network
  if(config$method == 'Delaunay'){
    if(!is.null(config$delaunay_method)){
      msg <- c(msg, checkTypeLengthValue(config, "delaunay_method", 
                                                 expect_type = "character", 
                                                 expect_len = 1, expect_range = "within",
                                                 expect_value = c("delaunayn_geometry", "deldir", "RTriangle")))
    } else {
      config$delaunay_method <- "delaunayn_geometry" # only method works for 3D
    }
    # distant cutoff = "auto" or non-negative number
    if(!is.null(config$maximum_distance_delaunay)){
      # check for character
      msg1 <- checkTypeLengthValue(config,"maximum_distance_delaunay", 
                                           expect_type = "character", expect_len = 1, 
                                           expect_range = "within", expect_value = "auto")
      if(length(msg1) > 0L){
        # not character, check for numeric values
        msg2 <- checkTypeLengthValue(config, "maximum_distance_delaunay",
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
      msg <- c(msg, checkTypeLengthValue(config, "options", expect_type = "character"))
    } else {
      config$options <- 'Pp'
    }
    
    if(!is.null(config$Y) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,checkTypeLengthValue(config,"Y", expect_type = "logical", 
                                                expect_len = 1))
    } else {
      config$Y <- TRUE
    }
    if(!is.null(config$j) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,checkTypeLengthValue(config,"j", expect_type = "logical", 
                                                expect_len = 1))
    } else {
      config$j <- TRUE
    }
    if(!is.null(config$S) & config$delaunay_method == "RTriangle"){
      msg <- c(msg,checkTypeLengthValue(config,"S", expect_type = c("numeric","integer"), 
                                                expect_len = 1, expect_range = "larger", 
                                                expect_value = (-1e-23)))
    } else {
      config$S <- 0
    }
    
  }
  
  
  if(length(msg) > 0L){
    stop( "Configuration Issues:\n" , paste( msg , collapse = "\n" ) )
  }
  return(config)
}


