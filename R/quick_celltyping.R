#' @title lldist
#' @description calculate the log-likelihood of cell belong to certain cell cluster given the reference profile using negative binomial model
#' @param x a vector of a reference cell cluster
#' @param mat a cells x genes matrix of expression levels in all cells
#' @param bg a vector for background level of each cell (default: 0.01)
#' @param size the parameters for dnbinom function (default: 10)
#' @param digits the number of digits for rounding
#'
#' @importFrom Matrix rowSums
#' @importFrom stats dnbinom
#' @return a cell x cell_type matrix of the log-likelihood
#' @export
lldist <- function(x, mat, bg = 0.01, size = 10, digits = 2) {
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat <- as(matrix(mat, nrow = 1), "dgCMatrix")
  } else if (is.matrix(mat)) {
    mat <- as(mat, "dgCMatrix")
  } else if (!is(mat, "dgCMatrix")) {
    errorMessage <-
      sprintf(
        "The `type` of parameter `mat` needs to be of one of dgCMatrix, vector, matrix, array, but is found to be of type %s",
        class(mat)
      )
    stop(errorMessage)
  }
  
  # accept a single value of bg if input by user:
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(mat))
  }
  
  # Check dimensions on bg and stop with informative error if not conformant
  if (is.vector(bg)) {
    if (!identical(length(bg), nrow(mat))) {
      errorMessage <- sprintf("Dimensions of count matrix and background are not conformant.\nCount matrix rows: %d, length of bg: %d",
                              nrow(mat), length(bg))
      stop(errorMessage)
    }
  }
  
  # calc scaling factor to put y on the scale of x:
  if (is.vector(bg) & nrow(mat) >1) {
    bgsub <- apply(mat, 2, function(xx) xx - bg)
    bgsub <- pmax(bgsub, 0)
  } else { 
    # bg is a matrix or mat only has 1 row
    bgsub <- pmax(mat - bg, 0)
  }
  
  sum_of_x <- sum(x)
  s <- Matrix::rowSums(bgsub) / sum_of_x
  # override it if s is negative:
  s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum_of_x
  
  # yhat is a n_cells x n_genes matrix of expected values 
  yhat <- s %*% t(x)

  # log-likelihood, cell x gene matrix
  res <- stats::dnbinom(x = as.matrix(mat), size = size, mu = yhat, log = TRUE)
  # get sum of loglike per cell under current cell cluster
  res <- rowSums(res)
  
  names(res) <- rownames(mat)
  return(round(res, digits))
}


#' Get number of cores for parallelized operations
#'
#' @return number of cores to use for mclapply
#' @export
numCores <- function() {
  num_cores <- 1
  if (.Platform$OS.type == "unix") {
    if (is.null(getOption("mc.cores"))) {
      num_cores <- parallel::detectCores() - 2
    } else {
      num_cores <- getOption("mc.cores")
    }
    
  }
  return(num_cores)
}

#' @title quick_celltype
#' @description Classify cells based on reference profiles given the maximum log-likelihood calculated via negative binomial model
#' @param x Counts matrix (or dgCMatrix), cells x genes.
#' @param bg a vector for background level of each cell (default = 0.01)
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters, genes x clusters.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @return A list, with the following elements:
#' \describe{
#'    \item{clust}{a vector given cells' cluster assignments, return NA for cells of zero counts. }
#'    \item{logliks}{a cells x clusters matrix of cells' log-likelihoods under each cluster, return -Inf for cells of zero counts. }
#'    \item{zeroCells}{a vector of cells of zero count, return NULL if none}
#' }
#' @export
quick_celltype <- function(x, bg = 0.01, reference_profiles, nb_size = 10, align_genes = TRUE) {
  
  if (any(rowSums(x) == 0)) {
    zeroCells <- rownames(x)[rowSums(x)==0]
    message(sprintf("%d cells with 0 counts are found. Return `clust = NA`, `logliks = -Inf` for those cells: `%s`.", 
                    length(zeroCells), paste0(zeroCells, collapse = "`, `")))
    
    # bg is a vector of same length as x
    if(is.vector(bg)){
      if(identical(length(bg), nrow(x))){
        bg <- bg[rowSums(x)!=0]
      }
      
    } else {
      # bg is a cell x gene matrix
      bg <- bg[rowSums(x)!=0, ]
    }
    
    # remove zero cells from x
    x <- x[rowSums(x)!=0, ]
    
  } else {
    zeroCells <- NULL
  }

  # accept a single value of bg if input by user:
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(x))
    names(bg) <- rownames(x)
  }
  
  
  # align genes:
  if (align_genes) {
    sharedgenes <- intersect(rownames(reference_profiles), colnames(x))
    lostgenes <- setdiff(colnames(x), rownames(reference_profiles))
    
    # subset:
    x <- x[, sharedgenes]
    reference_profiles <- reference_profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) && length(lostgenes) < 50) {
      message(paste0("The following genes in the count data are missing from reference_profiles and will be omitted from cell typing: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from reference_profiles and will be omitted from cell typing"))
    }
  }

  
  # get logliks
  logliks <- parallel::mclapply(asplit(reference_profiles, 2),
                                lldist,
                                mat = x,
                                bg = bg,
                                size = nb_size,
                                mc.cores = numCores())
  logliks <- do.call(cbind, logliks)
  
  
  # get remaining outputs
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  names(clust) <- rownames(logliks)
  
  out <- list(clust = clust,
              logliks = round(logliks, 4), 
              zeroCells = zeroCells)
  
  # add in zeroCell data
  if(!is.null(zeroCells)){
    clust <- rep(NA, length(zeroCells))
    names(clust) <- zeroCells
    
    logliks <- matrix(-Inf, nrow = length(zeroCells), ncol = ncol(out[['logliks']]), 
                      dimnames = list(zeroCells, colnames(out[['logliks']])))
    
    out[['clust']] <- c(out[['clust']], clust)
    out[['logliks']] <- rbind(out[['logliks']], logliks)
  }
  
  return(out)    
}
