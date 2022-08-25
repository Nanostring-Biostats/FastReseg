#' Calculate the likelihood of the vector x using the reference vector of y
#'
#' @param x a vector of a reference cell type
#' @param mat a matrix of expression levels in all cells
#' @param bg background level (default: 0.01)
#' @param size the parameters for dnbinom function (default: 10)
#' @param digits the number of digits for rounding
#'
#' @importFrom Matrix rowSums
#' @importFrom stats dnbinom
#'
#' @export
lldist <- function(x, mat, bg = 0.01, size = 10, digits = 2) {
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat <- as(matrix(mat, nrow = 1), "dgCMatrix")
  } else if (is.matrix(mat)){
    mat <- as(mat, "dgCMatrix")
  } else if(class(mat)!="dgCMatrix"){
    errorMessage <- sprintf( "The `type` of parameter `mat` needs to be of one of dgCMatrix, vector, matrix, array, but is found to be of type %s",class(mat))
    stop( errorMessage )
  }
  # Check dimensions on bg and stop with informative error if not
  #  conformant
  if ( is.vector( bg ) )
  {
    if ( !identical( length( bg ) , nrow( mat ) ) )
    {
      errorMessage <- sprintf( "Dimensions of count matrix and background are not conformant.\nCount matrix rows: %d, length of bg: %d" , nrow( mat ) , length( bg ) )
      stop( errorMessage )
    }
  }
  
  # calc scaling factor to put y on the scale of x:
  if ( is.vector( bg ) )
  {
    bgsub <- pmax( sweep( mat , 1 , bg , "-" ) , 0 )
  }
  else
  {
    bgsub <- pmax( mat - bg , 0 )
  }
  sum_of_x <- sum( x )
  s <- Matrix::rowSums( bgsub ) / sum_of_x
  # override it if s is negative:
  s[s <= 0] = Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum_of_x
  
  # expected counts:
  if ( is.vector( bg ) )
  {
    yhat <- as(sweep(Matrix::Matrix(s) %*% Matrix::t(x) , 1 , bg , "+" ), "dgCMatrix")
  }
  else
  {
    yhat <- as(Matrix::Matrix(s) %*% Matrix::t(x) + bg, "dgCMatrix")
  }
  # loglik:
  # lls <- stats::dnbinom(x = Matrix::as.matrix(mat), size = size, mu = yhat, log = TRUE)
  lls <- dnbinom_sparse(x = mat, mu = yhat, size_dnb = size)
  rownames(lls) <- rownames(mat)
  colnames(lls) <- colnames(mat)
  
  return(round(Matrix::rowSums(lls), digits))
}
