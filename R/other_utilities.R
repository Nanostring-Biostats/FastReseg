#' @title run_igraph_leiden
#' @description Run Leiden clustering with version-compatible resolution argument.
#' @param graph An igraph object representing the graph to cluster.
#' @param ... Additional arguments passed to `cluster_leiden()`, such as
#'   `objective_function`, `resolution_parameter`, `beta`, `initial_membership`,
#'   and `n_iterations`.
#' @return A clustering object returned by `cluster_leiden()`.
#' @details This function wraps `igraph::cluster_leiden()` and ensures compatibility
#' with both older (<2.1.0) and newer versions of the `igraph` package by renaming
#' the `resolution_parameter` argument to `resolution` if needed.
#' @importFrom igraph cluster_leiden
run_igraph_leiden <- function(graph, ...) {
  args <- list(...)
  
  # Rename resolution_parameter to resolution if igraph >= 2.1.0
  if (packageVersion("igraph") >= "2.1.0") {
    if ("resolution_parameter" %in% names(args)) {
      args$resolution <- args$resolution_parameter
      args$resolution_parameter <- NULL
    }
  }
  
  # Call cluster_leiden with modified arguments
  do.call(igraph::cluster_leiden, c(list(graph), args))
}


#' @title igraph_delete_edges
#' @description Delete edges from an igraph object with version compatibility
#' @param graph An igraph object from which edges will be deleted.
#' @param edges A vector of edge IDs or an edge selector to delete.
#'
#' @return An igraph object with the specified edges removed.
#' @details This function wraps `igraph::delete_edges()` and ensures compatibility
#' with older versions (<2.0.0) of `igraph` that used `delete.edges()`.
igraph_delete_edges <- function(graph, edges) {
  if (packageVersion("igraph") >= "2.0.0") {
    igraph::delete_edges(graph, edges)
  } else {
    igraph::delete.edges(graph, edges)
  }
}

