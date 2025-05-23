% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_spatial_networks.R
\name{createSpatialDelaunayNW_from_spatLocs}
\alias{createSpatialDelaunayNW_from_spatLocs}
\title{createSpatialDelaunayNW_from_spatLocs}
\usage{
createSpatialDelaunayNW_from_spatLocs(
  config_spatNW = list(name = "spatial_network", dimensions = "all", method = "Delaunay",
    minimum_k = 0, delaunay_method = "delaunayn_geometry", maximum_distance_delaunay =
    "auto", options = "Pp", Y = TRUE, j = TRUE, S = 0),
  spatLocs_df,
  ID_column = "cell_ID",
  spatLocs_column = c("sdimx", "sdimy", "sdimz")
)
}
\arguments{
\item{config_spatNW}{configuration list}

\item{spatLocs_df}{data.frame for spatial location of each entry for cell or transcript}

\item{ID_column}{column name for entry ID in \code{spatLocs_df}}

\item{spatLocs_column}{column name for 1st, 2nd, optional 3rd dimension of spatial coordinates in \code{spatLocs_df}}
}
\value{
delaunay_network_Obj, a spatial network object created by \code{GiottoClass} functions
}
\description{
generate delaunay network based on provided config and spatial location using \code{GiottoClass} functions
}
\details{
This function leverages \code{GiottoClass} package to create spatial networks from spatial coordinates. An example \code{config_spatNW} list is shown below with possible options on controlling the spatial network generation. For more details, see the manual for \code{GiottoClass::createSpatialNetwork}.
#' ' \describe{
\item{name}{spatial network name; default = 'spatial_network'}
\item{dimensions}{a vector for which spatial dimensions to use, default = 'all' to use all dimentions}
\item{method}{method name for creating a spatial network, default = 'Delaunay'}
\item{minimum_k}{minimum number of nearest neighbors if maximum_distance != NULL}
\item{delaunay_method}{Delaunay method to use, choose from c("delaunayn_geometry", "deldir", "RTriangle"), default = "delaunayn_geometry"}
\item{maximum_distance_delaunay}{distance cuttoff for nearest neighbors to consider for Delaunay network, default = "auto"}
\item{options}{(geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options; default = \code{Pp}, do not report precision problems)}
\item{Y}{(RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary}
\item{j}{(RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.}
\item{S}{(RTriangle) Specifies the maximum number of added Steiner points.}
}
}
