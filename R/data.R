#' Data for cell type specific reference profiles
#'
#' Semi-supervised cell typing is conducted on the example SMI RNA readout for FFPE melanoma tissue samples. The resulting cell typing outcomes are then used to estimate the mean profile of each identified cell type/cluster using \code{"estimate_MeanProfile"}. The corresponding average cell-type-specific profiles are then stored in \data{"refProfiles"} to be used as reference profiles. 
#'
#' @docType data
#'
#' @usage data(refProfiles)
#'
#' @format An object of class \code{"matrix"} with 960 genes in rows, 22 distinct cell types in columns.
#'
#' @keywords datasets
#'
#' @examples
#' data(refProfiles)
"refProfiles"



#' Transcript level data table from Example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. It contains information for 1375 cells with 756783 transcripts within a cropped region of one single FOV of 1 tissue section.
#'
#' @docType data
#'
#' @usage data(transcriptDT)
#'
#' @format An object of class \code{"data.table"} with 756783 transcripts in row and 9 variables
#' \describe{
#'    \item{transcript_id}{unique id for transcript}
#'    \item{cell_ID}{unique id for cell based on original cell segmentaion assignment}
#'    \item{x}{spatial coordinate of the transcript in x-axis of the given FOV, unit in micron}
#'    \item{y}{spatial coordinate of the transcript in y-axis of the given FOV, unit in micron}
#'    \item{z}{spatial coordinate of the transcript in z-axis of the given FOV, unit in micron}
#'    \item{target}{gene identity of given transcript}
#'    \item{slide}{slide ID of given transcript}
#'    \item{fov}{fov ID of given transcript}
#'    \item{CellId}{the cell label assignment of given transcript within given FOV based on original cell segmentation}
#' }
#'
#' @keywords datasets
"transcriptDT"



#' Cell level data table from Example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. It contains information for 1375 cells with 756783 transcripts within a cropped region of one single FOV of 1 tissue section.
#'
#' @docType data
#'
#' @usage data(ori_CellStatsDT)
#'
#' @format An object of class \code{"data.table"} with 1375 cells in row and 10 variables
#' \describe{
#'    \item{cell_ID}{unique id for cell based on original cell segmentaion assignment}
#'    \item{slide}{slide ID of given cell}
#'    \item{fov}{fov ID of given cell}
#'    \item{CellId}{the cell label assignment of given cell within given FOV based on original cell segmentation}
#'    \item{CenterX}{spatial coordinate of the cell centroid in x-axis of the given FOV, unit in micron}
#'    \item{CenterY}{spatial coordinate of the cell centroid in y-axis of the given FOV, unit in micron}
#'    \item{Width}{the bounding box length of given cell segment along x-axis, unit in micron}
#'    \item{Height}{the bounding box length of given cell segment along y-axis, unit in micron}
#'    \item{Area}{the area of given cell segment, unit in square micron}
#'    \item{AspectRatio}{the aspect ratio of boundingbox width over height for given cell segment}
#' }
#'
#' @keywords datasets
"ori_CellStatsDT"



#' Cell x Gene expression matrix from Example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. It contains the single cell gene expression profile for 1375 cells across 960 genes.
#'
#' @docType data
#'
#' @usage data(ori_RawExprs)
#'
#' @format An object of class \code{"matrix"} with 4619 cells in row and 960 genes in columns
#'
#' @keywords datasets
"ori_RawExprs"
