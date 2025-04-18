#' Data for cell type specific reference profiles
#'
#' Semi-supervised cell typing is conducted on the example SMI RNA readout for FFPE melanoma tissue samples that contains 2 tissue sections with 3 FOVs per section. The resulting cell typing outcomes are then used to estimate the mean profile of each identified cell type/cluster using \code{"estimate_MeanProfile"}. The corresponding average cell-type-specific profiles are then stored in \code{data("example_refProfiles")} to be used as reference profiles. 
#'
#' @docType data
#'
#' @usage data(example_refProfiles)
#'
#' @format An object of class \code{"matrix"} with 960 genes in rows, 22 distinct cell types in columns.
#'
#' @keywords datasets
#'
#' @examples
#' data(example_refProfiles)
"example_refProfiles"



#' Transcript level data frame from Example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. It contains information for 1375 cells with 756783 transcripts within a cropped region of one single FOV of 1 tissue section.
#'
#' @docType data
#'
#' @usage data(mini_transcriptDF)
#'
#' @format An object of class \code{"data.frame"} with 756783 transcripts in row and 9 variables
#' \describe{
#'    \item{UMI_transID}{unique id for transcript}
#'    \item{UMI_cellID}{unique id for cell based on original cell segmentaion assignment}
#'    \item{x}{spatial coordinate of the transcript in x-axis of the given FOV, unit in micron}
#'    \item{y}{spatial coordinate of the transcript in y-axis of the given FOV, unit in micron}
#'    \item{z}{spatial coordinate of the transcript in z-axis of the given FOV, unit in micron}
#'    \item{target}{gene identity of given transcript}
#'    \item{slide}{slide ID of given transcript}
#'    \item{fov}{fov ID of given transcript}
#'    \item{CellId}{the cell label assignment of given transcript within given FOV based on original cell segmentation}
#' }
#' @seealso [ori_CellStatsDF] for corresponding cell level dataset, [ori_RawExprs] for corresponding cell x gene expression matrix
#' @keywords datasets
"mini_transcriptDF"



#' Cell level data frame from Example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. It contains information for 1375 cells with 756783 transcripts within a cropped region of one single FOV of 1 tissue section.
#'
#' @docType data
#'
#' @usage data(ori_CellStatsDF)
#'
#' @format An object of class \code{"data.frame"} with 1375 cells in row and 10 variables
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
#' @seealso [mini_transcriptDF] for corresponding transcript level dataset, [ori_RawExprs] for corresponding cell x gene expression matrix
#' @keywords datasets
"ori_CellStatsDF"



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
#' @seealso [mini_transcriptDF] for corresponding transcript level dataset, [ori_CellStatsDF] for corresponding cell level dataset
#' 
#' @keywords datasets
"ori_RawExprs"


#' Percentile profile for single cell distribution of transcript number and average tLLR score from example dataset for spatial transcriptional profiling of tissue
#'
#' The example dataset is based on one run of SMI RNA readout for FFPE melanoma tissue samples. The single cell distribution was calculated for transcript number per cell and average transcript tLLR score per cell, and then grouped by the assigned cell types based on the original cell segmentation outcome to get the percentile profile for each cell type at 25% increment step. 
#'
#' @docType data
#'
#' @usage data(example_baselineCT)
#'
#' @format An object of class \code{"list"} with 2 elements
#' \describe{
#'    \item{span_score}{a matrix of average transcript tLLR score per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns}
#'    \item{span_transNum}{a matrix of transcript number per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns}
#'    }
#' @keywords datasets
"example_baselineCT"



#' per cell gene expression matrix based on raw transcript data files generated by SMI
#'
#' The example dataset is based on the cropped regions from two different FOVs of same SMI RNA readout run for FFPE melanoma tissue sample. 
#' The corresponding raw transcript data files are 
#' \itemize{
#'    \item{`extdata/Run4104_FOV001__complete_code_cell_target_call_coord.csv`}
#'    \item{`extdata/Run4104_FOV002__complete_code_cell_target_call_coord.csv`}
#' }
#' Each raw transcript data files have 1 transcript per row with 25 different meta information in columns. Some important columns are
#' \describe{
#'    \item{x}{spatial coordinate of the transcript in x-axis of the given FOV, unit in pixel, 0.18um per pixel}
#'    \item{y}{spatial coordinate of the transcript in y-axis of the given FOV, unit in pixel, 0.18um per pixel}
#'    \item{z}{spatial coordinate of the transcript in z-axis of the given FOV, unit in z-step, 0.8um per z-step}
#'    \item{target}{gene identity of given transcript}
#'    \item{CellId}{the cell label assignment of given transcript within given FOV based on original cell segmentation}
#' } 
#' The per cell gene expression matrix is calculated from the raw transcript data files to include only the 960 true target genes without including control probes.
#'
#' @docType data
#'
#' @usage data(example_CellGeneExpr)
#'
#' @format An object of class \code{"matrix"} with 754 cells in row and 960 genes in columns
#'
#' @seealso [example_clust] for corresponding cell cluster assignment
#' 
#' @keywords datasets
"example_CellGeneExpr"



#' Cluster assignment for cells based on per cell gene expression profile and reference profiles in example data set
#'
#' The example dataset is based on the cropped regions from two different FOVs of same SMI RNA readout run for FFPE melanoma tissue sample. 
#' The corresponding raw transcript data files are 
#' \itemize{
#'    \item{`extdata/Run4104_FOV001__complete_code_cell_target_call_coord.csv`}
#'    \item{`extdata/Run4104_FOV002__complete_code_cell_target_call_coord.csv`}
#' }
#' Each raw transcript data files have 1 transcript per row with 25 different meta information in columns. Some important columns are
#' \describe{
#'    \item{x}{spatial coordinate of the transcript in x-axis of the given FOV, unit in pixel, 0.18um per pixel}
#'    \item{y}{spatial coordinate of the transcript in y-axis of the given FOV, unit in pixel, 0.18um per pixel}
#'    \item{z}{spatial coordinate of the transcript in z-axis of the given FOV, unit in z-step, 0.8um per z-step}
#'    \item{target}{gene identity of given transcript}
#'    \item{CellId}{the cell label assignment of given transcript within given FOV based on original cell segmentation}
#' } 
#' The cluster assignment for each cell is calculated based on the corresponding per cell gene expression profiles, `data/example_CellGeneExpr.RData`, and reference profiles, `data/example_refProfiles.RData`. 
#'
#' @docType data
#'
#' @usage data(example_clust)
#'
#' @format An object of class \code{"character"} with 754 cells in same order as the cells in the row of `data/example_CellGeneExpr.RData`. 
#'
#' @seealso [example_CellGeneExpr] for corresponding cell x gene expression matrix, [example_refProfiles] for the reference profiles used for cluster assignment
#' 
#' @keywords datasets
"example_clust"


