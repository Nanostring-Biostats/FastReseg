#' @title plotSpatialScoreMultiCells
#' @description plot spatial plot of transcripts within chosen cells, colored by score, facet by cell_ID 
#' @param chosen_cells the cell_ID of chosen cells
#' @param cell_labels title labels for each cell's spatial plot
#' @param transcript_df the data.frame of transcript_ID, cell_ID, spatial locations
#' @param cellID_coln the column name of cell_ID in transcript_df
#' @param transID_coln the column name of transcript_ID in transcript_df
#' @param score_coln the column name of score in transcript_df 
#' @param spatLocs_colns column names for 1st, 2nd dimension of spatial coordinates in transcript_df 
#' @param point_size marker size for transcript spot in fig (default = 0.1)
#' @param plot_discrete flag to plot transcript score in discrete color 
#' @param title the text title for plot
#' @return ggplot
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient2 labs facet_wrap theme_classic theme element_text unit margin element_rect
#' @importFrom data.table .SD
#' @export
plotSpatialScoreMultiCells <- function(chosen_cells, cell_labels, 
                                       transcript_df, 
                                       cellID_coln = "CellId", 
                                       transID_coln = "transcript_id",
                                       score_coln = "score", 
                                       spatLocs_colns = c("x","y"),
                                       point_size = 0.1,
                                       plot_discrete = FALSE, 
                                       title = NULL){
  
  # check format of cell_labels
  if(length(chosen_cells) != length(cell_labels)){
    stop("chosen_cells and cell_labels are not in pair and have differenet length.")
  }
  names(cell_labels) <- chosen_cells
  
  if(length(spatLocs_colns) != 2){
    stop("spatLocs_colns must be the column names for 1st, 2nd dimension of spatial coordinates in transcript_df.")
  }
  
  # check format of transcript_df
  if(any(!c(cellID_coln, transID_coln, score_coln, spatLocs_colns) %in% colnames(transcript_df))){
    stop(sprintf("Not all necessary columns can be found in provided transcript_df, missing columns include `%s`.",
                 paste0(setdiff(c(cellID_coln, transID_coln, score_coln, spatLocs_colns), 
                                colnames(transcript_df)), collapse = "`, `")))
  }
  transcript_df <- data.table::as.data.table(transcript_df)
  transcript_df <- transcript_df[, .SD, .SDcols = c(cellID_coln, transID_coln,score_coln, spatLocs_colns)]
  transcript_df <- transcript_df[which(transcript_df[[cellID_coln]] %in% chosen_cells), ]
  
  transcript_df[['cell_labels']] <- cell_labels[transcript_df[[cellID_coln]]]
  
  if(plot_discrete){
    transcript_df[[score_coln]] <- as.factor(transcript_df[[score_coln]])
  }
  
  fig <- ggplot2::ggplot(transcript_df, ggplot2::aes(x = get(spatLocs_colns[1]), y = get(spatLocs_colns[2]), color = get(score_coln)))+
    ggplot2::geom_point(size = point_size, shape = 16)
  
  # spatial gradient for color if not discrete input
  if (! "factor" %in% class(transcript_df[[score_coln]])){
    # get mean score to for diverging color scheme
    score_mean <- mean(transcript_df[[score_coln]])
    
    fig <- fig +
      #scale_color_gradientn(name = score_coln, colours = topo.colors(10))+
      ggplot2::scale_color_gradient2(name = score_coln, midpoint = score_mean,
                                     low = "blue",mid = "gray90", high = "red")
  } 
  
  fig <- fig + ggplot2::labs(x = spatLocs_colns[1], y = spatLocs_colns[2], color = score_coln)+
    ggplot2::facet_wrap(~ cell_labels + get(cellID_coln), scales = "free")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = ggplot2::element_text(size = 6),
                   legend.position = "right",
                   legend.key.size = ggplot2::unit(16, "pt"),
                   legend.margin=ggplot2::margin(0,0,0,0),
                   legend.box.margin=ggplot2::margin(0,0,0,0),
                   strip.text = ggplot2::element_text(face = "bold", size = 6),
                   strip.background = ggplot2::element_rect(fill = "gray85"))
  
  if(!is.null(title)){
    fig <- fig + ggplot2::ggtitle(title)
  }
  
  return(fig)
}
