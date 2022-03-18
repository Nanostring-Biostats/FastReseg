#### pipeline to process multi-FOV multi-slide data
## The pipeline assumes the pre-existing of single cell typing results on original cell segmentation outcomes. 
## Firstly, calculate reference profiles based on the cell x gene expression matrix and cell typing results of entire data set ;
## Secondly, check the existence of fov offset position and transcript data.frame for each FOV in each slide; 
## Thirdly, loop through each FOV data.frame, reformat the data.frame to use unique transcript IDs and cell IDs and a global coordinate system, then score cells for segmentation errors and flag transcripts with low goodness-of-fit to current cell segment. 

library(ggplot2)
library(FastReseg)

main_output_dir <- '/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/melanoma/Run4104_melanoma_980plx/giotto_output'
sub_out_dir <- fs::path(main_output_dir, "fastReseg03")
if(!dir.exists(sub_out_dir)){
  dir.create(sub_out_dir, recursive = T)
}

setwd(sub_out_dir)



#### (1) common settings for detection system ----
# pixel size = 0.18um/pixel, z-step size = 0.8um/step
pixel_size = 0.18
zStep_size = 0.8

# cell type coln to use: `nb_clus` for nbclust outcomes, `leiden_clus` for leiden clustering outcomes
refClus_coln <- "nb_clus"

# exclude cells with low confident cell typing outcomes from reference profile estimation
# NULL, if include all cell typing outcomes
cellClus_to_exclude <- 'NotDet' 


# flag to remove FOVs with unpaired target call files and fov position information
# if FALSE, stop processing when missing target call files
removeUnpaired <- FALSE

# SMI TAP output folder
scAnalysis_dir <- "/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/melanoma/Run4104_melanoma_980plx/giotto_output/Run4104_cellpose_vs_oldDASH"
# get config files to get path to sample annotation file
source(fs::path(scAnalysis_dir, 'config_980plex.R'))



# path to giotto object containing cell typing and per cell info
path_to_SMIobject <- fs::path(config_main$resultspath, 'results/complete_giotto_object.RData')

rm(list = setdiff(ls(pattern = '^config_'), c('config_loading')))

# optional: exclude 33 high expressors from reference profile estimation in case of 980plx ISH panel
is_980plx <- TRUE
if(is_980plx){
  blacklist_genes <- fs::path("/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/melanoma/Run4104_melanoma_980plx/giotto_output/Run4104_redo_cellpose", 
                              "20201210_33HE_GeneName_in_1013plx.txt")
  blacklist_genes <- read.csv(blacklist_genes, header = TRUE, sep = '\t')
  blacklist_genes <- blacklist_genes[['GeneName']]
} else {
  blacklist_genes <- NULL
}

#### (2) prepare input for resegmentation from SMI object ----
smi_inputs <- prepSMI_for_fastReseg(path_to_SMIobject = path_to_SMIobject,  
                                    config_loading = config_loading, 
                                    refClus_coln = refClus_coln,
                                    cellClus_to_exclude = cellClus_to_exclude, 
                                    removeUnpaired = removeUnpaired,
                                    blacklist_genes = blacklist_genes,
                                    pixel_size = pixel_size)
# write `transDF_fov_fileInfo` into csv file
write.csv(smi_inputs[['transDF_fov_fileInfo']], 
          file = fs::path(sub_out_dir, 'transDF_fov_fileInfo.csv'))

#### (3) other parameters used in resegmentation workflow ----
# ## use default values for those parameters, change them as needed below and then pass to `fastReseg_internalRef` function
# # cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
# flagModel_TransNum_cutoff = 50 
# 
# # cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
# flagCell_lrtest_cutoff = 5
# 
# # cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
# svmClass_score_cutoff = -2 
# 
# # a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
# svm_args = list(kernel = "radial", 
#                 scale = FALSE, 
#                 gamma = 0.4)


#### (4) loop through each FOV to do resegmentation on pre-defined refProfiles, cutoffs ----
## open a log file to store all function outputs into disk
logFile <- file.path(sub_out_dir, "fastReseg_log.txt")
logHandle <- file(logFile, open = "wt")
sink(logHandle)
sink(logHandle, type = "message")

## run resegmentation
# # `findSegmentError_allFiles` function returns a list of 
# refProfiles: a genes * clusters matrix of cluster-specific reference profiles used in resegmenation pipeline.
# baselineData: a list of two matrice in cluster * percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.
# combined_modStats_ToFlagCells: a data.frame for spatial modeling statistics of each cell for all cells in the data set, output of `score_cell_segmentation_error` function.
# combined_flaggedCells: a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV.

reseg_outputs <- findSegmentError_allFiles(counts = smi_inputs[['counts']],
                                           clust = NULL,
                                           refProfiles = smi_inputs[['refProfiles']],
                                           transDF_fileInfo = smi_inputs[['transDF_fov_fileInfo']],
                                           filepath_coln = 'file_path',
                                           prefix_colns = c('slide','fov'),
                                           fovOffset_colns = c('offset_x', 'offset_y'), # match XY axes between stage and each FOV
                                           pixel_size = pixel_size, 
                                           zstep_size = zStep_size, 
                                           transcript_df = NULL,
                                           transID_coln = NULL, # row index as transcript_id
                                           transGene_coln = "target",
                                           cellID_coln = "CellId",
                                           spatLocs_colns = c("x","y","z"),
                                           extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data
                                           path_to_output = sub_out_dir, 
                                           combine_extra = TRUE) # if TRUE, extracellular and trimmed transcripts are included in the updated transcript data.frame


## revert output back to the console -- only then access the file!
sink(type = "message")
sink()

# # show the log file in pop up window
# file.show(logFile)

# print all lines in log file on console
cat(readLines(logFile), sep="\n")



### save resegmenation inputs, outcomes and new giotto object on file
save(smi_inputs, 
     reseg_outputs,
     file = fs::path(sub_out_dir, 'fastReseg03_inputs_outputs.RData'))

## check processing speed ----
## use file.info() to get file name, size and time created to see the timing of segmentation
# all updated_transDF.csv output
updated_transDF_list <- dir(path = sub_out_dir, full.names = TRUE, 
                            pattern = "[0-9]+_flagged_transDF.csv$")
fileMeta <- file.info(updated_transDF_list, extra_cols = FALSE)
fileMeta[['full_path']] <- rownames(fileMeta)
fileMeta[['file_name']] <- basename(fileMeta[['full_path']])
fileMeta[['idx']] <- as.numeric(sapply(strsplit(fileMeta[['file_name']], '_'), '[[', 1))
fileMeta <- fileMeta[order(fileMeta[['idx']]), ]
rownames(fileMeta) <- fileMeta[['idx']]
fileMeta[['endTime']] <- fileMeta[['mtime']]
fileMeta[['startTime']] <- c(fileMeta[['endTime']][1], fileMeta[['endTime']][1: nrow(fileMeta)-1])
fileMeta[['duration']] <- apply(fileMeta, 1, 
                                function(x) difftime(x[['endTime']], 
                                                     x[['startTime']], 
                                                     units = "mins"))

fileMeta[['fileSize']] <- utils:::format.object_size(fileMeta[['size']], units = "MB") 
fileMeta[['fileSize']] <- as.numeric(sapply(strsplit(fileMeta[['fileSize']], ' Mb'), '[[', 1))

plot_data <- fileMeta[2: nrow(fileMeta), c('fileSize','duration', 'idx')]
p <- ggplot(plot_data, aes(x = fileSize, y = duration, color = idx))+
  geom_point()+
  geom_smooth(method = lm)+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x = 'transcript file size in MB', y = 'resegment duration (min) of each FOV', 
       title = paste0('FastReseg speed: ', 
                      round(as.numeric(difftime(fileMeta$endTime[nrow(fileMeta)], 
                                                fileMeta$endTime[1], units = "hours")), 2),
                      ' hours for ', nrow(fileMeta)-1, ' FOVs of ', 
                      sum(fileMeta$fileSize) - fileMeta$fileSize[1], ' MB size'))+
  theme_linedraw()

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

colnames(plot_data) <- c('x','y')
p1 <- p+ annotate(geom = 'text', label = lm_eqn(plot_data), parse = TRUE,
                  x = -Inf, y = Inf, hjust = 0, vjust = 1, color = 'red', size = 5)

jpeg(filename = fs::path(sub_out_dir, "FastReseg03_process_speed_vs_fileSize.jpeg"), width = 600, height = 600)
print(p1)
dev.off()

# plot processing speed vs. processing order
plot_data <- fileMeta[2: nrow(fileMeta), c('fileSize','duration', 'idx')]
p <- ggplot(plot_data, aes(x = idx, y = duration, color = fileSize))+
  geom_point()+
  geom_smooth(method = lm)+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x = 'sequential order of processing', y = 'resegment duration (min) of each FOV', 
       title = paste0('FastReseg speed: ', 
                      round(as.numeric(difftime(fileMeta$endTime[nrow(fileMeta)], 
                                                fileMeta$endTime[1], units = "hours")), 2),
                      ' hours for ', nrow(fileMeta)-1, ' FOVs of ', 
                      sum(fileMeta$fileSize) - fileMeta$fileSize[1], ' MB size'))+
  theme_linedraw()

jpeg(filename = fs::path(sub_out_dir, "FastReseg03_process_speed_vs_order.jpeg"), width = 600, height = 600)
print(p)
dev.off()


### (5) visualization of segmenation outputs ----



