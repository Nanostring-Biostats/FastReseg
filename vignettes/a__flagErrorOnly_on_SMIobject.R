#### pipeline to process multi-FOV multi-slide data
## The pipeline assumes the pre-existing of single cell typing results on original cell segmentation outcomes. 
## Firstly, calculate reference profiles based on the cell x gene expression matrix and cell typing results of entire data set ;
## Secondly, check the existence of fov offset position and transcript data.frame for each FOV in each slide; 
## Thirdly, loop through each FOV data.frame, reformat the data.frame to use unique transcript IDs and cell IDs and a global coordinate system, then score cells for segmentation errors and flag transcripts with low goodness-of-fit to current cell segment. 
## Fourth, examples on how to change cutoff for flagging cells and transcript groups after initial processing

library(ggplot2)
library(FastReseg)

main_output_dir <- '/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/melanoma/Run4104_melanoma_980plx/giotto_output'
sub_out_dir <- fs::path(main_output_dir, "fastReseg03")
if(!dir.exists(sub_out_dir)){
  dir.create(sub_out_dir, recursive = T)
}

setwd(sub_out_dir)



#### (1) common settings for detection system ----
# flag to return per cell expression matrix after trimming all putative contaminating transctipts 
return_trimmed_perCell = TRUE 

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

# flag to include ctrl_genes in analysis
include_ctrlgenes <- FALSE

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
# # `fastReseg_flag_all_errors` function returns a list of 
# refProfiles: a genes * clusters matrix of cluster-specific reference profiles used in resegmenation pipeline.
# baselineData: a list of two matrice in cluster * percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.
# ctrl_genes: a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.
# combined_modStats_ToFlagCells: a data.frame for spatial modeling statistics of each cell for all cells in the data set, output of `score_cell_segmentation_error` function.
# combined_flaggedCells: a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV.

if(include_ctrlgenes){
  # if ctrl_genes = smi_inputs$ctrl_genes, the control genes like NegPrbs and FalseCodes would be included in analysis.
  ctrl_genes <- smi_inputs$ctrl_genes
} else {
  ctrl_genes <- NULL
}

reseg_outputs <- fastReseg_flag_all_errors(
  counts = smi_inputs[['counts']],
  clust = NULL,
  refProfiles = smi_inputs[['refProfiles']],
  
  transcript_df = NULL,
  transDF_fileInfo = smi_inputs[['transDF_fov_fileInfo']],
  filepath_coln = 'file_path',
  prefix_colns = c('slide','fov'),
  fovOffset_colns = c('offset_x', 'offset_y'), # match XY axes between stage and each FOV
  pixel_size = pixel_size, 
  zstep_size = zStep_size, 
  
  transID_coln = NULL, # row index as transcript_id
  transGene_coln = "target",
  cellID_coln = "CellId",
  spatLocs_colns = c("x","y","z"),
  extracellular_cellID = c(0), # CellId = 0 means extracelluar transcripts in raw data
  
  path_to_output = sub_out_dir, 
  combine_extra = TRUE,  # if TRUE, extracellular and trimmed transcripts are included in the updated transcript data.frame
  ctrl_genes = ctrl_genes,
  return_trimmed_perCell = return_trimmed_perCell
  )


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

#### (5) create new giotto object using the trimmed per cell count ----
if(return_trimmed_perCell){
  ## (5.1) prepare expression list for multi-slot giotto object
  # when use ctrl_genes = NULL, 
  # NegPrb and FalseCode genes are absence from refernece profiles and thus missing in the transcript data.frame with updated cell segmentaion. 
  # assume 0 values for both NegPrb and FalseCode genes for now
  if(is.null(ctrl_genes)){
    expr_lists <- list(rna = list(raw = reseg_outputs[['trimmed_perCellExprs']]))
    expr_lists[['negprobes']] <- list(raw = matrix(0, nrow = 5,
                                                   ncol = ncol(reseg_outputs[['trimmed_perCellExprs']]),
                                                   dimnames = list(paste0('NegPrb', seq_len(5)),
                                                                   colnames(reseg_outputs[['trimmed_perCellExprs']]))))
    expr_lists[['falsecode']] <- list(raw = matrix(0, nrow = 5,
                                                   ncol = ncol(reseg_outputs[['trimmed_perCellExprs']]),
                                                   dimnames = list(paste0('FalseCode', seq_len(5)),
                                                                   colnames(reseg_outputs[['trimmed_perCellExprs']]))))
  } else {
    tmpExp <- reseg_outputs[['trimmed_perCellExprs']]
    
    expr_lists <- list(rna = list(raw = tmpExp[!(rownames(tmpExp) %in% ctrl_genes), ]))
    
    expr_lists[['negprobes']] <- list(raw = tmpExp[rownames(tmpExp) %in% ctrl_genes[grepl('NegPrb', ctrl_genes)], ])
    expr_lists[['falsecode']] <- list(raw = tmpExp[rownames(tmpExp) %in% ctrl_genes[grepl('FalseCode', ctrl_genes)], ])
    
    rm(tmpExp)
    
  }
  
  
  ## (5.2) prepare cell annotation file, get it from original giotto object
  getOriginalCellMeta <- function(path_to_SMIobject, gname = "gem", extraColns){
    load(path_to_SMIobject)
    cell_annotDF <- merge(Giotto::pDataDT(get(gname)), 
                          Giotto:::select_spatial_locations(get(gname)), 
                          by = 'cell_ID')
    colns_to_keep <- intersect(colnames(cell_annotDF), 
                               c('cell_ID', 'sdimx', 'sdimy',
                                 'fov', 'Area', 'AspectRatio','Width','Height', 
                                 'slide_ID_numeric', extraColns))
    colns_to_keep <- c(colns_to_keep, grep('^Max|^Mean', colnames(cell_annotDF), value = T))
    cell_annotDF <- as.data.frame(cell_annotDF)[, unique(colns_to_keep)]
    return(cell_annotDF)
  }
  
  cell_annotDF <- getOriginalCellMeta(path_to_SMIobject = path_to_SMIobject, 
                                      extraColns = colnames(smi_inputs$sample_annot))
  gc()
  # rearrange cell order to be the same as expression matrix
  cell_annotDF <- cell_annotDF[match(colnames(reseg_outputs[['trimmed_perCellExprs']]), cell_annotDF[['cell_ID']]),]
  
  ## (5.3) create giotto object
  updated_SMIobj <- Giotto::createGiottoObject(expression= expr_lists,
                                               expression_feat= names(expr_lists),
                                               spatial_locs= cell_annotDF[, c('sdimx', 'sdimy')],
                                               cell_metadata= list(rna = cell_annotDF, 
                                                                   negprobes = cell_annotDF, 
                                                                   falsecode = cell_annotDF))
  
  ### save new giotto object on file
  saveRDS(updated_SMIobj, 
       file = fs::path(sub_out_dir, 'fastResegTrimmed_updated_SMIobj.RData'))
}



#### Below are example scripts for adjusting individual cutoffs ----
## no rerun of the whole pipeline, but only modify on a specific step.
## use only when you are looking for the impact of specific parameters on a few FOVs. 

#### (6) adjust the cutoff for flagging cells as needed ----
message(sprintf("%d cells, %.2f%% of all cells, are flagged for potential cell segemntation error in `reseg_output`. ", 
                length(unlist(reseg_outputs$combined_flaggedCells)), 
                length(unlist(reseg_outputs$combined_flaggedCells))/nrow(reseg_outputs$combined_modStats_ToFlagCells)))

# get spatial modeling statistics of each cell for all cells in the data set
combined_modStats_ToFlagCells <- reseg_outputs$combined_modStats_ToFlagCells


# cutoff of lrtest_-log10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile (default =5)
# lower values would flag more cells with potential segmentation error
flagCell_lrtest_cutoff = 3

# cells with potential segmentation errors, flagged by new cutoff
combined_flaggedCells <- combined_modStats_ToFlagCells[combined_modStats_ToFlagCells['lrtest_-log10P'] > flagCell_lrtest_cutoff, 'UMI_cellID']
message(sprintf("%d cells, %.2f%% of all cells, are flagged for potential cell segemntation error based on provided `flagCell_lrtest_cutoff` = %.2f. ", 
                length(combined_flaggedCells), 
                length(combined_flaggedCells)/nrow(reseg_outputs$combined_modStats_ToFlagCells), 
                flagCell_lrtest_cutoff))


#### (7) redo identification of low goodness-of-fit transcript groups as needed ----
## (7.1) cutoff for flagging transcript groups of low goodness-of-fit
# spatial dimension of provided data for evaluation
d2_or_d3 = 3

# cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
flagModel_TransNum_cutoff = 50

# cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
svmClass_score_cutoff = -2

# a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
svm_args = list(kernel = "radial",
                scale = FALSE,
                gamma = 0.4)

# transcript score matrix from reference profiles
transcript_loglik <- FastReseg::scoreGenesInRef(genes = rownames(reseg_outputs$refProfiles), ref_profiles = pmax(reseg_outputs$refProfiles, 1e-5))
tmp_max <- apply(transcript_loglik, 1, max)
score_GeneMatrix <- sweep(transcript_loglik, 1, tmp_max, '-')
rm(tmp_max, transcript_loglik)


# set tLLR score for control genes, same as `svmClass_score_cutoff`
if(!is.null(ctrl_genes)){
  message(sprintf("Include the following `ctrl_genes` in analysis: `%s`.\nIt's recommended to have total counts of those genes below 1%% of total counts of all genes in each cell.", 
                  paste0(ctrl_genes, collapse = "`, `")))
  
  if(any(ctrl_genes %in% rownames(score_GeneMatrix))){
    message(sprintf("Overwrite transcript score for %d `ctrl_genes` shared with `refProfiles`: `%s`.", 
                    sum(ctrl_genes %in% rownames(score_GeneMatrix)),
                    paste0(intersect(ctrl_genes, rownames(score_GeneMatrix)), collapse = "`, `")))
    
    score_GeneMatrix <- score_GeneMatrix[!(rownames(score_GeneMatrix) %in% ctrl_genes), ]
  }
  
  score_GeneMatrix <- rbind(score_GeneMatrix, 
                             matrix(svmClass_score_cutoff, 
                                    nrow = length(ctrl_genes), ncol = ncol(score_GeneMatrix),
                                    dimnames = list(ctrl_genes, colnames(score_GeneMatrix)))
  )
  
  
}

## (7.2) get file path to transcript data.frame
files_flagged_transDF <- dir(path = sub_out_dir, pattern = "^[0-9]+_flagged_transDF.csv", full.names = TRUE)

## (7.3) flag transcript groups of low goodness-of-fit based on the provided cutoffs and flagged cells
# output folder for new data
path_to_output <- fs::path(sub_out_dir, "newCutoff")
if(!dir.exists(path_to_output)){
  dir.create(path_to_output, recursive = T)
}


## define function for each file of transcript data.frame
myFun_flagTranscriptsSVM <- function(eachTransDF_path){
  # get idx from file name
  idx <- sapply(unlist(stringr::str_extract_all(basename(eachTransDF_path), "[:digit:]+_flagged_transDF")), function(x){
    as.numeric(gsub("_flagged_transDF", "", x))})
  
  message(sprintf("\n##############\nProcessing file `%d`: %s\n\n\n",
                  idx, eachTransDF_path))
  
  # load transcript data.frame
  eachTransDF <- read.csv(eachTransDF_path, header = TRUE)
  
  # subset to focus on flagged cells
  classDF_ToFlagTrans <- eachTransDF[eachTransDF[['UMI_cellID']] %in% combined_flaggedCells,]
  # remove original SVM class
  classDF_ToFlagTrans[['SVM_class']] <- NULL
  
  # perform SVM on flagged cells to identified transcript groups of low score
  tmp_df <- flag_bad_transcripts(chosen_cells = combined_flaggedCells,
                                score_GeneMatrix = score_GeneMatrix,
                                transcript_df = classDF_ToFlagTrans, 
                                cellID_coln = 'UMI_cellID', 
                                transID_coln = 'UMI_transID', 
                                score_coln = 'score_tLLRv2_maxCellType',
                                spatLocs_colns = c('x','y','z')[1:d2_or_d3], 
                                model_cutoff = flagModel_TransNum_cutoff, 
                                score_cutoff = svmClass_score_cutoff, 
                                svm_args = svm_args)
  
  # add in SVM results to flagged transcript, cells with all transcript score on same class are removed
  classDF_ToFlagTrans <- merge(classDF_ToFlagTrans, 
                               as.data.frame(tmp_df)[, c('UMI_transID','DecVal','SVM_class','SVM_cell_type')], 
                               by = 'UMI_transID')
  rm(tmp_df)
  
  # write into disk
  write.csv(classDF_ToFlagTrans, file = fs::path(path_to_output, paste0(idx, '_classDF_ToFlagTrans.csv')), row.names = FALSE)
  
  
  # flagged transcript ID, character vector
  flaggedSVM_transID3d <- classDF_ToFlagTrans[classDF_ToFlagTrans[['SVM_class']] ==0, 'UMI_transID']
  # assign SVM_class =0 for transcripts with low goodness-of-fit
  eachTransDF[['SVM_class']] <- 1- as.numeric(eachTransDF[['UMI_transID']] %in% flaggedSVM_transID3d)
  
  # save `updated_transDF` into csv file for each FOV 
  write.csv(eachTransDF, 
            file = fs::path(path_to_output, paste0(idx, "_flagged_transDF.csv")), 
            row.names = FALSE)
  
  # return only idx, file path as a data.frame
  res_to_return <- data.frame(file_idx = idx, 
                              classDF_ToFlagTrans = fs::path(path_to_output, paste0(idx, '_classDF_ToFlagTrans.csv')), 
                              flagged_transDF = fs::path(path_to_output, paste0(idx, "_flagged_transDF.csv")))
  
  return(res_to_return)
}


# processing each FOV in parallel
flagTrans_outputs <- parallel::mclapply(X = files_flagged_transDF, 
                                      mc.allow.recursive = TRUE,
                                      mc.cores = numCores(percentCores = 0.75),
                                      FUN = myFun_flagTranscriptsSVM)

flagTrans_outputs <- do.call(rbind, flagTrans_outputs)


## (7.4) get the per cell expression after trimming for the updated flagged transDF
trimmed_perCellExprs <- parallel::mclapply(
  X = flagTrans_outputs$flagged_transDF, 
  mc.allow.recursive = TRUE,
  mc.cores = numCores(percentCores = 0.75),
  FUN = function(flagTransDF_path){
    flagTransDF <- read.csv(flagTransDF_path)
    
    # get gene x cell count matrix from transcript data.frame
    res_to_return <- FastReseg::transDF_to_perCell_data(
      transcript_df = data.table::as.data.table(flagTransDF)[SVM_class == 1, ], 
      transGene_coln = 'target',
      cellID_coln = "UMI_cellID",
      spatLocs_colns = c('x','y','z'), 
      celltype_coln = "SVM_cell_type",
      return_cellMeta = FALSE)
    
    ## impute zero value for genes not in `res_to_return$perCell_expression` but in `score_GeneMatrix` 
    missingGenes <- setdiff(rownames(score_GeneMatrix), 
                            rownames(res_to_return$perCell_expression))
    if(length(missingGenes)>0){
      mockExprs <- matrix(0, nrow = length(missingGenes), 
                          ncol = ncol(res_to_return$perCell_expression), 
                          dimnames = list(missingGenes, 
                                          colnames(res_to_return$perCell_expression)))
      mockExprs <- Matrix::Matrix(mockExprs, sparse = TRUE)
      res_to_return$perCell_expression <- rbind(res_to_return$perCell_expression, 
                                                mockExprs)
    }
    return(res_to_return$perCell_expression)
  })

trimmed_perCellExprs <- do.call(cbind, trimmed_perCellExprs)

