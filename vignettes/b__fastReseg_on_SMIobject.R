#### pipeline to process multi-FOV multi-slide data
## The pipeline assumes the pre-existing of single cell typing results on original cell segmentation outcomes. 
## The pre-built object is expected to have same data structure as the ones built with `Giotto 2.0.0.9021`. 
## For pre-built object with different data structure, please refer to `vignettes\tutorial.Rmd` for preprocessing steps. 
## Firstly, calculate reference profiles and cutoffs for transcript number and transcript tLLR score, `cellular_distance_cutoff` on entire data set;
## Secondly, check the existence of fov offset position and transcript data.frame for each FOV in each slide; 
## Thirdly, estimate `molecular_distance_cutoff` that defining the neighborhood from 1st FOV of 1st slide. 
## Fourthly, loop through each FOV data.frame and perform re-segmentation using same cutoffs defined earlier. 

library(ggplot2)
library(FastReseg)

# install correct version of Giotto for SMITAP object
if(!requireNamespace("Giotto")){
  remotes::install_github("drieslab/Giotto/commit/0c8c2866e881b1c6b35ddc97c24dcb58b555c375")
} else if(packageVersion("Giotto") != "2.0.0.9021"){
  remotes::install_github("drieslab/Giotto/commit/0c8c2866e881b1c6b35ddc97c24dcb58b555c375")
}


# helper function for extracting data from pre-built Giotto object
source(fs::path(system.file("extdata", package = "FastReseg"), "SMITAP_interface.R"))

main_output_dir <- '/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/SMITAP_proj/SMI003_DavidTing_MGH'
sub_out_dir <- fs::path(main_output_dir, "fastReseg01")
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
removeUnpaired <- TRUE

# flag to include ctrl_genes in analysis
include_ctrlgenes <- FALSE

# SMI TAP output folder
scAnalysis_dir <- "/home/rstudio/smiqumulo/01 SMI TAP project/SMI-0003_DavidTing_MGH/6.4 Analysis combined"
# get config files to get path to sample annotation file
source(fs::path(scAnalysis_dir, 'config_v1.2.R'))

# path to giotto object containing cell typing and per cell info
path_to_SMIobject <- fs::path(config_main$resultspath, 'results/complete_giotto_object.RData')

rm(list = setdiff(ls(pattern = '^config_'), c('config_loading')))

# optional: exclude 33 high expressors from reference profile estimation in case of 980plx ISH panel
is_980plx <- FALSE
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

# save smi_inputs to disk
save(smi_inputs, file = fs::path(sub_out_dir, "smi_inputs.RData"))
gc()

#### (3) other parameters used in resegmentation workflow ----
# ## use default values for those parameters, change them as needed below and then pass to `fastReseg_internalRef` function
# # cutoff of transcript number to do spatial modeling for identification of wrongly segmented cells (default = 50)
# flagModel_TransNum_cutoff = 50 
# 
# # cutoff of lrtest_nlog10P to identify putative wrongly segemented cells with strong spatial dependency in transcript score profile
# flagCell_lrtest_cutoff = 5
# 
# # cutoff of transcript score to separate between high and low score transcripts in SVM (default = -2)
# svmClass_score_cutoff = -2 
# 
# # a list of arguments to pass to svm function for identifying low-score transcript groups in space, typically involve kernel, gamma, scale
# svm_args = list(kernel = "radial", 
#                 scale = FALSE, 
#                 gamma = 0.4)
# 
# # a list of configuration to pass to reticulate and `igraph::cluster_leiden` function, including objective_function, resolution_parameter, beta, n_iterations.
# leiden_args = list(objective_function = c("CPM", "modularity"),
#                    resolution_parameter = 1,
#                    beta = 0.01,
#                    n_iterations = 200) 
# 
# # minimal percentage of transcripts shared membership between query cell and neighbor cells in leiden clustering results for a valid merging event, default = 0.5 for 50% cutoff
# flagMerge_sharedLeiden_cutoff = 0.5
# 
# # percentage of cores used for parallel processing; try lower percentage if processing per FOV files with larger file size or running into memory issue
# percentCores = 0.75

#### (4) loop through each FOV to do resegmentation on pre-defined refProfiles, cutoffs ----
## open a log file to store all function outputs into disk
logFile <- file.path(sub_out_dir, "fastReseg_log.txt")
logHandle <- file(logFile, open = "wt")
sink(logHandle)
sink(logHandle, type = "message")

## run resegmentation
# # `fastReseg_internalRef` function returns a list of 
# refProfiles: a genes * clusters matrix of cluster-specific reference profiles used in resegmenation pipeline.
# baselineData: a list of two matrice in cluster * percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.
# cutoffs_list: a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`.
# ctrl_genes: a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not NULL.
# updated_perCellDT: a per cell data.table with mean spatial coordinates and new cell type after resegmentation, return when return_perCellData = TRUE.
# updated_perCellExprs: a gene x cell count sparse matrix for updated transcript data.frame after resegmentation, return when return_perCellData = TRUE.
# reseg_actions: a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, save when save_intermediates = TRUE.

# For automatic calculation of `molecular_distance_cutoff` based on 1st FOV, set `molecular_distance_cutoff` to NULL
# But for SMI dataset, it's recommended to use `molecular_distance_cutoff` = 2.7 um as hard cutoff given the footprint of fiducials

if(include_ctrlgenes){
  # if ctrl_genes = smi_inputs$ctrl_genes, the control genes like NegPrbs and FalseCodes would be included in analysis.
  ctrl_genes <- smi_inputs$ctrl_genes
} else {
  ctrl_genes <- NULL
}


reseg_outputs <- fastReseg_full_pipeline(
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
  
  molecular_distance_cutoff = 2.7, # recommend to use footprint of fiducials as cutoff in SMI dataset; if NULL, automatic calculation based on 1st FOV
  cellular_distance_cutoff = smi_inputs[['cellular_distance_cutoff']], # if NULL, automatic calculation based on 1st FOV
  score_baseline = smi_inputs[['score_baseline']],
  lowerCutoff_transNum = smi_inputs[['lowerCutoff_transNum']],
  higherCutoff_transNum= smi_inputs[['higherCutoff_transNum']],
  imputeFlag_missingCTs = TRUE,
  
  groupTranscripts_method = 'delaunay', # use 'delaunay' network to group low-score transcripts in space, option to use 'dbscan' instead
  
  path_to_output = sub_out_dir, 
  combine_extra = FALSE,  # if TRUE, extracellular and trimmed transcripts are included in the updated transcript data.frame
  ctrl_genes = ctrl_genes
  )


## revert output back to the console -- only then access the file!
sink(type = "message")
sink()

# # show the log file in pop up window
# file.show(logFile)

# print all lines in log file on console
cat(readLines(logFile), sep="\n")


#### (5) create new SMI giotto object with updated gene x cell count matrix, and per cell metadata file ----
## genes not presented in the reference profiles would not be included in the transcript data.frame with updated cell segmentaion.
## need a separate function to extract cell borders in space from updated cell ID assignment wihtin transcript data.frame,
## and then apply to the original transcript data.frame to assign cell_ID for those genes mising from reference profiles. 

## (5.1) prepare expression list for multi-slot giotto object
# when use ctrl_genes = NULL, 
# NegPrb and FalseCode genes are absence from refernece profiles and thus missing in the transcript data.frame with updated cell segmentaion. 
# assume 0 values for both NegPrb and FalseCode genes for now
if(is.null(ctrl_genes)){
  expr_lists <- list(rna = list(raw = reseg_outputs[['updated_perCellExprs']]))
  expr_lists[['negprobes']] <- list(raw = matrix(0, nrow = 5,
                                                 ncol = ncol(reseg_outputs[['updated_perCellExprs']]),
                                                 dimnames = list(paste0('NegPrb', seq_len(5)),
                                                                 colnames(reseg_outputs[['updated_perCellExprs']]))))
  expr_lists[['falsecode']] <- list(raw = matrix(0, nrow = 5,
                                                 ncol = ncol(reseg_outputs[['updated_perCellExprs']]),
                                                 dimnames = list(paste0('FalseCode', seq_len(5)),
                                                                 colnames(reseg_outputs[['updated_perCellExprs']]))))
} else {
  tmpExp <- reseg_outputs[['updated_perCellExprs']]
  
  expr_lists <- list(rna = list(raw = tmpExp[!(rownames(tmpExp) %in% ctrl_genes), ]))
  
  expr_lists[['negprobes']] <- list(raw = tmpExp[rownames(tmpExp) %in% ctrl_genes[grepl('NegPrb', ctrl_genes)], ])
  expr_lists[['falsecode']] <- list(raw = tmpExp[rownames(tmpExp) %in% ctrl_genes[grepl('FalseCode', ctrl_genes)], ])
  
  rm(tmpExp)
  
}


## (5.2) prepare cell annotation file
cell_annotDF <- as.data.frame(reseg_outputs[['updated_perCellDT']])
colnames(cell_annotDF) <- c('cell_ID', paste0('reSeg_', refClus_coln), 'CenterX', 'CenterY', 'CenterZ','reSeg_action')
cell_annotDF[['slide']] <- as.numeric(sapply(strsplit(cell_annotDF[['cell_ID']], '_'),'[[',2))
cell_annotDF[['fov']] <- as.numeric(sapply(strsplit(cell_annotDF[['cell_ID']], '_'),'[[',3))

## add in meta data for each slide that is stored in `config_loading$annotfile`
cell_annotDF <- merge(cell_annotDF, 
                      smi_inputs[['sample_annot']], 
                      by.x = 'slide', by.y = 'slide_ID_numeric', all.x = TRUE)

# rearrange cell order to be the same as expression matrix
cell_annotDF <- cell_annotDF[match(colnames(reseg_outputs[['updated_perCellExprs']]), cell_annotDF[['cell_ID']]),]

## (5.3) create giotto object

updated_SMIobj <- Giotto::createGiottoObject(expression= expr_lists,
                                             expression_feat= names(expr_lists),
                                             spatial_locs= cell_annotDF[, c('CenterX', 'CenterY')],
                                             cell_metadata= list(rna = cell_annotDF, 
                                                                 negprobes = cell_annotDF, 
                                                                 falsecode = cell_annotDF))

### save resegmenation inputs, outcomes and new giotto object on file
save(smi_inputs, 
     reseg_outputs,
     updated_SMIobj,
     file = fs::path(sub_out_dir, 'fastReseg01_outputs_SMIobj.RData'))

#### (6) visualization of segmenation outputs ----
