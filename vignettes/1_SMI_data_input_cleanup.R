#### pipeline to process multi-FOV multi-slot data
## part 1:
## (1.1) define config_dimension and config_spatNW2 used for resegmentation; 
## (1.2) reformat SMI data to have transcript data.frame using same cell_ID as SMITAP gem object
## (1.3) get mean cell type specific profile from existing SMITAP gem object

library(ggplot2)
library(FastReseg)

main_output_dir <- '/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/NSCLC/mock_TAP_1M_NSCLC'
sub_out_dir2 <- fs::path(main_output_dir, "FullSet_1Mcell_NSCLC_resegment")
if(!dir.exists(sub_out_dir2)){
  dir.create(sub_out_dir2, recursive = T)
}

setwd(sub_out_dir2)


#black list genes not present in ISH panel
if(TRUE){
  blacklist_genes <- fs::path("/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/melanoma/Run4104_melanoma_980plx/giotto_output/Run4104_redo_cellpose", 
                              "20201210_33HE_GeneName_in_1013plx.txt")
  blacklist_genes <- read.csv(blacklist_genes, header = TRUE, sep = '\t')
  blacklist_genes <- blacklist_genes[['GeneName']]
}


# spatial cutoff for the search, stitch at mm, then convert to um for analysis
config_dimension <- list(mm_per_pixel = 0.18/1000, 
                         um_per_pixel = 0.18)
# z step in mm
config_dimension[['zstep_mm']] = 0.8/1000
config_dimension[['zstep_um']] = 0.8
# distance between proximity transcripts, "auto" or absolute distance = 15 pixel = 2.7um
config_dimension[['TransDistance_cutoff']] = 2.7

# neighbor_distance = 100 pixel = 18um in xy, 4 z-step in z to select neighbor transcripts if search based on cell network
config_dimension[['CellNeighbor_xy']] = 18
config_dimension[['CellNeighbor_z']] = config_dimension[['zstep_um']]*4 
# when search transcript data.frame directly, use smaller search range for xy
config_dimension[['CellNeighbor_xy_in_transDF']] = config_dimension[['CellNeighbor_xy']]/3


# lrtest_-log10P cutoff for flagging wronlgy segmented cells
config_dimension[['flagCell_lrtestCutoff']] = 10 
# high and low tLLRv2 score cutoff for SVM
config_dimension[['tLLRv2_SVMcutoff']] = -2 

# config for low-score group identification
config_dimension[['svm_config']] <- list(kernel = "radial", 
                                         scale = TRUE, 
                                         gamma =0.8)

# config for leiden clustering for merging check in re-segmentation action
config_dimension[['leiden_config']] <- list(python_path = "/usr/bin/python3", 
                                            partition_type = "RBConfigurationVertexPartition",
                                            resolution =1,
                                            n_iterations = 1000,
                                            set_seed = T,
                                            seed_number = 1234)
config_dimension[['cutoff_sharedLeiden']] = 0.5



## record cell number and transcript number information
reseg_logInfo <- list()

#### load sample info and cell typing information ----
# SMI TAP output folder
scAnalysis_dir <- '/home/rstudio/smiqumulo/RnD/08. Mock TAP/02. NSCLC_1K/6.92 Analysis showcase full set - hades backup'
# get config files to get path to sample annotation file
source(fs::path(scAnalysis_dir, 'config_v1.2.R'))
annotfile <- read.csv(config_loading$annotfile)

# flag to remove black list genes for 980plex
blacklist_flag <- any(grepl('980', annotfile[[config_loading$panelColumn]]))

# load the giotto object containing cell typing and per cell info
load(fs::path(config_main$resultspath, 'results/complete_giotto_object.RData'))
cell_metadata <- Giotto::pDataDT(gem)
# cell_ID = c_[slide]_[FOV]_CellId
celltype_metadata <- data.frame(cell_ID = cell_metadata$cell_ID, 
                                FOV = cell_metadata$fov, 
                                slide = cell_metadata$slide_ID_numeric, 
                                slideName = cell_metadata[[config_loading$slidenameColumn]],
                                cell_type = cell_metadata$nb_clus)
celltype_metadata[['CellId']] <- sapply(strsplit(celltype_metadata$cell_ID, '_'),'[[',4)

# get rna target only raw expression matrix
exprs_tgrt <- Giotto:::select_expression_values(gem, feat_type = 'rna', values = 'raw')
targets <- rownames(exprs_tgrt)
if(blacklist_flag){
  targets <- setdiff(targets, blacklist_genes)
  exprs_tgrt <- exprs_tgrt[targets, ]
}

# get cell level spatial network
cell_delaunayNW_Obj <- Giotto:::select_spatialNetwork(gem, name = 'Delaunay_network', return_network_Obj = TRUE)
cellNetWorkDT <- cell_delaunayNW_Obj$networkDT                   

#### get fov position information ----
## slide_ID_numeric in cell_metadata and cell_ID is ordered by sample_annot file
slide_converter <- unique(celltype_metadata[, c('slide','slideName')])

## but UID in gem@parameters arrangeFOV step is ordered by config_loading$slidenameColumn
slide_converter <- slide_converter[order(slide_converter$slideName), ]
slide_converter[['slideID_position']] <- seq_len(nrow(slide_converter))


# get conversion formula from raw data to stitched complete gem
config_fovs <- gem@parameters[[which(grepl('_fov', names(gem@parameters)))[1]]]
fov_position <- data.frame(x = config_fovs$fov_position.x,
                           y = config_fovs$fov_position.y,
                           fov = config_fovs$fov_position.fov,
                           slideID_position = config_fovs$fov_position.slide,
                           UID = config_fovs$fov_position.UID)
# add in slide_ID_numeric and slideName
fov_position <- merge(fov_position, slide_converter, by = 'slideID_position')


# remove additional information to free some space
rm(list = setdiff(ls(pattern = '^config_'), c('config_loading','config_dimension','config_fovs')))
rm("gem")



#### config for spatial network for transcripts: ---------------------------
config_spatNW2 <- list(
  # name for spatial network (default = 'Delaunay_network' or 'kNN_network' for method = 'Delaunay' or 'kNN', respectively if NULL)
  name = 'transcript_delaunay_network',
  # which spatial dimensions to use (default = all)
  dimensions = "all",
  # which method to use to create a spatial network, choose from c("Delaunay", "kNN"). (default = Delaunay)
  method = 'Delaunay',
  # minimum number of nearest neighbors if maximum_distance != NULL; used by both Delaunay and kNN methods
  minimum_k = 0,
  
  #### Approach 1: create Delanuay network for HMRF module
  # Delaunay method to use, choose from c("deldir", "delaunayn_geometry", "RTriangle")
  delaunay_method = "delaunayn_geometry",
  # distance cuttoff for nearest neighbors to consider for Delaunay network. 
  maximum_distance_delaunay = "auto",
  
  ## only for delaunay_method = "delaunayn_geometry"
  # (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (http://www.qhull.org/html/qdelaun.htm) for the available options. (default = 'Pp', do not report precision problems)
  options = "Pp",
  
  ## only for delaunay_method = "RTriangle"
  # (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
  Y = TRUE,
  # (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
  j = TRUE,
  # (RTriangle) Specifies the maximum number of added Steiner points.
  S = 0,
  
  #### Approach 2: create kNN network for leiden clustering
  
  # method to create kNN network
  knn_method = "dbscan",
  # number of nearest neighbors based on physical distance
  k = 4,
  # distance cuttoff for nearest neighbors to consider for kNN network
  maximum_distance_knn = NULL
)


#### get refernce cell typing profiles from data, use expression matrix from gem direclty ----
targets <- rownames(exprs_tgrt)
select_cellmeta <- celltype_metadata
reseg_logInfo[['reference']] <- list(genes_included = targets, 
                                     feat_count = data.frame(cellNum = nrow(select_cellmeta)))

# Given cell assignments (or posterior probabilities), estimate the mean profile of each cluster.
# ignore background, use total count per cell as scaling factor
meanCelltype_profiles <- estimate_MeanProfile(counts = as.matrix(Matrix::t(exprs_tgrt[targets, select_cellmeta$cell_ID])),
                                              clust = as.character(select_cellmeta$cell_type), 
                                              s = Matrix::colSums(exprs_tgrt[targets, select_cellmeta$cell_ID]), 
                                              bg = rep(0,nrow(select_cellmeta)))
# save mean profiles for future usage on lung cancer
save(meanCelltype_profiles, file = fs::path(sub_out_dir2, "NSCLC_1Mcell_980plx_meanCelltype_profiles.RData"))

# remove exprs_tgrt and cell_metadata to save space
rm(exprs_tgrt, cell_metadata)

####similarity between clusters in reference profiles --------------------
# calculate correlation and p-values
res <- Hmisc::rcorr(meanCelltype_profiles, type = 'pearson')
# color scheme
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))

pdf(fs::path(sub_out_dir2,"nbclust6_ref_profiles_similarity.pdf"), 
    width = 6, height = 6)
par(mar = c(1,1,1,1))

corrplot::corrplot.mixed(res$r, upper = "color", lower = "number",
                         order = "hclust",lower.col = "black", 
                         number.cex = 0.5, tl.cex = 0.3,
                         p.mat = res$P, sig.level = 0.01, 
                         title = "Pearson correlation of reference profiles",
                         mar = c(0,0,1,0))

# no order, plot number
corrplot::corrplot(res$r, method = "number", type = "upper", 
                   number.cex = 0.5, tl.cex = 0.3, col = col1(20),
                   p.mat = res$P, sig.level = 0.01,
                   title = "Pearson correlation of reference profiles", 
                   mar = c(0,0,1,0))
dev.off()



#### score each transcript with respect to reference profiles, exclude "NotDet"  ---------------------------
# replace zero in mean profiles with 1E-5
meanCelltype_profiles <- pmax(meanCelltype_profiles, 1e-5)
# tLL score
transcript_loglik <- scoreGenesInRef(genes = targets, ref_profiles = meanCelltype_profiles)
# remove NotDet
transcript_loglik_cleaned <- transcript_loglik[, which(colnames(transcript_loglik) != 'NotDet')]

# tLLRv2 score, re-center on maximum per row/transcript
tmp_max <- apply(transcript_loglik_cleaned, 1, max)
tLLRv2_geneMatrix_cleaned <- sweep(transcript_loglik_cleaned, 1, tmp_max, '-')





#### break the data set into different files and processing only 1 at a time ----
#### the complete all_transDF is too big to create transcript level data.frame inside function
system.time(save(cell_delaunayNW_Obj, #cell delaunay network
                 celltype_metadata, # cell typing of gem
                 annotfile, # sample annotation file for loading
                 config_loading, # config for loading raw data for gem
                 slide_converter, # match slide name with slide ID
                 config_fovs, # stitch fovs
                 blacklist_flag, # flag for blacklisting
                 # unique data for re-segmentation
                 reseg_logInfo, # count information of raw vs. gem data
                 config_dimension, # spatial dimension of the search and images
                 config_spatNW2, # config for transcript spatial network
                 # all_transDF, # transcript metadata from raw
                 # all_cellstatsDF, # cell stats from raw
                 meanCelltype_profiles, # mean expression profiles
                 transcript_loglik_cleaned, # tLL Gscore, gene x cell type
                 tLLRv2_geneMatrix_cleaned, # recenter net Gscore based on maximum per row/transcript, gene x cell type
                 file = fs::path(sub_out_dir2, "NSCLC_1Mcell_980plx_allBlocks_data.RData")))
# user  system elapsed 
# 6.247   3.135  73.866


#### get raw data information for all FOVs, assigned new cell ID and new coordinates to be same as giotto object ----
## 8 slides; 239 total FOVs; 771,236 cells; 269,176,112 transcripts; 
## 36.4GB of all_transDF; 63.1MB of all_cellstatsDF; 1.2GB of exprs_tgrt; 474MB of cell_delaunayNW_Obj


#### seems too big to compile, continue to have problem when reading in data
#### break into 3 blocks for loading
## same order as annote file
slide_converter <- slide_converter[order(slide_converter$slide), ]

# others, lung9, lung5
lung9 <- which(grepl('Lung9', slide_converter$slideName))
lung5 <- which(grepl('Lung5', slide_converter$slideName))
others <- setdiff(seq_len(nrow(slide_converter)), c(lung9, lung5))

block_list <- list(others = slide_converter[others, ], 
                   lung9 = slide_converter[lung9, ], 
                   lung5 = slide_converter[lung5, ])

# look over block to load data and save
for(blockID in names(block_list)){
  sub_out_dir3 <- fs::path(sub_out_dir2, paste0("Block_", blockID))
  if(!dir.exists(sub_out_dir3)){
    dir.create(sub_out_dir3, recursive = T)
  }
  
  currentBlock <- block_list[[blockID]]
  
  all_transDF <- list()
  all_cellstatsDF <- list()
  ## notice some additional Cell_Stats file in folder, need to search just the root folder 
  ## fov_position in gem@parameters have more FOVs than fov in gem@cell_metadata
  ## seems voted file missing for some fovs
  for (idx in currentBlock$slide){
    # for fov position only
    slide_position <- slide_converter$slideID_position[which(slide_converter$slideName == annotfile[[config_loading$slidenameColumn]][idx])]
    # for correct cell_ID name matching gem
    slide_Number <- slide_converter$slide[which(slide_converter$slideName == annotfile[[config_loading$slidenameColumn]][idx])]
    
    datafolder <- fs::path(annotfile[[config_loading$folderpathColumn]][idx], 
                           annotfile[[config_loading$slidefoldersColumn]][idx], 
                           annotfile[[config_loading$votedfoldersColumn]][idx])
    exprs_files_list <- dir(path = datafolder, full.names = TRUE, recursive = TRUE, 
                            pattern = "^Run[0-9]+_FOV[0-9]+__complete_code_cell_target_call_coord.csv")
    exprs_FOV <- sapply(unlist(stringr::str_extract_all(basename(exprs_files_list), "FOV[:digit:]+")), 
                        function(x){as.numeric(gsub("FOV", "", x))})
    
    cellstatsfolder <- fs::path(annotfile[[config_loading$folderpathColumn]][idx], 
                                annotfile[[config_loading$slidefoldersColumn]][idx], 
                                'CellStatsDir')
    # only search root folder 
    cellstats_files_list <- dir(path = cellstatsfolder, full.names = TRUE, recursive = FALSE, 
                                pattern = "^Run[0-9]+_[0-9]+_[0-9]+[_]*.*_Cell_Stats_F[0-9]+.csv")
    annot_FOV <- sapply(unlist(stringr::str_extract_all(basename(cellstats_files_list), "F[:digit:]+")), 
                        function(x){as.numeric(gsub("F", "", x))})
    
    # pair expression files with cell stats files
    exprs_files_list <- exprs_files_list[order(exprs_FOV)]
    exprs_FOV <- sort(exprs_FOV)
    cellstats_files_list <- cellstats_files_list[order(annot_FOV)]
    annot_FOV <- sort(annot_FOV)
    
    if(!identical(unname(exprs_FOV), unname(annot_FOV))){
      printout_msg <- sprintf("For slide %s:\nCell stats files are missing for %s.\nExpression files are missing for %s. ", 
                              annotfile[[config_loading$slidenameColumn]][idx], 
                              paste0(sprintf("FOV%s", setdiff(exprs_FOV, annot_FOV)), collapse = ","), 
                              paste0(sprintf("FOV%s", setdiff(annot_FOV, exprs_FOV)), collapse = ","))
      message(printout_msg)
      
      # remove unpaired data
      if ( length(intersect(exprs_FOV, annot_FOV))==0 ){
        stop("No paired expression files and anntation files are found.")
      } else {
        exprs_files_list <- exprs_files_list[which(exprs_FOV %in% annot_FOV)]
        cellstats_files_list <- cellstats_files_list[which(annot_FOV %in% exprs_FOV)]
        exprs_FOV <- exprs_FOV[which(exprs_FOV %in% annot_FOV)]
        annot_FOV <- annot_FOV[which(annot_FOV %in% exprs_FOV)]
      }
    }
    
    # read in individual FOVs and get new cell_ID
    perSlide_TransDF <- list()
    perSlide_CellsDF <- list()
    for (i in seq_len(length(exprs_FOV))){
      # current FOV locations
      abs_locs <- unlist(fov_position[which(fov_position$UID==paste0(slide_position, '_', exprs_FOV[i])), c("x","y")])
      # skip the files if no fov_position information in gem 
      if(length(abs_locs) == 0){
        message(sprintf("No fov postioin found in gem, skip loading for UID = `%s`. ", 
                        paste0(slide_position, '_', exprs_FOV[i])))
      }else{
        # get just needed columns and add FOV, slideID columns; not include CellComp since not present in oldDASH data
        tmp_df <- read.csv(exprs_files_list[[i]])
        tmp_df <- data.table::setDT(tmp_df[, c('fov','x','y','z','target', 'CellId', 'CellComp')])
        tmp_df[, cell_ID := paste0('c_', slide_Number,'_', fov,'_', CellId)]
        tmp_df[, transcript_id := paste0('t_', slide_Number,'_', fov,'_', .I)]
        tmp_df[, slide := slide_Number]
        
        # update coordinates to be the same as gem
        raw_locs <- tmp_df[, .SD, .SDcols = c('x','y','z')]
        # flip y coordinates
        raw_locs[['y']] <- -raw_locs[['y']]
        raw_locs[, 1:2] <- sweep(raw_locs[, 1:2] * config_fovs$pixel_size, 2, abs_locs,"+")
        raw_locs[, 3] <- raw_locs[, 3]*config_dimension$zstep_mm
        # convert to um
        tmp_df[, match(colnames(raw_locs), colnames(tmp_df))] <- raw_locs*1000
        perSlide_TransDF[[exprs_FOV[i]]] <- data.table::copy(tmp_df)
        
        # get cell stats info
        tmp_df <- data.table::setDT(read.csv(cellstats_files_list[[i]]))
        tmp_df[, cell_ID := paste0('c_', slide_Number,'_', annot_FOV[i],'_', CellId)]
        
        raw_locs <- tmp_df[, .SD, .SDcols = c('CenterX','CenterY')]
        raw_locs[['CenterY']] <- -raw_locs[['CenterY']]
        raw_locs[, 1:2] <- sweep(raw_locs[, 1:2] * config_fovs$pixel_size, 2, abs_locs,"+")
        # convert to um
        tmp_df[, match(colnames(raw_locs), colnames(tmp_df))] <- raw_locs *1000
        tmp_df[, slide := slide_Number]
        tmp_df[, fov := annot_FOV[i]]
        # no change to the Area, Width and Height of cell stats file
        perSlide_CellsDF[[annot_FOV[i]]] <- data.table::copy(tmp_df)
      }
      
    }
    all_transDF[[idx]] <- do.call(rbind, perSlide_TransDF)
    all_cellstatsDF[[idx]] <- do.call(rbind, perSlide_CellsDF)
  }
  
  
  all_transDF <- do.call(rbind, all_transDF)
  # fill NA for missing columns
  all_cellstatsDF <- do.call(plyr::rbind.fill, all_cellstatsDF)
  reseg_logInfo[['raw_data']] <- list(feat_count = data.frame(cellNum = nrow(all_cellstatsDF), 
                                                              transNum = nrow(all_transDF)))
  
  ## check out the fov number remainin for each slide
  print(block_list[[blockID]])
  print(aggregate(all_cellstatsDF$fov, by = list(all_cellstatsDF$slide), max))
  
  # cells in current block
  common_cells <- intersect(celltype_metadata$cell_ID, unique(all_transDF$cell_ID))
  select_cellmeta <- celltype_metadata[which(celltype_metadata$cell_ID %in% common_cells), ]
  
  block_cellNetWorkDT <- cellNetWorkDT[which(cellNetWorkDT$from %in% common_cells & cellNetWorkDT$to %in% common_cells), ]
  
  # 39.24min to save 1 block of 3 slides
  system.time(save(currentBlock, # slide, slide name, slide position of current block
                   block_cellNetWorkDT, #cell delaunay network data table only
                   celltype_metadata, # cell typing of gem
                   select_cellmeta, # cell meta data for currentblock
                   annotfile, # sample annotation file for loading
                   config_loading, # config for loading raw data for gem
                   slide_converter, # match slide name with slide ID
                   config_fovs, # stitch fovs
                   blacklist_flag, # flag for blacklisting
                   # unique data for re-segmentation
                   reseg_logInfo, # count information of raw vs. gem data
                   config_dimension, # spatial dimension of the search and images
                   config_spatNW2, # config for transcript spatial network
                   all_transDF, # transcript metadata from raw
                   all_cellstatsDF, # cell stats from raw
                   meanCelltype_profiles, # mean expression profiles
                   transcript_loglik_cleaned, # tLL Gscore, gene x cell type
                   tLLRv2_geneMatrix_cleaned, # recenter net Gscore based on maximum per row/transcript, gene x cell type
                   file = fs::path(sub_out_dir3, paste0("NSCLC_1Mcell_980plx_Block_", blockID,"_data.RData"))))
  # user   system  elapsed 
  # 206.850  100.190 2354.862 
  
}


rm(perSlide_CellsDF, perSlide_TransDF, raw_locs)


