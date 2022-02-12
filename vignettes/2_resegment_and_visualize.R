#### pipeline to process multi-FOV multi-slot data
## part 2:
## (2.0) load exisitng RData containing data.frame for all transcripts and cells and the reference cell type profiles
## (2.1) flag cells based on linear regression of transcript score using lrtest_-log10P cutoff
## (2.2) use SVM~hyperplane to identify the connected transcripts group based on transcript score
## (2.3) do network analysis on flagged transcript to split flagged transcript groups in space
## (2.4) re-segmentation of each flagged transcript groups based on their neighborhood

library(ggplot2)
library(FastReseg)


main_output_dir <- '/home/rstudio/NAS_data/lwu/testRun/SpatialTest/giotto_test/NSCLC/mock_TAP_1M_NSCLC'
sub_out_dir2 <- fs::path(main_output_dir, "FullSet_1Mcell_NSCLC_resegment")
setwd(sub_out_dir2)

# load previous block data
blockID = 'lung9'
sub_out_dir3 <- fs::path(sub_out_dir2, paste0("Block_", blockID))
load(fs::path(sub_out_dir3, paste0("NSCLC_1Mcell_980plx_Block_", blockID,"_data.RData")))

# use cell-to-cell distance to find neighborhood, search range = 25um in xy
config_dimension[['CellNeighbor_xy']] = 25
config_dimension[['CellNeighbor_xy_in_transDF']] = 25

# change flaging cutoff to 5
config_dimension[['flagCell_lrtestCutoff']] = 5

# change svm configuraiton
config_dimension[['svm_config']] <- list(kernel = "radial", 
                                         scale = FALSE, 
                                         gamma =0.4)

# create output folder for new configuration
sub_out_dir2 <- fs::path(main_output_dir, "NSCLC1m_lr5_svm0p4")
sub_out_dir3 <- fs::path(sub_out_dir2, paste0("Block_", blockID))
if(!dir.exists(sub_out_dir3)){
  dir.create(sub_out_dir3, recursive = T)
}
setwd(sub_out_dir2)


# extract info from saved RData
targets <- rownames(meanCelltype_profiles)
# cells in current block
common_cells <- intersect(celltype_metadata$cell_ID, unique(all_transDF$cell_ID))


#### start to process current block ----

# filter for the common cells and genes present in cell typing data and transcript data
all_genes <- unique(all_transDF$target)
negprobes <- all_genes[grepl('NegPrb',all_genes)]
falsecodes <- all_genes[grepl('FalseCode',all_genes)]

# 993 true genes, but only 960 genes in expression matrix for cell typing
targets <- intersect(targets, setdiff(all_genes, c(negprobes, falsecodes)))

reseg_logInfo[['common_data']] <- list(common_cells = common_cells, 
                                       common_tgrts = targets, 
                                       genes_omited = list(negprobes = negprobes, 
                                                           falsecodes = falsecodes), 
                                       feat_count = data.frame(cellNum_common = length(common_cells), 
                                                               transNum_tgrts = sum(all_transDF$target %in% targets)))

#### 
# > reseg_logInfo$raw_data$feat_count
# cellNum transNum
# 1  242476 60880312
# > reseg_logInfo$common_data$feat_count
# cellNum_common transNum_tgrts
# 1         227110       59447174

message(sprintf("%d common cells,  %.3f have cell type info; %d transcripts, %.3f are true targets. ", 
                reseg_logInfo$common_data$feat_count$cellNum_common, 
                reseg_logInfo$raw_data$feat_count$cellNum/reseg_logInfo$common_data$feat_count$cellNum_common, 
                reseg_logInfo$raw_data$feat_count$transNum, 
                reseg_logInfo$common_data$feat_count$transNum_tgrts/reseg_logInfo$raw_data$feat_count$transNum))
# 227110 common cells,  1.068 have cell type info; 
# 60880312 transcripts, 0.976 are true targets. 


#### block2 has 227,110 cells, 60,880,312 transcripts, 2 slides, 65 FOVs ----

# 1.6min to get new cell type based on maxium score
system.time(tmp_df <- getCellType_maxScore(score_GeneMatrix = tLLRv2_geneMatrix_cleaned, 
                                           transcript_df = all_transDF, 
                                           transID_coln = 'transcript_id',
                                           transGene_coln = "target",
                                           cellID_coln = "cell_ID", 
                                           return_transMatrix = FALSE))

# user  system elapsed 
# 431.748  61.822  96.555 

colnames(tmp_df[['cellType_DF']]) <- c('cell_ID','cleaned_tLLRv2_maxCellType')

## cells with no prior cell typing info in gem would be removed, due to failed QC 
select_cellmeta <- merge(select_cellmeta, tmp_df[['cellType_DF']], by = 'cell_ID')

# get score on assigned cell types at transcript level
## option 1: not used
#### keep all cells in all_transDF, update the cell_type of missing cells based on cleaned_tLLRv2_maxCellType
#### remove c_0 from data.frame
## option 2: only keep the transcripts have cell_type information in gem
all_transDF <- merge(all_transDF, 
                     select_cellmeta[, c('cell_ID','cell_type','cleaned_tLLRv2_maxCellType')], 
                     by = 'cell_ID')

# need to have NotDet in reference if calculate tLLRv2 score for orignla cell type
# only calculate tLLRv2 score based on new cell type, excluding 'NotDet'
# some transcripts (not targets) are missing tLLRv2 score and take long time to merge with transcript_id
tmp_df <- getScoreCellType_gene(score_GeneMatrix = tLLRv2_geneMatrix_cleaned, 
                                transcript_df = all_transDF, 
                                transID_coln = "transcript_id",
                                transGene_coln = "target",
                                celltype_coln = 'cleaned_tLLRv2_maxCellType')
all_transDF <- merge(all_transDF, tmp_df, by = 'transcript_id')

reseg_logInfo[['common_data']][['feat_count']][['cellNum']] <-  nrow(select_cellmeta)
reseg_logInfo[['common_data']][['feat_count']][['transNum']] <- nrow(all_transDF)

## print % cells and transcripts being removed from analysis
message(sprintf("Raw data has %d cells, %d transcripts; only %.4f cells have original cell type info and %.4f transcripts are targets in common cells.",
                reseg_logInfo[['raw_data']]$feat_count$cellNum, reseg_logInfo[['raw_data']]$feat_count$transNum,
                reseg_logInfo[['common_data']]$feat_count$cellNum/reseg_logInfo[['raw_data']]$feat_count$cellNum,
                reseg_logInfo[['common_data']]$feat_count$transNum/reseg_logInfo[['raw_data']]$feat_count$transNum))

# Raw data has 242476 cells, 60880312 transcripts; 
# only 0.9366 cells have original cell type info and 
# 0.7543 transcripts are targets in common cells.

#### linear regression on tLLRv2 to get wrongly segmented cells ----
# check p-value range and the corresponding score overlay segmentation
# 34043 cells = 17.6% cells with <50 transcripts; proceed with 193067 cells
# user   system  elapsed 
# 4528.534 3272.732 1590.532

if(TRUE){
  system.time(tmp_df <- score_cell_segmentation_error(chosen_cells = common_cells, 
                                                      transcript_df = all_transDF, 
                                                      cellID_coln = 'cell_ID', 
                                                      transID_coln = 'transcript_id', 
                                                      score_coln = 'score_cleaned_tLLRv2_maxCellType',
                                                      spatLocs_colns = c('x','y','z'), 
                                                      model_cutoff = 50))
  #-log10(P)
  tmp_df[['lrtest_-log10P']] <- -log10(tmp_df[['lrtest_Pr']])
  modStats_cleaned_tLLRv2_3D <- merge(tmp_df, select_cellmeta[, c('cell_ID','cleaned_tLLRv2_maxCellType','slide')], by = 'cell_ID')
  
  # visualize extreme cells, 2D plots
  if(TRUE){
    tmp_df <- data.table::copy(modStats_cleaned_tLLRv2_3D)
    tmp_df <- tmp_df[order(tmp_df[['lrtest_-log10P']]),]
    tmp_df[['labels']] <- paste0(tmp_df[['slide']],'_', tmp_df[['cleaned_tLLRv2_maxCellType']], 
                                 ', -log10P=', round(tmp_df[['lrtest_-log10P']],2))
    chosen_cells <- c(tmp_df$cell_ID[1:9], tmp_df$cell_ID[(nrow(tmp_df)-9):nrow(tmp_df)])
    fig1 <- plotSpatialScoreMultiCells(chosen_cells = chosen_cells, 
                                       cell_labels = tmp_df[match(chosen_cells, tmp_df$cell_ID), 'labels'], 
                                       transcript_df = all_transDF,
                                       cellID_coln = "cell_ID", 
                                       transID_coln = "transcript_id",
                                       score_coln = "score_cleaned_tLLRv2_maxCellType", 
                                       spatLocs_colns = c("x","y"))  
    
    pdf(fs::path(sub_out_dir3, paste0(blockID,"_SpatialPlot_SpatModel2_tLLRv2_-log10P_cleanNBclust6_flagExtreme.pdf")), 
        width = 8.5, height = 6)
    print(fig1)
    dev.off()
    
    # histogram for linear regression -log10P values
    fig <- ggplot(tmp_df, aes(x = get('lrtest_-log10P'))) + 
      geom_histogram(aes(y=..density..), fill = 'blue',color = 'black')+
      geom_vline(xintercept= quantile(tmp_df[['lrtest_-log10P']], 0.9),
                 linetype="dashed", color = 'red')+
      labs(y = 'density', x = 'lrtest_-log10P', 
           title = paste0(nrow(tmp_df),' cells above model_cutoff, skip ', length(common_cells) - nrow(tmp_df), ' cells'))
    ggsave(plot = fig, filename = fs::path(sub_out_dir3, paste0(blockID,"_histogram_lrtest_-log10P_allCells.jpeg")))
    
    
    
  }
}

####re-segmentation based on tLLRv2 score (1) flag cells, identify transcript groups ----
## (1) flag cells based on linear regression of tLLRv2, lrtest_-log10P
#5640 cells, 0.0292 of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > 5.0.
flagged_cells_cleaned <- modStats_cleaned_tLLRv2_3D[['cell_ID']][which(modStats_cleaned_tLLRv2_3D[['lrtest_-log10P']] > config_dimension[['flagCell_lrtestCutoff']])]
message(sprintf("%d cells, %.4f of all evaluated cells, are flagged for resegmentation with lrtest_-log10P > %.1f.", 
                length(flagged_cells_cleaned), length(flagged_cells_cleaned)/nrow(modStats_cleaned_tLLRv2_3D),config_dimension[['flagCell_lrtestCutoff']]))


## (2) use SVM~hyperplane to identify the connected transcripts group based on tLLRv2 score ----
# it turns out SVM can separate continuous low score transcript from the rest.
# but observed flagged cells with no flagged transcripts or multiple groups of flagged transcripts
flag_tLLRv2_cutoff = config_dimension[['tLLRv2_SVMcutoff']]
flagged_transDF3d_cleaned <- all_transDF[which(all_transDF[['cell_ID']] %in% flagged_cells_cleaned),]


## SVM in 3D, vectorized operation
# take ~1.5min to run 5640 cells, 45920416 transcripts for spatial 3D SVM, remove 2 cell
# coordinate in um when doing SVM
reseg_logInfo[['flagging_cleaned_SVM']] <- list(svm_config = config_dimension[['svm_config']])

system.time(tmp_df <- flagTranscripts_SVM(chosen_cells = flagged_cells_cleaned,
                                          score_GeneMatrix = transcript_loglik_cleaned,
                                          transcript_df = flagged_transDF3d_cleaned, 
                                          cellID_coln = 'cell_ID', 
                                          transID_coln = 'transcript_id', 
                                          score_coln = 'score_cleaned_tLLRv2_maxCellType',
                                          spatLocs_colns = c('x','y','z'), 
                                          model_cutoff = 50, 
                                          score_cutoff = flag_tLLRv2_cutoff, 
                                          svm_args = reseg_logInfo[['flagging_cleaned_SVM']][['svm_config']]))

# user  system elapsed 
# 314.317  12.442  80.250 

# add in SVM results to flagged transcript, cells with all transcript score on same class are removed
flagged_transDF_SVM3 <- tmp_df[, c('transcript_id','DecVal','SVM_class','SVM_cell_type')]
flagged_transDF_SVM3 <- merge(flagged_transDF_SVM3, flagged_transDF3d_cleaned, by = 'transcript_id')

message(sprintf("Remove %d cells with raw transcript score all in same class based on cutoff %.2f when running spatial SVM model.", 
                length(flagged_cells_cleaned) - length(unique(flagged_transDF_SVM3[['cell_ID']])), config_dimension[['tLLRv2_SVMcutoff']]))


# flagged transcript ID, character vector
flaggedSVM_transID3d <- flagged_transDF_SVM3[flagged_transDF_SVM3[['SVM_class']] ==0, get('transcript_id')]

reseg_logInfo[['flagging_cleaned_SVM']][['flagged_cells_cleaned']] <- flagged_cells_cleaned
reseg_logInfo[['flagging_cleaned_SVM']][['flagged_transID']] <- flaggedSVM_transID3d


# visualize the spatial plot of score
# choose 500 cells to plot, space across -log10P values
tmp_df <- modStats_cleaned_tLLRv2_3D[modStats_cleaned_tLLRv2_3D[['cell_ID']] %in% flagged_cells_cleaned, ]
tmp_df <- tmp_df[order(tmp_df[['lrtest_-log10P']]),]
cells_for_plots <- tmp_df[['cell_ID']][seq(1, length(flagged_cells_cleaned), by = round(length(flagged_cells_cleaned)/500))]
cells_for_plots <- unique(cells_for_plots)

if(TRUE){
  tmp_df <- flagged_transDF_SVM3
  
  for(score_coln in c('score_cleaned_tLLRv2_maxCellType','DecVal')){
    # visualize the flagged transcripts in chosen_cells
    if(score_coln == 'DecVal') {
      score_mid = 0
      score_range <- range(tmp_df[[score_coln]])
    } else {
      score_mid = flag_tLLRv2_cutoff
      score_range <- c(-5, 0)
    }
    n_perPage <- 25
    
    pdf(fs::path(sub_out_dir3, paste0(blockID,"_SpatialPlot_cleaned-SVM-3D_tLLRv2_", flag_tLLRv2_cutoff,"_flagged_transcripts_",score_coln,".pdf")))
    for(idx in seq(1, length(cells_for_plots), by = n_perPage)){
      if(idx +n_perPage-1 < length(cells_for_plots)){
        plot_data <- tmp_df[which(tmp_df[['cell_ID']] %in% cells_for_plots[idx : (idx+n_perPage-1)]), ]
      }else{
        plot_data <- tmp_df[which(tmp_df[['cell_ID']] %in% cells_for_plots[idx : length(cells_for_plots)]), ]
      }
      
      fig <- ggplot(plot_data, aes(x = x, y = y, group = as.factor(SVM_class), shape = as.factor(SVM_class)))+
        geom_point(colour = 'green', aes(size = as.factor(SVM_class)))+
        scale_size_manual(values = c(1, 0.1))+
        geom_point(size = 0.7, aes(color = get(score_coln)))+
        scale_color_gradientn(colours = c("blue","gray90","red"), 
                              values = scales::rescale(c(score_range[1], score_mid,score_range[2])), 
                              limits = score_range,
                              name = score_coln)+
        labs(size = "SVM_class", shape = "SVM_class")+
        facet_wrap(~cell_ID, scales = "free")+
        theme_classic()+
        theme(title = element_text(face = "bold"),
              axis.title = element_blank(), 
              axis.text = element_text(size = 6),
              legend.position = "top",
              legend.box = "horizontal",
              legend.key.size = unit(16, "pt"),
              legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(0,0,0,0),
              strip.text = element_text(face = "bold"),
              strip.background = element_rect(fill = "gray85"))
      print(fig)
    }
    dev.off()
  }
}




## (3) do network analysis on flagged transcript in vectorized operation to see if more than 1 connected group ----
## try to do denaulay on flagged transcript only and identity groups in network
# https://bookdown.org/markhoff/social_network_analysis/finding-groups-in-networks.html
# take ~1.5min to run 65651 transcripts in 5640 cells (451 cells = 8.00% with same class based on SVM)

system.time(flaggedSVM_transGroupDF3d <- groupTranscripts_Delanuay(chosen_transcripts = flaggedSVM_transID3d, 
                                                                   config_spatNW_transcript = config_spatNW2, 
                                                                   distance_cutoff = config_dimension[['TransDistance_cutoff']],
                                                                   transcript_df = flagged_transDF3d_cleaned, 
                                                                   cellID_coln = "cell_ID", 
                                                                   transID_coln = "transcript_id",
                                                                   transSpatLocs_coln = c('x','y','z')))
# user  system elapsed 
# 135.371   6.814  99.740

reseg_logInfo[['flagging_cleaned_SVM']][['cellsWflagged_transGroup']] <- unique(flaggedSVM_transGroupDF3d[['cell_ID']][which(flaggedSVM_transGroupDF3d[['transcript_group']] >0)])

reseg_logInfo[['flagging_cleaned_SVM']][['feat_count']] <- data.frame(cellNum = length(unique(flaggedSVM_transGroupDF3d[['cell_ID']])), 
                                                                      transNum = sum(flaggedSVM_transGroupDF3d[['transcript_group']] >0), 
                                                                      SVM_cellNum = length(unique(flagged_transDF_SVM3[['cell_ID']])))
message(sprintf("SVM spatial model further identified %d cells with transcript score all in same class, exclude from transcript group analysis.", 
                reseg_logInfo[['flagging_cleaned_SVM']][['feat_count']]$SVM_cellNum - reseg_logInfo[['flagging_cleaned_SVM']][['feat_count']]$cellNum))




# add in transcript group based on connectivity
flagged_transDF_SVM3[['connect_group']] <- 1- as.numeric(as.character(flagged_transDF_SVM3[['SVM_class']]))
group_converter <- flaggedSVM_transGroupDF3d[['transcript_group']] 
names(group_converter) <- flaggedSVM_transGroupDF3d[['transcript_id']]

tmp_idx <- which(flagged_transDF_SVM3[['transcript_id']] %in% flaggedSVM_transGroupDF3d[['transcript_id']])
flagged_transDF_SVM3[['connect_group']][tmp_idx] <- group_converter[flagged_transDF_SVM3[['transcript_id']][tmp_idx]]

# get new cell_id and cell type for each group
flagged_transDF_SVM3[, tmp_cellID := ifelse(connect_group == 0, cell_ID, paste0(cell_ID,'_g', connect_group))]

# get new cell type of each group based on maximum
tmp_df <- getCellType_maxScore(score_GeneMatrix = tLLRv2_geneMatrix_cleaned, 
                               transcript_df = flagged_transDF_SVM3, 
                               transID_coln = 'transcript_id',
                               transGene_coln = "target",
                               cellID_coln = "tmp_cellID", 
                               return_transMatrix = FALSE)
colnames(tmp_df[['cellType_DF']]) <- c('tmp_cellID','group_maxCellType')
flagged_transDF_SVM3 <- merge(flagged_transDF_SVM3, tmp_df[['cellType_DF']], by = 'tmp_cellID', all.x = TRUE)


#### (4) re-segmentation in neighborhood ----
## (4.1) get baseline cutoff from old assigned cell type in original segmentation ----
if(TRUE){
  # get quantile for per transcript tLLRv2 score of old assigned cell type in original segmentation
  # remove the cells with 'NotDet' in original cell type
  tmp_idx <- which(all_transDF[['cell_type']] != 'NotDet')
  span_tLLRv2_CellType <- tapply(all_transDF[['score_cleaned_tLLRv2_maxCellType']][tmp_idx], 
                                 all_transDF[['cleaned_tLLRv2_maxCellType']][tmp_idx], 
                                 function(x) quantile(x, probs = seq(0, 1, 0.25)))
  span_tLLRv2_CellType <- do.call(rbind, span_tLLRv2_CellType)
  # span of transcript number 
  tmp_df <- aggregate(all_transDF[['transcript_id']][tmp_idx], 
                      by = list(all_transDF[['CellId']][tmp_idx], 
                                all_transDF[['cleaned_tLLRv2_maxCellType']][tmp_idx]), 
                      FUN = length)
  colnames(tmp_df) <- c('cell_ID','cleaned_tLLRv2_maxCellType','transcript_num')
  span_transNum_CellType <- tapply(tmp_df[['transcript_num']][tmp_idx], 
                                   tmp_df[['cleaned_tLLRv2_maxCellType']][tmp_idx], 
                                   function(x) quantile(x, probs = seq(0, 1, 0.25)))
  span_transNum_CellType <- do.call(rbind, span_transNum_CellType)
  
  reseg_logInfo[['baseline_celltype']] <- list(span_score = span_tLLRv2_CellType, 
                                               span_transNum = span_transNum_CellType)
}


# pause to save data
if(FALSE){
  # 1min to save before evaluate neighborhood
  system.time(save(currentBlock, # slide, slide name, slide position of current block
                   block_cellNetWorkDT, #cell delaunay network data table only
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
                   all_transDF, # transcript metadata from raw
                   all_cellstatsDF, # cell stats from raw
                   select_cellmeta, # cell typing info for all cells, exclude extracellular
                   meanCelltype_profiles, # mean expression profiles
                   transcript_loglik_cleaned, # tLL Gscore, gene x cell type
                   tLLRv2_geneMatrix_cleaned, # recenter net Gscore based on maximum per row/transcript, gene x cell type
                   modStats_cleaned_tLLRv2_3D, # 3D spatial fitting to linear regression hyperplane using tLLRv2 under old cell type
                   # for resegmentation in 3D
                   flagged_cells_cleaned, # cells flagged for resegmenation
                   flagged_transDF_SVM3, # 3D SVM spatial model to flag transcripts in flagged cells
                   flaggedSVM_transID3d, # flagged transcript ID
                   flaggedSVM_transGroupDF3d, # group the flagged transcripts based on delanuay spatial network
                   # neighborReSeg_df_cleanSVM, # neighborhood info for flagged transcript groups
                   # post_reseg_results_cleanSVM_leiden, # post resegmentation results
                   # # for visualization of resegmentation impact (no LDA)
                   # altered_cells_cleanSVM_leiden, # cells changed during resegmentation
                   # neighborTransDF_cleanSVM_leiden, # neighborhood transcript data for cells receiving merge
                   # alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D, # linear regression on altered cells after resegmentation
                   file = fs::path(sub_out_dir3, paste0(blockID,"_NSCLC_1Mcell_980plx_reseg_cleanNBclust6_SVMleiden_log.RData"))))
  # user  system elapsed 
  # 44.150   7.367  65.584
  
}


## (4.2) prepare all_transDF for re-segmentation ----
# update the transcript_df with flagged transcript_group
# take ~2min for merge 39916952 rows with 147948 rows
reseg_transcript_df <- merge(all_transDF, 
                             flagged_transDF_SVM3[, c('transcript_id', 'connect_group','tmp_cellID','group_maxCellType')], 
                             by = 'transcript_id', all.x = TRUE)
# fill in the missing values for unflagged cells
# take more than 2min to get index for na; 115MB size for tmp_idx; faster operation with idx than TRUE/FALSE
tmp_idx <- which(is.na(reseg_transcript_df[['connect_group']]))
reseg_transcript_df[['connect_group']][tmp_idx]<-rep(0, length(tmp_idx))
reseg_transcript_df[['tmp_cellID']][tmp_idx] <- reseg_transcript_df[['cell_ID']][tmp_idx]
reseg_transcript_df[['group_maxCellType']][tmp_idx] <- reseg_transcript_df[['cleaned_tLLRv2_maxCellType']][tmp_idx]

## (4.3) evaluate the neighborhood of each group for re-segmentation ----
cells_to_use <- unique(flagged_transDF_SVM3[which(flagged_transDF_SVM3[['connect_group']]!=0),][['tmp_cellID']])
score_baseline <- reseg_logInfo[['baseline_celltype']][['span_score']][,"25%"]
transNum_lowCutoff  <- reseg_logInfo[['baseline_celltype']][['span_transNum']][,"25%"]
transNum_highCutoff  <- reseg_logInfo[['baseline_celltype']][['span_transNum']][,"50%"]

reseg_logInfo[['resegmentation_leiden']] <- list(score_baseline = score_baseline, 
                                                 transNum_lowCutoff = transNum_lowCutoff, 
                                                 transNum_highCutoff = transNum_highCutoff)

# did not consider extra cellular transcripts for neighbor identification. 

### search within absolute distance, 25um in xy for cell level search, consider 15 pixel = 2.7um to be direct neighbor at transcript level.
# new function ussing spatstat to locate neighbor cells and rank them by minimal molecular distance to query cell
# 13038 transcript groups flagged among 240148 common groups, done with in 5.67hr in 2D, 4.82hr in 3D
system.time(neighborReSeg_df_cleanSVM <- neighborhood_for_resegment_spatstat(chosen_cells = cells_to_use,
                                                                             score_GeneMatrix = tLLRv2_geneMatrix_cleaned,
                                                                             score_baseline = score_baseline,
                                                                             neighbor_distance_xy = config_dimension[['CellNeighbor_xy_in_transDF']],
                                                                             distance_cutoff = config_dimension[['TransDistance_cutoff']],
                                                                             transcript_df = reseg_transcript_df,
                                                                             cellID_coln = "tmp_cellID",
                                                                             celltype_coln = "group_maxCellType",
                                                                             transID_coln = "transcript_id",
                                                                             transGene_coln = "target",
                                                                             transSpatLocs_coln = c('x','y','z')))


# ### when distance cutoff evaluation done in 2D
# user    system   elapsed 
# 18708.929  5745.984 20422.588 

### when distance cutoff evaluation done in 3D
# user   system  elapsed 
# 13636.29  5294.41 17352.50

# pause to save data
if(FALSE){
  # 42.14min to save after neighborhood evaluation before resegmentation action
  system.time(save(currentBlock, # slide, slide name, slide position of current block
                   block_cellNetWorkDT, #cell delaunay network data table only
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
                   all_transDF, # transcript metadata from raw
                   all_cellstatsDF, # cell stats from raw
                   select_cellmeta, # cell typing info for all cells, exclude extracellular
                   meanCelltype_profiles, # mean expression profiles
                   transcript_loglik_cleaned, # tLL Gscore, gene x cell type
                   tLLRv2_geneMatrix_cleaned, # recenter net Gscore based on maximum per row/transcript, gene x cell type
                   modStats_cleaned_tLLRv2_3D, # 3D spatial fitting to linear regression hyperplane using tLLRv2 under old cell type
                   # for resegmentation in 3D
                   flagged_cells_cleaned, # cells flagged for resegmenation
                   flagged_transDF_SVM3, # 3D SVM spatial model to flag transcripts in flagged cells
                   flaggedSVM_transID3d, # flagged transcript ID
                   flaggedSVM_transGroupDF3d, # group the flagged transcripts based on delanuay spatial network
                   neighborReSeg_df_cleanSVM, # neighborhood info for flagged transcript groups
                   # post_reseg_results_cleanSVM_leiden, # post resegmentation results
                   # # for visualization of resegmentation impact (no LDA)
                   # altered_cells_cleanSVM_leiden, # cells changed during resegmentation
                   # neighborTransDF_cleanSVM_leiden, # neighborhood transcript data for cells receiving merge
                   # alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D, # linear regression on altered cells after resegmentation
                   file = fs::path(sub_out_dir3, paste0(blockID,"_NSCLC_1Mcell_980plx_reseg_cleanNBclust6_SVMleiden_log.RData"))))
  # user   system  elapsed 
  # 304.348  127.800 2528.149
  
}
#### (4.4) decide resegmentation operation: merge, new cell, or discard ----
reseg_logInfo[['resegmentation_leiden']][['leiden_config']] <- config_dimension[['leiden_config']]
reseg_logInfo[['resegmentation_leiden']][['cutoff_sharedLeiden']] = config_dimension[['cutoff_sharedLeiden']]

# it takes 22.22min to make resegemntation decision if use leiden clustering to judge on 2243 merging events
system.time(reseg_actions <- decide_ReSegment_Operations_leidenCut(neighborhood_df = neighborReSeg_df_cleanSVM, 
                                                                   score_baseline = score_baseline, 
                                                                   lowerCutoff_transNum = transNum_lowCutoff, 
                                                                   higherCutoff_transNum = transNum_highCutoff,
                                                                   config_spatNW_transcript = config_spatNW2,
                                                                   transcript_df = reseg_transcript_df, 
                                                                   cellID_coln = "tmp_cellID",
                                                                   transID_coln = "transcript_id",
                                                                   transSpatLocs_coln = c('x','y','z'), 
                                                                   leiden_config = reseg_logInfo[['resegmentation_leiden']][['leiden_config']], 
                                                                   cutoff_sharedLeiden = reseg_logInfo[['resegmentation_leiden']][['cutoff_sharedLeiden']]))
# user   system  elapsed 
# 1440.335  112.845 1333.174  

reseg_logInfo[['resegmentation_leiden']][['reseg_actions']] <- reseg_actions
reseg_logInfo[['resegmentation_leiden']][['feat_count']] <- data.frame(cellNum_discard = length(reseg_actions$cells_to_discard), 
                                                                       cellNum_update = length(reseg_actions$cells_to_update), 
                                                                       cellNum_keep = length(reseg_actions$cells_to_keep), 
                                                                       cellNum_total = nrow(neighborReSeg_df_cleanSVM)) 





#### (4.5) update with new cells_to_update ----
# update neighborReSeg_df_cleanSVM with resegmentaion actions
neighborReSeg_df_cleanSVM[['corrected_CellId']] <- reseg_actions[['reseg_full_converter']][neighborReSeg_df_cleanSVM[['CellId']]]

# # update the resegmentation data frame, remove cells_to_discard, take ~6.3min to update and get count matrix
# > dim(reseg_transcript_df)
# [1] 45920416       16
# > dim(all_transDF)
# [1] 45920416       13

system.time(post_reseg_results_cleanSVM_leiden <- update_transDF_ResegActions(transcript_df = reseg_transcript_df, 
                                                                              reseg_full_converter = reseg_actions$reseg_full_converter, 
                                                                              score_GeneMatrix = tLLRv2_geneMatrix_cleaned, 
                                                                              transGene_coln = 'target', 
                                                                              cellID_coln = 'tmp_cellID', 
                                                                              celltype_coln = 'group_maxCellType', 
                                                                              spatLocs_colns = c('x','y','z'), 
                                                                              return_perCellDF = TRUE))

# # return gene x cell count table
# user  system elapsed 
# 501.915  94.484 378.297 
# > dim(post_reseg_results_cleanSVM_leiden$updated_transDF)
# [1] 45867194       18
# > dim(post_reseg_results_cleanSVM_leiden$perCell_DT)
# [1] 227181      5

# get tLLRv2 score under updated cell type, 1min for 45867194 transcripts
system.time(tmp_df <- getScoreCellType_gene(score_GeneMatrix = tLLRv2_geneMatrix_cleaned, 
                                            transcript_df = post_reseg_results_cleanSVM_leiden$updated_transDF, 
                                            transID_coln = "transcript_id",
                                            transGene_coln = "target",
                                            celltype_coln = 'updated_celltype'))
# user  system elapsed 
# 34.244  18.445  51.620    

post_reseg_results_cleanSVM_leiden$updated_transDF <- merge(post_reseg_results_cleanSVM_leiden$updated_transDF, tmp_df, by = 'transcript_id')


#### visualize the resegmentation outcomes ----
# altered cells
altered_cells_cleanSVM_leiden <- list(
  # original cell_ID with transcripts got discarded
  oriCells_trimmed = unique(sapply(strsplit(reseg_actions$cells_to_discard,"_g"),"[[",1)), 
  # original cell_ID with transcripts got split, the split groups were kept as new cells
  oriCells_split = unique(sapply(strsplit(reseg_actions$cells_to_keep,"_g"),"[[",1)), 
  # new cell_ID that received merge cells
  updatedCells_merged = unique(reseg_actions$cells_to_update), 
  # new cell_ID that got split and kept as separate new cells
  updatedCells_kept = unique(reseg_actions$cells_to_keep))

# get neighborhood transcript data.frame for cells receiving merge
# 25um in xy search range for neighbor cells
# 29min to get info for 1238 cells among 227181 all cells
system.time(neighborTransDF_cleanSVM_leiden <- getNeighbors_transDF_spatstat(chosen_cells = altered_cells_cleanSVM_leiden$updatedCells_merged, 
                                                                             neighbor_distance_xy = config_dimension[['CellNeighbor_xy_in_transDF']],
                                                                             transcript_df = post_reseg_results_cleanSVM_leiden$updated_transDF,
                                                                             cellID_coln = "updated_cellID",
                                                                             transID_coln = "transcript_id",
                                                                             transSpatLocs_coln = c('x','y','z')))
# user   system  elapsed 
# 2526.725  677.435 1739.795 

# add in cell type and cell_ID information
neighborTransDF_cleanSVM_leiden <- merge(neighborTransDF_cleanSVM_leiden, 
                                         post_reseg_results_cleanSVM_leiden$updated_transDF[, .SD, .SDcols = c('transcript_id','updated_celltype', 
                                                                                                               'cell_ID','cleaned_tLLRv2_maxCellType','tmp_cellID','group_maxCellType')], 
                                         by = 'transcript_id')
data.table::setkeyv(neighborTransDF_cleanSVM_leiden, c('query_CellId','updated_cellID', 'transcript_id'))


## for plotting
tmp_df <- neighborTransDF_cleanSVM_leiden
# cell type as color
color_code <- Giotto::getDistinctColors(ncol(meanCelltype_profiles))
names(color_code) <- colnames(meanCelltype_profiles)

# black border via cell size
size_code <- c(0.5, 0.01)
names(size_code) <- c("TRUE","FALSE")

# number of query cell per page of plot
n_perPage <- 16

# new cell type and new  cell id
## for plotting, plot 500 cells
tmp_cellID <- unique(tmp_df[['query_CellId']])
if(length(tmp_cellID) >510){
  tmp_cellID <- tmp_cellID[seq(1, length(tmp_cellID), by = round(length(tmp_cellID)/500))]
  tmp_cellID <- unique(tmp_cellID)
}


score_coln <- 'updated_celltype'
label_coln <- 'updated_cellID'

if(TRUE){
  pdf(fs::path(sub_out_dir3, paste0(blockID,"_Neighborhood_SVM3d", flag_tLLRv2_cutoff,"_resegmented_",score_coln,'_', label_coln,".pdf")), 
      width = 16, height = 15)
  for(idx in seq(1, length(tmp_cellID), by = n_perPage)){
    if(idx +n_perPage-1 < length(tmp_cellID)){
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : (idx+n_perPage-1)])
    }else{
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : length(tmp_cellID)])
    }
    plot_data <- tmp_df[tmp_idx,]
    plotDT_perCell <- plot_data[, lapply(.SD, mean), by = .(query_CellId, get(label_coln)), .SDcols = c('x','y')]
    colnames(plotDT_perCell) <- c('query_CellId',label_coln, 'x','y')
    
    fig <- ggplot()+
      geom_point(data = plot_data, color = 'black', aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                        size = as.factor(in_query_cell)))+
      scale_size_manual(values = size_code)+
      geom_point(data = plot_data,size = 0.2, alpha = 0.8, aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                               color = as.factor(get(score_coln))))+
      scale_color_manual(values = color_code)+
      geom_label(data = plotDT_perCell, aes(x = x, y = y, label = as.factor(get(label_coln))), 
                 size = 2, label.size = NA, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), 
                 fontface = "bold", fill = alpha(c("white"),0.5))+
      labs(color = score_coln, size = 'in_query_cell')+
      #guides(color=guide_legend(ncol=5), size=guide_legend(ncol=2))+
      facet_wrap(~query_CellId, scales = "free")+
      theme_classic()+
      theme(title = element_text(face = "bold"),
            axis.title = element_blank(), 
            axis.text = element_text(size = 6),
            # legend.position = "top",
            # legend.box = "horizontal",
            legend.key.size = unit(16, "pt"),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6, face = "bold"), 
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = "gray85"))
    
    print(fig)
  }
  dev.off()
  
}



# old cell type and old cell id
score_coln <- 'cleaned_tLLRv2_maxCellType'
label_coln <- 'cell_ID'

if(TRUE){
  pdf(fs::path(sub_out_dir3, paste0(blockID,"_Neighborhood_SVM3d", flag_tLLRv2_cutoff,"_resegmented_",score_coln,'_', label_coln,".pdf")), 
      width = 16, height = 15)
  for(idx in seq(1, length(tmp_cellID), by = n_perPage)){
    if(idx +n_perPage-1 < length(tmp_cellID)){
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : (idx+n_perPage-1)])
    }else{
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : length(tmp_cellID)])
    }
    plot_data <- tmp_df[tmp_idx,]
    plotDT_perCell <- plot_data[, lapply(.SD, mean), by = .(query_CellId, get(label_coln)), .SDcols = c('x','y')]
    colnames(plotDT_perCell) <- c('query_CellId',label_coln, 'x','y')
    
    fig <- ggplot()+
      geom_point(data = plot_data, color = 'black', aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                        size = as.factor(in_query_cell)))+
      scale_size_manual(values = size_code)+
      geom_point(data = plot_data,size = 0.2, alpha = 0.8, aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                               color = as.factor(get(score_coln))))+
      scale_color_manual(values = color_code)+
      geom_label(data = plotDT_perCell, aes(x = x, y = y, label = as.factor(get(label_coln))), 
                 size = 2, label.size = NA, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), 
                 fontface = "bold", fill = alpha(c("white"),0.5))+
      labs(color = score_coln, size = 'in_query_cell')+
      #guides(color=guide_legend(ncol=5), size=guide_legend(ncol=2))+
      facet_wrap(~query_CellId, scales = "free")+
      theme_classic()+
      theme(title = element_text(face = "bold"),
            axis.title = element_blank(), 
            axis.text = element_text(size = 6),
            # legend.position = "top",
            # legend.box = "horizontal",
            legend.key.size = unit(16, "pt"),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6, face = "bold"), 
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = "gray85"))
    
    print(fig)
  }
  dev.off()
  
}


# per group cell_ID and cell type
score_coln <- 'group_maxCellType'
label_coln <- 'tmp_cellID'

if(TRUE){
  pdf(fs::path(sub_out_dir3, paste0(blockID,"_Neighborhood_SVM3d", flag_tLLRv2_cutoff,"_resegmented_",score_coln,'_', label_coln,".pdf")), 
      width = 16, height = 15)
  for(idx in seq(1, length(tmp_cellID), by = n_perPage)){
    if(idx +n_perPage-1 < length(tmp_cellID)){
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : (idx+n_perPage-1)])
    }else{
      tmp_idx <- which(tmp_df[['query_CellId']] %in% tmp_cellID[idx : length(tmp_cellID)])
    }
    plot_data <- tmp_df[tmp_idx,]
    plotDT_perCell <- plot_data[, lapply(.SD, mean), by = .(query_CellId, get(label_coln)), .SDcols = c('x','y')]
    colnames(plotDT_perCell) <- c('query_CellId',label_coln, 'x','y')
    
    fig <- ggplot()+
      geom_point(data = plot_data, color = 'black', aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                        size = as.factor(in_query_cell)))+
      scale_size_manual(values = size_code)+
      geom_point(data = plot_data,size = 0.2, alpha = 0.8, aes(x = x, y = y, group = as.factor(get(score_coln)), 
                                                               color = as.factor(get(score_coln))))+
      scale_color_manual(values = color_code)+
      geom_label(data = plotDT_perCell, aes(x = x, y = y, label = as.factor(get(label_coln))), 
                 size = 2, label.size = NA, label.padding = unit(0.1, "lines"), label.r = unit(0, "lines"), 
                 fontface = "bold", fill = alpha(c("white"),0.5))+
      labs(color = score_coln, size = 'in_query_cell')+
      #guides(color=guide_legend(ncol=5), size=guide_legend(ncol=2))+
      facet_wrap(~query_CellId, scales = "free")+
      theme_classic()+
      theme(title = element_text(face = "bold"),
            axis.title = element_blank(), 
            axis.text = element_text(size = 6),
            # legend.position = "top",
            # legend.box = "horizontal",
            legend.key.size = unit(16, "pt"),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6, face = "bold"), 
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = "gray85"))
    
    print(fig)
  }
  dev.off()
  
}




#### do linear regression on altered cells ----
# check p-value range, skip 631 cells, only 2076 altered cells being evlauated 
if(TRUE){
  tmp_cellID <- unique(c(altered_cells_cleanSVM_leiden$updatedCells_merged, altered_cells_cleanSVM_leiden$updatedCells_kept))
  tmp_df <- score_cell_segmentation_error(chosen_cells = tmp_cellID, 
                                          transcript_df = post_reseg_results_cleanSVM_leiden$updated_transDF, 
                                          cellID_coln = 'updated_cellID', 
                                          transID_coln = 'transcript_id', 
                                          score_coln = 'score_updated_celltype',
                                          spatLocs_colns = c('x','y','z'), 
                                          model_cutoff = 50)
  
  #-log10(P)
  tmp_df[['lrtest_-log10P']] <- -log10(tmp_df[['lrtest_Pr']])
  alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D <- merge(tmp_df, 
                                                               post_reseg_results_cleanSVM_leiden$perCell_DT[, .SD, .SDcols = c('updated_cellID','updated_celltype')], 
                                                               by.x = 'cell_ID', by.y = 'updated_cellID')
  colnames(alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D)[1] <- 'updated_cellID'
  
  
  # visualize extreme cells, 2D plots
  if(TRUE){
    tmp_df <- data.table::copy(alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D)
    colnames(tmp_df)[1] <- 'cell_ID'
    tmp_df <- tmp_df[order(tmp_df[['lrtest_-log10P']]),]
    tmp_df[['labels']] <- paste0(tmp_df[['slide']],'_', tmp_df[['updated_celltype']], 
                                 ', -log10P=', round(tmp_df[['lrtest_-log10P']],2))
    chosen_cells <- c(tmp_df$cell_ID[1:9], tmp_df$cell_ID[(nrow(tmp_df)-9):nrow(tmp_df)])
    fig1 <- plotSpatialScoreMultiCells(chosen_cells = chosen_cells, 
                                       cell_labels = tmp_df[match(chosen_cells, tmp_df$cell_ID), 'labels'], 
                                       transcript_df = post_reseg_results_cleanSVM_leiden$updated_transDF,
                                       cellID_coln = "cell_ID", 
                                       transID_coln = "transcript_id",
                                       score_coln = "score_updated_celltype", 
                                       spatLocs_colns = c("x","y"))  
    
    pdf(fs::path(sub_out_dir3, paste0(blockID,"_SpatialPlot_SpatModel2_tLLRv2_-log10P_resegmented_AlteredOnly_cleanNBclust6_SVMleiden.pdf")), 
        width = 8.5, height = 6)
    print(fig1)
    dev.off()
    
    # histogram for linear regression -log10P values
    fig <- ggplot(tmp_df, aes(x = get('lrtest_-log10P'))) + 
      geom_histogram(aes(y=..density..), fill = 'blue',color = 'black')+
      geom_vline(xintercept= quantile(tmp_df[['lrtest_-log10P']], 0.9),
                 linetype="dashed", color = 'red')+
      labs(y = 'density', x = 'lrtest_-log10P', 
           title = paste0(nrow(tmp_df),' cells above model_cutoff, skip ', length(tmp_cellID) - nrow(tmp_df), ' cells'))
    ggsave(plot = fig, filename = fs::path(sub_out_dir3, paste0(blockID,"_histogram_lrtest_-log10P_resegmented_AlteredOnly_cleanNBclust6_SVMleiden.jpeg")))
    
  }
}



# compare stats across slides
if(TRUE){
  # total cells
  tmp_df <- all_cellstatsDF 
  tmp2_df <- tmp_df[, c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'rawCells')
  stats_df <- tmp2_df
  
  # cells being evaluated for segmentation error
  tmp_df <- modStats_cleaned_tLLRv2_3D  
  tmp2_df <- tmp_df[, c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'lrEvaluatedCells')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  
  # flagged cell distribution
  # histogram for linear regression -log10P values
  tmp_df <- modStats_cleaned_tLLRv2_3D 
  tmp_cellID <- reseg_logInfo$flagging_cleaned_SVM$flagged_cells_cleaned
  tmp2_df <- tmp_df[which(tmp_df$cell_ID %in% tmp_cellID), c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'flaggedCells')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  fig <- ggplot(tmp_df, aes(x = get('lrtest_-log10P'), group = as.factor(slide), 
                            color = as.factor(slide), fill = as.factor(slide))) + 
    geom_histogram(position="dodge2", alpha = 0.5)+
    geom_vline(xintercept= 10,
               linetype="dashed", color = 'black')+
    labs(x = 'lrtest_-log10P', color = "slide", fill = "slide",
         title = paste0("slide ", paste0(tmp2_df$slide, collapse = ":"), "w/ flagged cell #, ", paste0(tmp2_df$flaggedCells, collapse = ":")))+
    theme_linedraw()
  ggsave(plot = fig, filename = fs::path(sub_out_dir3, paste0(blockID,"_histogram_lrtest_-log10P_slide-comparison_cleanNBclust6_SVMleiden.jpeg")))
  
  
  
  # resegmented cell distribution
  # cells with transcript groups identified 
  tmp_df <- flaggedSVM_transGroupDF3d[flaggedSVM_transGroupDF3d$transcript_group !=0, ]
  tmp_df[['slide']] <- sapply(strsplit(tmp_df$cell_ID, "_"),"[[", 2)
  tmp_df <- unique(tmp_df[, c('cell_ID','slide')])
  tmp2_df <- tmp_df[, c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'cellNum_TransGroup')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  
  # cells with transcript groups got discarded
  tmp_cellID <- altered_cells_cleanSVM_leiden$oriCells_trimmed
  tmp_df <- modStats_cleaned_tLLRv2_3D 
  tmp2_df <- tmp_df[which(tmp_df$cell_ID %in% tmp_cellID), c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'cellNum_GroupDiscarded')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  # cells got split and kept the split groups as new cells
  tmp_cellID <- altered_cells_cleanSVM_leiden$oriCells_split
  tmp_df <- modStats_cleaned_tLLRv2_3D 
  tmp2_df <- tmp_df[which(tmp_df$cell_ID %in% tmp_cellID), c('cell_ID', 'slide')]
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'cellNum_GroupKept')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  # new cells being updated during resementation 
  tmp_cellID <- sapply(strsplit(altered_cells_cleanSVM_leiden$updatedCells_merged,"_g"),"[[",1)
  tmp_df <- data.frame(cell_ID = tmp_cellID)
  tmp_df[['slide']] <- sapply(strsplit(tmp_df$cell_ID, "_"),"[[", 2)
  tmp2_df <- tmp_df
  tmp2_df <- aggregate(tmp2_df$cell_ID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'mergeNum')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  # cell number after resegmentation
  tmp_df <- post_reseg_results_cleanSVM_leiden$perCell_DT
  tmp_df[['slide']] <- sapply(strsplit(tmp_df$updated_cellID, "_"),"[[", 2)
  tmp2_df <- tmp_df
  tmp2_df <- aggregate(tmp2_df$updated_cellID, by = list(tmp2_df$slide), length)
  colnames(tmp2_df) <- c('slide', 'cellNum_postReseg')
  stats_df <- merge(stats_df, tmp2_df, by = 'slide')
  
  # percentage to rawCells
  reseg_logInfo[['featcount_cleanSVM_leiden']] <- list(cellNum = stats_df, 
                                                       normalized_cellNum = stats_df[, 2:ncol(stats_df)]/stats_df$rawCells)
  plot_data <- reshape2::melt(stats_df, id.var = "slide")
  fig <- ggplot(plot_data, aes(x = slide, fill = as.factor(slide)))+
    geom_bar(aes(y = value), stat = "identity")+
    geom_text(y = 0, aes(label = value), vjust = 0)+
    facet_wrap(~variable)+
    labs(x = 'slide', y = 'cell number')+
    theme_linedraw()+
    theme(legend.position = "none", 
          strip.text = element_text(face = "bold", color = "black"),
          strip.background = element_rect(fill = "gray85"))
  ggsave(plot = fig, filename = fs::path(sub_out_dir3, paste0(blockID,"_histogram_slide-comparison_resegmentation-process_cleanNBclust6_SVMleiden.jpeg")))
  
  # plot score cell type and flagged for FOV1 only
  # visualize the spatial plot of score
  # choose 500 cells to plot, space across -log10P values
  tmp_df <- modStats_cleaned_tLLRv2_3D[modStats_cleaned_tLLRv2_3D[['cell_ID']] %in% flagged_cells_cleaned, ]
  tmp_df[['slide']] <- sapply(strsplit(tmp_df$cell_ID, "_"),"[[", 2)
  tmp_df[['fov']] <- sapply(strsplit(tmp_df$cell_ID, "_"),"[[", 3)
  tmp_df <- tmp_df[which(tmp_df$fov ==1), ]
  tmp_df <- tmp_df[order(tmp_df[['lrtest_-log10P']]),]
  cells_for_plots <- tmp_df[['cell_ID']]
  
  if(length(cells_for_plots) >510){
    cells_for_plots <- cells_for_plots[1:510]
  }
  
  if(TRUE){
    tmp_df <- flagged_transDF_SVM3
    
    for(score_coln in c('score_cleaned_tLLRv2_maxCellType','DecVal')){
      # visualize the flagged transcripts in chosen_cells
      if(score_coln == 'DecVal') {
        score_mid = 0
        score_range <- range(tmp_df[[score_coln]])
      } else {
        score_mid = flag_tLLRv2_cutoff
        score_range <- c(-5, 0)
      }
      n_perPage <- 25
      
      pdf(fs::path(sub_out_dir3, paste0(blockID,"_FOV1_SpatialPlot_cleanNBclust6_SVMleiden_tLLRv2_", flag_tLLRv2_cutoff,"_flagged_transcripts_",score_coln,".pdf")))
      for(idx in seq(1, length(cells_for_plots), by = n_perPage)){
        if(idx +n_perPage-1 < length(cells_for_plots)){
          plot_data <- tmp_df[which(tmp_df[['cell_ID']] %in% cells_for_plots[idx : (idx+n_perPage-1)]), ]
        }else{
          plot_data <- tmp_df[which(tmp_df[['cell_ID']] %in% cells_for_plots[idx : length(cells_for_plots)]), ]
        }
        
        fig <- ggplot(plot_data, aes(x = x, y = y, group = as.factor(SVM_class), shape = as.factor(SVM_class)))+
          geom_point(colour = 'green', aes(size = as.factor(SVM_class)))+
          scale_size_manual(values = c(1, 0.1))+
          geom_point(size = 0.7, aes(color = get(score_coln)))+
          scale_color_gradientn(colours = c("blue","gray90","red"), 
                                values = scales::rescale(c(score_range[1], score_mid,score_range[2])), 
                                limits = score_range,
                                name = score_coln)+
          labs(size = "SVM_class", shape = "SVM_class")+
          facet_wrap(~cell_ID, scales = "free")+
          theme_classic()+
          theme(title = element_text(face = "bold"),
                axis.title = element_blank(), 
                axis.text = element_text(size = 6),
                legend.position = "top",
                legend.box = "horizontal",
                legend.key.size = unit(16, "pt"),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(0,0,0,0),
                strip.text = element_text(face = "bold"),
                strip.background = element_rect(fill = "gray85"))
        print(fig)
      }
      dev.off()
    }
  }
  
  # save reseg log as separate data for comparison across runs
  save(reseg_logInfo, file = fs::path(sub_out_dir3, paste0(blockID,"_NSCLC_1Mcell_980plx_reseg_cleanNBclust6_SVMleiden_3dSegRecords.RData")))
  
  
  
}


#### save outputs --------
# 3min to save after resegmentation
system.time(save(currentBlock, # slide, slide name, slide position of current block
                 block_cellNetWorkDT, #cell delaunay network data table only
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
                 all_transDF, # transcript metadata from raw
                 all_cellstatsDF, # cell stats from raw
                 select_cellmeta, # cell typing info for all cells, exclude extracellular
                 meanCelltype_profiles, # mean expression profiles
                 transcript_loglik_cleaned, # tLL Gscore, gene x cell type
                 tLLRv2_geneMatrix_cleaned, # recenter net Gscore based on maximum per row/transcript, gene x cell type
                 modStats_cleaned_tLLRv2_3D, # 3D spatial fitting to linear regression hyperplane using tLLRv2 under old cell type
                 # for resegmentation in 3D
                 flagged_cells_cleaned, # cells flagged for resegmenation
                 flagged_transDF_SVM3, # 3D SVM spatial model to flag transcripts in flagged cells
                 flaggedSVM_transID3d, # flagged transcript ID
                 flaggedSVM_transGroupDF3d, # group the flagged transcripts based on delanuay spatial network
                 neighborReSeg_df_cleanSVM, # neighborhood info for flagged transcript groups
                 post_reseg_results_cleanSVM_leiden, # post resegmentation results
                 # for visualization of resegmentation impact (no LDA)
                 altered_cells_cleanSVM_leiden, # cells changed during resegmentation
                 neighborTransDF_cleanSVM_leiden, # neighborhood transcript data for cells receiving merge
                 alteredOnly_modStats_resegcleanSVM_leiden_tLLRv2_3D, # linear regression on altered cells after resegmentation
                 file = fs::path(sub_out_dir3, paste0(blockID,"_NSCLC_1Mcell_980plx_reseg_cleanNBclust6_SVMleiden_log.RData"))))

# user  system elapsed 
# 96.553  17.391 155.516


