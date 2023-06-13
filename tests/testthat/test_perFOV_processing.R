library(FastReseg)
library(testthat)

# example data
data("mini_transcriptDF")
data("example_refProfiles")
data("example_baselineCT")

extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID']

score_baseline = example_baselineCT[["span_score"]][, "25%"]
lowerCutoff_transNum = example_baselineCT[["span_transNum"]][, "25%"]
higherCutoff_transNum = example_baselineCT[["span_transNum"]][, "50%"]

# calculate log-likelihood of each gene under each cell type and center the score matrix on per gene basis
score_GeneMatrix <- scoreGenesInRef(genes = rownames(example_refProfiles), 
                                    ref_profiles = pmax(example_refProfiles, 1e-5))

# set cutoffs
flagCell_lrtest_cutoff <- 5
svmClass_score_cutoff <- -2
molecular_distance_cutoff <- 2.7
cellular_distance_cutoff <- 20

### run perFOV processing pipeline wrapper 
# case 1: return all intermediates and perCell data
res1 <- fastReseg_perFOV_full_process(score_GeneMatrix = score_GeneMatrix, 
                                      transcript_df = mini_transcriptDF, 
                                      transID_coln = 'UMI_transID',
                                      transGene_coln = "target",
                                      cellID_coln = 'UMI_cellID', 
                                      spatLocs_colns = c('x','y','z'), 
                                      extracellular_cellID = extracellular_cellID, 
                                      flagModel_TransNum_cutoff = 50, 
                                      flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
                                      svmClass_score_cutoff = svmClass_score_cutoff, 
                                      molecular_distance_cutoff = molecular_distance_cutoff,
                                      cellular_distance_cutoff = cellular_distance_cutoff,
                                      score_baseline = score_baseline, 
                                      lowerCutoff_transNum = lowerCutoff_transNum, 
                                      higherCutoff_transNum = higherCutoff_transNum,
                                      
                                      groupTranscripts_method = "dbscan",
                                      spatialMergeCheck_method = "leidenCut", 
                                      cutoff_spatialMerge = 0.5,
                                      return_intermediates = TRUE,
                                      return_perCellData = TRUE, 
                                      includeAllRefGenes = TRUE, 
                                      seed_process = 123)

test_that("fastReseg_perFOV_full_process(), per FOV pipeline processing, returns the expected output", {
  expect_type(res1, "list")
  expect_true("modStats_ToFlagCells" %in% names(res1))
  expect_true("groupDF_ToFlagTrans" %in% names(res1))
  expect_true("neighborhoodDF_ToReseg" %in% names(res1))
  expect_true("reseg_actions" %in% names(res1))
  expect_true("updated_transDF" %in% names(res1))
  expect_true("updated_perCellDT" %in% names(res1))
  expect_true("updated_perCellExprs" %in% names(res1))
  
  # Test the class and structure of the data frames and matrices
  expect_true(is.data.frame(res1$modStats_ToFlagCells))
  expect_true(is.data.frame(res1$groupDF_ToFlagTrans))
  expect_true(is.data.frame(res1$neighborhoodDF_ToReseg))
  expect_type(res1$reseg_actions, "list")
  expect_true(is.data.frame(res1$updated_transDF))
  expect_true("data.table" %in% class(res1$updated_perCellDT))
  expect_true("dgCMatrix" %in% class(res1$updated_perCellExprs))
  
  # Test the presence of specific columns in the data frames or elements in list
  expect_true(all(c('UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'updated_cellID', 'updated_celltype') %in% colnames(res1$updated_transDF)))
  expect_true(all(c('UMI_cellID', 'tLLR_maxCellType', 'transcript_num', 'lrtest_Pr', 'lrtest_-log10P', 'flagged') %in% colnames(res1$modStats_ToFlagCells)))
  expect_true(all(c('tmp_cellID', 'UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'DecVal', 'SVM_class', 'connect_group', 'group_maxCellType') %in% colnames(res1$groupDF_ToFlagTrans)))
  expect_true(all(c('CellId', 'cell_type', 'transcript_num', 'self_celltype', 'score_under_self', 'neighbor_CellId', 'neighbor_celltype', 'score_under_neighbor', 'corrected_CellId') %in% colnames(res1$neighborhoodDF_ToReseg)))
  expect_true(all(c('cells_to_discard', 'cells_to_update', 'cells_to_keep', 'reseg_full_converter') %in% names(res1$reseg_actions)))
  
  # all genes included 
  expect_true(all(rownames(score_GeneMatrix) %in% rownames(res1$updated_perCellExprs)))

})


### run modular wrappers for individual task 

## run segmentation evaluation 
outs <- runSegErrorEvaluation(
  score_GeneMatrix= score_GeneMatrix, 
  transcript_df = mini_transcriptDF, # intracellular transcripts only
  cellID_coln = 'UMI_cellID', 
  transID_coln = 'UMI_transID',
  transGene_coln = 'target',
  spatLocs_colns = c('x','y','z'),
  flagModel_TransNum_cutoff = 50) 

modStats_ToFlagCells <- outs [['modStats_ToFlagCells']]

modStats_ToFlagCells[['flagged']] <- (modStats_ToFlagCells[['lrtest_-log10P']] > flagCell_lrtest_cutoff )
flagged_cells <- modStats_ToFlagCells[['UMI_cellID']][modStats_ToFlagCells[['flagged']]]

test_that("runSegErrorEvaluation(), cell-level segmentation error evaluation, returns the expected outupts", {
  expect_type(outs, "list")
  expect_true(all(c("modStats_ToFlagCells", "transcript_df") %in% names(outs)))
  expect_s3_class(outs$modStats_ToFlagCells, "data.frame")
  expect_s3_class(outs$transcript_df, "data.frame")
  
  expect_true(all(c('UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'tLLR_maxCellType', 'score_tLLR_maxCellType') %in% colnames(res1$updated_transDF)))
  
  # same results as run in perFOV pipeline wrapper
  expect_true(identical(modStats_ToFlagCells, res1$modStats_ToFlagCells))
  
})

## run transcript level evaluation
# some randomness in spatial network generation and SVM! 
groupDF_ToFlagTrans <- runTranscriptErrorDetection(chosen_cells = flagged_cells,
                                                   score_GeneMatrix = score_GeneMatrix, 
                                                   transcript_df = outs[['transcript_df']], 
                                                   cellID_coln = "UMI_cellID", 
                                                   transID_coln = "UMI_transID", 
                                                   # column for transcript score in current cell segment
                                                   score_coln = 'score_tLLR_maxCellType',
                                                   spatLocs_colns = c("x","y","z"),
                                                   model_cutoff = 50, 
                                                   score_cutoff = svmClass_score_cutoff, 
                                                   distance_cutoff = molecular_distance_cutoff, 
                                                   groupTranscripts_method = "dbscan",
                                                   seed_transError = 123)

test_that("runTranscriptErrorDetection(), transcript-level segmentation error evaluation, returns the expected outupts", {
  expect_s3_class(groupDF_ToFlagTrans, "data.frame")
  expect_true(all(c('tmp_cellID', 'UMI_transID', 'UMI_cellID', 'x', 'y', 'z', 'target', 'DecVal', 'SVM_class', 'connect_group', 'group_maxCellType') %in% colnames(groupDF_ToFlagTrans)))
  
  # same results as run in perFOV pipeline wrapper
  expect_true(identical(groupDF_ToFlagTrans[, c('tmp_cellID', 'UMI_transID', 'UMI_cellID', 'x', 'y', 'z', 'target', 'DecVal', 'SVM_class', 'connect_group', 'group_maxCellType')], 
                        res1$groupDF_ToFlagTrans[, c('tmp_cellID', 'UMI_transID', 'UMI_cellID', 'x', 'y', 'z', 'target', 'DecVal', 'SVM_class', 'connect_group', 'group_maxCellType')]))
})

## case 2: run transcript grouping with delaunay 
# input transcript df contains the edge cases of solo transcripts, collinear transcripts
groupDF_ToFlagTrans2 <- runTranscriptErrorDetection(chosen_cells = flagged_cells,
                                                   score_GeneMatrix = score_GeneMatrix, 
                                                   transcript_df = outs[['transcript_df']], 
                                                   cellID_coln = "UMI_cellID", 
                                                   transID_coln = "UMI_transID", 
                                                   # column for transcript score in current cell segment
                                                   score_coln = 'score_tLLR_maxCellType',
                                                   spatLocs_colns = c("x","y","z"),
                                                   model_cutoff = 50, 
                                                   score_cutoff = svmClass_score_cutoff, 
                                                   distance_cutoff = molecular_distance_cutoff, 
                                                   groupTranscripts_method = "delaunay",
                                                   seed_transError = 123)

tmp_idx <- match(groupDF_ToFlagTrans$UMI_transID, groupDF_ToFlagTrans2$UMI_transID)

test_that("runTranscriptErrorDetection() returns the expected outupts with delaunay network for transcript grouping", {
  expect_s3_class(groupDF_ToFlagTrans2, "data.frame")
  expect_true(all(c('tmp_cellID', 'UMI_transID', 'UMI_cellID', 'x', 'y', 'z', 'target', 'DecVal', 'SVM_class', 'connect_group', 'group_maxCellType') %in% colnames(groupDF_ToFlagTrans2)))
  
  # same SVM class, but different grouping 
  expect_true(identical(groupDF_ToFlagTrans2[tmp_idx, 'SVM_class'], groupDF_ToFlagTrans[, 'SVM_class']))
  expect_false(identical(groupDF_ToFlagTrans2[tmp_idx, 'connect_group'], groupDF_ToFlagTrans[, 'connect_group']))
})

rm(tmp_idx, groupDF_ToFlagTrans2)


## get ready to seg refinement
# update the transcript_df with flagged transcript_group
reseg_transcript_df <- merge(outs[['transcript_df']], 
                             groupDF_ToFlagTrans[, c("UMI_transID", 'connect_group','tmp_cellID','group_maxCellType')], 
                             by = "UMI_transID", all.x = TRUE)
# fill in the missing values for unflagged cells
tmp_idx <- which(is.na(reseg_transcript_df[['connect_group']]))
reseg_transcript_df[['connect_group']][tmp_idx] <- rep(0, length(tmp_idx))
reseg_transcript_df[['tmp_cellID']][tmp_idx] <- reseg_transcript_df[["UMI_cellID"]][tmp_idx]
reseg_transcript_df[['group_maxCellType']][tmp_idx] <- reseg_transcript_df[['tLLR_maxCellType']][tmp_idx]
rm(tmp_idx)

# cells or group IDs for neighborhood evaluation 
groups_to_reseg <- unique(groupDF_ToFlagTrans[which(groupDF_ToFlagTrans[['connect_group']]!=0),][['tmp_cellID']])

## run segmentation refinement
finalRes <- runSegRefinement(
  score_GeneMatrix = score_GeneMatrix,  
  chosen_cells = groups_to_reseg, 
  reseg_transcript_df = reseg_transcript_df, 
  reseg_cellID_coln = "tmp_cellID", 
  reseg_celltype_coln = "group_maxCellType", 
  transID_coln = "UMI_transID",
  transGene_coln = "target", 
  transSpatLocs_coln = c('x','y','z'),
  score_baseline = score_baseline, 
  lowerCutoff_transNum = lowerCutoff_transNum, 
  higherCutoff_transNum= higherCutoff_transNum, 
  neighbor_distance_xy = cellular_distance_cutoff,
  distance_cutoff = molecular_distance_cutoff,
  spatialMergeCheck_method = "leidenCut", 
  cutoff_spatialMerge = 0.5,
  return_intermediates = TRUE,
  return_perCellData = TRUE, 
  includeAllRefGenes = TRUE, 
  seed_segRefine = 123 
)

test_that("runSegRefinement() returns the expected outupts", {
  expect_type(finalRes, "list")
  expect_true("neighborhoodDF_ToReseg" %in% names(finalRes))
  expect_true("reseg_actions" %in% names(finalRes))
  expect_true("updated_transDF" %in% names(finalRes))
  expect_true("updated_perCellDT" %in% names(finalRes))
  expect_true("updated_perCellExprs" %in% names(finalRes))
  
  # Test the class and structure of the data frames and matrices
  expect_true(is.data.frame(finalRes$neighborhoodDF_ToReseg))
  expect_type(finalRes$reseg_actions, "list")
  expect_true(is.data.frame(finalRes$updated_transDF))
  expect_true("data.table" %in% class(finalRes$updated_perCellDT))
  expect_true("dgCMatrix" %in% class(finalRes$updated_perCellExprs))
  
  # Test the presence of specific columns in the data frames or elements in list
  expect_true(all(c('UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'updated_cellID', 'updated_celltype') %in% colnames(finalRes$updated_transDF)))
  expect_true(all(c('CellId', 'cell_type', 'transcript_num', 'self_celltype', 'score_under_self', 'neighbor_CellId', 'neighbor_celltype', 'score_under_neighbor', 'corrected_CellId') %in% colnames(finalRes$neighborhoodDF_ToReseg)))
  expect_true(all(c('cells_to_discard', 'cells_to_update', 'cells_to_keep', 'reseg_full_converter') %in% names(finalRes$reseg_actions)))
  
  # all genes included 
  expect_true(all(rownames(score_GeneMatrix) %in% rownames(finalRes$updated_perCellExprs)))
  
  # same results as run in perFOV pipeline wrapper
  expect_true(identical(finalRes$neighborhoodDF_ToReseg, res1$neighborhoodDF_ToReseg))
  expect_true(identical(finalRes$reseg_actions, res1$reseg_actions))
  expect_true(identical(finalRes$updated_perCellDT, res1$updated_perCellDT))
  expect_true(identical(finalRes$updated_perCellExprs, res1$updated_perCellExprs))
  expect_true(identical(finalRes$updated_transDF[, c("UMI_transID", "UMI_cellID", "updated_cellID", "updated_celltype")], 
                        res1$updated_transDF[, c("UMI_transID", "UMI_cellID", "updated_cellID", "updated_celltype")]))
})

rm(res1, outs, modStats_ToFlagCells, groupDF_ToFlagTrans, reseg_transcript_df, groups_to_reseg, finalRes)

