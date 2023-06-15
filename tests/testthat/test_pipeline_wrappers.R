library(testthat)
library(FastReseg)

# example data
data("example_CellGeneExpr")
data("example_clust")
data("example_refProfiles")

dataDir = gsub("tests/testthat", "data", getwd())

# create `transDF_fileInfo` for multiple per FOV transcript data.frame 
# coordinates for each FOV, `stage_x` and `stage_y`, should have units in micron.
transDF_fileInfo <- data.frame(file_path = fs::path(dataDir, 
                                                    c("Run4104_FOV001__complete_code_cell_target_call_coord.csv",
                                                      "Run4104_FOV002__complete_code_cell_target_call_coord.csv")),
                               slide = c(1, 1),
                               fov = c(1,2),
                               stage_X = 1000*c(5.13, -2.701),
                               stage_Y = 1000*c(-0.452, 0.081))

# set cutoffs
flagCell_lrtest_cutoff <- 5
svmClass_score_cutoff <- -2

## flag error in multiple files
# case 1: known `counts` and `clust`
outDir1 <- "res1f_multiFiles"
errorFlagRes1 <- fastReseg_flag_all_errors(counts = example_CellGeneExpr,
                                           clust = example_clust,
                                           refProfiles = NULL,
                                           transcript_df = NULL,
                                           transDF_fileInfo = transDF_fileInfo,
                                           filepath_coln = 'file_path',
                                           prefix_colns = c('slide','fov'),
                                           fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
                                           pixel_size = 0.18, 
                                           zstep_size = 0.8,
                                           transID_coln = NULL, # row index as transcript_id
                                           transGene_coln = "target",
                                           cellID_coln = "CellId",
                                           spatLocs_colns = c("x","y","z"),
                                           extracellular_cellID = c(0), 
                                           
                                           flagCell_lrtest_cutoff = flagCell_lrtest_cutoff, 
                                           svmClass_score_cutoff = svmClass_score_cutoff, 
                                           path_to_output = outDir1,
                                           return_trimmed_perCell = TRUE)

outFiles <- list.files(path = outDir1)
flagged_transDF <- read.csv(fs::path(outDir1, "1_flagged_transDF.csv"))

test_that("fastReseg_flag_all_errors() returns the expected output", {
  expect_type(errorFlagRes1, "list")
  expect_true("refProfiles" %in% names(errorFlagRes1))
  expect_true("baselineData" %in% names(errorFlagRes1))
  expect_true("combined_modStats_ToFlagCells" %in% names(errorFlagRes1))
  expect_true("combined_flaggedCells" %in% names(errorFlagRes1))
  expect_true("trimmed_perCellExprs" %in% names(errorFlagRes1))
  
  # Test the class and structure of the matrices, data frames, and vectors
  expect_true(is.matrix(errorFlagRes1$refProfiles))
  expect_type(errorFlagRes1$baselineData, "list")
  expect_true(is.data.frame(errorFlagRes1$combined_modStats_ToFlagCells))
  expect_type(errorFlagRes1$combined_flaggedCells, "list")
  expect_true("dgCMatrix" %in% class(errorFlagRes1$trimmed_perCellExprs))

  # Test the presence of specific elements within the errorFlagRes1
  expect_true(all(c("span_score", "span_transNum")  %in%  names(errorFlagRes1$baselineData)))
  expect_length(errorFlagRes1$combined_flaggedCells, nrow(transDF_fileInfo))
  expect_true(all(c('UMI_cellID', 'tLLR_maxCellType', 'transcript_num', 'lrtest_Pr', 'lrtest_-log10P', 'flagged') %in% colnames(errorFlagRes1$combined_modStats_ToFlagCells)))
  
  expect_true(identical(unique(errorFlagRes1$combined_modStats_ToFlagCells$file_idx), 
                        seq_len(nrow(transDF_fileInfo))))
  
  # Test if the per FOV output files are saved in the specified directory
  expect_true(file.exists(outDir1))
  expect_true(all(paste0(seq_len(nrow(transDF_fileInfo)), "_modStats_ToFlagCells.csv") %in% outFiles))
  expect_true(all(paste0(seq_len(nrow(transDF_fileInfo)), "_classDF_ToFlagTrans.csv") %in% outFiles))
  expect_true(all(paste0(seq_len(nrow(transDF_fileInfo)), "_flagged_transDF.csv") %in% outFiles))
  
  # check on 1st transDF
  expect_true(all(c('UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'tLLR_maxCellType', 'score_tLLR_maxCellType', 'SVM_class') %in% colnames(flagged_transDF)))
  expect_true(all(flagged_transDF[flagged_transDF$SVM_class ==0, 'UMI_cellID'] %in% errorFlagRes1$combined_flaggedCells[[1]]))
})


## seg refinement in multiple files
# case 2: known `counts` and `clust`, unknown cutoffs
outDir2 <- "res2_multiFiles"
res2 <- fastReseg_full_pipeline(counts = example_CellGeneExpr,
                                clust = example_clust,
                                refProfiles = NULL,

                                transcript_df = NULL,
                                transDF_fileInfo = transDF_fileInfo,
                                filepath_coln = 'file_path',
                                prefix_colns = c('slide','fov'),
                                fovOffset_colns = c('stage_Y','stage_X'),
                                pixel_size = 0.18,
                                zstep_size = 0.8,
                                transID_coln = NULL,
                                transGene_coln = "target",
                                cellID_coln = "CellId",
                                spatLocs_colns = c("x","y","z"),
                                extracellular_cellID = c(0),
                                
                                molecular_distance_cutoff = NULL,
                                cellular_distance_cutoff = NULL,
                                score_baseline = NULL,
                                lowerCutoff_transNum = NULL,
                                higherCutoff_transNum= NULL,
                                imputeFlag_missingCTs = TRUE,

                                flagCell_lrtest_cutoff = 5,
                                svmClass_score_cutoff = -2,
                                groupTranscripts_method = "dbscan",
                                spatialMergeCheck_method = "leidenCut", 
                                cutoff_spatialMerge = 0.5, 
                                
                                path_to_output = outDir2,
                                save_intermediates = TRUE, 
                                return_perCellData = TRUE,
                                combine_extra = TRUE)
outFiles <- list.files(path = outDir2)
updated_transDF <- read.csv(fs::path(outDir2, "1_updated_transDF.csv"))

test_that("fastReseg_full_pipeline() returns the expected output", {
  expect_type(res2, "list")
  expect_true("refProfiles" %in% names(res2))
  expect_true("baselineData" %in% names(res2))
  expect_true("cutoffs_list" %in% names(res2))
  expect_true("reseg_actions" %in% names(res2))
  expect_true("updated_perCellDT" %in% names(res2))
  expect_true("updated_perCellExprs" %in% names(res2))
  
  # Test the class and structure of the matrices, data frames, and vectors
  expect_true(is.matrix(res2$refProfiles))
  expect_type(res2$baselineData, "list")
  expect_type(res2$cutoffs_list, "list")
  expect_type(res2$reseg_actions, "list")
  expect_true("data.table" %in% class(res2$updated_perCellDT))
  expect_true("dgCMatrix" %in% class(res2$updated_perCellExprs))
  
  # Test the presence of specific elements within the res2
  expect_true(all(c("span_score", "span_transNum")  %in%  names(res2$baselineData)))
  expect_true(all(c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum", "cellular_distance_cutoff", "molecular_distance_cutoff") %in% names(res2$cutoffs_list)))
  expect_true(all(c('cells_to_discard', 'cells_to_update', 'cells_to_keep', 'reseg_full_converter') %in% names(res2$reseg_actions)))
  
  # same baseline and flagging outcomes as the flag only wrapper 
  expect_true(identical(res2$refProfiles, errorFlagRes1$refProfiles))
  expect_true(identical(res2$baselineData, errorFlagRes1$baselineData))
  
  # Test if the per FOV output files are saved in the specified directory
  expect_true(file.exists(outDir2))
  expect_true(all(paste0(seq_len(nrow(transDF_fileInfo)), "_updated_transDF.csv") %in% outFiles))
  expect_true(all(paste0(seq_len(nrow(transDF_fileInfo)), "_each_segRes.RData") %in% outFiles))
  expect_true("combined_updated_perCellDT_perCellExprs.RData" %in% outFiles)
  
  # check on 1st transDF
  expect_true(all(c('UMI_transID', 'UMI_cellID', 'target', 'x', 'y', 'z', 'updated_cellID', 'updated_celltype') %in% colnames(updated_transDF)))
  expect_true(all(c('intraC','trimmed', 'extraC') %in% unique(updated_transDF[["transComp"]])))

})

rm(errorFlagRes1, flagged_transDF, res2, updated_transDF, outFiles)

# option to delete the output directory when finished
if(FALSE){
  unlink(outDir1,recursive=TRUE)
  unlink(outDir2,recursive=TRUE)
}
