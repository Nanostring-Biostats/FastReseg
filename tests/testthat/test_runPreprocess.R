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

## run preprocessing 
# case 1: known `counts` and `clust`, calculate all cutoffs 
prep_res1 <- runPreprocess(
  counts = example_CellGeneExpr, 
  clust = example_clust, 
  refProfiles = NULL,
  
  score_baseline = NULL, 
  lowerCutoff_transNum = NULL, 
  higherCutoff_transNum= NULL, 
  imputeFlag_missingCTs = FALSE, 
  ctrl_genes = NULL,
  svmClass_score_cutoff = -2,
  molecular_distance_cutoff = NULL,
  cellular_distance_cutoff = NULL,
  
  transcript_df = NULL, 
  transDF_fileInfo = transDF_fileInfo, 
  filepath_coln = 'file_path', 
  prefix_colns = c('slide','fov'), 
  fovOffset_colns = c('stage_X','stage_Y'), 
  pixel_size = 0.18, 
  zstep_size = 0.8, 
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = 'CellId', 
  spatLocs_colns = c('x','y','z'), 
  extracellular_cellID = 0 
)


test_that("runPreprocess() returns the expected output", {
  expect_type(prep_res1, "list")
  expect_true("clust" %in% names(prep_res1))
  expect_true("refProfiles" %in% names(prep_res1))
  expect_true("baselineData" %in% names(prep_res1))
  expect_true("cutoffs_list" %in% names(prep_res1))
  expect_true("score_GeneMatrix" %in% names(prep_res1))

  # Test the structure and dimensions of the nested lists and matrices
  expect_type(prep_res1$clust, "character")
  expect_true(is.matrix(prep_res1$refProfiles))
  expect_type(prep_res1$baselineData, "list")
  expect_type(prep_res1$cutoffs_list, "list")
  expect_true(is.matrix(prep_res1$score_GeneMatrix))
  
  expect_true(all(sapply(prep_res1$baselineData, is.matrix)))
  expect_true(all(sapply(prep_res1$cutoffs_list, is.numeric)))

  # Test the presence of specific elements within the nested lists and matrices
  expect_true(all(c("span_score", "span_transNum")  %in% names(prep_res1$baselineData)))
  expect_true(all(c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum", "cellular_distance_cutoff", "molecular_distance_cutoff") %in% names(prep_res1$cutoffs_list)))
})

test_that("runPreprocess() returns `processed_1st_transDF` in correct format when either distance cutoff is set to NULL", {
  expect_true("processed_1st_transDF" %in% names(prep_res1))
  expect_type(prep_res1$processed_1st_transDF, "list")
  expect_true(all(c("intraC", "extraC") %in% names(prep_res1$processed_1st_transDF)))
  expect_true(all(sapply(prep_res1$processed_1st_transDF, is.data.frame)))
  
  expect_true(all(sapply(prep_res1$processed_1st_transDF, 
                         function(x) all(c("UMI_cellID","UMI_transID","target","x" ,"y", "z" )  %in% colnames(x)))))
})


# case 2: known `counts`, `refProfiles`, distance cutoffs but unknown `clust` 
prep_res2 <- runPreprocess(
  counts = example_CellGeneExpr, 
  clust = NULL, 
  refProfiles = example_refProfiles,
  
  score_baseline = NULL, 
  lowerCutoff_transNum = NULL, 
  higherCutoff_transNum= NULL, 
  imputeFlag_missingCTs = TRUE, 
  ctrl_genes = NULL,
  svmClass_score_cutoff = -2,
  molecular_distance_cutoff = prep_res1$cutoffs_list$molecular_distance_cutoff,
  cellular_distance_cutoff = prep_res1$cutoffs_list$cellular_distance_cutoff,
  
  transcript_df = NULL, 
  transDF_fileInfo = transDF_fileInfo, 
  filepath_coln = 'file_path', 
  prefix_colns = c('slide','fov'), 
  fovOffset_colns = c('stage_X','stage_Y'), 
  pixel_size = 0.18, 
  zstep_size = 0.8, 
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = 'CellId', 
  spatLocs_colns = c('x','y','z'), 
  extracellular_cellID = 0 
)


test_that("runPreprocess() skips transcript data.frame loading when both distance cutoffs are provided", {
  expect_true(!"processed_1st_transDF" %in% names(prep_res2))
})


test_that("runPreprocess() returns the expected outputs when `clust` is not provided", {
  expect_type(prep_res2, "list")
  expect_true("clust" %in% names(prep_res2))
  expect_true("refProfiles" %in% names(prep_res2))
  expect_true("baselineData" %in% names(prep_res2))
  expect_true("cutoffs_list" %in% names(prep_res2))
  expect_true("score_GeneMatrix" %in% names(prep_res2))
  
  # Test the structure and dimensions of the nested lists and matrices
  expect_type(prep_res2$clust, "character")
  expect_true(is.matrix(prep_res2$refProfiles))
  expect_type(prep_res2$baselineData, "list")
  expect_type(prep_res2$cutoffs_list, "list")
  expect_true(is.matrix(prep_res2$score_GeneMatrix))
  
  expect_true(all(sapply(prep_res2$baselineData, is.matrix)))
  expect_true(all(sapply(prep_res2$cutoffs_list, is.numeric)))
  
  # Test the presence of specific elements within the nested lists and matrices
  expect_true(all(c("span_score", "span_transNum")  %in% names(prep_res2$baselineData)))
  expect_true(all(c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum", "cellular_distance_cutoff", "molecular_distance_cutoff") %in% names(prep_res2$cutoffs_list)))

  # all cell types included
  expect_true(all(colnames(example_refProfiles) %in% colnames(prep_res2$score_GeneMatrix)))
  
  # same distance cutoffs as provided
  expect_true(identical(prep_res2$cutoffs_list[4:5], prep_res1$cutoffs_list[4:5]))
  
})

# case 3: input transcript data.frame directly, unknown `clust` and distance cutoffs, include extra `ctrl_genes` 
data("mini_transcriptDF")

prep_res3 <- runPreprocess(
  counts = example_CellGeneExpr, 
  clust = NULL, 
  refProfiles = example_refProfiles,
  
  score_baseline = NULL, 
  lowerCutoff_transNum = NULL, 
  higherCutoff_transNum= NULL, 
  imputeFlag_missingCTs = TRUE, 
  ctrl_genes = c('NegPrb1', 'NegPrb2'),
  
  svmClass_score_cutoff = -2,
  molecular_distance_cutoff = NULL,
  cellular_distance_cutoff = NULL,
  
  transcript_df = mini_transcriptDF, 
  transDF_fileInfo = NULL, 
  pixel_size = 0.18, 
  zstep_size = 0.8, 
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = 'CellId', 
  spatLocs_colns = c('x','y','z'), 
  extracellular_cellID = 0
)

test_that("runPreprocess() returns expected score matrix with `ctrl_genes`.", {
  expect_true(all(c('NegPrb1', 'NegPrb2') %in% rownames(prep_res3$score_GeneMatrix)))
  
  # same score for other genes, given same counts and refProfiles used 
  expect_true(identical(prep_res3$score_GeneMatrix[rownames(prep_res2$score_GeneMatrix), ], 
                        prep_res2$score_GeneMatrix))
})

test_that("runPreprocess() returns expected outputs when direct `transcript_df` is provided", {
  # outputs based on overall data
  expect_true(identical(prep_res2$clust, prep_res3$clust))
  expect_true(identical(prep_res2$refProfiles, prep_res3$refProfiles))
  expect_true(identical(prep_res2$baselineData, prep_res3$baselineData))
  expect_true(identical(prep_res2$cutoffs_list[c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")], 
                        prep_res3$cutoffs_list[c("score_baseline", "lowerCutoff_transNum", "higherCutoff_transNum")]))
  
  # outputs based on transcript data.frame
  expect_true("processed_1st_transDF" %in% names(prep_res3))
  expect_type(prep_res3$processed_1st_transDF, "list")
  expect_true(all(names(prep_res3$processed_1st_transDF) %in% c("intraC", "extraC")))
  
  # no extracurricular transcript in this example df input
  expect_true(is.data.frame(prep_res3$processed_1st_transDF$intraC))
  expect_true(is.null(prep_res3$processed_1st_transDF$extraC))
  
  expect_true(all(c("UMI_cellID","UMI_transID","target","x" ,"y", "z" ) %in% colnames(prep_res3$processed_1st_transDF$intraC)))
  
  # expect different distance cutoffs given different input transcript df
  expect_false(identical(prep_res3$cutoffs_list[c("cellular_distance_cutoff", "molecular_distance_cutoff")], 
                         prep_res1$cutoffs_list[c("cellular_distance_cutoff", "molecular_distance_cutoff")]))
  
})

rm(prep_res1, prep_res2, prep_res3)

