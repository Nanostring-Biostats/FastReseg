#### Specs for fastReseg_flag_all_errors:

-   Returns a genes x clusters matrix of cluster-specific reference profiles used in resegmenation pipeline. -- test: test_pipeline_wrappers.R#L54
-   Returns a list of two matrices in cluster x percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell. -- test: test_pipeline_wrappers.R#L55 -- test: test_pipeline_wrappers.R#L68
-   Returns a `data.frame` for spatial modeling statistics of each cell for all cells in the data set. -- test: test_pipeline_wrappers.R#L56
-   Returns a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV. -- test: test_pipeline_wrappers.R#L57 -- test: test_pipeline_wrappers.R#L70
-   If run with `return_trimmed_perCell = TRUE`,
    -   Returns a gene x cell count sparse matrix where all putative contaminating transcripts are trimmed. -- test: test_pipeline_wrappers.R#L58 -- test: test_pipeline_wrappers.R#L65
-   If run with `transDF_export_option =1`,
    -   Saves each of the per FOV outputs as individual files in `path_to_output` directory; `flagged_transDF`, `modStats_ToFlagCells` and `classDF_ToFlagTrans` would be saved as csv file, respectively. -- test: test_pipeline_wrappers.R#L72-79
    -   The saved `flagged_transDF` csv file should be a transcript `data.frame` with columns for unique IDs of transcripts `UMI_transID` and cells `UMI_cellID`, for global coordinate system `x`, `y`, `z`, and for the goodness-of-fit in original cell segment `SMI_class`. -- test: test_pipeline_wrappers.R#L82

#### Specs for fastReseg_full_pipeline:

-   Returns a genes x clusters matrix of cluster-specific reference profiles used in resegmenation pipeline -- test: test_pipeline_wrappers.R#L129
-   Returns a list of two matrices in cluster x percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell. -- test: test_pipeline_wrappers.R#L130 -- test: test_pipeline_wrappers.R#L145
-   Returns a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`. -- test: test_pipeline_wrappers.R#L131 -- test: test_pipeline_wrappers.R#L146
-   If run with `return_perCellData = TRUE`,
    -   Returns a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation. -- test: test_pipeline_wrappers.R#L133 -- test: test_pipeline_wrappers.R#L141
    -   Returns a gene x cell count sparse matrix for updated transcript data.frame after resegmentation. -- test: test_pipeline_wrappers.R#L134 -- test: test_pipeline_wrappers.R#L142
-   If run with `save_intermediates = TRUE`,
    -   Returns a list of 4 elements describing how the resegmenation would be performed on original transcript `data.frame`. -- test: test_pipeline_wrappers.R#L132 -- test: test_pipeline_wrappers.R#L147
    -   Saves all intermediate files and resegmenation outputs of each FOV as single `.RData` object per FOV. -- test: test_pipeline_wrappers.R#L154-157
-   If run with `save_intermediates = TRUE` and `transDF_export_option =1`,
    -   Saves each of the per FOV outputs as individual files with columns for resegmented outcomes, e.g. `updated_cellID` and `updated_celltype`. -- test: test_pipeline_wrappers.R#L155 -- test: test_pipeline_wrappers.R#L160

#### Specs for runPreprocess:

-   Runs without error when using `counts`, `refProfiles` as input but no `clust`. -- test: test_runPreprocess.R#L22-49
-   Runs without error when using `counts`, `clust` as input but no `refProfiles`. -- test: test_runPreprocess.R#L86-113
-   Runs without error when using transcript `data.frame` directly as input along with `counts`, `refProfiles` but no `clust`. -- test: test_runPreprocess.R#L151-178
-   For all input combinations,
    -   Returns a vector of cluster assignments for each cell. -- test: test_runPreprocess.R#L54 -- test: test_runPreprocess.R#L123 -- test: test_runPreprocess.R#L190
    -   Returns a genes x clusters matrix of cluster-specific reference profiles used in resegmenation pipeline. -- test: test_runPreprocess.R#L55 -- test: test_runPreprocess.R#L124 -- test: test_runPreprocess.R#L144 -- test: test_runPreprocess.R#L191
    -   Returns a list of two matrices in cluster x percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell. -- test: test_runPreprocess.R#L56 -- test: test_runPreprocess.R#L71 -- test: test_runPreprocess.R#L125 -- test: test_runPreprocess.R#L140 -- test: test_runPreprocess.R#L192
    -   Returns a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`. -- test: test_runPreprocess.R#L57 -- test: test_runPreprocess.R#L72 -- test: test_runPreprocess.R#L126 -- test: test_runPreprocess.R#L141 -- test: test_runPreprocess.R#L193-194 -- test: test_runPreprocess.R#L208-209
    -   Returns a gene x cell-type score matrix to use in resegmenation pipeline. -- test: test_runPreprocess.R#L58 -- test: test_runPreprocess.R#L127
    -   Returns a list of 2 elements for the intracellular and extracellular transcript `data.frame` of the processed outcomes of 1st transcrip file. -- test: test_runPreprocess.R#L75-83 -- test: test_runPreprocess.R#L196-203

#### Specs for fastReseg_perFOV_full_process:

-   Returns a updated transcript `data.farme` with `updated_cellID` and `updated_celltype` columns. -- test: test_perFOV_processing.R#L57 -- test: test_perFOV_processing.R#L71
-   If run with `return_perCellData = TRUE`,
    -   Returns a per cell `data.table` with mean spatial coordinates, new cell type and resegmentation action after resegmentation. -- test: test_perFOV_processing.R#L58 -- test: test_perFOV_processing.R#L67
    -   Returns a gene x cell count sparse matrix for the updated transcript data.frame after resegmentation. -- test: test_perFOV_processing.R#L59 -- test: test_perFOV_processing.R#L68
-   If run with `return_intermediates = TRUE`,
    -   Returns a `data.frame` for spatial modeling statistics of each cell. -- test: test_perFOV_processing.R#L53 -- test: test_perFOV_processing.R#L72
    -   Returns a `data.frame` for the group assignment of transcripts within putative wrongly segmented cells. -- test: test_perFOV_processing.R#L54 -- test: test_perFOV_processing.R#L73
    -   Returns a `data.frame` for neighborhood environment of low-score transcript groups. -- test: test_perFOV_processing.R#L55 -- test: test_perFOV_processing.R#L74
    -   Returns a list of 4 elements describing how the resegmenation would be performed on original transcript `data.frame`. -- test: test_perFOV_processing.R#L56 -- test: test_perFOV_processing.R#L75

#### Specs for runSegErrorEvaluation:

-   Returns a `data.frame` contains evaluation model statistics in columns for each cell's potential to have segmentation error. -- test: test_perFOV_processing.R#L102-103 -- test: test_perFOV_processing.R#L109
-   Returns a transcript `data.frame` with 2 additional columns: `tLLR_maxCellType` for cell types of maxmium transcript score under current segments and `score_tLLR_maxCellType` for the corresponding transcript score for each transcript. -- test: test_perFOV_processing.R#L102-106

#### Specs for runTranscriptErrorDetection:

-   Runs without error when using `dbscan` as transcript grouping method. -- test: test_perFOV_processing.R#L115-136
-   Runs without error when using `delaunay` as transcript grouping method. -- test: test_perFOV_processing.R#L140-163
-   For both transcript grouping methods,
-   Returns a transcript `data.frame` containing information for transcript score classifications and spatial group assignments as well as new cell/group ID for downstream resegmentation. -- test: test_perFOV_processing.R#L129-136 -- test: test_perFOV_processing.R#L156-163

#### Specs for runSegRefinement:

-   Returns the updated transcript data.frame after resegmentation with `updated_cellID` and `updated_celltype` columns. -- test: test_perFOV_processing.R#L202 -- test: test_perFOV_processing.R#L214 -- test: test_perFOV_processing.R#L226
-   If run with `return_intermediates = TRUE`,
    -   Returns a `data.frame` for neighborhood environment of low-score transcript groups. -- test: test_perFOV_processing.R#L200 -- test: test_perFOV_processing.R#L215 -- test: test_perFOV_processing.R#L222
    -   Returns a list of 4 elements describing how the resegmenation would be performed on original transcript `data.frame`. -- test: test_perFOV_processing.R#L201 -- test: test_perFOV_processing.R#L216 -- test: test_perFOV_processing.R#L223
-   If run with `return_perCellData = TRUE`,
    -   Returns a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation. -- test: test_perFOV_processing.R#L203 -- test: test_perFOV_processing.R#L210 -- test: test_perFOV_processing.R#L224
    -   Returns a gene x cell count sparse matrix for updated transcript data.frame after resegmentation. -- test: test_perFOV_processing.R#L204 -- test: test_perFOV_processing.R#L211 -- test: test_perFOV_processing.R#L219 -- test: test_perFOV_processing.R#L225
