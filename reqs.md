#### Reqs for fastReseg_flag_all_errors:

`fastReseg_flag_all_errors()` is a wrapper function to process multiple files of one dataset for segmentation error detection in transcript level. The function reformats the individual transcript `data.frame` to have unique IDs and a global coordinate system and save into disk, then scores each cell for segmentation error and flags transcripts that have low goodness-of-fit to current cells.

##### Inputs:

-   Counts matrix for entire data set, cells X genes.
-   Either a vector of cluster assignments for each cell, or a matrix of genes X clusters reference profiles that could be used for internal cluster assignment.
-   Either a `data.frame` of transcript level information with unique cell ids, or a `data.frame` with each row for each individual file of per FOV transcript `data.frame` within which the coordinates and cell ids are unique, columns include the file path of per FOV transcript `data.frame` file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set.
-   additional arguments describing input data structures and output data format.
-   additional arguments for finer control on error detection and flagging.

##### Outputs:

A list, with the following elements:

-   refProfiles: a genes x clusters matrix of cluster-specific reference profiles used in resegmenation pipeline
-   baselineData: a list of two matrices in cluster x percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.
-   ctrl_genes: a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not `NULL`.
-   combined_modStats_ToFlagCells: a `data.frame` for spatial modeling statistics of each cell for all cells in the data set, output of `score_cell_segmentation_error` function
-   combined_flaggedCells: a list with each element to be a vector of `UMI_cellID` for cells flagged for potential cell segmentation errors within each FOV
-   trimmed_perCellExprs: a gene x cell count sparse matrix where all putative contaminating transcripts are trimmed, return when `return_trimmed_perCell = TRUE`
-   flagged_transDF_list: a list of per-FOV transcript `data.frame` with flagging information in `SVM_class` column, return when `transDF_export_option = 2`

When `transDF_export_option =1`, the function would save the each per FOV output as individual file in `path_to_output` directory; `flagged_transDF`, `modStats_ToFlagCells` and `classDF_ToFlagTrans` would be saved as csv file, respectively.

-   flagged_transDF: a transcript `data.frame` for each FOV, with columns for unique IDs of transcripts `UMI_transID` and cells `UMI_cellID`, for global coordinate system `x`, `y`, `z`, and for the goodness-of-fit in original cell segment `SMI_class`; the original per FOV cell ID and pixel/index-based coordinates systems are saved under columns, `CellId`, `pixel_x`, `pixel_y`, `idx_z`.
-   modStats_ToFlagCells: a `data.frame` for spatial modeling statistics of each cell, output of `score_cell_segmentation_error()` function.
-   classDF_ToFlagTrans: `data.frame` for the class assignment of transcripts within putative wrongly segmented cells, output of `flag_bad_transcripts()` functions.

#### Reqs for fastReseg_full_pipeline:

`fastReseg_full_pipeline()` is a wrapper for full resegmentation pipeline using internal reference profiles and cutoffs. This function first estimates proper reference profiles and cutoffs from the provided data and then use `fastReseg_perFOV_full_process()` function to process each transcript `data.frame`. For each transcript `data.frame`, the pipeline would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including triming to extracellular space, merging to neighbor cell or labeling as new cell.

##### Inputs:

-   Counts matrix for entire data set, cells X genes.
-   Either a vector of cluster assignments for each cell, or a matrix of genes X clusters reference profiles that could be used for internal cluster assignment.
-   Either a `data.frame` of transcript level information with unique cell ids, or a `data.frame` with each row for each individual file of per FOV transcript `data.frame` within which the coordinates and cell ids are unique, columns include the file path of per FOV transcript `data.frame` file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set.
-   additional arguments describing input data structures and output data format.
-   additional arguments for finer control on error detection, flagging and correction.

##### Outputs:

A list, with the following elements:

-   refProfiles: a genes X clusters matrix of cluster-specific reference profiles used in resegmenation pipeline.
-   baselineData: a list of two matrice in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.
-   cutoffs_list: a list of cutoffs used in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`.
-   ctrl_genes: a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not `NULL`.
-   updated_perCellDT: a per cell `data.table` with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData = TRUE`.
-   updated_perCellExprs: a gene x cell count sparse matrix for updated transcript `data.frame` after resegmentation, return when `return_perCellData = TRUE`.
-   reseg_actions: a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations()` function, return when `save_intermediates = TRUE`.
-   updated_transDF_list: a list of per-FOV transcript `data.frame` with updated cell segmenation in `updated_cellID` and `updated_celltype` columns, return when `transDF_export_option = 2`.

The pipeline function would save the each per FOV output as individual file in `path_to_output` directory;`updated_transDF`would be saved as csv file. When`save_intermediates = TRUE`, all intermediate files and resegmenation outputs of each FOV would be saved as single `.rds` object in 1 list containing the following elements:

-   modStats_ToFlagCells: a `data.frame` for spatial modeling statistics of each cell, output of `score_cell_segmentation_error()` function, save when `save_intermediates = TRUE`.
-   groupDF_ToFlagTrans: `data.frame` for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts()` and `groupTranscripts_Delaunay()` or `groupTranscripts_dbscan()` functions, save when `save_intermediates = TRUE`.
-   neighborhoodDF_ToReseg: a `data.frame` for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content()` function, save when `save_intermediates = TRUE`.
-   reseg_actions: a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations()` function, save when `save_intermediates = TRUE`.
-   updated_transDF: the updated transcript_df with `updated_cellID` and ``` updated_celltyp`` column based on ```reseg_full_converter`, write to disk when`transDF_export_option =1\`.
-   updated_perCellDT: a per cell `data.table` with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData = TRUE`
-   updated_perCellExprs: a gene x cell count sparse matrix for updated transcript `data.frame` after resegmentation, return when `return_perCellData = TRUE`.

The pipeline would also combine per cell data for all FOVs and return the combined data when `return_perCellData = TRUE`; `updated_perCellDT` and `updated_perCellExprs` would also be saved in a list as single `.rds` object in `path_to_output` directory when `transDF_export_option = 1`.

-   updated_perCellDT: a per cell `data.table` with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData = TRUE`.
-   updated_perCellExprs: a gene x cell count sparse matrix for updated transcript `data.frame` after resegmentation, return when `return_perCellData = TRUE`.

#### Reqs for fastReseg_perFOV_full_process:

`fastReseg_perFOV_full_process()` function is the core wrapper for resegmentation pipeline using transcript score matrix derived from external reference profiles and preset cutoffs. This function would score each transcript based on the provided cell type-specific reference profiles, evaluate the goodness-of-fit of each transcript within original cell segment, identify the low-score transcript groups within cells that has strong spatial dependency in transcript score profile, evaluate the neighborhood environment of low-score transcript groups and perform resegmentation actions including trimming to extracellular space, merging to neighbor cell or labeling as new cell.

##### Inputs:

-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates.
-   additional arguments describing input data structures and output data format.
-   additional arguments for finer control on error detection, flagging and correction.

##### Outputs:

A list, with the following elements:

-   modStats_ToFlagCells: a `data.frame` for spatial modeling statistics of each cell, output of `score_cell_segmentation_error()` function, return when `return_intermediates` = TRUE}
-   groupDF_ToFlagTrans: `data.frame` for the group assignment of transcripts within putative wrongly segmented cells, merged output of `flag_bad_transcripts()` and `groupTranscripts_Delaunay()` or `groupTranscripts_dbscan()` functions, return when `return_intermediates = TRUE`.
-   neighborhoodDF_ToReseg: a `data.frame` for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content` function, return when `return_intermediates = TRUE`.
-   reseg_actions: a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations_leidenCut` function, return when `return_intermediates = TRUE`.
-   updated_transDF: the updated transcript_df with `updated_cellID` and `updated_celltype` column based on reseg_full_converter}
-   updated_perCellDT: a per cell data.table with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData = TRUE`.
-   updated_perCellExprs: a gene x cell count sparse matrix for updated transcript `data.frame` after resegmentation, return when `return_perCellData = TRUE`.

#### Reqs for runPreprocess:

`runPreprocess()` function is a modular wrapper to get baseline data and cutoffs from entire dataset.

##### Inputs:

-   Counts matrix for entire data set, cells X genes.
-   Either a vector of cluster assignments for each cell, or a matrix of genes X clusters reference profiles that could be used for internal cluster assignment.
-   Either a `data.frame` of transcript level information with unique cell ids, or a `data.frame` with each row for each individual file of per FOV transcript `data.frame` within which the coordinates and cell ids are unique, columns include the file path of per FOV transcript `data.frame` file, annotation columns like slide and fov to be used as prefix when creating unique cell_ID across entire data set.
-   additional arguments describing input data structures and output data format
-   additional arguments for finer control on error detection, flagging and correction.
-   optional inputs of external baseline and cutoffs for transcript scores and transcript number to skip calculation based on the provided dataset.

##### Outputs:

A list, with the following elements:

-   clust: vector of cluster assignments for each cell in `counts`, used in caculating `baselineData`.
-   refProfiles: a genes X clusters matrix of cluster-specific reference profiles to use in resegmenation pipeline.
-   baselineData: a list of two matrices in cluster X percentile format for the cluster-specific percentile distribution of per cell value; `span_score` is for the average per molecule transcript tLLR score of each cell, `span_transNum` is for the transcript number of each cell.\
-   cutoffs_list: a list of cutoffs to use in resegmentation pipeline, including, `score_baseline`, `lowerCutoff_transNum`, `higherCutoff_transNum`, `cellular_distance_cutoff`, `molecular_distance_cutoff`.
-   ctrl_genes: a vector of control genes whose transcript scores are set to fixed value for all cell types, return when `ctrl_genes` is not `NULL`.
-   score_GeneMatrix: a gene x cell-type score matrix to use in resegmenation pipeline, the scores for `ctrl_genes` are set to be the same as `svmClass_score_cutoff`.
-   processed_1st_transDF: a list of 2 elements for the intracellular and extracellular transcript `data.frame` of the processed outcomes of 1st transcript file.

The `cutoffs_list` is a list containing

-   score_baseline: a named vector of score baseline under each cell type listed in `refProfiles` such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence.
-   lowerCutoff_transNum: a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is.
-   higherCutoff_transNum: a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
-   cellular_distance_cutoff: maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, unit in micron. - molecular_distance_cutoff: maximum molecule-to-molecule distance within connected transcript group, unit in micron.

#### Reqs for runSegErrorEvaluation:

`runSegErrorEvaluation()` function is a modular wrapper to flag cell segmentation error.

##### Inputs:

-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` for each transcript with columns for transcript_id, target or gene name, cell_id, spatial coordinates.
-   additional arguments describing input data structures.
-   additional arguments for finer control on error detection, and flagging at cell level.

##### Outputs:

A list, with the following elements:

-   modStats_ToFlagCells: a `data.frame` contains evaluation model statistics in columns for each cell's potential to have segmentation error.
-   transcript_df: transcript `data.frame` with 2 additional columns: `tLLR_maxCellType` for cell types of maxmium transcript score under current segments and `score_tLLR_maxCellType` for the corresponding transcript score for each transcript.

#### Reqs for runTranscriptErrorDetection:

`runTranscriptErrorDetection()` function is a modular wrapper to identify transcript groups of poor fit to current cell segments in space.

##### Inputs:

-   a vector of cell ID for cells of interest, typically the cells flagged by `runSegErrorEvaluation()` function.
-   a gene x cell-type matrix of log-like score of gene in each cell type
-   a `data.frame` for each transcript with columns for transcript_id, cell_id, spatial coordinates and transcript score, typically the outputs of `runSegErrorEvaluation()` function.
-   a string indicating how to group transcripts in space, use either "dbscan" or "delaunay" method.
-   additional arguments describing input data structures.
-   additional arguments for finer control on error detection at transcript level.

##### Outputs:

a `data.frame` for transcripts in cells of interest only, containing information for transcript score classifications and spatial group assignments as well as new cell/group ID for downstream resegmentation.

#### Reqs for prepResegDF:

`prepResegDF()` function is a supporting function for `fastReseg_perFOV_full_process()` function, combine `runTranscriptErrorDetection()` output with transcript `data.frame` to prep for `runSegRefinement()`

##### Inputs:

-   a `data.frame` for each transcript with columns for transcript_id, target or gene name, original cell_id, spatial coordinates and cell type under which each transcript group gives the maximum transcript score, typically the outputs of `runSegErrorEvaluation()` function.
-   a `data frame` for transcripts in cells of interest only, with columns for `connect_group`,`tmp_cellID`,`group_maxCellType`, output of `runTranscriptErrorDetection()` function

##### Outputs:

A list, with the following elements: - reseg_transcript_df: `data.frame` with transcript_id, target or gene name, x, y, cell_id for all transcript groups in `tmp_cellID` column and the cell type of maximum transcript scores for each transcript group in `group_maxCellType` column. - groups_to_reseg: vector of chosen transcript groups need to be evaluate for re-segmentation.

#### Reqs for runSegRefinement:

`runSegRefinement()` function is a modular wrapper to evaluate transcript groups in neighborhood, decide resegmentation operations and execute.

##### Inputs:

-   a vector of transcript group ids for cells of interest, typically the transcript groups in the cells flagged by `runTranscriptErrorDetection()` function.
-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` for each transcript with columns for transcript_id, target or gene name, spatial coordinates, original cell id, group id for all transcript groups and the cell type of maximum transcript scores for each transcript group, typically the outputs of `prepResegDF()` function.
-   a named vector of score baseline for all cell type such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence, typically the output of `runPreprocess()` function.
-   a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is, typically the output of `runPreprocess()` function.
-   a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type, typically the output of `runPreprocess()` function.
-   additional arguments describing input data structures and output data format.
-   additional arguments for finer control on error correction.

##### Outputs:

A list, with the following elements:

-   updated_transDF: the updated `transcript_df` with `updated_cellID` and `updated_celltype` column based on `reseg_full_converter`.
-   neighborhoodDF_ToReseg: a `data.frame` for neighborhood enviornment of low-score transcript groups, output of `get_neighborhood_content()` function, return when `return_intermediates = TRUE`.
-   reseg_actions: a list of 4 elements describing how the resegmenation would be performed on original `transcript_df` by the group assignment of transcripts listed in `groupDF_ToFlagTrans`, output of `decide_ReSegment_Operations()` function, return when `return_intermediates = TRUE`.
-   updated_perCellDT: a per cell `data.table` with mean spatial coordinates, new cell type and resegmentation action after resegmentation, return when `return_perCellData = TRUE`.
-   updated_perCellExprs: a gene x cell count sparse matrix for updated transcript `data.frame` after resegmentation, return when `return_perCellData = TRUE`.

#### Reqs for score_cell_segmentation_error:

`score_cell_segmentation_error()` function scores each cell for how much their transcripts change their goodness-of-fit over space. It is a supporting function for `runSegErrorEvaluation()` modular wrapper function.

##### Inputs:

-   a `data.frame` of transcript_id, cell_id, transcript score, spatial coordinates.
-   a cutoff of transcript number to do spatial modeling
-   additional arguments describing input data structures.

##### Outputs:

A `data.frame` with columns for cell_id, number of transcripts in given cell, r.squared value of alternative model where transcript score has significant spatial dependency, lrtest chi-squared value and p-value for lrtest probability larger than chi-squared value.

#### Reqs for flag_bad_transcripts:

`flag_bad_transcripts()` function finds out the spatially connected transcripts among chosen_transcripts based on SVM spatial model which scores each cell for how much their transcripts change their goodness-of-fit over space. It is a supporting function for `runTranscriptErrorDetection()` modular wrapper function.

##### Inputs:

-   a vector of cell ids for cells of interest, typically the cells flagged by `runSegErrorEvaluation()` function.
-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` for each transcript with columns for transcript_id, transcript score, spatial coordinates, cell_id.
-   additional arguments describing input data structures.
-   additional arguments for finer control on spatial modeling and error detection.

##### Outputs:

A `data.frame` with columns for transcript_id, target or gene name, cell_id, spatial coordinates, transcript score, SVM class with 0 for below cutoff and 1 for above cutoff, decision values of svm model output, new cell type for each transcript groups within each cells.

#### Reqs for groupTranscripts_Delaunay:

`groupTranscripts_Delaunay()` function groups the flagged transcript within each cell based on spatial connectivity of their transcript delaunay network. It's a supporting function for `runTranscriptErrorDetection()` modular wrapper function.

##### Inputs:

-   a vector of transcript ids of interest, typically the cells flagged by `flag_bad_transcripts()` function.
-   a transcript `data.frame` with transcript_id, target or gene Name, spatial coordinates and cell_id
-   a configure list on controlling the spatial network generation.For more details, see manual for `createSpatialDelaunayNW_from_spatLocs()` function for more details.
-   the maximum distance allowed within connected transcript group

##### Outputs:

A `data.frame` for transcripts of interest only, containing columns for cell ids, transcript ids, spatial coordinates, group id for spatially connected transcripts.

#### Reqs for createSpatialDelaunayNW_from_spatLocs:

`createSpatialDelaunayNW_from_spatLocs()` function generates delaunay network based on provided config and spatial location. It is a supporting function for `groupTranscripts_Delaunay()`.

##### Inputs:

-   a `data.frame` for spatial location of each entry for cell or transcript
-   a configure list on controlling the spatial network generation.For more details, see the manual for `GiottoClass::createSpatialNetwork`.

##### Outputs:

a `delaunay_network_Obj`, a spatial network object created by `GiottoClass` functions. For more details, see the manual for `GiottoClass::createSpatialNetwork`.

#### Reqs for groupTranscripts_dbscan:

`groupTranscripts_dbscan()` functiongroups the flagged transcript within each cell based on spatial clustering using `dbscan`. It's a supporting function for `runTranscriptErrorDetection()` modular wrapper function.

##### Inputs:

-   a vector of transcript ids of interest, typically the transcripts of the cells flagged by `flag_bad_transcripts()` function.
-   a transcript `data.frame` with transcript_id, target or gene Name, spatial coordinates and cell_id
-   the maximum distance allowed within connected transcript group

##### Outputs:

A `data.frame` for transcripts of interest only, containing columns for cell ids, transcript ids, spatial coordinates, group id for spatially connected transcripts.

#### Reqs for get_neighborhood_content:

`get_neighborhood_content()` function finds neighbor cells with transcripts that are direct neighbor of chosen cells, check transcript score under neighbor cell type, return neighborhood information. It is a supporting function for `runSegRefinement()` modular wrapper function.

##### Inputs:

-   a vector of transcript group ids for cells of interest, typically the transcript groups in the cells flagged by `runTranscriptErrorDetection()` function.
-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a named vector of transcript score baseline for all cell types.
-   a `data.frame` for each transcript with columns for transcript_id, spatial coordinates, group ids for all transcript groups including the original cell ids for cells not being flagged.
-   additional arguments describing input data structures.
-   additional arguments for finer control on defining neighborhood.

##### Outputs:

A neighborhood information `data.frame` for transcript groups of interest, containing the following columns:

-   CellId: original cell or transcript group id of chosen cells.
-   cell_type: original cell type of chosen cells.
-   transcript_num: number of transcripts in chosen cells.
-   self_celltype: cell type give maximum score for query cell only.
-   score_under_self: score in query cell under its own maximum celltype.
-   neighbor_CellId: cell id of neighbor cell whose cell type gives maximum score in query cell among all neighbors, not including query cell itself.
-   neighbor_celltype: cell type that gives maximum score in query cell among all non-self neighbor cells. -score_under_neighbor: score in query cell under neighbor_celltype.

#### Reqs for decide_ReSegment_Operations:

`decide_ReSegment_Operations()` function evaluates neighborhood information against score and transcript number cutoff to decide the resegmetation operations. Use either leiden clustering or geometry statistics to determine whether a merge event is allowed. It is a supporting function for `runSegRefinement()` modular wrapper function.

##### Inputs:

-   a neighborhood information `data.frame` for transcript groups of interest, typically the output of `get_neighborhood_content()` function.
-   a string indicating use either "leidenCut" (in 2D or 3D) or "geometryDiff" (in 2D only) method to determine whether a cell pair merging event is allowed in space.
-   a named vector of score baseline for all cell type such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence, typically the output of `runPreprocess()` function.
-   a named vector of transcript number cutoff under each cell type such that higher than the cutoff is required to keep query cell as it is, typically the output of `runPreprocess()` function.
-   a named vector of transcript number cutoff under each cell type such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type, typically the output of `runPreprocess()` function.
-   a cutoff on spatial constraint on a valid merging event between two source transcript groups.
-   additional arguments describing input data structures.

##### Outputs:

A list, with the following elements:

-   cells_to_discard: a vector of cell ID that should be discarded during resegmentation.
-   cells_to_update: a named vector of cell ID where the cell ID in name would be replaced with cell_ID in value.
-   cells_to_keep: a vector of cell ID that should be kept as it is.
-   reseg_full_converter: a single named vector of cell ID to update the original cell ID, assign `NA` for `cells_to_discard`.

#### Reqs for update_transDF_ResegActions:

`update_transDF_ResegActions()` function updates transcript `data.frame` based on resegmentation action, calculates the new cell type and mean per cell spatial coordinates. It is a supporting function for `runSegRefinement()` modular wrapper function.

##### Inputs:

-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` for each transcript to be updated, containing columns for transcript_id, target or gene name, spatial coordinates, original cell id, group id for all transcript groups and the cell type of maximum transcript scores for each transcript group, typically the outputs of `prepResegDF()` function.
-   a single named vector of cell ID to update the original cell ID, typicallly outputs of `decide_ReSegment_Operations()` function.
-   additional arguments describing input data structures and output data format.

##### Outputs:

A list, with the following elements:

-   updated_transDF: the updated transcript_df with `updated_cellID` and `updated_celltype` column based on `reseg_full_converter`.
-   perCell_DT: a per cell `data.table` with mean spatial coordinates and new cell type when `return_perCellDF = TRUE`.
-   perCell_expression: a gene x cell count sparse matrix for updated transcript `data.frame` when `return_perCellDF = TRUE`.

#### Reqs for prepare_perFOV_transDF:

`prepare_perFOV_transDF()` function convert per FOV unique IDs and spatial coordinates of cells and transcripts to the ones unique for the entire dataset. It also converts pixel coordinates to um. It is supporting function for `runPreprocess()` modular wrapper function, `fastReseg_flag_all_errors()` and `fastReseg_full_pipeline()` pipeline wrapper functions.

##### Inputs:

-   a `data.frame` for per FOV transcript level information within which the coordinates and cell ids are unique.
-   a named vector of fov 2D coordinates.
-   additional arguments describing input data structures.
-   additional arguments for finer control on how to stitch per FOV data into entire dataset of multiple FOVs.

##### Outputs:

A list, with the following elements:

-   intraC: a `data.frame` for intracellular transcript, `UMI_transID` and `UMI_cellID` as column names for unique transcript_id and cell_id, `target` as column name for target gene name.
-   extraC: a `data.frame` for extracellular transcript, same structure as the `intraC` data frame in returned list.

#### Reqs for get_baselineCT:

`get_baselineCT()` function gets cluster-specific quantile distribution of transcript number and per cell per molecule transcript score in the provided cell x gene expression matrix based on the reference profiles and cell cluster assignment. The function would also recommend the cutoff for transcript score and transcript number to be used in re-segmentation pipeline based on the calculated quantile distribution. It is supporting function for `runPreprocess()` modular wrapper function.

##### Inputs:

-   a matrix of genes X clusters reference profiles that could be used for internal cluster assignment.
-   Counts matrix for entire data set, cells X genes.
-   optional input of external cluster assignments for each cell to skip nternal cluster assignment based on the provided reference profiles.

##### Outputs:

A list, with the following elements:

-   span_score: a matrix of average transcript tLLR score per molecule per cell for 22 distinct cell types in rows, percentile at (0%, 25%, 50%, 75%, 100%) in columns.
-   span_transNum: a matrix of transcript number per cell for each distinct cell types in row, percentile at (0%, 25%, 50%, 75%, 100%) in columns.
-   score_baseline: a named vector of 25% quantile of cluster-specific per cell transcript score, to be used as score baseline such that per cell transcript score higher than the baseline is required to call a cell type of high enough confidence,
-   lowerCutoff_transNum: a named vector of 25% quantile of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that higher than the cutoff is required to keep query cell as it is.
-   higherCutoff_transNum: a named vector of median value of cluster-specific per molecule per cell transcript number, to be used as transcript number cutoff such that lower than the cutoff is required to keep query cell as it is when there is neighbor cell of consistent cell type.
-   clust_used: a named vector of cluster assignments for each cell used in baseline calculation, cell_ID in `counts` as name.

#### Reqs for choose_distance_cutoff:

`choose_distance_cutoff()` function chooses appropriate cellular distance cutoff and molecular distance cutoff based on input transcript `data.frame` for downstream resegmentation; cellular distance cutoff is defined as the search radius of direct neighbor cell, while molecular distance cutoff is defined as the maximum distance between two neighbor transcripts from same source cells. It is supporting function for `runPreprocess()` modular wrapper function.

##### Inputs:

-   a `data.frame` for per FOV transcript level information within which the coordinates and cell ids are unique.
-   additional arguments describing input data structures.
-   additional arguments for finer control on sub-sampling the input data for molecular distance cutoff estimation.

##### Outputs:

A list, with the following elements:

-   cellular_distance_cutoff: maximum cell-to-cell distance in x, y between the center of query cells to the center of neighbor cells with direct contact, same unit as input spatial coordinate.
-   perCell_coordDT: a `data.table` with cell in row, spatial XY coordinates of centroid and dimensions of bounding box in column.
-   molecular_distance_cutoff: maximum molecule-to-molecule distance within connected transcript group, same unit as input spatial coordinate; return if `run_molecularDist = TRUE`.
-   distance_profile: a named vector for the quantile profile of minimal molecular distance between transcripts belong to different cells at step size of 10% quantile; return if `run_molecularDist = TRUE`,

#### Reqs for scoreGenesInRef:

`scoreGenesInRef()` function calculates log-likilhood score of each gene based on reference expression profiles and returns the centered score matrix. It is utility function used by other wrapper functions.

##### Inputs:

-   a vector of gene name to score
-   a gene X cell_type expression matrix for reference profiles
-   flag to center the score matrix per gene before return

##### Outputs:

-   a gene X cell type matrix of loglik score for each gene under each cell type.

#### Reqs for getCellType_maxScore:

`getCellType_maxScore()` function gets the cell type give maximum transcript score. It is utility function used by other wrapper functions.

##### Inputs:

-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` of transcript level information with unique cell and transcript ids.
-   additional arguments describing input data structures.

##### Outputs:

a named vector with cell type in values and cell ID in names.

#### Reqs for getScoreCellType_gene:

`getScoreCellType_gene()` function gets each transcript's score based on score matrix and chosen cell-type. It is utility function used by other wrapper functions.

##### Inputs:

-   a gene x cell-type matrix of log-like score of gene in each cell type.
-   a `data.frame` of transcript level information with unique cell and transcript ids, and cell type.
-   additional arguments describing input data structures.

##### Outputs:

a named vector with score of given cell type in values and transcript_id in names

#### Reqs for estimate_MeanProfile:

`estimate_MeanProfile()` function estimates the mean profile of each cluster, given the input cell type assignments. It is utility function used by other wrapper functions.

##### Inputs:

-   Counts matrix for entire data set, cells X genes.
-   a vector of cluster assignments for each cell.
-   a vector of scaling factors for each cell.
-   a numeric value for expected background in count matrix

##### Outputs:

A matrix of genes X clusters profiles observed in dataset.
