# FastReseg

# FastReseg 1.0.0

* Addresses compatibility with new latest.fovs that has additional column for acquisition order
* Fixes error associated with 1-cell-per-fov data and non-standard-formatted transcript data.frame
* Add in requirements and specifications
* Switch to `GiottoClass` package for creating Delaunay spatial networks  
* Move data used in `vignettes` to `inst\extdata` folder, including example raw transcript files and helper functions for SMITAP objects
* Save outputs as `.rds` instead of `.RData` in `fastReseg_full_pipeline()` wrapper

# FastReseg 0.1.1

* First release version used with nanopipeline
* Includes parallel computation for multi-FOV processing with pipeline wrapper functions
