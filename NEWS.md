# FastReseg

# FastReseg 1.1.0

* Add fallback for spatial network building of colinear input points as identified in issue #45
* Tag release version `v1.1.0`

# FastReseg 1.0.4

* Ensure probe order when combining results across FOVs as identified in issue #55.

# FastReseg 1.0.3

* Include an option to turn off Y axis inversion in pipeline functions.  
* Add support for gzip compressed transcript data file. 
* Expand version compatibility to igraph 2.1.0 and above.

# FastReseg 1.0.2

* Fix 2d data failure with `runTranscriptErrorDetection()` as identified in issue #39. 

# FastReseg 1.0.1

* Update DESCRIPTION with github remotes and include installation instructions for R 4.3.x

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
