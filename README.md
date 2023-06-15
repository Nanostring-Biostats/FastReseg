# FastReseg
 An R package for detection and correction of cell segmentation error based on spatial profile of transcripts
 
### dev notes 
- `fastReseg_flag_all_errors()` function flags segmentation error in transcript data.frame for all the input file path.
- `fastReseg_perFOV_full_process()` function performs error detection and segmentation refinement on 1 input transcript data.frame.
- `fastReseg_full_pipeline()` function first extracts reference profiles based on input `counts` and `clust` of the whole dataset and then process individual FOVs for segmentation error detection and correction for all files provided. 

### System requirements
- R (>= 3.5.0)
- UNIX, Mac or Windows
- see DESCRIPTION for full dependencies

### Demo
See the "vignettes" folder. 
- `tutorial.Rmd` and `tutorial.html` for example usages of streamline pipeline wrappers and modular functions for individual task.
- `a__flagErrorOnly_on_SMIobject.R` for flagging segmentation errors without correcting, interfacing `FastReseg` with SMI TAP pipeline (`Giotto`).
- `b__fastReseg_on_SMIobject.R` for runing entire resegmentation workfkow on a given dataset, example dataset, interfacing `FastReseg` with SMI TAP pipeline (`Giotto`).



#### workflow:
![image](vignettes/FastReseg_diagram.png)


### Installation
```
git clone https://github.com/Nanostring-Biostats/FastReseg.git
devtools::install()
```
