# FastReseg
 An R package for detection and correction of cell segmentation error based on spatial profile of transcripts
 
### dev notes 
- `findSegmentError_allFiles` function flags segmentation error in transcript data.frame for all the input file path
- `fastReseg_core_externalRef` function performs error detection and segmentation refinement on 1 input transcript data.frame
- `fastReseg_internalRef` function first extracts reference profiles based on input `counts` and `clust` of the whole dataset and then process individual FOVs for segmentation error detection and correction for all files provided. 

### System requirements
- R (>= 3.5.0)
- UNIX, Mac or Windows
- see DESCRIPTION for full dependencies

### Demo
See the "vignettes" folder. 
- `0_flagErrorOnly_on_SMIobject.R` for flagging segmentation errors without correcting. 
- `1_fastReseg_on_SMIobject.R` for runing entire resegmentation workfkow on a given dataset, example dataset.



#### workflow:
![image](https://user-images.githubusercontent.com/62775692/234708514-991badac-a6d4-41cc-9c53-0b2bc92cc56d.png)


### Installation
```
git clone https://github.com/Nanostring-Biostats/FastReseg.git
devtools::install()

```