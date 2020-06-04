### Integrated pharmaco-proteogenomics defines two subgroups in isocitrate dehydrogenase wild-type glioblastoma with prognostic and therapeutic opportunities
This repository provides R scripts to carry out the analyses presented in 'Integrated pharmaco-proteogenomics inform prognostic and therapeutic opportunities in isocitrate dehydrogenase wild-type glioblastoma'.

### Contents
The `code` folder contains three R scripts to reproduce our results.
1) `nature_comm2020_R_functions.r` : R functions for basic data processing, statistical test and plotting.

2) `nature_comm2020_R_sub.r` : R functions to make figures in our paper.

3) `nature_comm2020_R_script.r` : R script showing overall workflow of plotting. This script requires 1) and 2).

The `output` folder contains expected results and figures.

### Dependencies
The main analysis workflow is run in

  `R version 3.5.3`

Additionally it requires the following R libraries.
```
matrixStats  
biomaRt  
pheatmap  
data.table  
tidyr  
plyr  
dplyr  
ConsensusClusterPlus  
riverplot  
digest  
wordspace  
DNAcopy  
grid  
ggplot2  
trackViewer  
gridExtra  
circlize  
graphics  
readxl  
rjson  
RCurl  
stringr  
survival  
grDevices  
ggsignif
 ```
