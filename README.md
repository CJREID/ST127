# ST127_project
Scripts to run analysis associated with 2021 ST127 manuscript

## Overview
The following repository allows conscientious readers of the ST127 manuscript to reproduce all the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`data`__ (the contents of which should be self explanatory) and generates a __`data/outputs`__ folder with __`figures`__ and __`data`__ subdirectories.


## Installation
### Software requirements
These scripts are currently functional on mac OS Big Sur 11.5.2 using RStudio 1.4.1106 and R version 4.0.5. We cannot guarantee they will work on other distributions of R or RStudio. Your OS should not be an issue provided you use these versions of R and RStudio though.

### Packages required
You will need to install the following packages and versions to work with the scripts:
- data.table_1.14.0 
- magrittr_2.0.1    
- paletteer_1.4.0   
- gridExtra_2.3     
- abricateR_0.1.0  
- RColorBrewer_1.1-2
- ggtree_3.1.0      
- pheatmap_1.0.12   
- reshape2_1.4.4    
- ggpubr_0.4.0      
- forcats_0.5.1     
- stringr_1.4.0     
- dplyr_1.0.7       
- purrr_0.3.4       
- readr_2.0.1       
- tidyr_1.1.4       
- tibble_3.1.5      
- ggplot2_3.3.5     
- tidyverse_1.3.1   

### Issues with ggtree
There have recently been some issues with ggtree in the way it interacts with dplyr. The solution is to install the latest version of ggtree directly from github instead of via BiocManager. You can do this in the console on RStudio with the __`remotes`__ package like so:
```
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
```

## Usage
Clone this repository
```
git clone https://github.com/CJREID/ST127_project.git
cd ST127_project
pwd
```
Open the ST127_analyis.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 16 to the output of __`pwd`__ above and save the script.

Run the ST127_analyis.R script and your __`outputs`__ folder will populate with figures and tables.

## Outputs
### Figures
1. Figure 1. Summary of geography and sources of collection
2. Figure 2. Core gene maximum-likelihood phylogenetic tree with metadata (Legends manually edited for publication)
3. Figure 3. Summary of clusters, sources and pUTI89-like plasmid carriage.
4. Figure 4. Heatmap of low SNP distances between human and companion animal isolates
5. Figure 5. pUTI89 heatmap visualisation (Exported in 2 parts; manually edited for publication)
6. Figure 6. Comparison of antimicrobial resistance gene carriage in relation to pUTI89 and intI1 carriage

### Supplementary
#### Tables
1. Table S1. Metadata, accession numbers and gene screening results for 299 ST127 isolates used in this study
2. Table S2. Integron co-carriage data; summary of assembly contigs containing the intI1 gene in conjunction with ARGs

#### Figures
1. Fig S1. Presence/absence of ARGs mapped to core gene phylogeny
2. Fig S2. Presence/absence of VAGs mapped to core gene phylogeny
3. Fig S3. Presence/absence of plasmid-associated genes mapped to core gene phylogeny

