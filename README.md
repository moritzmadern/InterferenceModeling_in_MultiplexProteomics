# InterferenceModeling_in_MultiplexProteomics
An R implementation for interference modeling and subsequent interference-correction in MS2-based multiplex proteomics. Further contains a demo dataset + userguide to get familiar with the workflow.



## Contents:

The main R script that performs the entire computational workflow on the basis of specified input parameters:
-  **IM.Rmd**


Functions to be sourced in IM.Rmd:
-  **functions_IM.R**


A demo dataset + userguide:
- **Demo**
Check out Userguide.pdf in the Demo folder for detailed instructions and explanations on the demo and the workflow in general.



## Required Data Input:

- A PSM table. Currently supported are MaxQuants **msms.txt** and Fragpipe's **psm.tsv**. Other formats might require minor adjustments to the script.

- Corresponding **Thermo raw files** of the search, located in a separate Folder.

- The **rawStallion** Windows command line application to read relevant data from Thermo raw files and write to tsv files. Download [here](https://github.com/fstanek/rawStallion).

- An isotopic impurity matrix for isotopic impurity correction. Details on the required format are described in the parameter section of IM.Rmd.



## Data Output:

- a modified PSM table named **modified_PSM.txt**. This PSM table contains normalized reporter ion intensities (suffix `_norm`), normalized interference-corrected reporter ion intensities (suffix `_norm__interference_corrected`), as well as calculated PSM-wise metrics such as Estimated Interference Level (EIL), Precursor Purity Fraction (PPF), and more.



## Session Info

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] msqrob2_1.2.0               QFeatures_1.4.0             MultiAssayExperiment_1.20.0 DESeq2_1.34.0              
 [5] SummarizedExperiment_1.24.0 MatrixGenerics_1.6.0        matrixStats_0.61.0          GenomicRanges_1.46.1       
 [9] GenomeInfoDb_1.30.0         IRanges_2.28.0              limma_3.50.0                MSnbase_2.20.0             
[13] ProtGenerics_1.26.0         S4Vectors_0.32.3            mzR_2.28.0                  Rcpp_1.0.7                 
[17] Biobase_2.54.0              BiocGenerics_0.40.0         cowplot_1.1.1               fields_13.3                
[21] viridis_0.6.2               viridisLite_0.4.0           spam_2.8-0                  doParallel_1.0.17          
[25] iterators_1.0.13            foreach_1.5.2               rlist_0.4.6.2               gridExtra_2.3              
[29] MASS_7.3-54                 plot3D_1.4                  pracma_2.3.8                forcats_0.5.1              
[33] stringr_1.4.0               dplyr_1.0.7                 purrr_0.3.4                 readr_2.1.1                
[37] tidyr_1.1.4                 tibble_3.1.6                ggplot2_3.3.6               tidyverse_1.3.1 
```


