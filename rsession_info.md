R Session Info
================

All analyses are implemented in **R (4.4.1)** with a R-studio docker
container available
[here](https://hub.docker.com/r/bakeronit/rstudio_hpc_r4.4).
Alternatively, you can reproduce the R environment with the
[renv.lock](renv.lock) file:

``` r
install.packages("renv")
renv::restore()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggrepel_0.9.6                     ggalluvial_0.12.5                
    ##  [3] ggparty_1.0.0                     partykit_1.2-22                  
    ##  [5] mvtnorm_1.2-5                     libcoin_1.0-10                   
    ##  [7] patchwork_1.2.0                   colorspace_2.1-1                 
    ##  [9] ggmosaic_0.3.3                    waffle_1.0.2                     
    ## [11] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.72.0                  
    ## [13] rtracklayer_1.64.0                BiocIO_1.14.0                    
    ## [15] Biostrings_2.72.1                 XVector_0.44.0                   
    ## [17] GenomicRanges_1.56.1              GenomeInfoDb_1.40.1              
    ## [19] IRanges_2.38.1                    S4Vectors_0.42.1                 
    ## [21] BiocGenerics_0.50.0               CHORD_2.03                       
    ## [23] ggh4x_0.2.8                       signature.tools.lib_2.4.4        
    ## [25] googlesheets4_1.1.1               RColorBrewer_1.1-3               
    ## [27] gt_0.11.1                         vcfR_1.15.0                      
    ## [29] data.table_1.15.4                 lubridate_1.9.3                  
    ## [31] forcats_1.0.0                     stringr_1.5.1                    
    ## [33] dplyr_1.1.4                       purrr_1.0.2                      
    ## [35] readr_2.1.5                       tidyr_1.3.1                      
    ## [37] tibble_3.2.1                      ggplot2_3.5.1                    
    ## [39] tidyverse_2.0.0                   yaml_2.3.10                      
    ## [41] here_1.0.1                       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.16.0           jsonlite_1.8.8             
    ##  [3] magrittr_2.0.3              rmarkdown_2.27             
    ##  [5] fs_1.6.4                    zlibbioc_1.50.0            
    ##  [7] vctrs_0.6.5                 Rsamtools_2.20.0           
    ##  [9] RCurl_1.98-1.16             htmltools_0.5.8.1          
    ## [11] S4Arrays_1.4.1              curl_5.2.1                 
    ## [13] cellranger_1.1.0            Formula_1.2-5              
    ## [15] SparseArray_1.4.8           htmlwidgets_1.6.4          
    ## [17] plyr_1.8.9                  plotly_4.10.4              
    ## [19] GenomicAlignments_1.40.0    lifecycle_1.0.4            
    ## [21] iterators_1.0.14            pkgconfig_2.0.3            
    ## [23] Matrix_1.7-0                R6_2.5.1                   
    ## [25] fastmap_1.2.0               GenomeInfoDbData_1.2.12    
    ## [27] MatrixGenerics_1.16.0       digest_0.6.36              
    ## [29] rprojroot_2.0.4             vegan_2.6-6.1              
    ## [31] fansi_1.0.6                 timechange_0.3.0           
    ## [33] httr_1.4.7                  abind_1.4-5                
    ## [35] mgcv_1.9-1                  compiler_4.4.1             
    ## [37] gargle_1.5.2                withr_3.0.1                
    ## [39] backports_1.5.0             BiocParallel_1.38.0        
    ## [41] Rttf2pt1_1.3.12             MASS_7.3-60.2              
    ## [43] DelayedArray_0.30.1         rjson_0.2.21               
    ## [45] permute_0.9-7               tools_4.4.1                
    ## [47] googledrive_2.1.1           ape_5.8                    
    ## [49] extrafontdb_1.0             glue_1.8.0                 
    ## [51] restfulr_0.0.15             inum_1.0-5                 
    ## [53] nlme_3.1-164                checkmate_2.3.2            
    ## [55] cluster_2.1.6               generics_0.1.3             
    ## [57] gtable_0.3.5                tzdb_0.4.0                 
    ## [59] pinfsc50_1.3.0              hms_1.1.3                  
    ## [61] xml2_1.3.6                  utf8_1.2.4                 
    ## [63] foreach_1.5.2               pillar_1.9.0               
    ## [65] splines_4.4.1               lattice_0.22-6             
    ## [67] survival_3.6-4              tidyselect_1.2.1           
    ## [69] knitr_1.47                  gridExtra_2.3              
    ## [71] SummarizedExperiment_1.34.0 xfun_0.44                  
    ## [73] Biobase_2.64.0              matrixStats_1.3.0          
    ## [75] DT_0.33                     stringi_1.8.4              
    ## [77] UCSC.utils_1.0.0            lazyeval_0.2.2             
    ## [79] evaluate_0.24.0             codetools_0.2-20           
    ## [81] extrafont_0.19              cli_3.6.3                  
    ## [83] rpart_4.1.23                munsell_0.5.1              
    ## [85] Rcpp_1.0.13                 XML_3.99-0.17              
    ## [87] parallel_4.4.1              bitops_1.0-8               
    ## [89] viridisLite_0.4.2           scales_1.3.0               
    ## [91] crayon_1.5.3                rlang_1.1.4
