# Paedomorphosis at the origin of marsupial mammals

__Authors:__
[Heather E White](mailto:heather.white.17@ucl.ac.uk), 
Anjali Goswami

__To cite the paper__: 

>White HW, Tucker AS, Fernandez V, Portela Miguez R, Hautier L, Herrel A, Urban DJ, Sears KE, Goswami A. Paedomorphosis at the origin of marsupial mammals. Current Biology. 2023.

Available at: https://github.com/HeatherEWhite/mammal_cranial_development

If using any of this code or data please cite the paper above and this repo

__To cite this repo__: 

> White HW, Goswami A. Paedomorphosis at the origin of marsupial mammals. Current Biology. 2023. Github: https://github.com/HeatherEWhite/mammal_cranial_development plus the Zenodo DOI: https://doi.org/10.5281/zenodo.7850303

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7850303.svg)](https://doi.org/10.5281/zenodo.7850303)



## Data

In this folder you will find data stored in the `Data` folder used for running the scripts stored in the main folder.

All .csv documents have single tabs 

The data used here are provided in the `Data` folder
1. `LMs` folder contains the raw .pts files for all 165 specimens used within this analysis
2. `plys_ascii` folder contains surface meshes for two scans used to visualise landmark placement. All other surface meshes are available on MorphoSource.
3. `Centroid_size` - centroid size as calculated from Procrusted landmarks for each specimen (.csv)
4. `LHS_absent_Procrusted_skull_LMs_slid` - LHS landmarks only, following mirroring, Procrustes superimposition, deletion of RHS, and sliding (.Rdata)
5. `LHS_only_Procrusted_LMs_mirrored` - LHS landmarks only, following mirroring, Procrustes superimposition, and deletion of RHS (.Rdata)
6. `LM_list` - cranial landmark numbers and cranial bone associations (.csv)
7. `Mirrored_LMs` - array containing all LHS mirrored to RHS landmarks, prior to Procrustes superimposition (.Rdata)
8. `PCAresults_all` - results of PCA for all specimens, including PC axis variation (.Rdata)
9. `PC_scores_all_specimens` - Principal component scores for each PC for each specimen (.Rdata)
10. `Precocial_altricial_proportion_of_sig_diffs_allometry` - species' altriciality and proportion of allometric trajectory significant differences (.csv)
11. `Specimen_info` - full dataset specimen details, including: continuous and discrete age, altriciality, species, centroid size (.csv)
12. `Specimen_info_adults` - details for the adult only specimens (n=22) as above (.csv)
13. `Specimen_names` - a list containing all 165 specimen names (.csv)
14. `Upham_NDexp_MCC` - mammal phylogeny as published in Upham et al. (2019) (.tre.txt)
15. `absent_LMs` - variably present bones for specimens across the dataset, listed as the first landmark defining the cranial bone (.csv)
16. `adult_Procrusted_LMs` - array containing landmarks subjected to Procrustes superimposition for adult specimens only (n=22) (.Rdata)
17. `mirrored_Procrusted_LMs` - array containing mirrored landmarks (LHS and RHS) subjected to Procrustes superimposition (.Rdata)
18. `my_mammal_tree` - phylogeny containing all 22 species within the dataset, trimmed from Upham et al. (2019) (.nexus)
19. `my_mammal_tree_placentals` - phylogeny containing all placental mammal species within the dataset (n=15), trimmed from Upham et al. (2019) (.nexus)
20. `tree_taxa_names` - taxa names matching those within the phylogeny for the full dataset (.csv)
21. `tree_taxa_names_placentals` - taxa names matching those within the placental mammal phylogeny (.csv)

## Analysis
In this repository you will find raw data (.csv, .Rdata, .nexus) and scripts for analyses (scripts supplied as .R files)

 :file_folder:
* **Data**

All data is outlined above in the 'Data' section. The `Data` folder includes .csv files, .Rdata files, .nexus files for analysis, as well as additional folders, including `LMs` containing all original .pts landmark files for the dataset and `plys_ascii` containing a couple of surface meshes for visualisation purposes only. All other surface meshes are available to download from MorphoSource.

 :file_folder:
* **Code_for_analyses**

`MANOVAs.R`

`PCA_3D.R` This script contains the relevant code to produce the 3D PCA which is captured in Figure 4 and can be visualised in 3D here.

`PCAs.R`

`PC_extreme_shapes.R`

`PlottingValues_Function.R`

`allometry.R`

`ancestral_states.R`

`mirroring_LMs.R`

`ontogenetic_trajectories.R`

`phylogeny.R`

`phylomorphospace.R`

`variably_present.R`



## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/HeatherEWhite/suture_complexity_metrics/blob/master/LICENSE) file for details

## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication. 

```{r}
─ Session info ──────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.1 (2019-07-05)
 os       macOS Catalina 10.15.7      
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_GB.UTF-8                 
 ctype    en_GB.UTF-8                 
 tz       Europe/London               
 date     2023-01-11                  

─ Packages ──────────────────────────────────────────────────────────────────────────────────────
 package     * version date       lib source        
 assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.6.0)
 cachem        1.0.5   2021-05-15 [1] CRAN (R 3.6.2)
 callr         3.7.0   2021-04-20 [1] CRAN (R 3.6.2)
 cli           3.0.1   2021-07-17 [1] CRAN (R 3.6.2)
 crayon        1.4.1   2021-02-08 [1] CRAN (R 3.6.1)
 desc          1.2.0   2018-05-01 [1] CRAN (R 3.6.0)
 devtools      2.3.2   2020-09-18 [1] CRAN (R 3.6.2)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 3.6.2)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 3.6.2)
 fs            1.5.0   2020-07-31 [1] CRAN (R 3.6.2)
 glue          1.4.2   2020-08-27 [1] CRAN (R 3.6.2)
 lifecycle     1.0.1   2021-09-24 [1] CRAN (R 3.6.2)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 3.6.1)
 memoise       2.0.0   2021-01-26 [1] CRAN (R 3.6.2)
 pkgbuild      1.2.0   2020-12-15 [1] CRAN (R 3.6.2)
 pkgload       1.1.0   2020-05-29 [1] CRAN (R 3.6.2)
 prettyunits   1.1.1   2020-01-24 [1] CRAN (R 3.6.0)
 processx      3.5.2   2021-04-30 [1] CRAN (R 3.6.2)
 ps            1.6.0   2021-02-28 [1] CRAN (R 3.6.2)
 purrr         0.3.4   2020-04-17 [1] CRAN (R 3.6.2)
 R6            2.5.1   2021-08-19 [1] CRAN (R 3.6.2)
 remotes       2.4.0   2021-06-02 [1] CRAN (R 3.6.2)
 rlang         1.0.6   2022-09-24 [1] CRAN (R 3.6.1)
 rprojroot     2.0.2   2020-11-15 [1] CRAN (R 3.6.2)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 3.6.2)
 sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.6.0)
 testthat      3.0.1   2020-12-17 [1] CRAN (R 3.6.2)
 usethis       2.0.0   2020-12-10 [1] CRAN (R 3.6.2)
 withr         2.4.2   2021-04-18 [1] CRAN (R 3.6.2)
```
