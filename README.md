# exprsex
## E Flynn and Annie Chang
R package for sex labeling expression data. Still beta version.

To install, in R
 1. Install three Bioconductor packages for data mapping:
   `BiocManager::install(c("org.Hs.eg.db", "org.Rn.eg.db", "org.Mm.eg.db"))`
 2. Make sure you have devtools installed and loaded, then run: `devtools::install_github(erflynn/exprsex)`. 
  This will also prompt you to install the required dependencies.


To run:
 1. Download an expression dataset using `getPrepGSE(gse)`
 2. Prepare the input by running:
      `reorderRank(expr)`
 3. Then, either:
     * Train a sex labeler: 
         `fit <- trainSexLab(expr, sex_lab)`
     * Or run sex labeling:
         `predSexLab(fit, expr)`

Sex labeling genes are from ISEXs (Bongen et al., in preparation). 


