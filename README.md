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
         the default fit can be accessed with `data(default_fit)`

For a full example, see `vignettes/sex-label-vignette.Rmd`.

We also provide code to map a wide range of GEO platforms to entrez IDs in `gpl_map()` functions.
This is done as part of the `getPrepGSE` function. To map an individual platform:
 `ref_tab <- parse_entrez_from_gpl(platform_id)`

