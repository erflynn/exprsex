# exprsex
## E Flynn
R package for sex labeling expression data. Still beta version.

To run:
 1. Download an expression dataset
 2. Extract expression matrix and probe to gene mappings ("keys")
 3. Prepare the input by running:
      `prepInput(expr, keys)`
 4. Then, either:
     * Train a sex labeler: 
         `fit <- trainSexLab(expr, sex_lab)`
     * Or run sex labeling:
         `predSexLab(fit, expr)`

Sex labeling genes are from ISEXs (Bongen et al., in preparation). 


