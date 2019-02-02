require('pROC')
require('psych')


filterGenesByVar <- function(dat, genes, cut.frac="1st Qu."){
  # filter genes based on their standard deviations
  #  removes genes with variance below the quartile cut
  genes <- intersect(genes, rownames(dat))
  sds <- apply(dat[genes, colnames(dat)], 1, sd, na.rm=TRUE)
  cut <- summary(sds)[cut.frac]
  genes2 <- genes[sds > cut]
  genes2 <- genes2[!is.na(genes2)]
  return(genes2)
}


trainSexLab <- function(train_dat, train_lab, female.genes=NULL, male.genes=NULL,
                        cut.frac="1st Qu."){
  # trains a model, saves it
  # read in the reference female, male genes if not provided - using ISEXs
  if (is.null(female.genes)){
    female.genes <- read.csv("data/list_f_genes_i.csv", stringsAsFactors = FALSE)
    female.genes <- female.genes[,1]
  } 
  if (is.null(male.genes)){
    male.genes <- read.csv("data/list_m_genes_i.csv", stringsAsFactors=FALSE)
    male.genes <- male.genes[,1]
  } 
  
  # filter probes based on variance --> remove low variance probes
  if (cut.frac == "Min."){
    m.genes <- male.genes
    f.genes <- female.genes
  } else {
    m.genes <- filterGenesByVar(train_dat, male.genes, cut.frac)
    f.genes <- filterGenesByVar(train_dat, female.genes, cut.frac)    
  }
  
  
  # calculate scores for the training data
  preds <- geomMeanScore(train_dat,f.genes, m.genes)
  preds2 <- preds[!is.nan(preds)]
  labels <- train_lab[colnames(train_dat)][!is.nan(preds)]
  
  # compute the threshold 
  roc_plot_train <- roc(labels, preds2, plot=FALSE) 
  threshold_gpl <- coords(roc_plot_train, "best", best.method="closest.topleft", ret=c("threshold"))
  
  # fit object
  return(list("f"=f.genes, "m"=m.genes, "threshold"=threshold_gpl))
}

predSexLab <- function(fit, test_dat, numeric_lab=FALSE){
  f.genes <- intersect(fit$f, rownames(test_dat))
  m.genes <- intersect(fit$m, rownames(test_dat))
  threshold_gpl <- fit$threshold
  preds_test <- geomMeanScore(test_dat, f.genes, m.genes)
  if (numeric_lab){
    sex_lab <- ifelse(preds_test > threshold_gpl, 1, 0)
    
  } else{
    sex_lab <- ifelse(preds_test > threshold_gpl, "male", "female")
  }
  return(sex_lab)
}

expDataToRanks <- function(gse.obj, list.genes=NULL){
  # convert the expression data in a GSE object to ranks 
  # uses the list of genes as a the full set of genes
  if (is.null(list.genes)){
    list.genes <- read.csv("data/list_all_genes.csv")[,1]
  }
  exp.genes <- convertToGenes(gse.obj, list.genes)
  rank.dat <- apply(exp.genes, 2, rank, na.last="keep")
  return(rank.dat)
} 


geomMeanScore <- function(dat, female.probes, male.probes){
  # calculate the score
  # geometric mean of male probes - female probes
  mean.m <- sapply(1:ncol(dat), function(x){
    m.dat <- unlist(dat[male.probes,x])
    return(geometric.mean(m.dat, na.rm=TRUE)    )
  })
  
  mean.f <- sapply(1:ncol(dat), function(x){
    f.dat <- unlist(dat[female.probes,x])
    return(geometric.mean(f.dat, na.rm=TRUE)    )
  })
  mean.diff <- mean.m-mean.f
  return(mean.diff)
}