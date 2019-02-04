


#' Checks that an expression matrix has the correct format
#'
#' @param expr_mat
#'
.checkExprMatFormat <- function(expr_mat){
}

#' Checks that expression matrix and sex labels have the correct structure
#'
#' @param expr_mat the expression matrix, rows should be gene names, columns are samples
#' @param sex_lab list of sex labels, names of list should match column names
.checkInputFormat <- function(expr_mat, sex_lab){
  .checkExprMatFormat(expr_mat)
  assertthat::are_equal(dim(expr_mat)[2], length(sex_lab))
  assertthat::are_equal(colnames(expr_mat), names(sex_lab))
  # TODO assert that sex_lab are 0/1 or male/female
}

#' Checks that a fit object has the correct structure
#'
#' @param fit
.checkFitObject <- function(fit){
  assertthat::has_name(fit, "f")
  assertthat::has_name(fit, "m")
  assertthat::has_name(fit, "threshold") #TODO check
  assertthat::is.number(fit$threshold)
  # TODO error messages
  # TODO check the format of the output
}


# --------- Main functions --------- #


#' Run sex labeling using default settings.
#'
#' @param expr_mat expression data matrix - rows are genes, columns are samples
#' @param platform the microarray platform (e.g. GPL570), if specified will use a model specific to that platform
#' @param is_rank TRUE if the data is already in rank form, otherwise will perform the conversion
#'
#' @return sex_lab - a list of sex labels (0/1 or "female"/"male")
runSexLab <- function(expr_mat, platform=NULL, is_rank=FALSE, numeric_lab=FALSE){
  # TODO
  #   assert that input is correct
  #   assert that the ranked data is actually ranked if it says it is

  # if it is not ranked, ranked the data
  if (!is_rank){

  }

  if (is.null(platform)){
    # use the standard fit
  } else {
    # read in the platform list

    # is it in the platform list? -- if not, throw and error

    # load the fit
  }

  # sex label
  sex_lab <- predSexLab(fit, test_dat, numeric_lab)
  return(sex_lab)
}


#' Train a sex labeling classifier using input training data.
#'
#' @param train_dat training data frame - rows are genes (HGNC_symbol), columns are samples, values are ranks
#' @param train_lab named list of training labels - 0 for female, 1 for male
#' @param female_genes a list of female-specific genes (HGNC_symbol), if not provided - default genes will be used
#' @param male_genes a list of male-specific genes (HGNC_symbol), if not provided, default genes will be used
#' @param cut_frac
#'
#' @return fit
trainSexLab <- function(train_dat, train_lab, female_genes = NULL, male_genes = NULL,
                        cut_frac = "1st Qu.") {
  # TODO:
  #  make sure input is correct
  # assert that rownames are HGNC_symbol
  # assert that column names == training labels
  # assert that label values are 1 or 0
  # assert that female_genes and male_genes are lists and their values are in the rownames
  # get rid of cut_frac
  require('pROC')
  # trains a model, saves it read in the reference female, male genes if not provided

  if (is.null(female_genes)) {
    female_genes <- sex_lab_genes$f
  }
  if (is.null(male_genes)) {
    male_genes <- sex_lab_genes$m
  }

  # filter probes based on variance --> remove low variance probes
  if (cut_frac == "Min.") {
    m_genes <- male_genes
    f_genes <- female_genes
  } else {
    m_genes <- .filterGenesByVar(train_dat, male_genes, cut_frac)
    f_genes <- .filterGenesByVar(train_dat, female_genes, cut_frac)
  }

  # calculate scores for the training data
  preds <- .geomMeanScore(train_dat, f_genes, m_genes)
  preds2 <- preds[!is.nan(preds)]
  labels <- train_lab[colnames(train_dat)][!is.nan(preds)]

  # compute the threshold
  roc_plot_train <- roc(labels, preds2, plot = FALSE)
  threshold_gpl <- coords(roc_plot_train, "best", best.method = "closest.topleft", ret = c("threshold"))

  # fit object
  return(list(f = f_genes, m = m_genes, threshold = threshold_gpl))
}

#' Use a particular sex labeling fit to label an input dataset
#' @param fit a sex labeling fit object
#' @param expr_mat data frame to label - rows are genes (HGNC_symbol), columns are samples, values are ranks
#' @param numeric_lab if TRUE output is numeric (0/1), otherwise it is ("female"/"male")
#'
#' @return sex_lab - a list of sex labels (0/1 or "female"/"male")
predSexLab <- function(fit, expr_mat, numeric_lab = FALSE) {
  # TODO
  # assert correct formatting
  require('pROC')
  # find the
  f_genes <- intersect(fit$f, rownames(expr_mat))
  m_genes <- intersect(fit$m, rownames(expr_mat))
  # TODO assert that there are some genes

  threshold_gpl <- fit$threshold
  preds_test <- .geomMeanScore(expr_mat, f_genes, m_genes) # calculate the score
  if (numeric_lab) {
    sex_lab <- ifelse(preds_test > threshold_gpl, 1, 0)
  } else {
    sex_lab <- ifelse(preds_test > threshold_gpl, "male", "female")
  }
  # TODO - add names to the sex lab if not already there
  #   output should be a named list
  return(sex_lab)
}

# --------- Helper functions --------- #

#' Helper function to filter genes based on their standard deviations.
#'
#' Removes sex labeling genes with variance below the provided quartile cut. This can
#' be helpful for removing genes on a platform that have low variance,
#' and such are not informative for sex labeling.
#'
#' @param expr_mat expr_mat - matrix of expression genes
#' @param genes list of sex labeling genes
#' @param cut_frac the fraction
#'
#' @return a filtered set of sex labeling genes
.filterGenesByVar <- function(expr_mat, genes, cut_frac = "1st Qu.") {
  # TODO
  #  genes is confusing - just
  genes <- intersect(genes, rownames(expr_mat))
  sds <- apply(expr_mat[genes, colnames(expr_mat)], 1, sd, na.rm = TRUE)
  cut <- summary(sds)[cut_frac]
  genes2 <- genes[sds > cut]
  genes2 <- genes2[!is.na(genes2)]
  return(genes2)
}

#' Convert an expression dataset in probe form to ranks.
#'
#' Briefly, takes an expression matrix (rows are probes) and mapping to genes,
#' and converts the matrix to genes, and then ranks each column, 1 to n (number of genes)
#' in order of decreasing expression data. Missing data is ranked last.
#'
#' @param probe_mat an expression matrix with probes as rows and columns as samples
#' @param probe_map list mapping from probes to genes, names are probes, values are genes
#' @param list_genes list of all genes to extract, if not provided will use default list
#' @return rank_dat ranked dataset with rows as genes
expDataToRanks <- function(probe_mat, probe_map, gene_list = list_genes) {

  # TODO
  #   should work if already genes
  #   decreasing or increasing?

  expr_mat <- convertToGenes(probe_mat, probe_map, gene_list)
  rank_dat <- apply(expr_mat, 2, rank, na.last = "keep")
  return(rank_dat)
}

#' Calculate the score geometric mean of male genes - female genes
#'
#' Relies on the psych package to calculate geometric means for each set of genes.
#' The end result is a score for the difference in means which is used for classification.
#'
#' @param expr_mat expression dataset: rows are genes, columns are samples
#' @param female_genes list of genes with female-specific expression
#' @param male_genes list of genes with male-specific expression
#'
#' @return the difference in means
.geomMeanScore <- function(expr_mat, female_genes, male_genes) {
  # TODO
  #  change so that the nomenclature for dat is always the same
  require('psych')
  mean_m <- sapply(1:ncol(expr_mat), function(x) {
    m_dat <- unlist(expr_mat[male_genes, x])
    return(geometric.mean(m_dat, na.rm = TRUE))
  })

  mean_f <- sapply(1:ncol(expr_mat), function(x) {
    f_dat <- unlist(expr_mat[female_genes, x])
    return(geometric.mean(f_dat, na.rm = TRUE))
  })
  mean_diff <- mean_m - mean_f
  return(mean_diff)
}

#' Convert an expression matrix with probes as rows to genes.
#'
#' Briefly, objects downloaded by MetaIntegrator or GEOquery contain an expression matrix and a list of
#' key mapping from the matrix rows (probes) to genes. This function uses a list of genes
#' and takes the average of the values of all probes pointing to a particular gene.
#' If no probes map to that gene, the gene value is NA.
#'
#' @param probe_mat an expression matrix with probes as rows and columns as samples
#' @param probe_map list mapping from probes to genes, names are probes, values are genes
#' @param list_genes list of genes for the rows
#' @return expr_mat - an expression matrix with probes as genes
convertToGenes <- function(probe_mat, probe_map, list_genes){
  # TODO
  #  optimize, this is slow!! <-- possibly switch to MetaIntegrator function
  #  make sure the object contains expression data, keys

  expData <- probe_mat
  keys <- probe_map
  list.keys <- keys[keys %in% list_genes]
  key.df <- data.frame(list.keys, names(list.keys))
  colnames(key.df) <- c("gene", "probes")
  key.df$gene <- as.character(key.df$gene)
  key.df$probes <- as.character(key.df$probes)

  gene.to.probe <- split(key.df$probes,  key.df$gene) # this is slow... mb store for each platform
  expData2 <- do.call(cbind, lapply(1:length(gene.to.probe), function(x) {
    # get the gene and the probe
    g <- names(gene.to.probe)[x]
    p <- unlist(gene.to.probe[g])
    if (length(p)>1){
      expD <- expData[p,]
      df <- (data.frame(colMeans(expD, na.rm=TRUE)))
      return(df)
    }
    else {
      df <- data.frame(expData[p,])
      return(df)
    }})) ### ALSO SLOW...

  colnames(expData2) <- names(gene.to.probe)
  expData2.2 <- data.frame(t(expData2)) # columns are samples, rows are genes

  # create a data fram of NAs for missing genes
  missing.genes <- setdiff(list_genes, list.keys)
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)

  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expr_mat <- expDataPlusMiss[list_genes,] # REORDER so it matches other data

  return(expr_mat)
}
