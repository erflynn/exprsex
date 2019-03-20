#' Run sex labeling using default settings.
#'
#' @param expr_mat expression data matrix - rows are genes, columns are samples
#' @param platform the microarray platform (e.g. GPL570), if specified will use a model specific to that platform
#' @param is_rank TRUE if the data is already in rank form, otherwise will perform the conversion
#'
#' @return sex_lab - a list of sex labels (0/1 or "female"/"male")
runSexLab <- function(expr_mat, platform=NULL, is_rank=FALSE, numeric_lab=TRUE){
  # TODO
  #   assert that input is correct
  #   assert that the ranked data is actually ranked if it says it is
  error("This function is not implemented yet, please use predSexLab instead")
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
  require('pROC')

  # TODO:
  # get rid of cut_frac

  .checkTrainInputFormat(train_dat, train_lab)

  # loads genes if not specified
  if (is.null(female_genes)) {
    female_genes <- sex_lab_genes$f
  }
  if (is.null(male_genes)) {
    male_genes <- sex_lab_genes$m
  }

  # make sure the genes are in the rownames
  .checkGeneFormat(train_dat, female_genes)
  .checkGeneFormat(train_dat, female_genes)

  # remove low variance probes
  m_genes <- .filterGenesByVar(train_dat, male_genes, cut_frac)
  f_genes <- .filterGenesByVar(train_dat, female_genes, cut_frac)

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
#' @param numeric_lab if FALSE, the output is "female/male", otherwise it is numeric (0/1)
#'
#' @return sex_lab - a list of sex labels (0/1 or "female"/"male")
predSexLab <- function(fit, expr_mat, numeric_lab = TRUE) {
  require('pROC')

  # check input
  .checkFitFormat(fit)
  .checkExprMatFormat(expr_mat)

  # find the fit format
  f_genes <- intersect(fit$f, rownames(expr_mat))
  m_genes <- intersect(fit$m, rownames(expr_mat))

  # assert that there are some sex labeling genes
  .my_assert("The expresion matrix does not contain any of the f genes", !is.null(f_genes))
  .my_assert("The expresion matrix does not contain any of the m genes", !is.null(m_genes))

  threshold_gpl <- fit$threshold

  # caclulate the score
  preds_test <- .geomMeanScore(expr_mat, f_genes, m_genes) # calculate the score
  if (numeric_lab) {
    sex_lab <- ifelse(preds_test > threshold_gpl, 1, 0)
  } else {
    sex_lab <- ifelse(preds_test > threshold_gpl, "male", "female")
  }

  # Add names to the list of sex labels
  names(sex_lab) <- colnames(expr_mat)
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

  # if there is no cutoff, just return the list
  if (cut_frac == "Min.") {
    return(genes)
  }

  # otherwise, find the subset of genes, and filter
  genes <- intersect(genes, rownames(expr_mat))
  sds <- apply(expr_mat[genes, colnames(expr_mat)], 1, sd, na.rm = TRUE)
  cut <- summary(sds)[cut_frac]
  genes2 <- genes[sds > cut]
  genes2 <- genes2[!is.na(genes2)]
  return(genes2)
}



#' Helper function to calculate geometric mean of a set of genes across all samples in a matrix.
#'
#' Computes the geometric mean of the expression of a particular set of genes for every sample
#' (column) in a matrix. NAs are ignored during this calculation. The result is a vector of
#' geometric means, each corresponding to one sample (column) in the original matrix.
#'
#' @param expr_mat expression dataset, rows are genes, columns are samples
#' @param genes list of genes to calculate the geometric mean of
#'
#' @return a vector of geometric means
.geomMeanAcrossGenes <- function(expr_mat, genes){
  require('psych')
  sapply(1:ncol(expr_mat), function(x) {
    dat <- unlist(expr_mat[genes, x])
    return(geometric.mean(dat, na.rm = TRUE))
  })
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
#' @return a vector of mean difference scores
.geomMeanScore <- function(expr_mat, female_genes, male_genes) {

  mean_m <- .geomMeanAcrossGenes(expr_mat, male_genes)
  mean_f <- .geomMeanAcrossGenes(expr_mat, female_genes)

  mean_diff <- mean_m - mean_f
  return(mean_diff)
}

