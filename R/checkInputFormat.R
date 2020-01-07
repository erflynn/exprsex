
#' Custom assertion function, prints formatted error messages when they fail.
#'
#' Prints an error message if the assertion fails, otherwise does nothing.
#' Error message is of the format (and does not include the function call):
#'    "Error: <error_msg>"
#'
#' @param error_msg error message to print on failure
#' @param test the value to test
.my_assert <- function(error_msg, test){
  if (!test){
    stop(error_msg, call.=FALSE)
  }
}


#' Checks that an input expression matrix has the correct format
#'
#' Makes sure the object is a matrix, and the values in it are numeric.
#'
#' @param expr_mat the expression matrix
.checkExprMatFormat <- function(expr_mat){
  .my_assert("the expr_mat must be a matrix", class(expr_mat)=="matrix")
  # print an error if it is a data frame
  .my_assert("the expr_mat must contain only numeric values", is.numeric(expr_mat))
  #   update to work for a data frame, which is.numeric will not evaluate to TRUE
}



#' Checks that input sex labels have the correct format
#'
#' Make sure the sex labels are all 0/1 (NAs ok, but not all)
#'
#' @param sex_lab list of sex labels formatted as 0s/1s
.checkSexLabFormat <- function(sex_lab){
  .my_assert("sex labels should be formatted as 0 (f), 1 (m).", all(sex_lab %in% c(0, 1, NA)))
  .my_assert("sex label vector is all NAs.", !all(is.na(sex_lab)))
  .my_assert("the names of the sex label vector should correspond to the samples.", !is.null(names(sex_lab)))
}

#' Checks that sex labels and expression matrix match
#'
#' @param sex_lab list of sex labels formatted as 0s/1s
.checkMatching <- function(expr_mat, sex_lab){
  .my_assert("the columns of the expr_mat do not match the sex label names",
             all(colnames(expr_mat)==names(sex_lab)))
}


#' Checks that expression matrix and sex labels have the correct structure for training
#'
#' @param expr_mat the expression matrix, rows should be gene names, columns are samples
#' @param sex_lab list of sex labels, names of list should match column names
.checkTrainInputFormat <- function(expr_mat, sex_lab){
  .checkExprMatFormat(expr_mat)
  .checkSexLabFormat(sex_lab)
  .checkMatching(expr_mat, sex_lab)
}

#' Checks that the expression matrix rownames contain the genes of interest.
#'
#' @param mat the expression matrix, rows should be gene names
#' @param genes list of genes
.checkGeneFormat <- function(mat, genes){
  .my_assert("the expr_mat rows do not contain the genes",
         length(intersect(genes, rownames(mat))) > 0)
}


#' Checks that a fit object has the correct structure
#'
#' @param fit a sex labeling fit with attributes "f" (list of f genes), "m" (list of m genes), and "threshold" (cutoff score)
.checkFitFormat <- function(fit){
  .my_assert("the fit object is missing f genes", assertthat::has_name(fit, "f"))
  .my_assert("the fit object is missing m genes", assertthat::has_name(fit, "m"))
  .my_assert("the fit object is missing a threshold", assertthat::has_name(fit, "threshold"))
  .my_assert("the threshold must be numeric", assertthat::is.number(fit$threshold))
  .my_assert("the fit object is missing f genes", length(fit$f > 0) & !is.null(fit$f))
  .my_assert("the fit object is missing m genes", length(fit$m > 0) & !is.null(fit$m))
}

#' Checks that a probe mapping has the correct format.
#'
#' Briefly, the probe mapping should be formatted as a list with names as probes and values as genes.
#' The names of the probe map should match the rownames of the matrix and the values should match
#' the gene list.
#' This function looks at the intersections, and throws an error if there is no intersection
#' and prints a warning if the intersection is small (10% of either).
#'
#' @param probe_mat an expression matrix with probes as rows and columns as samples
#' @param probe_map list mapping from probes to genes, names are probes, values are genes
#' @param gene_list list of all genes to extract
.checkProbeMapping <- function(probe_mat, probe_map, gene_list){

  # // TODO - update this!

  # throw an error if there is no intersection
  intersect_w_mat <- length(intersect(rownames(probe_mat), names(probe_map)))
  intersect_w_gene_list <- length(intersect(probe_map, gene_list))
  .my_assert("The probe map names do not match the matrix rows", intersect_w_mat > 0)
  .my_assert("The probe map values do not match the gene list", intersect_w_gene_list > 0)

  # throw a warning
  if (intersect_w_mat < 0.1*nrow(probe_mat) | intersect_w_gene_list < 0.1*length(gene_list)){
    warning("There is little overlap between probe mapping and matrix rows or gene list", call.=FALSE)
  }
}
