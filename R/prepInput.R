




#' Prepare the input data for sex labeling by converting it to genes
#' and getting the ranks.
#'
#' @param gse.obj the object
#' @param gene_list the list of genes to use
#' @return the gse.obj prepared for sex labeling
prepInput <- function(gse.obj, gene_list=list_genes){
  .checkExprMatFormat(gse.obj$expr)
  .checkProbeMapping(gse.obj$expr, gse.obj$keys, gene_list)

  if (!exists(gse.obj$exp_mat)){
    gse.obj$exp_mat <- .convertToGenes(gse.obj, gene_list)
  }
  if (!exists(gse.obj$rank)){
    gse.obj$rank <- .expDataToRanks(gse.obj$exp_mat)
  }
  return(gse.obj)
}

#' Take a GEOQuery object and convert it to a MetaIntegrator gse.obj
#'
#' @param geoQueryObj that is the output of
#' @return gse.obj with expression data and keys extracted
createGSEObj <- function(geoQueryObj){
 gse.obj <- list()
 gse.obj$expr <- exprs(geoQueryObj)
 gse.obj$keys <- .parseKeysFromFData(fData(geoQueryObj))
 return(gse.obj)
}


#' #' Prepare the input data for sex labeling.
#' #'
#' #' Alternate method for
#' #'
#' #' @param probe_mat an expression matrix with probes as rows and columns as samples
#' #' @param probe_map list mapping from probes to genes, names are probes, values are genes
#' #' @param gene_list list of all genes to extract, if not provided will use default list
#' #' @return rank_dat ranked dataset with rows as genes
#' prepInputProbe <- function(probe_mat, probe_map, gene_list=list_genes){
#'   .checkExprMatFormat(probe_mat)
#'   .checkProbeMapping(probe_mat, probe_map, gene_list)
#'
#'   expr_mat <- convertToGenes(probe_mat, probe_map, gene_list)
#'   rank_dat <- .expDataToRanks(expr_mat)
#' }

