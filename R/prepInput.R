
#' Prepare the input data for sex labeling.
#'
#' Briefly, takes an expression matrix with rows
#'
#' @param probe_mat an expression matrix with probes as rows and columns as samples
#' @param probe_map list mapping from probes to genes, names are probes, values are genes
#' @param gene_list list of all genes to extract, if not provided will use default list
#' @return rank_dat ranked dataset with rows as genes
prepInput <- function(probe_mat, probe_map, gene_list=list_genes){
  .checkExprMatFormat(probe_mat)
  .checkProbeMapping(probe_mat, probe_map, gene_list)

  expr_mat <- convertToGenes(probe_mat, probe_map, gene_list)
  rank_dat <- .expDataToRanks(expr_mat)
}

#' Convert an expression dataset to ranks.
#'
#' Briefly, takes an expression matrix, and then ranks each column, 1 to n (number of genes)
#' in order of increasing expression values. Missing data is ranked last.
#'
#' @param expr_mat an expression matrix with genes as rows and columns as samples
#' @param gene_list list of all genes to extract, if not provided will use default list
#' @return rank_dat ranked dataset
.expDataToRanks <- function(expr_mat, gene_list = list_genes) {
  rank_dat <- apply(expr_mat, 2, rank, na.last = "keep")
  return(rank_dat)
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
#' @param gene_list list of genes for the rows
#' @return expr_mat - an expression matrix with probes as genes
convertToGenes <- function(probe_mat, probe_map, gene_list){
  # TODO
  #  optimize, this is slow!! <-- possibly switch to MetaIntegrator function

  expData <- probe_mat
  keys <- probe_map
  list.keys <- keys[keys %in% gene_list]
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
  missing.genes <- setdiff(gene_list, list.keys)
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)

  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expr_mat <- expDataPlusMiss[gene_list,] # REORDER so it matches other data

  return(expr_mat)
}
