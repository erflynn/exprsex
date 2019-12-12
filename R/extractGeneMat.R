
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


#' Get the list of genes for a GPL
#'
#' @param gpl.name the name of the platform
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return probe/gene data frame
.getGenes <- function(gpl.name, ref_dir=NULL){
  print(sprintf("Extracting genes from %s", gpl.name))
  tryCatch({
    return(parse_entrez_from_gpl(gpl.name, ref_dir))
  }, error = function(err){
    print(sprintf("error loading %s", gpl.name))
    return(list())
  })
}




#' Take a list with probe to gene mapping and formate it as gene to probe mapping
#'  In the process, separate out to allow for multi-mapping
#'
#' @param keys a list probe to gene
#' @return list with gene to probe mapping, genes are names, probes
.formatGeneProbe <- function(keys){
  keys <- keys[!is.null(keys)]
  key.df <- data.frame(cbind("gene"=keys, "probe"=names(keys)))
  key.df <- key.df[key.df$gene!="NULL",]
  gene.to.probe <- key.df$probe # // TODO may cause errors for multi-map
  names(gene.to.probe) <- key.df$gene
  return(gene.to.probe)
}



#' Grab the gene to probe mapping for a given object
#'
#' @param platform the platform to download
#' @param gse_keys the keys from the gse_obj
#' @param platform_dir option to set platform dir  (defaults to tempdir)
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return gene.to.probe a mapping from genes to probes
.getGeneToProbe <- function(platform, gse_keys=NULL, platform_dir=NULL, ref_dir=NULL){
  if (is.null(platform_dir)){
    platform_dir <- tempdir()
  }
  platform.path <- sprintf("%s/%s.RData", platform_dir, platform)
  if (file.exists(platform.path)){
    miceadds::load.Rdata(platform.path, "gene.to.probe")
  }  else {

    # extract the keys and then convert it to a table
    keys <- gse_keys
    if (is.null(gse_keys) | (length(keys[!is.na(gse_keys)]) < 100)){
      keys <- .getGenes(platform, ref_dir)
    }
    gene.to.probe <- .formatGeneProbe(keys)
    save(gene.to.probe, file=platform.path)
  }
  return(gene.to.probe)
}


#' Convert an expression matrix with probes as rows to genes.
#'
#' Briefly, objects downloaded by MetaIntegrator or GEOquery contain an expression matrix and a list of
#' key mapping from the matrix rows (probes) to genes. This function uses a list of genes
#' and takes the average of the values of all probes pointing to a particular gene.
#' If no probes map to that gene, the gene value is NA.
#'
#' @param gse.obj the gene object
#' @param gene_list list of genes for the rows
#' @return expr_mat - an expression matrix with probes as genes
.convertToGenes <- function(gse.obj, gene_list){
  expData <- gse.obj$expr
  gene.to.probe <- .getGeneToProbe(gse.obj$platform, gse.obj$keys)

  # filter to remove hugely multi-mapping??
  gene.to.probe <- gene.to.probe[gene.to.probe %in% rownames(expData)]
  gene.to.probe2 <- gene.to.probe[!is.na(names(gene.to.probe))]

  # TODO - if gene_list is null?
  if (!is.null(gene_list)){
    # filter for the genes in the gene_list
    gene.to.probe <- gene.to.probe[intersect(gene_list, names(gene.to.probe))]
  }

  expData2 <- do.call(cbind, lapply(1:length(gene.to.probe), function(x) {
    # get the gene and the probe
    g <- names(gene.to.probe)[x]
    p <- unlist(gene.to.probe[g])

    # take the average
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
  missing.genes <- setdiff(gene_list, names(gene.to.probe))
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)

  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expr_mat <- expDataPlusMiss[gene_list,] # REORDER so it matches other data
  return(expr_mat)
}
