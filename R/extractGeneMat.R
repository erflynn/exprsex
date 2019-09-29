
#'


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



#' Function to parse genes from fData
#' Can be used directly if we have output from GEOquery
#'
#' @param fData feature/probe data table
#' @return list of genes with names as probes
.parseKeysFromFData <- function(fData){
  require("MetaIntegrator")
  gene_parse <- MetaIntegrator:::.GEO_fData_key_parser(fData)
  gene <- gene_parse$keys
  gene <- gene[!is.na(gene)]
  if (length(gene) < 100){
    print("re-extracting probe mapping with alternate parsing")
    gene <- .parseGenesAlt(fData)
    gene <- gene[!is.na(gene)]
  }
  return(gene)
}

#' Get the list of genes for a GPL
#'
#' @param gpl.name
#' @return list of genes with names as probes
.getGenes <- function(gpl.name){
  require('GEOquery')
  print(sprintf("Extracting genes from %s", gpl.name))
  tryCatch({
    gpl.obj <- getGEO(gpl.name, AnnotGPL=TRUE)
    fData <- gpl.obj@dataTable@table
    return(.parseKeysFromFData(fData))
  }, error = function(err){
    print(sprintf("error loading %s", gpl.name))
    return(list())
  })
}


#' Reformat keys into a data frame
#'
#' @param keys list with probe to gene mapping
#' @return gene dataframe with columns gene and probes
.getGeneDf <- function(keys){
  require('dplyr')
  require('tidyr')
  require('stringr')
  if (length(keys)==0){
    print("error - no keys available")
    return(data.frame("gene"=character()))
  }
  key.df <- data.frame(keys, names(keys))
  colnames(key.df) <- c("gene", "probes")
  key.df2 <- key.df %>%
    mutate(gene = as.character(gene),
           probes = as.character(probes)) %>%
    separate_rows(gene, sep=",")  %>%
    mutate(gene=str_trim(gene))
  return(key.df2)
}


#' Take a list with probe to gene mapping and formate it as gene to probe mapping
#'  In the process, separate out to allow for multi-mapping
#'
#' @param keys list with probe to gene mapping
#' @return list with gene to probe mapping, genes are names, probes
.formatGeneProbe <- function(keys){
  key.df <- .getGeneDf(keys)
  gene.to.probe <- split(key.df$probes,  key.df$gene)
  return(gene.to.probe)
}



#' Grab the gene to probe mapping for a given object
#'
#' @param platform the platform to download
#' @param gse_keys the keys from the gse_obj
#' @param platform_dir option to set platform  (defaults to tempdir)
#' @return gene.to.probe a mapping from genes to probes
.getGeneToProbe <- function(platform, gse_keys=NULL, platform_dir=NULL){
  require('miceadds')
  if (is.null(platform_dir)){
    platform_dir <- tempdir()
  }
  platform.path <- sprintf("%s/%s.RData", platform_dir, platform)
  if (file.exists(platform.path)){
    load.Rdata(platform.path, "gene.to.probe")
  }  else {

    # extract the keys and then convert it to a table
    keys <- gse_keys
    if (is.null(gse_keys) | (length(keys[!is.na(keys)]) < 100)){
      keys <- .getGenes(platform)
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
#' @param gse.obj
#' @param gene_list list of genes for the rows
#' @return expr_mat - an expression matrix with probes as genes
.convertToGenes <- function(gse.obj, gene_list){
  expData <- gse.obj$expr
  gene.to.probe <- .getGeneToProbe(gse.obj$platform, gse.obj$keys)

  # TODO - if gene_list is null?
  if (!is.null(gene_list)){
    # filter for the genes in the gene_list
    gene.to.probe <- gene.to.probe[gene_list]
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

