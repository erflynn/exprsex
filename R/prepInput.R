

#' Get and prepare a gse file for sex labeling by converting it
#' to genes. Optionally this will also reorder by a list of genes and
#' convert it to ranks.
#'
#' @param gse the ID of the GSE to download
#' @param to.ranks whether to convert to ranks, default is false
#' @param gse.dir directory that contains gse files, if empty, will download the GSE
#' @param ref.dir directory where GPL ref data is location, empty defaults to tempdir()
#' @param out.dir directory to write output, if empty will just be returned and not be written out
#' @param gene_list the list of genes to use, this is helpful if you want to compare across multiple studies
#'
#' @return gene expression matrix, with rows as genes and columns as samples
#'         values are gene expression levels unless converted to
getPrepGSE <- function(gse, to.ranks=FALSE, gse.dir=NULL,
                       ref.dir=NULL,
                       out.dir=NULL,
                       gene_list=NULL){

  if (!is.null(gse.dir)){
    series.mat.f <- sprintf("%s/%s_series_matrix.txt.gz", gse.dir, gse)
    if (file.exists(series.mat.f)){
      geo.obj <- MetaIntegrator::getGEOData(gse, filename = series.mat.f)
    } else {
      geo.obj <- MetaIntegrator::getGEOData(gse, destdir = gse.dir)
    }
  } else {
    geo.obj <- MetaIntegrator::getGEOData(gse)
  }

  # // TODO: check that the object downloaded

  # grab platform and platform mapping
  platform.id <- unique(geo.obj$originalData[[1]]$platform)
  if (!is.null(ref.dir)){
    gpl.f <- sprintf("%s/%s_map.RData", ref.dir, platform.id)
    if (!file.exists(gpl.f)){
      map_mult_gpl(c(platform.id), ref.dir)
    }
    load(gpl.f) # --> ref_tab
  } else {
    ref_tab <- parse_entrez_from_gpl(platform.id)
  }

  # remove missing data from ref tab
  ref_tab2 <- dplyr::filter(ref_tab,
                            !is.na(gene) &
                              !is.na(probe) &
                              gene != "" &
                              probe != "" &
                              gene != "NA" &
                              probe != "NA")


  # map to genes
  exp_mat <- geo.obj$originalData[[1]]$expr
  gene_mat <- .convertGenes(exp_mat, ref_tab2)

  # convert to ranks
  if (to.ranks){
    # // TODO
  }

  # write out for sex labeling
  if (!is.null(out.dir)){
    write.table(gene_mat,
                file=sprintf("%s/%s_gene.txt", out.dir, gse),
                row.names=TRUE, quote=FALSE)
  }
  geo.obj2 <- geo.obj$originalData[[1]] # // TODO - does this work if mult platforms?
  geo.obj2$expr_mat <- gene_mat
  return(geo.obj2)
}

#' Convert the expression matrix to genes.
#' This takes the average of probes to map them to a gene.
#'
#' @param expData the expression matrix with rows as probes and columns as conditions
#' @param probe_gene a mapping frame with a probe and gene column
#' @return an expression matrix with genes as rows
.convertGenes <- function(expData, probe_gene){
  .checkExprMatFormat(expData) # other input checks?

  gene.to.probe <- gene.to.probe <- split(sapply(probe_gene$probe, as.character),
                                          sapply(probe_gene$gene, as.character))
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
    }})) ### SLOW...

  colnames(expData2) <- names(gene.to.probe)
  expData2.2 <- data.frame(t(expData2)) # columns are samples, rows are genes

  # // TODO: reordering?
  return(expData2.2)
}


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

