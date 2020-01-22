

#' Get and prepare a gse file for sex labeling by converting it to genes.
#' Optionally this will also reorder by a list of genes and convert it to ranks.
#'
#' @param gse the ID of the GSE to download
#' @param gse.dir directory that contains gse files, if empty, will download the GSE
#' @param gpl.dir directory where GPL ref data is location, empty defaults to tempdir()
#' @param out.dir directory to write output, if empty will just be returned and not be written out
#'
#' @return list of gene expression matrices, with rows as genes and columns as samples
#'         values are gene expression levels
getPrepGSE <- function(gse, gse.dir=NULL,gpl.dir=NULL, out.dir=NULL){

  geo.res <- NULL
    # if there is a prespecified directory
    if (!is.null(gse.dir)){
      # look for files in the directory
      series.f <- list.files(gse.dir, pattern=sprintf("%s[-|_]+", gse))
      if (length(series.f) > 0){
        geo.res <- lapply(series.f, function(f.name){
          geo.plat <- GEOquery::getGEO(file=sprintf("%s/%s", gse.dir, f.name), getGPL=FALSE)
        })
        names(geo.res) <- series.f
      }
      # if they aren't there, then load the study but save it to the directory
      else {
        tryCatch({
          geo.res <- GEOquery::getGEO(gse, destdir=gse.dir, getGPL=FALSE)
        }, err = function(e){
          print(sprintf("%s error loading", gse))
          return(NULL)
        })
      }
    } else {
      tryCatch({
        geo.res <- GEOquery::getGEO(gse, getGPL=FALSE)
      }, err = function(e){
        print(sprintf("%s error loading", gse))
        return(NULL)
      })
    }


  if (is.null(geo.res)){
    print(sprintf("No data for %s", gse))
    return(NULL)
  }

  # go through the geo obj (multiple if there is more than one platform for a study)
  geo.obj.list <- lapply(geo.res, function(geo.plat) {
    geo.obj <- list("expr"=Biobase::exprs(geo.plat), "pheno"=Biobase::pData(geo.plat),
                    "platform"=unique(geo.plat$platform), "note"="")

    # check that the object downloaded
    if ((is.null(dim(geo.obj$expr))) | (nrow(geo.obj$expr) < 1)){
      print(sprintf("expression data is missing for %s", gse))
      return(NA)
    }

    # log-transform if needed
    # // TODO: we rely on a hidden MetaIntegrator function to do this, update
    if (! MetaIntegrator:::.GEM_log_check(geo.obj)){
      min_value <- min(geo.obj$expr, na.rm = T)

      # adjust up if the minimum value is less than zero
      if (min_value < 0) {
        geo.obj$expr <- geo.obj$expr - min_value + 1
        geo.obj$note <- sprintf("%s; min val adjusted", geo.obj$note)
      }

      # log-transform if needed
      geo.obj$expr <- log2(geo.obj$expr)
      geo.obj$note <- sprintf("%s; log2 transformed", geo.obj$note)
    }

    # grab platform and platform mapping
    platform.id <- sapply(unique(geo.obj$platform), as.character)
    if (!is.null(gpl.dir)){
      gpl.f <- sprintf("%s/%s_map.RData", gpl.dir, platform.id)
      if (!file.exists(gpl.f)){
        map_mult_gpl(c(platform.id), gpl.dir)
      }
      if (file.exists(gpl.f)){
        miceadds::load.Rdata(gpl.f, "ref_tab")
      } else {
        ref_tab <- NA
      }
    } else {
      ref_tab <- parse_entrez_from_gpl(platform.id, verbose=FALSE)
    }

    # return NA if the mapping does not exist
    if (is.null(dim(ref_tab))){
      print(sprintf("probe mapping data is missing for %s so we cannot map %s_%s",
                    platform.id, gse, platform.id))
      return(NA)
    }

    # map to genes
    exp_mat <- geo.obj$expr
    gene_mat <- .convertGenes(exp_mat, ref_tab)

    # write out for sex labeling
    # // TODO - should we include this?
    if (!is.null(out.dir)){
      write.table(gene_mat,
                  file=sprintf("%s/%s_gene.txt", out.dir, gse),
                  row.names=TRUE, quote=FALSE)
    }
    geo.obj2 <- geo.obj
    geo.obj2$gene_mat <- gene_mat
    return(geo.obj2)
  })

  return(geo.obj.list)

}

#' Prepare data from an existing expression matrix
#'
#' @param expr.mat expression matrix
#' @param mapping.mat the matrix of, defaults to NULL
#' @param to.genes whether the data needs to be converted to genes, defaults to FALSE
#'        if this is TRUE, either the mapping matrix or platform.id need to be provided
#' @param platform.id the platform ID, defaults to NULL
#' @param gpl.dir directory where GPL ref data is location, empty defaults to tempdir()
#' @param pheno pheno data frame to add if desired, defaults to NULL
#'
#' @return list of gene expression matrices, with rows as genes and columns as samples
#'         values are gene expression levels
prepFromExpr <- function(expr.mat, mapping.mat=NULL, to.genes=FALSE, platform.id=NULL, gpl.dir=NULL, pheno=NULL){

  # // TODO: make sure the expression matrix is the correct format

  geo.obj <- list()

  # map to genes if this is required
  if (!to.genes){
    # warn if provide mapping data
    if (!is.null(mapping.mat) | !is.null(platform.id)){
      warning("You have provided a mapping matrix and/or platform id but specified not to convert the data to genes, so the data will remain as is.",
              call.=FALSE)
    }
    gene_mat <- expr.mat


  } else {
    geo.obj$expr <- expr.mat

    # raise an error - we need something to map
    .my_assert("Neither a platform id or mapping matrix is provided to map to genes",
               is.null(mapping.mat) & is.null(platform.id))

    if (!is.null(mapping.mat)){
      # // TODO: make sure the mapping matrix is in the correct format
      if (!is.null(platform.id)){
        warning("You have provided a mapping matrix and a platform id, we will use the matrix.",
                call.=FALSE)
      }
      ref_tab <- mapping.mat
    } else {
      # map based on the platform id

      if (!is.null(gpl.dir)){
        gpl.f <- sprintf("%s/%s_map.RData", gpl.dir, platform.id)
        if (!file.exists(gpl.f)){
          map_mult_gpl(c(platform.id), gpl.dir)
        }
        miceadds::load.Rdata(gpl.f, "ref_tab")
      } else {
        ref_tab <- parse_entrez_from_gpl(platform.id, verbose=FALSE)
      }
      # return NA if the mapping does not exist
      if (is.null(dim(ref_tab))){
        print("probe mapping data is missing so we cannot map")
        return(NA)
      }
    }
    gene_mat <- .convertGenes(expr.mat, ref_tab)

  }
  geo.obj$gene_mat <- gene_mat
  geo.obj$platform <- platform.id
  geo.obj$pheno <- pheno

  return(list(geo.obj))
}


#' Create a consensus gene list from a set of expression matrices.
#' The goal of this function is to get a unified list of genes across studies.
#'
#' This consensus gene list can then be used for reordering and ranking a dataset.
#'
#' @param list.exprs the list of dataset objects from prepFromGSE or prepFromExpr
#' @param min.fraction the fraction number of studies that we allow, defaults to 0.6
#' @param min.genes the minimum number of genes that the studies need to contain to be considered
#'   (if you use the package GPL mapping, it will already filter to a minimum threshold,
#'   but use this to adjust further)
#' @param na.max maximum NA fraction to keep a row, defaults to 0.3
#'
#' @return list of genes present in min.fraction of the studies
getConsensusGenes <- function(list.dats, min.fraction=0.6, min.genes=8000, na.max=0.3){

  # // TODO: add input checks

  gene.names <- sapply(list.dats, function(ds)
    sapply(ds, function(d) {
      if (length(d)==1){
        return(NA)
      }
      expr <- d$gene_mat;
      if (is.null(dim(expr)) | nrow(expr) < min.genes){
        return(NA)
      }
      # remove rows with large numbers of NAs
      na.counts <- apply(expr, 1, function(x) sum(is.na(x)))
      expr2 <- expr[floor(na.counts/ncol(expr)) <= na.max,]
      return(rownames(expr2))
    }))

  gene.names2 <- gene.names[!is.na(gene.names)]
  num.studies <- length(gene.names2)
  gene_counts <- table(unlist(gene.names2))
  gene.list <- names(gene_counts)[gene_counts > floor(min.fraction*num.studies)]

  return(gene.list)
}

#' Reorder and rank a gene expression matrix. We then use this for sex labeling and comparison.
#'
#' @param gene_mat the expression matrix to reorder and rank
#' @param gene_list the list of genes to use, this is helpful if you want to compare across multiple studies
#' @param to.ranks whether to convert to ranks, default is TRUE
#'
#' @return expression matrix reordered based on the gene list, with expression converted to ranks
reorderRank <- function(gene_mat, gene_list=list_genes, to.ranks=TRUE){
  # re-order based on gene list
  gene_mat2 <- .reorderGeneMat(gene_mat, gene_list)

  # convert to ranks
  if (to.ranks){
    gene_mat2 <- .expDataToRanks(gene_mat2)
  }
  return(gene_mat2)
}

#' Convert the expression matrix to genes.
#' This takes the average of probes to map them to a gene.
#'
#' @param expData the expression matrix with rows as probes and columns as conditions
#' @param probe_gene a mapping frame with a probe and gene column
#' @return an expression matrix with genes as rows
.convertGenes <- function(expData, probe_gene){
  .checkExprMatFormat(expData) # other input checks?
  probe_gene <- probe_gene[probe_gene$gene !="" & probe_gene$probe !="",]
  probe_gene$gene <- sapply(probe_gene$gene, as.character)
  probe_gene$probe <- sapply(probe_gene$probe, as.character)

  probe_gene2 <- dplyr::filter(probe_gene, probe %in% rownames(expData))

  gene.to.probe <- split(probe_gene2$probe, probe_gene2$gene)

  # filter to remove hugely multi-mapping??


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

  return(expData2.2)
}



#' Reorder the expression matrix so it contains the pre-specified list
#' of genes and fill in missing data with NAs.
#'
#' @param gene_mat the expression matrix to reorder
#' @param gene_list the list of genes to use
#'
#' @return a reordered expression mat with NAs for missing genes
.reorderGeneMat <- function(gene_mat, gene_list=list_genes){
  # create a data frame of NAs for missing genes
  missing.genes <- setdiff(gene_list, rownames(gene_mat))
  missing.vec <- rep(NA, ncol(gene_mat))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(gene_mat)

  # put together and reorder
  expDataPlusMiss <- rbind(gene_mat, missing.df )
  gene_mat2 <- expDataPlusMiss[gene_list,]
  return(gene_mat2)
}

#' Convert an expression dataset to ranks.
#'
#' Briefly, takes an expression matrix, and then ranks each column, 1 to n (number of genes)
#' in order of increasing expression values. Missing data is ranked last.
#'
#' @param gene_mat an expression matrix with genes as rows and columns as samples
#'
#' @return rank_dat ranked dataset
.expDataToRanks <- function(gene_mat) {
  rank_dat <- apply(gene_mat, 2, rank, na.last="keep") # keeps NAs as rank NA
  return(rank_dat)
}



