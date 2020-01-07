
#' Load the desired reference data mapping table.
#' The reference tables are used later for mapping gpl data to entrez ids.
#' This first checks if the mapping exists in the ref_dir and then if not,
#' calls a helper function to extract it.
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @param type type of data desired - one of "ensembl", "refseq", "hgnc", "unigene", "genbank", or "entrezids"
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @param save_dat whether to save the data in the reference directory, defaults to TRUE
#' @return the desired mapping data
.load_ref <- function(organism, type, ref_dir=NULL, save_dat=TRUE){
  if (is.null(ref_dir)){
    ref_dir <- tempdir()
  }
  if (type %in% c("ensembl", "refseq", "hgnc", "gene_map")){
    type <- "gene_map"
  }
  # checks first for the data - if it exists, load it
  dat.path <- sprintf("%s/%s_%s.RData", ref_dir, organism, type)
  if (file.exists(dat.path)){
    miceadds::load.Rdata(dat.path, "ref_dat")
    return(ref_dat)
  } else {
    # if the data does not exist - generate it and save it in the ref_dir
    ref_dat <- .get_ref_data(organism, type, ref_dir)
    if (save_dat){
      save(ref_dat, file=dat.path)
    }
    return(ref_dat)
  }
}

#' Generate the desired mapping data for that organism
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @param type type of data desired
#' @param ref_dir the reference directory to pull from, this is either provided or a tempdir
#' @return the desired mapping data
.get_ref_data <- function(organism, type, ref_dir){
  if (type == "gene_map"){
    ref_dat <- .extract_biomart_ref(organism)
  } else if (type == "entrezids") {
    ref_dat <- .extract_entrez(organism, ref_dir)
  } else if (type == "genbank"){
    ref_dat <- .get_genbank(organism)
  } else if (type =="unigene"){
    ref_dat <- .get_unigene(organism, ref_dir)
  } else if (type == "xy_genes"){
    ref_dat <- .get_xy_genes(organism, ref_dir)
  } else {
    .my_assert("error type not implemented", 1==2)
  }
  return(ref_dat)
}


#' Gets a mapping of unigene ids to entrez ids.
#' This pulls the list from the entrezids, which is either in the ref_dir or is generated.
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @param ref_dir the reference directory to pull from, this is either provided or a tempdir
#' @return table with unigene and entrez ids
.get_unigene <- function(organism, ref_dir){
  library(mygene)
  entrezids <- .load_ref(organism, "entrezids", ref_dir)

  df_unigene <- queryMany(entrezids,
                          scopes = "entrezgene",
                          fields = "unigene",
                          pecies = organism,
                          returnall = TRUE)

  table(sapply(df_unigene$response[,"unigene"], is.null))
  df_unigene$response$unigene <- sapply(df_unigene$response$unigene,
                                        function(x) paste(x, collapse=";"))

  # remove null entries
  my_df <- df_unigene$response[df_unigene$response$unigene!="",]

  unigene <- data.frame(my_df) %>%
    dplyr::select("query", "unigene") %>%
    dplyr::rename(entrezgene_id=query) %>%
    tidyr::separate_rows(unigene, sep=";")
  return(unigene)
}

#' Gets a mapping of genbank to entrez.
#' Note - these tables are very large.
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @return a table with genbank and entrez ids
.get_genbank <- function(organism){

  if (organism=="human"){
    library(org.Hs.eg.db)
    x <- org.Hs.egACCNUM
  } else if (organism == "mouse"){
    library(org.Mm.eg.db)
    x <- org.Mm.egACCNUM
  } else if (organism == "rat"){
    library(org.Rn.eg.db)
    x <- org.Rn.egACCNUM
  } else {
    print("Only implemented for human, rat, and mouse")
    return(NA)
  }

  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])

  # convert to a dataframe with one mapping per row
  collapsed <- lapply(xx, function(y) paste(unique(y), collapse=" /// "))
  df <- data.frame(cbind("entrezgene_id"= names(xx), "genbank"=collapsed)  )
  rownames(df) <- NULL
  df2 <- tidyr::separate_rows(df, genbank, sep=" /// ")
  df2$entrezgene_id <- sapply(df2$entrezgene_id, unlist)

  return(df2)
}

#' Gets a list of all entrez ids for an organism.
#' This pulls the list from the gene map, which is either in the ref_dir or is generated.
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @param ref_dir the reference directory to pull from, this is either provided or a tempdir
#' @return a list of entrez ids
.extract_entrez <- function(organism, ref_dir){
  library('dplyr') # for #%>%
  dat.map <- .load_ref(organism, "gene_map", ref_dir)

  entrez_df <- dat.map %>% dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::select(entrezgene_id) %>% unique()
  return(entrez_df)
}

#' Extract a gene mapping table from biomaRt, this contains ensembl, refseq, entrez ids
#' and chromsome names. HGNC symbols are included for human and mouse.
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @return mapping table for use in converting between IDs
.extract_biomart_ref <- function(organism){

  ORG.STR.LIST <- list("human"="hsapiens", "rat"="rnorvegicus", "mouse"="mmusculus")
  ensembl <- biomaRt::useMart("ensembl", dataset=sprintf("%s_gene_ensembl",
                                                ORG.STR.LIST[[organism]]))
  if (organism %in% c("human", "mouse")){
    attribute_list <- c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "entrezgene_id",  "chromosome_name")
  } else {
    attribute_list <- c("ensembl_gene_id",  "refseq_mrna", "entrezgene_id",  "chromosome_name")
  }
  gene_map <- biomaRt::getBM(attributes=attribute_list, mart=ensembl)

  # clean up the data
  gene_map <- .cleanMapping(gene_map)
  # note: we are retaining data with odd chromosome mapping

  return(gene_map)
}



#' Gets a table of all the X,Y genes for that organism with the entrez ids
#'
#' @param organism the organism to grab data for: rat, mouse, or human
#' @param ref_dir the reference directory to pull from, this is either provided or a tempdir
#' @return a table with two columns "gene", the entrez gene id, and "chr",
#'        the chromosome name, which will be one of X and Y
.get_xy_genes <- function(organism, ref_dir){
  dat.map <- .load_ref(organism, "gene_map", ref_dir)
  dat.map %>% dplyr::filter(chromosome_name %in% c("X", "Y")) %>%
    dplyr::select(entrezgene_id, chromosome_name) %>%
    dplyr::arrange(chromosome_name) %>%
    dplyr::rename("gene"=entrezgene_id, "chr"=chromosome_name)
}

#' Code to generate all reference tables for rat, mouse, and human in a reference directory.
#' This is particularly useful if you need to load data for a ton of different GPLs.
#' Note that this produces about 400 MB of data - so make sure you have the space.
#'
#' @export
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @param organisms specify which organisms you want, can be any of "mouse", "rat"
#'                  or "human" (defaults to all three)
generate_all_ref <- function(ref_dir=NULL, organisms=c("mouse", "rat", "human")){
  # check organisms are in mouse, rat, or human
  .my_assert("Organisms must be one of mouse, rat, or human",
             length(setdiff(organisms, c("mouse", "rat", "human")))==0)

  if (is.null(ref_dir)){
    ref_dir <- tempdir()
  } else { # make sure the ref_dir exists

  }
  ref.types=c("gene_map", "entrezids", "xy_genes", "genbank", "unigene")
  lapply(organisms, function(organism){
    lapply(ref.types, function(type){
      tryCatch({
        res <- .load_ref(organism, type, ref_dir, save_dat=TRUE)
        rm(res)
        return(1)
      }, error = function(err){
        print(err)
        print(sprintf("error loading %s, %s", organism, type))
        return(0)
      })
    })
  })
}


# clean mappings
.cleanMapping <- function(df){
  # entrezids should be converted to characters
  df <- .convertEntrez(df)
  # missing data should be set up as "NA"s
  df <- .convertNullChar(df)
  # the data should be expanded out so that it is tidy
  #df <- .checkExpandCol(df)
  return(df)
}

# convert Entrez to character
.convertEntrez <- function(df){
  dplyr::mutate(df, entrezgene_id=as.character(entrezgene_id))
}

# convert alternate null characters to NAs
.convertNullChar <- function(df){
  df[df==""] <- NA
  df[df=="NA"] <- NA
  return(df)
}

# .checkExpandCol <- function(df){
#
#   check_slash <- apply(df, 2, function(col) any(str_detect(col, "//"), na.rm=TRUE))
#   check_semi <- apply(df, 2, function(col) any(str_detect(col, ";"), na.rm=TRUE))
#   if (any(check_slash)){
#
#   }
#   if (any(check_semi)){
#     # expand the column
#     which(check_semi)
#   }
#   return(df)
# }



