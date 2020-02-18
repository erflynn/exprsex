

#' Wrapper to parse and map data from multiple gpls.
#'
#' @param gpl.list list of platforms to map
#' @param ref_dir optional reference directory, defaults to a temporary directory
#' @param parallelize optionally parallelize the mapping
map_mult_gpl <- function(gpl.list, ref_dir=NULL, parallelize=FALSE){
  if (is.null(ref_dir)){
    ref_dir <- tempdir()
  }
  if (parallelize){
    library('parallel')
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    parLapply(cl, gpl.list, function(gpl){
      ref_tab <- parse_entrez_from_gpl(gpl, ref_dir)
      if (is.data.frame(ref_tab)){
        save(ref_tab, file=sprintf("%s/%s_map.RData", ref_dir, gpl))
      }
    })
    stopCluster(cl)

  } else {
    lapply(gpl.list, function(gpl){
      print(gpl)
      ref_tab <- parse_entrez_from_gpl(gpl, ref_dir)
      if (is.data.frame(ref_tab)){
        save(ref_tab, file=sprintf("%s/%s_map.RData", ref_dir, gpl))
      }
    })
  }
}


#' Parse and map entrez data from a GPL annotation.
#'
#' @param gpl.name the name gpl to load
#' @param ref_dir optional reference directory, defaults to a temporary directory
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default is 8000
#' @param verbose print addition logging data if it fails, defaults to TRUE
#' @return a data frame with "probe" and "gene" columns, where gene is the entrezid
parse_entrez_from_gpl <- function(gpl.name, ref_dir=NULL, MIN.OVERLAP=8000, verbose=TRUE){
  if (is.null(ref_dir)){
    ref_dir <- tempdir()
  }

  # read in GPL data
  tryCatch({
    gpl <- GEOquery::getGEO(gpl.name)

  }, error = function(e) {
    print(sprintf("Error loading %s:", gpl.name))
    print(e)
    return(NA)
  })

  # check the organism
  organism <- paste(gpl@header$taxid, collapse=";")
  if (stringr::str_detect(organism, "9606")){
    org.name <- "human"
  } else if (stringr::str_detect(organism, "10090")){
    org.name <- "mouse"
  } else if (stringr::str_detect(organism, "10116")){
    org.name <- "rat"
  } else {
    print(sprintf("Error for %s, gpl key extraction is only implemented for human, rat, and mouse.",
                  gpl.name))
    return(NA)
  }

  gpl.df <- gpl@dataTable@table
  probe.ids <- gpl.df[,1]

  log_str <- sprintf("LOG %s:",gpl.name)

  # ------ Find columns labeled ENTREZ IDs ------ #
  entrez.col <- .detect_entrez_ids_name(gpl@dataTable@columns)
  if (length(entrez.col) != 0){
    log_str <- sprintf("%s\n\tFound entrez column by name %s", log_str,  paste(entrez.col, collapse=";"))
    # look for overlap
    overlap.lengths <- .check_entrez_overlap(gpl.df, entrez.col,
                                             org.name, ref_dir)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      print(sprintf("parsed %s from entrez", gpl.name))
      return(.reform_entrez_df(df))
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, max(overlap.lengths))

  }

  # ------ Find any all integer columns, check for overlap w entrez ------ #
  entrez.col1 <- .detect_int_cols(gpl.df, MIN.OVERLAP)

  if (length(entrez.col1) != 0){
    log_str <- sprintf("%s\n\tFound int only column %s", log_str,  paste(entrez.col1, collapse=";"))
    overlap.lengths <- .check_entrez_overlap(gpl.df, entrez.col1, org.name, ref_dir)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      print(sprintf("parsed %s from entrez", gpl.name))
      return(.reform_entrez_df(df))
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, max(overlap.lengths))

  }

  # ------ Look for integers *within* a column, check for overlap w entrez ----- #
  entrez.col2 <- .find_entrez_w_in_block(gpl.df, MIN.OVERLAP)
  entrez.col3 <- setdiff(entrez.col2, entrez.col1) # remove already examined columns
  if (length(entrez.col3) != 0){
    log_str <- sprintf("%s\n\tFound int within column %s", log_str,  paste(entrez.col3, collapse=";"))
    overlap.lengths <- .check_entrez_overlap(gpl.df, entrez.col3, org.name, ref_dir)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      ref <- .load_ref(org.name, "entrezids", ref_dir)
      mapped <- .find_col_loc(gpl.df, entrez.col, "[0-9]+", check.overlap=TRUE, ref)
      print(sprintf("parsed %s from entrez", gpl.name))
      return(dplyr::rename(mapped, gene=gene_col))
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, max(overlap.lengths))
  }

  # ------ MAP FROM GenBank ------ #
  genbank_col <- .detect_genbank_cols_name(gpl@dataTable@columns)
  if (length(genbank_col) != 0){
    # TODO - do this differently if this is a genbank column...

    mappings <- lapply(genbank_col, function(my_col){
        log_str <- sprintf("%s\n\tFound genbank column by name %s", log_str, paste(genbank_col, collapse=";"))
        mapped <- .map_from_genbank(gpl.df, my_col, org.name, ref_dir, col.given=TRUE)
        return(mapped)
      })
    if (length(mappings) > 1){
        mapped <- mappings[which.max(lapply(mappings, nrow))]
      } else {
        mapped <- mappings[[1]]
      }

    if (length(unique(mapped$gene)) > MIN.OVERLAP ){
      print(sprintf("parsed %s from GenBank", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))
  }
  genbank_col <- .detect_genbank_cols(gpl.df, org.name)
  if (length(genbank_col) != 0){
    # TODO - do this differently if this is a genbank column...

    mappings <- lapply(genbank_col, function(my_col){
      log_str <- sprintf("%s\n\tFound genbank column by name %s", log_str, paste(genbank_col, collapse=";"))
      mapped <- .map_from_genbank(gpl.df, my_col, org.name, ref_dir, col.given=FALSE)
      return(mapped)
    })
    if (length(mappings) > 1){
      mapped <- mappings[which.max(lapply(mappings, nrow))]
    } else {
      mapped <- mappings[[1]]
    }

    if (length(unique(mapped$gene)) > MIN.OVERLAP ){
      print(sprintf("parsed %s from GenBank", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))
  }
  # <-- genbank is done first because often contains REFSEQ cols ---> #

  # ------ MAP FROM REFSEQ ------ #
  refseq_col <- .detect_refseq_cols(gpl.df, MIN.OVERLAP)

  if (length(refseq_col) != 0){
    log_str <- sprintf("%s\n\tFound refseq ids in a column %s", log_str,  paste(refseq_col, collapse=";"))
    mapped <- .map_from_refseq(gpl.df, refseq_col, org.name, ref_dir)
    if (length(unique(mapped$gene)) > MIN.OVERLAP){
      print(sprintf("parsed %s from refseq", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))
  }

  # ------ MAP FROM ENSEMBL ------ #
  ensembl_col <- .detect_ensembl_cols(gpl.df, org.name, MIN.OVERLAP)
  if (length(ensembl_col) != 0){
    log_str <- sprintf("%s\n\tFound ensembl ids in a column %s", log_str,  paste(ensembl_col, collapse=";"))

    mapped <- .map_from_ensembl(gpl.df, ensembl_col, org.name, ref_dir)
    if (length(unique(mapped$gene)) > MIN.OVERLAP){
      print(sprintf("parsed %s from ensembl", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))

  }

  # ------  MAP FROM HGNC  ------ #
  hgnc_col <- .detect_hgnc_cols(gpl.df)
  if (length(hgnc_col) != 0){
    log_str <- sprintf("%s\n\tFound hgnc ids in a column %s", log_str,  paste(hgnc_col, collapse=";"))
    mapped <- .map_from_hgnc(gpl.df, hgnc_col, org.name, ref_dir)
    if (length(unique(mapped$gene)) > MIN.OVERLAP){
      print(sprintf("parsed %s from HGNC", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))

  }

  # ------ MAP FROM UNIGENE ----- #
  unigene_col <- .detect_unigene_cols(gpl.df, MIN.OVERLAP)
  if (length(unigene_col) != 0){
    log_str <- sprintf("%s\n\tFound unigene ids in a column %s", log_str,  paste(unigene_col, collapse=";"))

    mapped <- .map_from_unigene(gpl.df, unigene_col, org.name, ref_dir)
    if (length(unique(mapped$gene)) > MIN.OVERLAP){
      print(sprintf("parsed %s from unigene", gpl.name))
      return(mapped)
    }
    log_str <- sprintf("%s\n\t\tAfter extraction, insufficient ids (n=%s)", log_str, length(unique(mapped$gene)))

  }

  print(sprintf("No mapping for %s", gpl.name))
  if (verbose){print(log_str)}
  return(NA)
}


#' Helper function to reformat entrez gene/probe data frame.
#' Expands out the gene field so that there is one gene per row.
#'
#' @param df the data frame to expand
#' @param entrez.col.name the name of the column containing entrez ids, defaults to "gene"
#' @return data reformatted with one entrez id per row
.reform_entrez_df <- function(df, entrez.col.name="gene"){
  df[,entrez.col.name] <- sapply(df[,entrez.col.name], function(x)
    paste(stringr::str_extract_all(x, "[0-9]+")[[1]], collapse=";"))
  .clean_mapping(tidyr::separate_rows(df, entrez.col.name, sep=";"))
}



# ------------------- MAPPING ------------------------ #

#' Helper function to extract and map refseq to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the refseq data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_refseq <- function(gpl.df, my.col, organism, ref_dir=NULL){
  # check that it matches the whole column -- if not, parse
  library(dplyr) # for pipe
  probe.gene <- .find_col_loc(gpl.df, my.col, "N[R|M][_][\\d]+[_.-]*[\\w\\d]*")
  if (nrow(probe.gene)==0){
    return(probe.gene)
  }
  colnames(probe.gene) <- c("probe", "refseq_mrna")

  gene_map <- .load_ref(organism, "refseq", ref_dir)
  ref_entr <- gene_map %>%
    dplyr::select(refseq_mrna, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::filter(refseq_mrna != "") %>%
    unique()

  INTERNAL.SEP <-"///|//|,|\\|"
  SUFFIX.SEP <- "\\."
  probe.gene$refseq_mrna <- sapply(probe.gene$refseq_mrna, function(x)
    unlist(stringr::str_split(x, SUFFIX.SEP)[[1]][[1]]))
  mapping.df <- dplyr::inner_join(probe.gene, ref_entr) %>%
    dplyr::rename(gene=entrezgene_id)
  mapping.df$probe <- sapply(mapping.df$probe, unlist)
  return(.clean_mapping(mapping.df))
}

#' Helper function to extract and map ensembl to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the ensembl data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_ensembl <- function(gpl.df, my.col, organism, ref_dir=NULL){
  library('dplyr') # for pipe

  ORG.ENSEMBL.SUFF <- list("human"="G", "mouse"="MUSG", "rat"="RNOG")
  org.ensembl.id <- sprintf("ENS%s", ORG.ENSEMBL.SUFF[[organism]])

  probe.gene <- .parse_multi_col(gpl.df, my.col,
                                sprintf("%s[\\d]+[.-]*[\\w\\d]*",
                                        org.ensembl.id))
  colnames(probe.gene) <- c("probe", "ensembl_gene_id")

  INTERNAL.SEP <-"///|//|,|\\|"
  SUFFIX.SEP <- "\\."
  probe.gene$ensembl_gene_id <- sapply(probe.gene$ensembl_gene_id, function(x)
    unlist(stringr::str_split(x, SUFFIX.SEP)[[1]][[1]]))
  gene_map <- .load_ref(organism, "ensembl", ref_dir)

  ens_entr <- gene_map %>%
    dplyr::select(ensembl_gene_id, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
   probe.gene %>% dplyr::filter(ensembl_gene_id != "") %>%
    dplyr::inner_join(ens_entr) %>%
    dplyr::rename(gene=entrezgene_id) %>%
    .clean_mapping()
}

#' Helper function to extract and map unigene to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the unigene data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_unigene <- function(gpl.df, my.col, organism, ref_dir=NULL){
  library('dplyr') # for pipe
  ORG.UNIGENE.PRE <- list("human"="Hs", "mouse"="Mm", "rat"="Rn")
  probe.gene <- .parse_multi_col(gpl.df, my.col, sprintf("%s\\.[0-9]+",ORG.UNIGENE.PRE[[organism]]))
  colnames(probe.gene) <- c("probe", "unigene")
  unigene <- .load_ref(organism, "unigene", ref_dir)

  dplyr::inner_join(probe.gene, unigene) %>%
    dplyr::rename(gene=entrezgene_id) %>%
    .clean_mapping()
}

#' Helper function to extract and map HGNC to entrez
#' Note - this does not work for rat.
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the HGNC data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_hgnc <- function(gpl.df, my.col, organism, ref_dir=NULL){
  library('dplyr') # for pipe

  if (organism=="rat"){
    print("error no HGNC symbol mapping for rat")
    return(data.frame("gene"=c(), "probe"=c()))
  }
  probe.gene <- .find_col_loc(gpl.df, my.col, "GAPDH|gapdh|Gapdh")
  if (nrow(probe.gene)==0){
    return(data.frame("gene"=c(), "probe"=c()))
  }
  colnames(probe.gene) <- c("probe", "hgnc_symbol")
  gene_map <- .load_ref(organism, "hgnc", ref_dir)

  hgnc_entr <- gene_map %>%
    dplyr::select(hgnc_symbol, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, hgnc_entr) %>%
    dplyr::rename(gene=entrezgene_id) %>%
    .clean_mapping()
}

#' Helper function to extract and map genbank to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the genbank data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @param col.given whether we know the column by name --> don't find loc (defaults to FALSE)
#' @return a parsed probe/gene data frame
.map_from_genbank <- function(gpl.df, my.col, organism, ref_dir=NULL, col.given=FALSE){
  library('dplyr') # for pipe

  LIST.GENBANK.STR <-
    list("rat"="AAA40814|AAA41193|AAA41795|AAB19105|AAH59110",
         "mouse" ="AAA37659|AAH82592|AAH83065|AAH83079|AAH83080",
         "human"="AAA52496|AAA52518|AAA52519|AAA53191|AAA86283")
  genbank_str <- LIST.GENBANK.STR[[organism]]

  if (!col.given) {
    probe.gene <- .find_col_loc(gpl.df, my.col, genbank_str)
    colnames(probe.gene) <- c("probe", "genbank")

  } else {
    probe.gene <- gpl.df[,c(1,my.col)]
    colnames(probe.gene) <- c("probe", "genbank")

  }

  INTERNAL.SEP <-"///|//|,|\\|"
  SUFFIX.SEP <- "\\."
  probe.gene <- tidyr::separate_rows(probe.gene, genbank, sep=INTERNAL.SEP)
  probe.gene <- dplyr::mutate(probe.gene, genbank=stringr::str_trim(genbank))
  probe.gene$genbank <- sapply(probe.gene$genbank, function(x)
    unlist(stringr::str_split(x, SUFFIX.SEP)[[1]][[1]]))

  genbank <- .load_ref(organism, "genbank", ref_dir)

  res <- dplyr::inner_join(probe.gene, genbank)
  dplyr::rename(res, gene=entrezgene_id) %>%
    .clean_mapping()
}


#' Clean up the mapping output
#'
#' @param mapped a data frame of mapped data
#' @return cleaned mapping with NAs/empty strings removed
.clean_mapping <- function(mapped){
  mapped2 <- unique(dplyr::select(mapped, gene, probe))
  return(dplyr::filter(mapped2, !is.na(gene) & !is.na(probe)))
}


