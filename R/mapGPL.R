
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
      if (!is.na(ref_tab)){
        write.table(ref_tab, file=sprintf("%s/%s_map.txt", ref_dir, gpl), sep="\t", quote=FALSE, row.names=FALSE)
      }
    })
    stopCluster(cl)

  } else {
    lapply(gpl.list, function(gpl){
      print(gpl)
      ref_tab <- parse_entrez_from_gpl(gpl, ref_dir)
      if (!is.na(ref_tab)){
        write.table(ref_tab, file=sprintf("%s/%s_map.txt", ref_dir, gpl), sep="\t", quote=FALSE, row.names=FALSE)
      }
      })
  }
}

#' Parse and map entrez data from a GPL annotation.
#'
#' @param gpl.name the name gpl to load
#' @param ref_dir optional reference directory, defaults to a temporary directory
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default is 8000
#' @return a data frame with "probe" and "gene" columns, where gene is the entrezid
parse_entrez_from_gpl <- function(gpl.name, ref_dir=NULL, MIN.OVERLAP=8000){
  if (is.null(ref_dir)){
    ref_dir <- tempdir()
  }

  # read in GPL data
  tryCatch({
    gpl <- GEOquery::getGEO(gpl.name)

  }, error = function(e) {
    print(sprintf("Error loading %s:", gpl.name))
    print(error)
    return(NA)
  })

  # check the organism
  organism <- gpl@header$taxid
  if (stringr::str_detect(organism, "9606")){
    org.name <- "human"
  } else if (stringr::str_detect(organism, "10090")){
    org.name <- "mouse"
  } else if (stringr::str_detect(organism, "10116")){
    org.name <- "rat"
  } else {
    print(sprintf("Error for %s, gpl key extraction is only implemented for human, rat, and mouse.",
                  gpl_name))
    return(NA)
  }

  # TODO - what if the probe column is the row name
  gpl.df <- gpl@dataTable@table
  probe.ids <- gpl.df[,1]

  # ------ Find columns labeled ENTREZ IDs ------ #
  entrez.col1 <- which(grepl("entrez",
                             gpl@dataTable@columns$Column, ignore.case = TRUE))
  entrez.col2 <- which(grepl("entrez",
                             gpl@dataTable@columns$Description, ignore.case = TRUE))
  entrez.col <- unique(entrez.col1, entrez.col2)
  if (length(entrez.col) != 0){
    # look for overlap
    overlap.lengths <- .check_entrez_overlap(gpl.df, entrez.col, org.name, ref_dir)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      return(.reform_entrez_df(df))
    }
  }

  # ------ Find any all integer columns, check for overlap ------ #
  int_only_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "^[0-9]+$"), na.rm=TRUE)
  })
  int_only_cols <- (int_only_cols > MIN.OVERLAP)
  col_idx <- which(int_only_cols)
  if (length(col_idx) != 0){
    overlap.lengths <- .check_entrez_overlap(gpl.df, col_idx, org.name, ref_dir)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      return(.reform_entrez_df(df))
    }
  }

  # ------ MAP FROM REFSEQ ------ #
  refseq_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "NM_[0-9]+"), na.rm=TRUE)
  })
  if (max(refseq_cols) > MIN.OVERLAP){
    refseq_col <- which.max(refseq_cols)[[1]]
    print("parsed from refseq")
    return(.map_from_refseq(gpl.df, refseq_col, org.name, ref_dir))
  }

  # ------ MAP FROM ENSEMBL ------ #
  ORG.ENSEMBL.SUFF <- list("human"="G", "mouse"="MUSG", "rat"="RNOG")
  org.ensembl.id <- sprintf("ENS%s", ORG.ENSEMBL.SUFF[[organism]])
  ensembl_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], sprintf("%s[0-9]+",org.ensembl.id)), na.rm=TRUE)
  })
  if (max(ensembl_cols) > MIN.OVERLAP){
    ensembl_col <- which.max(ensembl_cols)
    print("parsed from ensembl")
    return(.map_from_ensembl(gpl.df, ensembl_col, org.name, org.ensembl.id, ref_dir))
  }

  # ------ MAP FROM UNIGENE ----- #
  unigene_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "Hs.|Mm.|Rn."), na.rm=TRUE)
  })
  if (max(unigene_cols) > MIN.OVERLAP){
    unigene_col <- which.max(unigene_cols)
    print("parsed from unigene")
    return(.map_from_unigene(gpl.df, unigene_col, org.name, ref_dir))
  }

  # ------  MAP FROM HGNC  ------ #
  hgnc_cols <- sapply(1:ncol(gpl.df), function(i){
    grepl("GAPDH", gpl.df[,i], ignore.case=TRUE)
  })
  if (length(hgnc_cols) > 0){
    print("parsed from HGNC")
    return(.map_from_hgnc(gpl.df, hgnc_cols[[1]], org.name, ref_dir))

  }

  # ------ MAP FROM GenBank ------ #
  #rat_genbank <- paste((genbank %>% filter(entrezgene_id=="24383"))[1:5,"genbank"], collapse="|")
  #mouse_genbank <- paste((genbank %>% filter(entrezgene_id=="14433"))[1:5,"genbank"], collapse="|")
  #human_genbank <- paste((genbank %>% filter(entrezgene_id=="2597"))[1:5,"genbank"], collapse="|")
  LIST.GENBANK.STR <-
    list("rat"="AAA40814|AAA41193|AAA41795|AAB19105|AAH59110",
         "mouse" ="AAA37659|AAH82592|AAH83065|AAH83079|AAH83080",
         "human"="AAA52496|AAA52518|AAA52519|AAA53191|AAA86283")
  genbank_str <- LIST.GENBANK.STR[[org.name]]
  genbank_cols <- sapply(1:ncol(gpl.df), function(i){
    grepl(genbank_str, gpl.df[,i])
  })
  if (length(genbank_cols) > 0){
    print("parsed from GenBank")
    return(.map_from_genbank(gpl.df, genbank_cols[[1]], org.name, genbank_str, ref_dir))
  }
  print("No mapping")
  return(NA)
}


#' Helper function to reformat entrez gene/probe data frame.
#' Expands out the gene field so that there is one gene per row.a
#'
#' @param df the data frame to expand
#' @param entrez.col.name the name of the column containing entrez ids, defaults to "gene"
#' @return data reformatted with one entrez id per row
.reform_entrez_df <- function(df, entrez.col.name="gene"){
  df[,entrez.col.name] <- sapply(df[,entrez.col.name], function(x)
    paste(stringr::str_extract_all(x, "[0-9]+")[[1]], collapse=";"))
  tidyr::separate_rows(df, entrez.col.name, sep=";")
}


#' Heper function to find the location within a column of an ID.
#' Used for columns with many fields separated by double and triple slashes
#'
#' @param df the data.frame with feature data
#' @param my.col the column to examine
#' @param ex.str an example ID of that type
#' @return an expanded data frame with probe/gene data
.find_col_loc <- function(df, my.col, ex.str){

  # find example row
  ex_rows <- which(stringr::str_detect(df[,my.col], ex.str))
  my_str <- df[ex_rows[[1]], my.col]

  # split into fields
  mult_fields <- stringr::str_split(my_str, " /// ")[[1]]
  lst_fields <- stringr::str_split(mult_fields, " // ")
  idces <- lapply(lst_fields, function(x)
    which(sapply(x, function(y)
      stringr::str_detect(y, sprintf("^%s", ex.str)))))

  # find the location
  idx <- as.numeric(unique(unlist(idces))[[1]])

  # grab IDs from that location
  gene_vals <- lapply(df[,my.col], function(x){
    mult_fields <- stringr::str_split(x, " /// ")[[1]]
    lst_fields <- stringr::str_split(mult_fields, " // ")
    genes <- unique(lapply(lst_fields, function(x) x[[idx]]))
    genes <- genes[genes!="---"]
    paste(genes, collapse=" /// ")
  })

  # format into a data frame
  df2 <- data.frame(cbind("probe"=df[,1], "gene_col"=gene_vals))
  df3 <- tidyr::separate_rows(df2, gene_col, sep=" /// ")
  return(df3)
}

#' Helper function to parse a pattern out of a multi-column.
#' Used for columns with many fields contained within and a consistent
#' defined pattern. This is an alternate to `.find_col_loc()` for cases
#' where all strings match a pattern (e.g. "ENSG") and there may be
#' any type of delimiter within the column.
#'
#' @param df the data.frame with feature data
#' @param my.col the column to examine
#' @param pattern the pattern to extract
#' @return an expanded data frame with probe/gene data
.parse_multi_col <- function(df, my.col, pattern){

  # check whether it matches the whole pattern
  if (any(stringr::str_detect(df[,my.col], sprintf("^%s$", pattern)))){
    df2 <- df[,c(1, my.col)]
    colnames(df2) <- c("probe", "gene_col")
    df3 <- tidyr::separate_rows(df2, gene, sep=" /// ")
    return(df3)
  }

  # extract all of the data from that column
  gene_vals <- str_extract_all(df[,my.col], pattern)
  gene_vals_col <- lapply(gene_vals, function(x)
    paste(x, collapse=" /// "))
  probe.ids <- df[,1]
  df2 <- data.frame(cbind("probe"=probe.ids, "gene_col"=gene_vals_col))
  df3 <- tidyr::separate_rows(df2, gene_col, sep=" /// ")
  return(df3)
}

#' Count the length of overlap with entrez ids for a sequence of cols.
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param cols the columns you think contain entrez ids to check for overlap with
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return the length of overlap for all the columns, where names are the column indices
.check_entrez_overlap <- function(gpl.df, cols, organism, ref_dir=NULL){
  # grab the entrezids associated with the organism
  entrezids <- .load_ref(organism, "entrezids", ref_dir)

  # check the length of overlap
  overlap.lengths <- sapply(1:length(cols), function(i){
    length(intersect(entrezids, stringr::str_extract_all(gpl.df[,cols[i]], "[0-9]+")))
  })
  names(overlap.lengths) <- cols
  return(overlap.lengths)
}

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
  colnames(probe.gene) <- c("probe", "refseq_mrna")
  gene_map <- .load_ref(organism, "refseq", ref_dir)
  ref_entr <- gene_map %>%
    dplyr::select(refseq_mrna, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::filter(refseq_mrna != "") %>%
    unique()
  probe.gene$refseq_mrna <- sapply(probe.gene$refseq_mrna, function(x)
    unlist(stringr::str_split(x, "\\.")[[1]][[1]]))
  mapping.df <- dplyr::inner_join(probe.gene, ref_entr) %>%
    dplyr::rename(gene=entrezgene_id)
  mapping.df$probe <- sapply(mapping.df$probe, unlist)
  return(mapping.df)
}

#' Helper function to extract and map ensembl to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the ensembl data
#' @param organism the organism: rat, mouse, or human
#' @param org.ensembl.id the prefix for the organism (e.g. "ENSG")
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_ensembl <- function(gpl.df, my.col, organism, org.ensembl.id, ref_dir=NULL){
  library('dplyr') # for pipe
  probe.gene <- .parse_multi_col(gpl.df, my.col,
                                sprintf("[%s][\\d]+[.-]*[\\w\\d]*",
                                        org.ensembl.id))
  colnames(probe.gene) <- c("probe", "ensembl_gene_id")
  gene_map <- .load_ref(organism, "ensembl", ref_dir)

  ens_entr <- gene_map %>%
    dplyr::select(ensembl_gene_id, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, ens_entr) %>%
    dplyr::rename(gene=entrezgene_id)
}

#' Helper function to extract and map unigene to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the unigene data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_unigene <- function(gpl.df, my.col, organism, ref_dir=NULL){
  ORG.UNIGENE.PRE <- list("human"="Hs", "mouse"="Mm", "rat"="Rn")
  probe.gene <- .parse_multi_col(gpl.df, my.col, sprintf("%s.[0-9]+",ORG.UNIGENE.PRE[[organism]]))
  colnames(probe.gene) <- c("probe", "unigene")
  unigene <- .load_ref(organism, "unigene", ref_dir)

  dplyr::inner_join(probe.gene, unigene) %>%
    dplyr::rename(gene=entrezgene_id)
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
  if (organism=="rat"){
    print("error no HGNC symbol mapping for rat")
    return(NA)
  }
  probe.gene <- .find_col_loc(gpl.df, my.col, "GAPDH|gapdh|Gapdh")
  colnames(probe.gene) <- c("probe", "hgnc_symbol")
  gene_map <- .load_ref(organism, "hgnc", ref_dir)

  hgnc_entr <- gene_map %>%
    dplyr::select(hgnc_symbol, entrezgene_id) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, ref_entr) %>%
    dplyr::rename(gene=entrezgene_id)
}

#' Helper function to extract and map genbank to entrez
#'
#' @param gpl.df the feature data from the GPL data frame
#' @param my.col the column containing the genbank data
#' @param organism the organism: rat, mouse, or human
#' @param ref_dir option to set reference directory  (defaults to tempdir)
#' @return a parsed probe/gene data frame
.map_from_genbank <- function(gpl.df, my.col, organism, genbank_str, ref_dir=NULL){
  probe.gene <- .find_col_loc(gpl.df, my.col, genbank_str)
  colnames(probe.gene) <- c("probe", "genbank")

  genbank <- .load_ref(organism, "genbank", ref_dir)

  res <- dplyr::inner_join(probe.gene, genbank)
  dplyr::rename(res, gene=entrezgene_id)
}


