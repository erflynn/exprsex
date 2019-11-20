
#' Looks for entrez ids by trying to find a column name named "Entrez"
#'
#' @param gpl_cols the column data from a GPL
#' @return a list of columns containing these locations, an empty list if no columns
.detect_entrez_ids_name <- function(gpl_cols){
  entrez.col1 <- which(grepl("entrez",
                             gpl_cols$Column, ignore.case = TRUE))
  entrez.col2 <- which(grepl("entrez",
                             gpl_cols$Description, ignore.case = TRUE))
  entrez.col <- unique(entrez.col1, entrez.col2)
  return(entrez.col)
}

#' Looks for genbank cols by trying to find a column name named "GB_ACC"
#' or genbank
#'
#' @param gpl_cols the column data from a GPL
#' @return a list of columns containing these locations, an empty list if no columns
.detect_genbank_cols_name <- function(gpl_cols){
  genbank.col1 <- which(grepl("GB[_|\\.][ACC|LIST]",
                             gpl_cols$Column, ignore.case = TRUE))
  genbank.col2 <- which(grepl("Genbank",
                             gpl_cols$Column, ignore.case = TRUE))
  genbank.col <- unique(genbank.col1, genbank.col2)
  return(genbank.col)
}

#' Identify integer only only columns, hopefully these are Entrez.
#' (that will be tested later)
#' @param gpl.df the gpl data frame
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default=8000
#' @return a list of the indices of the columns containing these ids
.detect_int_cols <- function(gpl.df, MIN.OVERLAP=8000){
  int_only_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "^[0-9]+$"), na.rm=TRUE)
  })
  col_idx <- which(int_only_cols > MIN.OVERLAP)
  return(col_idx)
}

#' Identify columns that might have entrez in them. This would be
#' a block column where the entrez column is separated by "///" or ","
#'
#' @param gpl.df the gpl data frame
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default=8000
#' @return a list of the indices of the columns containing these ids
.find_entrez_w_in_block <- function(gpl.df, MIN.OVERLAP=8000){
  int_cols <- sapply(1:ncol(gpl.df), function(i){
    ids <- stringr::str_replace_all(unlist(stringr::str_split(gpl.df[,i],
                                                              "///|//|,")), " ", "")
    sum(stringr::str_detect(ids, "^[0-9]+$"), na.rm=TRUE)
  })
  return(which(int_cols > MIN.OVERLAP))
}

#' Identify the refseq columns and return the one with the most overlap.
#' (Note - this looks for "NM_" ids, which is not comprehensive.)
#'
#' @param gpl.df the gpl data frame
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default=8000
#' @return the index of the column containing these ids
.detect_refseq_cols <- function(gpl.df, MIN.OVERLAP=8000 ){
  refseq_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "NM_[0-9]+"), na.rm=TRUE)
  })
  if (max(refseq_cols) > MIN.OVERLAP){
    refseq_col <- which.max(refseq_cols)[[1]]
    return(refseq_col)
  }
  return(c())
}

#' Identify the ensembl columns and return the one with the max overlap.
#'
#' @param gpl.df the gpl data frame
#' @param organism the organism we are looking at
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default=8000
#' @return the index of the column containing these ids
.detect_ensembl_cols <- function(gpl.df, organism, MIN.OVERLAP=8000){
  ORG.ENSEMBL.SUFF <- list("human"="G", "mouse"="MUSG", "rat"="RNOG")
  org.ensembl.id <- sprintf("ENS%s", ORG.ENSEMBL.SUFF[[organism]])
  ensembl_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], sprintf("%s[0-9]+",org.ensembl.id)), na.rm=TRUE)
  })
  if (max(ensembl_cols) > MIN.OVERLAP){
    ensembl_col <- which.max(ensembl_cols)[[1]]
    return(ensembl_col)
  }
  return(c())
}

#' Identify the unigene columns and return the one with the most overalp
#'
#' @param gpl.df the gpl data frame
#' @param MIN.OVERLAP the minimum number of mapped genes allowed, default=8000
#' @return the index of the column containing these ids
.detect_unigene_cols <- function(gpl.df, MIN.OVERLAP=8000){
  unigene_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "(Hs\\.[0-9]+)|(Rn\\.[0-9]+)|(Mm\\.[0-9]+)"), na.rm=TRUE)
  })
  if (max(unigene_cols) > MIN.OVERLAP){
    unigene_col <- which.max(unigene_cols)[[1]]
    return(unigene_col)
  }
  return(c())
}

#' Identify the HGNC columns and return the first
#' (Note - this is by looking for Gapdh which may not be the best metric...)
#'
#' @param gpl.df the gpl data frame
#' @return the index of the column containing these ids
.detect_hgnc_cols <- function(gpl.df){
  hgnc_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(grepl("GAPDH", gpl.df[,i], ignore.case=TRUE), na.rm=TRUE)
  })
  hgnc_idces <- which(hgnc_cols > 0)
  if (length(hgnc_idces) > 0){
    return(which.max(hgnc_cols)[[1]])
  }
  return(c())
}


#' Identify the genbank columns and return the first
#' (Note - this is by looking for a handful of ids...)
#'
#' @param gpl.df the gpl data frame
#' @param organism the organism we are looking at
#' @return the index of the column containing these ids
.detect_genbank_cols <- function(gpl.df, organism){
  #rat_genbank <- paste((genbank %>% filter(entrezgene_id=="24383"))[1:5,"genbank"], collapse="|")
  #mouse_genbank <- paste((genbank %>% filter(entrezgene_id=="14433"))[1:5,"genbank"], collapse="|")
  #human_genbank <- paste((genbank %>% filter(entrezgene_id=="2597"))[1:5,"genbank"], collapse="|")
  LIST.GENBANK.STR <-
    list("rat"="AAA40814|AAA41193|AAA41795|AAB19105|AAH59110",
         "mouse" ="AAA37659|AAH82592|AAH83065|AAH83079|AAH83080",
         "human"="AAA52496|AAA52518|AAA52519|AAA53191|AAA86283")
  genbank_str <- LIST.GENBANK.STR[[organism]]
  genbank_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(grepl(genbank_str, gpl.df[,i]), na.rm=TRUE)
  })
  genbank_idces <- which(genbank_cols > 0)
  if (length(genbank_idces) > 0){
    return(which.max(genbank_cols)[[1]])
  }
  return(c())
}



#' Heper function to find the location within a column of an ID.
#' Used for columns with many fields separated by double and triple slashes
#'
#' @param df the data.frame with feature data
#' @param my.col the column to examine
#' @param ex.str an example ID of that type
#' @param check.overlap defaults to false
#' @param id.list defaults to null
#' @return an expanded data frame with probe/gene data
.find_col_loc <- function(df, my.col, ex.str,
                          check.overlap=FALSE,
                          id.list=NULL){

  # find example row
  ex_rows <- which(stringr::str_detect(df[,my.col], ex.str))

  my.row <- ex_rows[[1]]

  if (check.overlap & !is.null(id.list)){
    i <- 1
    while(i < length(ex_rows)){
      my.row <- ex_rows[[i]]
      tokens <- stringr::str_extract_all(df[my.row, my.col], ex.str)[[1]]
      tokens2 <- tokens[stringr::str_length(tokens) >4]
      tokens3 <- intersect(tokens2, id.list)
      if (length(tokens3) > 0){
        ex.str <- tokens3[[1]]
        i <- length(ex_rows) # break out of the loop because we have found an example!
      }
      i <- i+1
    }
  }
  my_str <- df[my.row, my.col]

  # split into fields
  mult_fields <- stringr::str_trim(stringr::str_split(my_str, "///")[[1]])
  lst_fields <- stringr::str_split(mult_fields, "//")
  idces <- lapply(lst_fields, function(x)
    which(sapply(x, function(y)
      stringr::str_detect(stringr::str_trim(y), sprintf("^%s", ex.str)))))

  # find the location
  idx <- as.numeric(unique(unlist(idces))[[1]])

  # grab IDs from that location
  gene_vals <- lapply(df[,my.col], function(x){
    mult_fields <- stringr::str_trim(stringr::str_split(x, "///")[[1]])
    lst_fields <- stringr::str_split(mult_fields, "//")
    genes <- unique(lapply(lst_fields[sapply(lst_fields, length) >= idx],
                           function(x) stringr::str_trim(x[[idx]])))
    genes <- genes[genes!="---"]

    paste(genes, collapse="///")
  })

  # format into a data frame
  df2 <- data.frame(cbind("probe"=df[,1], "gene_col"=gene_vals))
  df3 <- tidyr::separate_rows(df2, gene_col, sep="///")
  return(dplyr::filter(df3, gene_col != ""))
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
    df3 <- tidyr::separate_rows(df2, gene_col, sep="///|,")
    return(df3)
  }

  # extract all of the data from that column
  gene_vals <- stringr::str_extract_all(df[,my.col], pattern)
  gene_vals_col <- lapply(gene_vals, function(x)
    paste(x, collapse="///"))
  probe.ids <- df[,1]
  df2 <- data.frame(cbind("probe"=probe.ids, "gene_col"=gene_vals_col))
  df3 <- tidyr::separate_rows(df2, gene_col, sep="///")
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
    length(intersect(entrezids,
                     unlist(stringr::str_extract_all(gpl.df[,cols[i]], "[0-9]+"))))
  })
  names(overlap.lengths) <- cols
  return(overlap.lengths)
}
