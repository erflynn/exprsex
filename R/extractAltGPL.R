
# <---- DEAL WITH GPL PARSING ERRORS -----> #
# Note - all of this is adapted from MetaIntegrator to handle edge cases,
# this is not code I have written from scratch
# https://cran.r-project.org/web/packages/MetaIntegrator/index.html

#' Alternate GPL parsing
#'
.parseGenesAlt <- function(fData_table){
  # slight edits to MetaIntegrator to allow for these patterns
  symbol_pos <- which(sapply(1:ncol(fData_table),
                             function(i)
                               length(grep("^GAPDH(S)*$|^ICAM1$|^Gapdh(s)*$",fData_table[,i]))>0));
  if (length(symbol_pos) < 1){
    symbol_pos <- which(sapply(1:ncol(fData_table),
                               function(i)#bug fix-> create pattern match for correctly mapping REFSEQ
                                 length(grep("^NM_002046.*$|^NM_198317.*$",fData_table[,i]))>0))
    if(length(symbol_pos)> 1){symbol_pos <- symbol_pos[1]}
    if(length(symbol_pos)==1){

      #convert genes into gene names by matching accession numbers
      #bug fix -> removed versioning so that it is now possible to match stuff
      ucsc_table  <- MetaIntegrator:::.ucsc_refseq_table()

      #if REFSEQ IDs have _at at the end add it to the table
      if(length(grep("_at",fData_table[1,symbol_pos]))>0){
        ucsc_table$name <- gsub("$","_at",ucsc_table$name)
      }

      probe2table <- match(gsub("\\..$","",fData_table[,symbol_pos]),
                           ucsc_table$name)

      #store data
      keys        <- ucsc_table$name2[probe2table]
      names(keys) <- as.character(fData_table[[1]])

      comment     <- "Human REFSEQ format"
    }
  }
  keys <- as.character(fData_table[,symbol_pos])
  names(keys) <- as.character(fData_table[[1]])
  return(keys)
}

# ---- parse genes non-human ---- #
