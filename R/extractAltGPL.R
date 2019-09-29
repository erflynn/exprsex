
# <---- DEAL WITH GPL PARSING ERRORS -----> #

#' Alternate GPL parsing
#'
.parseGenesAlt <- function(fData_table){
  # adapted from MetaIntegrator
  symbol_pos <- which(sapply(1:ncol(fData_table),
                             function(i)
                               length(grep("^GAPDH(S)*$|^ICAM1$|^Gapdh(s)*$",fData_table[,i]))>0));
  keys <- as.character(fData_table[,symbol_pos])
  names(keys) <- as.character(fData_table[[1]])
  return(keys)
}
