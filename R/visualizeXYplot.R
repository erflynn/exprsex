
#'This function graphs the XIST vs. RPS4Y1 probes for a given study.
#'
#' @param expr_matrix expression data matrix - rows are genes, columns are samples
#' @param keys mapping between genes and probes
#' @param sex_lab sex labels (0 - female, 1 - male) for each of the samples
#'
#' @return a graph of RPS4Y1 vs. XIST probes corresponding to the sex labels of studies

visualizeXYplot <- function (expr_matrix, keys, sex_lab, x_axis = "XIST", y_axis = "RPS4Y1") {
  require('dplyr')

  keys <- as.data.frame(keys)
  keys$probe <- rownames(keys)
  names(keys) <- c("gene", "probe")

  if ('XIST' %in% keys$gene) {
    xist.probes <- filter(keys, gene == "XIST")
    xist.vals <- as.data.frame(expr_matrix[xist.probes$probe, ])
    xist.vals$samp_ids <- names(expr_matrix[xist.probes$probe, ])
    colnames(xist.vals)[ncol(xist.vals)] <- "samp_ids"
  }

  else {
    return("Error: XIST gene missing.")
  }

  #   obtain RPS4Y1 values from expression data
  if ('RPS4Y1' %in% keys$gene) {
    rps4y1.probes <- filter(keys, gene == "RPS4Y1")
    rps4y1.vals <- as.data.frame(expr_matrix[rps4y1.probes$probe, ])
    rps4y1.vals$samp_ids <- names(expr_matrix[rps4y1.probes$probe, ])
    colnames(rps4y1.vals) <- c(rps4y1.probes$probe, "samp_ids")
  }

    else {
    return("Error: RPS4Y1 gene missing.")
  }

  probe_map <- unlist(c(rep("XIST", length(xist.probes)), rep("RPS4Y1", length(rps4y1.probes))))
  names(probe_map) <- unlist(c(xist.probes, rps4y1.probes))
  #   plot the ranked XIST and RPS4Y1 data
  rank_data <- convertToGenes(expr_matrix, probe_map, c("XIST", "RPS4Y1"))
  xy_vals <- rank_data[c('XIST', 'RPS4Y1'),]
  sex_lab <- as.data.frame(sex_lab)
  names(sex_lab) <- c("sex_lab")
  my.colors <- sapply(sex_lab$sex_lab, function(x) { ifelse(x == 1, "blue", "red")})
  plot(t(xy_vals), xlab = 'XIST', ylab = "RPS4Y1", col=my.colors)
}
