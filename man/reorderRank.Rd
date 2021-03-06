% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepInput.R
\name{reorderRank}
\alias{reorderRank}
\title{Reorder and rank a gene expression matrix. We then use this for sex labeling and comparison.}
\usage{
reorderRank(gene_mat, gene_list = list_genes, to.ranks = TRUE)
}
\arguments{
\item{gene_mat}{the expression matrix to reorder and rank}

\item{gene_list}{the list of genes to use, this is helpful if you want to compare across multiple studies}

\item{to.ranks}{whether to convert to ranks, default is TRUE}
}
\value{
expression matrix reordered based on the gene list, with expression converted to ranks
}
\description{
Reorder and rank a gene expression matrix. We then use this for sex labeling and comparison.
}
