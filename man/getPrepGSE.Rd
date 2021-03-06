% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepInput.R
\name{getPrepGSE}
\alias{getPrepGSE}
\title{Get and prepare a gse file for sex labeling by converting it to genes.
Optionally this will also reorder by a list of genes and convert it to ranks.}
\usage{
getPrepGSE(gse, gse.dir = NULL, gpl.dir = NULL, out.dir = NULL)
}
\arguments{
\item{gse}{the ID of the GSE to download}

\item{gse.dir}{directory that contains gse files, if empty, will download the GSE}

\item{gpl.dir}{directory where GPL ref data is location, empty defaults to tempdir()}

\item{out.dir}{directory to write output, if empty will just be returned and not be written out}
}
\value{
list of gene expression matrices, with rows as genes and columns as samples
        values are gene expression levels
}
\description{
Get and prepare a gse file for sex labeling by converting it to genes.
Optionally this will also reorder by a list of genes and convert it to ranks.
}
