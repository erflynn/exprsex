% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplMap.R
\name{.map_from_refseq}
\alias{.map_from_refseq}
\title{Helper function to extract and map refseq to entrez}
\usage{
.map_from_refseq(gpl.df, my.col, organism, ref_dir = NULL)
}
\arguments{
\item{gpl.df}{the feature data from the GPL data frame}

\item{my.col}{the column containing the refseq data}

\item{organism}{the organism: rat, mouse, or human}

\item{ref_dir}{option to set reference directory  (defaults to tempdir)}
}
\value{
a parsed probe/gene data frame
}
\description{
Helper function to extract and map refseq to entrez
}
