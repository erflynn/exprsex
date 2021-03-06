% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplMapDetectGeneID.R
\name{.detect_refseq_cols}
\alias{.detect_refseq_cols}
\title{Identify the refseq columns and return the one with the most overlap.
(Note - this looks for "NM_" ids, which is not comprehensive.)}
\usage{
.detect_refseq_cols(gpl.df, MIN.OVERLAP = 8000)
}
\arguments{
\item{gpl.df}{the gpl data frame}

\item{MIN.OVERLAP}{the minimum number of mapped genes allowed, default=8000}
}
\value{
the index of the column containing these ids
}
\description{
Identify the refseq columns and return the one with the most overlap.
(Note - this looks for "NM_" ids, which is not comprehensive.)
}
