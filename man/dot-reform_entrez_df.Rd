% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplMap.R
\name{.reform_entrez_df}
\alias{.reform_entrez_df}
\title{Helper function to reformat entrez gene/probe data frame.
Expands out the gene field so that there is one gene per row.}
\usage{
.reform_entrez_df(df, entrez.col.name = "gene")
}
\arguments{
\item{df}{the data frame to expand}

\item{entrez.col.name}{the name of the column containing entrez ids, defaults to "gene"}
}
\value{
data reformatted with one entrez id per row
}
\description{
Helper function to reformat entrez gene/probe data frame.
Expands out the gene field so that there is one gene per row.
}
