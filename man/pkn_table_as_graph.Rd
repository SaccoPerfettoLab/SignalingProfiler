% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shortest_paths_custom.R
\name{pkn_table_as_graph}
\alias{pkn_table_as_graph}
\title{Convert a Prior Knowledge Network Table to an igraph Object}
\usage{
pkn_table_as_graph(PKN_table)
}
\arguments{
\item{PKN_table}{Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.}
}
\value{
An `igraph` object representing the Prior Knowledge Network.
}
\description{
Transforms a tabular Prior Knowledge Network (`PKN_table`) into an
`igraph` object, creating nodes and edges with SIGNOR-compliant attributes
}
