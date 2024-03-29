% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carnival_analysis.R
\name{union_of_graphs}
\alias{union_of_graphs}
\title{union_of_graphs}
\usage{
union_of_graphs(
  graph_1,
  graph_2,
  proteins_df,
  files,
  path_sif = "./union_graph.sif",
  path_rds = "./union_graph.rds"
)
}
\arguments{
\item{graph_1}{run1 igraph object of run_carnival_and_create_graph}

\item{graph_2}{run2 igraph object of run_carnival_and_create_graph}

\item{proteins_df}{tibble of inferred proteins in the naive network}

\item{files}{boolean value, TRUE if you want output files}

\item{path_rds}{}
}
\description{
union_of_graphs
}
