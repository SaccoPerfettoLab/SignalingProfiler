% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carnival_analysis.R
\name{run_carnival_and_create_graph}
\alias{run_carnival_and_create_graph}
\title{run_carnival_and_create_graph}
\usage{
run_carnival_and_create_graph(
  source_df,
  target_df,
  naive_network,
  carnival_options,
  proteins_df,
  organism,
  topbottom = FALSE,
  with_atlas = TRUE,
  direct = FALSE,
  files = TRUE,
  path_sif = "./optimized_network.sif",
  path_rds = "./optimized_SP_oject.RDS"
)
}
\arguments{
\item{source_df}{if NULL, inverse CARNIVAL,
otherwise, tibble with source nodes discretized among 1 and -1}

\item{target_df}{tibble with target nodes in a continuous range of activity}

\item{naive_network}{tibble with naive network in SIF format}

\item{carnival_options}{list of options returned by default_CARNIVAL_options}

\item{proteins_df}{all inferred proteins df}

\item{organism}{string, organism you are using}

\item{topbottom}{Boolean value, default FALSE, if TRUE optimization priorities sources}

\item{with_atlas}{Boolean value, default TRUE, if FALSE excludes Kinome Altas derived regulons}

\item{direct}{Boolean value, default FALSE, if TRUE uses only direct interactions}

\item{files}{boolean value, TRUE if you want output files}

\item{path_sif}{path of the sif output file of network}

\item{path_rds}{path of the rds output file of network}
}
\value{
list with igraph object, nodes df and edges df
}
\description{
run_carnival_and_create_graph
}
