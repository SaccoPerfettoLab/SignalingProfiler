% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naive_network_construction.R
\name{two_layer_naive_network}
\alias{two_layer_naive_network}
\title{two_layer_naive_network}
\usage{
two_layer_naive_network(
  starts_gn,
  intermediate_gn,
  targets_gn,
  PKN_table,
  max_length_1,
  max_length_2,
  rds_path = "two_layer_naive.RDS",
  sif_path = "two_layer_naive.sif",
  connect_all = FALSE
)
}
\arguments{
\item{starts_gn}{gene names of starting nodes}

\item{intermediate_gn}{gene names of intermediate nodes (e.g. kins and phos)}

\item{targets_gn}{gene names of target nodes (e.g. transcription factors)}

\item{PKN_table}{tibble of all causal interactions}

\item{max_length_1}{max_length of shortest path from start to intermediates}

\item{max_length_2}{max_length of shortest path from intermediates to targets}

\item{rds_path}{path of network rds file}

\item{sif_path}{path of network sif file}

\item{connect_all}{Boolean, if TRUE connect intermediate nodes}
}
\value{
naive network
}
\description{
two_layer_naive_network
}
