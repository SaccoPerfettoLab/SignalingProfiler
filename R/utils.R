#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Convert an igraph Object to a SIF (Simple Interaction Format) File
#'
#' This function exports an `igraph` object into a SIF (Simple Interaction Format) file,
#' commonly used for representing network interactions.
#'
#' @param inGraph An `igraph` object representing the network.
#' @param outfile A string specifying the output file path for the SIF file; default: `"output.sif"`
#' @param edgeLabel A string specifying the edge attribute to use as an interaction label in the SIF file; default: `"label"`
#'
#' @return `NULL`. The function creates a SIF file named as defined by `outfile`.
#'
#' @details
#' The SIF format represents network interactions using a three-column structure:
#' \tabular{lll}{
#'   **Node1** \tab **Interaction Type** \tab **Node2** \cr
#' }
#' Each row represents an interaction between two nodes, with the interaction type
#' taken from the specified edge attribute (`edgeLabel`).
#' If a file with the same name as `outfile` already exists, it will be overwritten.
#'
#' @export
#'
#' @examples
#'
#' library(igraph)
#' g <- graph_from_edgelist(matrix(c("A", "B", "B", "C"), byrow=TRUE, ncol=2))
#' E(g)$label <- c(1, -1)
#' igraphToSif(g, outfile="network.sif", edgeLabel="label")
#'
igraphToSif <- function(inGraph, outfile="output.sif", edgeLabel="label") {

  if(file.exists(outfile)){file.remove(outfile)}

  file_conn <- file(outfile, open="wt")

  edges <- igraph::as_edgelist(inGraph)
  attributes <- igraph::edge_attr(inGraph, edgeLabel)

  for (i in 1:nrow(edges)) {
    edge <- edges[i,]
    attribute <- attributes[i]
    cat(paste0(edge[1], '\t', attribute, '\t', edge[2],  '\n'), file = file_conn)
  }

  close(file_conn)
}
