#' Split and Sort Elements in a String
#'
#' This function splits a string by underscores (`_`), sorts the elements
#' alphabetically, removes duplicates, and then concatenates them back
#' using an underscore separator.
#'
#' @param x A character vector where each element consists of substrings
#'   separated by underscores (`_`).
#'
#' @return A character vector where each element is sorted and unique,
#'   joined by underscores (`_`).
split_and_sort <- function(x) {
  x <- strsplit(x, "_")
  lapply(x, function(y) paste(sort(unique(y)), collapse = "_"))
}

#' Extract Nodes from an iGraph Object
#'
#' This function extracts node (vertex) information from an `igraph` network
#' and formats it into a tibble with relevant attributes.
#'
#' @param graph An `igraph` object representing a network.
#'
#' @return A tibble with columns:
#'   \item{gene_name}{Node names in the network (renamed from `name`).}
#'   \item{UNIPROT}{UniProt accession IDs (if available in the graph).}
#'   \item{carnival_activity}{CARNIVAL activity scores (if available in the graph).}
#'
#' @export
#'
#' @examples
#' library(igraph)
#' data('toy_opt_network')
#' extract_carnival_nodes(toy_opt_network$igraph_network)
#' 
extract_carnival_nodes <- function(graph) {
  igraph::as_data_frame(graph, what = "vertices") %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_name = name) %>%
    dplyr::select(gene_name, UNIPROT, carnival_activity)
}

#' Extract Edges from an iGraph Object
#'
#' This function extracts edge information from an `igraph` network
#' and returns it as a tibble.
#'
#' @param graph An `igraph` object representing a network.
#'
#' @return A tibble where each row represents an edge, containing at least
#'   the `source` and `target` nodes.
#'
#' @export
#'
#' @examples
#' library(igraph)
#' data('toy_opt_network')
#' extract_carnival_edges(toy_opt_network$igraph_network)
extract_carnival_edges <- function(graph) {
  igraph::as_data_frame(graph, what = "edges") %>%
    tibble::as_tibble()
}
