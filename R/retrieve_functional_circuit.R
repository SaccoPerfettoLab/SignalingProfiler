#' Extract QueryNodes-to-Phenotype Circuit from SignalingProfiler Network
#'
#' This function extracts a subnetwork linking specified start nodes
#' to phenotypes within a given path length (`k`). It generates a
#' directed subgraph where only paths from phenotypes to start nodes
#' are retained.
#'
#' @param SP_object A list containing the SignalingProfiler network (igraph object), nodes table, and edges table.
#' @param start_nodes A character vector of nodes representing the start points in the model.
#' @param phenotypes A character vector of phenotype names to be linked to the start nodes.
#' @param k An integer specifying the maximum path length between phenotypes and start nodes.
#' @param start_to_top Logical, whether to remove incoming edges to start nodes (`TRUE` or `FALSE`, default: `FALSE`).
#'
#' @return An `igraph` object representing the extracted subnetwork with paths from phenotypes to start nodes.
#'
#' @details
#' The function:
#' - Identifies shortest paths from each phenotype to the defined start nodes.
#' - Retains paths that are at most `k` steps long.
#' - Extracts a subgraph containing the identified paths.
#' - (Optional) Removes incoming edges to start nodes if `start_to_top = TRUE`.
#'
#' If no valid paths are found, a warning is issued.
#' @export
#'
#' @examples
#' data('toy_sp_output')
#' pheno_to_start_circuit(SP_object = toy_sp_output, start_nodes = 'FLT3', phenotypes = c('PROLIFERATION', 'APOPTOSIS'), k = 4)
#'
pheno_to_start_circuit <- function(SP_object, start_nodes, phenotypes, k, start_to_top = FALSE) {

  SP_graph <- SP_object$igraph_network

  final_nodes <- c()
  for (i in c(1:length(phenotypes))) {
    phenotype <- phenotypes[i]

    # Convert phenotype name to a valid format for matching
    formatted_phenotype <- stringr::str_replace_all(phenotype, "[[:space:]\\\\/]", "_")

    # Identify start nodes present in the graph
    valid_start_nodes <- igraph::V(SP_graph)$name[igraph::V(SP_graph)$name %in% start_nodes]

    # Extract paths from phenotype to start nodes
    all_paths <- igraph::all_simple_paths(
      graph = SP_graph,
      from = formatted_phenotype,
      to = valid_start_nodes,
      mode = "in",
      cutoff = k
    )

    # Extract unique nodes in the identified paths
    all_paths_nodes <- unique(names(unlist(all_paths)))

    # Warning if no paths are found for the given phenotype
    if (length(all_paths_nodes) == 0) {
      warning(paste0("No path of length ", k, " have been found for ",
                     phenotype))
    }

    # Store the identified nodes
    final_nodes <- c(final_nodes, all_paths_nodes)
  }

  # Extract subgraph containing the selected nodes
  pheno_circuit <- igraph::induced_subgraph(SP_graph, unique(final_nodes))

  # Extract edges from the subgraph
  pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                               what = c("edges"))

  # Remove incoming edges to start nodes if specified
  if(start_to_top == TRUE){
    pheno_circuit_edges <- pheno_circuit_edges %>%
      dplyr::filter(!to %in% start_nodes)
  }

  # Create final subgraph with adjusted edges
  pheno_circuit_no_incoming <-igraph::graph_from_data_frame(
    d = pheno_circuit_edges,
    vertices = igraph::as_data_frame(pheno_circuit,
                                     what = c("vertices")))

  return(pheno_circuit_no_incoming)
}
