#' Convert CARNIVAL Output Nodes into the Next Input Format
#'
#' This function extracts and formats node attributes from a CARNIVAL optimization
#' result (`run_carnival_and_create_graph()` output) to be used as input
#' for subsequent analyses with CARNIVAL.
#'
#' @param carnival_result A list output from `run_carnival_and_create_graph()`,
#'        containing optimized network information.
#'
#' @return A tibble containing the formatted node attributes with:
#' \itemize{
#'   \item `UNIPROT` - UniProt ID of the node.
#'   \item `gene_name` - Gene symbol of the node.
#'   \item `final_score` - The inferred activity score of the node (CARNIVAL activity).
#'   \item `mf` - Molecular function category (e.g., TF, kinase).
#'   \item `method` - Method used for inference (e.g., `VIPER`, `CARNIVAL`).
#' }
#'
#' @export
#'
#' @examples
#' # Example CARNIVAL output
#' data('toy_opt_network')
#' # Convert nodes into the next input format
#' formatted_nodes <- convert_output_nodes_in_next_input(toy_opt_network)
#' print(formatted_nodes)
convert_output_nodes_in_next_input <- function(carnival_result){
  if (!"nodes_df" %in% names(carnival_result)) {
    stop("Invalid input: `carnival_result` must contain a `nodes_df` element.")
  }

  nodes <- carnival_result$nodes_df

  formatted_nodes <- nodes %>%
    dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf, method)

  return(formatted_nodes)
}


#' Merge Two CARNIVAL Graphs into a Unified Network
#'
#' This function takes two CARNIVAL-generated graphs, merges their nodes and edges,
#' and computes the union of the networks while resolving duplicate nodes and edges.
#'
#' @param graph_1 An `igraph` object from `run_carnival_and_create_graph()`, representing the first network.
#' @param graph_2 An `igraph` object from `run_carnival_and_create_graph()`, representing the second network.
#' @param proteins_df A tibble containing inferred proteins from the naive network.
#' @param files Logical, if `TRUE`, saves the resulting network as a `.sif` file and an `.rds` object.
#' @param path_sif Character, the file path for saving the network in `.sif` format. Default is `'./union_graph.sif'`.
#' @param path_rds Character, the file path for saving the network as an `.rds` object. Default is `'./union_graph.rds'`.
#'
#' @return *SP_object*, a list  containing:
#' \itemize{
#'   \item `igraph_network` - The merged `igraph` network.
#'   \item `nodes_df` - A tibble with node attributes.
#'   \item `edges_df` - A tibble with edge attributes.
#' }
#'
#' @export
#'
#' @examples
#' # Example CARNIVAL outputs (simplified)
#' graph_1 <- toy_opt_network$igraph_network
#' graph_2 <-  toy_opt_network$igraph_network
#'
#' proteins_df <- tibble::tibble(
#'   gene_name = c("A", "B", "C", "D", "E", "F", "G"),
#'   UNIPROT = c("P12345", "P67890", "Q12345", "Q67890", "R12345", "R67890", "S12345"),
#'   mf = c("tf", "kin", "tf", "kin", "tf", "kin", "tf"),
#'   final_score = c(1.1, -0.4, 0.9, -0.6, 0.8, -1.0, 0.3),
#'   method = c("CARNIVAL", "CARNIVAL", "CARNIVAL", "CARNIVAL", "CARNIVAL", "CARNIVAL", "CARNIVAL")
#' )
#'
#' # Merge graphs
#' merged_network <- union_of_graphs(graph_1, graph_2, proteins_df, files = FALSE)
union_of_graphs <- function(graph_1,
                            graph_2,
                            proteins_df,
                            files,
                            path_sif = './union_graph.sif',
                            path_rds = './union_graph.rds'){

  # Ensure the proteins dataframe does not contain redundant UNIPROT column
  if("UNIPROT" %in% colnames(proteins_df)) {
    proteins_df <- proteins_df %>% dplyr::select(-UNIPROT)
  }

  nodes_rec_kin <- extract_carnival_nodes(graph_1)
  nodes_kin_tf <- extract_carnival_nodes(graph_2)

  # Combine node attributes and remove duplicates
  nodes <- dplyr::bind_rows(nodes_rec_kin, nodes_kin_tf) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene_name)

  # Merge with protein metadata
  nodes <- dplyr::left_join(nodes, proteins_df, by = "gene_name") %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score, method) %>%
    dplyr::relocate(gene_name)

  # Resolve duplicate nodes by averaging activity scores
  nodes_avg_activity <- nodes %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(carnival_activity = mean(carnival_activity, na.rm = TRUE))

  nodes <- dplyr::distinct(nodes, gene_name, UNIPROT, mf, final_score, method)
  nodes <- dplyr::left_join(nodes_avg_activity, nodes, by = "gene_name") %>%
    dplyr::relocate(gene_name)

  # Flag discordant nodes where CARNIVAL activity conflicts with final score
  nodes <- nodes %>%
    dplyr::mutate(discordant = carnival_activity * final_score < 0)

  # Extract edges from both graphs and combine them
  edges <- dplyr::bind_rows(extract_carnival_edges(graph_1), extract_carnival_edges(graph_2))

  # Create unified igraph object
  union_graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

  # Package output
  sp_object <- list(igraph_network = union_graph,
                    nodes_df = nodes,
                    edges_df = edges)

  # Save output if required
  if (files) {
    igraphToSif(union_graph, path_sif, 'sign')
    saveRDS(sp_object, path_rds)
  }

  return(sp_object)
}
