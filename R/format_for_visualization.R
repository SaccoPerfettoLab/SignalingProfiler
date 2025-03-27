#' Format Network for Cytoscape Visualization
#'
#' This function prepares a **CARNIVAL-optimized** network for visualization in **Cytoscape**
#' by formatting node and edge attributes, including **phosphorylation site labeling, fold-change direction,**
#' and **highlighting significant interactions**.
#'
#' This function should be used after [expand_and_map_edges()], which integrates phosphoproteomics data into the network.
#'
#' @param sp_object A list containing:
#'   \itemize{
#'     \item `igraph_network` - An igraph object representing the expanded CARNIVAL network.
#'     \item `nodes_df` - A dataframe with node attributes, including regulatory activity scores.
#'     \item `edges_df` - A dataframe with edge attributes, including phosphorylation site details.
#'   }
#'
#' @return A modified list (*SP_object*) with attributes optimized for **Cytoscape visualization**:
#' \item{igraph_network}{An igraph object containing formatted node and edge attributes.}
#' \item{nodes_df}{A dataframe of nodes with **highlighting for phosphorylation-regulated proteins**.}
#' \item{edges_df}{A dataframe of edges with attributes for Cytoscape visualization, including:
#'   \itemize{
#'     \item `label` - The **one-letter phosphosite notation** (S, T, Y) for visualization.
#'     \item `FC_sign` - The **direction of fold-change** (e.g., `"up"` for activation, `"down"` for inhibition).
#'     \item `highlight` - `"Y"` for significant phosphorylation events, `"N"` otherwise.
#'   }
#' }
#'
#' @details
#' - The function converts **phosphorylation site amino acid labels** (`Ser`, `Thr`, `Tyr`)
#'   into **one-letter notation** (`S`, `T`, `Y`) for easier visualization in Cytoscape.
#' - Significantly modulated phosphosites are labeled with a `*` (e.g., `"S326*"` for significant Serine-326).
#' - **Nodes associated with significant phosphosites** are flagged for easy identification.
#' - This function should be applied **after** [expand_and_map_edges()] to ensure the integration of phosphoproteomics data.
#'
#' @seealso [expand_and_map_edges()], [run_carnival_and_create_graph()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load an expanded CARNIVAL network with phosphoproteomics integration
#' data('toy_opt_network')
#' cytoscape_ready_network <- format_for_visualization(toy_opt_network)
#' }
#'
format_for_visualization <- function(sp_object) {

  # Ensure required columns exist
  if (!"aminoacid" %in% colnames(sp_object$edges_df)) {
    stop("Error: `aminoacid` column is missing in edges_df.")
  }

  edges_table <- sp_object$edges_df

  # Convert amino acids to one-letter notation
  edges_table <- edges_table %>%
    dplyr::mutate(
      label = stringr::str_replace_all(aminoacid, c("Ser" = "S", "Tyr" = "Y", "Thr" = "T"))
    )

  # Remove asterisks (*) from labels and handle multiple labels correctly
  edges_table_new <- edges_table %>%
    tidyr::separate_rows(label, sep = ";") %>%
    dplyr::mutate(label = ifelse(grepl("\\*", label), label, "")) %>%
    dplyr::group_by(across(-label)) %>%
    dplyr::summarize(label = paste(unique(label), collapse = ";"), .groups = "drop")

  # Add flag for fold-change direction
  edges_table_new <- edges_table_new %>%
    dplyr::mutate(
      label = stringr::str_remove_all(label, "\\*"),
      FC_sign = ifelse(grepl("-", FC), "down", "up")
    )

  # Create an igraph object with the updated network
  SP_graph <- igraph::graph_from_data_frame(edges_table_new, vertices = sp_object$nodes_df, directed = TRUE)

  # Return the updated structured object
  return(list(
    igraph_network = SP_graph,
    nodes_df = sp_object$nodes_df,
    edges_df = edges_table_new
  ))
}
