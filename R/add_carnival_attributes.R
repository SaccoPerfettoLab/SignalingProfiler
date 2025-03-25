
#' Add Attributes to Edges from CARNIVAL Optimization Output
#'
#' This function extracts and annotates edge attributes from the output of
#' a CARNIVAL optimization run. It creates a data frame containing source and
#' target nodes, the sign of the interaction, and the corresponding CARNIVAL weight.
#'
#' @param carnival_result A list returned from `run_carnival_and_create_graph()`,
#' typically including the `weightedSIF` component with columns \code{Node1},
#' \code{Node2}, \code{Sign}, and \code{Weight}.
#'
#' @return A \code{tibble} with the following columns:
#' \describe{
#'   \item{source}{The Gene Symbol of source node of the interaction.}
#'   \item{target}{The Gene Symbol of target node of the interaction.}
#'   \item{sign}{The sign of the interaction (1 for activation, -1 for inhibition).}
#'   \item{carnival_weight}{The weight assigned by CARNIVAL to the edge.}
#' }
#'
#' @examples
#' \dontrun{
#' edges_df <- add_output_carnival_edges_attributes(carnival_result)
#' }
#'
add_output_carnival_edges_attributes <- function(carnival_result){
  optimal_edges <- tibble::as_tibble(carnival_result$weightedSIF)

  edges_df <- tibble::tibble(
    source = optimal_edges$Node1,
    target = optimal_edges$Node2,
    sign = optimal_edges$Sign,
    carnival_weight = as.numeric(optimal_edges$Weight))
  return(edges_df)
}

#' Add Attributes to Nodes from CARNIVAL Optimization Output
#'
#' This function processes the nodes inferred by CARNIVAL and enriches them
#' with molecular attributes such as UNIPROT identifiers, molecular function scores,
#' and final activity scores. It also computes a `discordant` flag indicating
#' disagreement between inferred activity and experimental evidence.
#'
#' Depending on the selected organism (human or mouse) and the source of prior knowledge
#' (e.g., with or without kinome atlas), it joins CARNIVAL output with relevant
#' protein databases and experimental data.
#'
#' @param carnival_result A list returned from `run_carnival_and_create_graph()`,
#' \code{$nodesAttributes} and \code{$weightedSIF}.
#'
#' @param proteins_df A data frame containing inferred protein information, including
#' columns such as \code{gene_name}, \code{UNIPROT}, \code{mf} (molecular function),
#' \code{final_score}, and \code{method}.
#' @param organism A character string, either \code{"mouse"} or \code{"human"},
#' indicating the species for which the prior knowledge network (PKN) should be used.
#' @param with_atlas Logical, default \code{TRUE}. If \code{FALSE}, kinome atlas
#' interactions are excluded.
#' @param direct Logical, default \code{FALSE}. If \code{TRUE}, uses only direct
#' interactions in the network.
#'
#' @return A \code{tibble} containing nodes and enriched attributes, including:
#' \describe{
#'   \item{gene_name}{Gene symbol of the node.}
#'   \item{carnival_activity}{Inferred activity score from CARNIVAL.}
#'   \item{UNIPROT}{UNIPROT identifier.}
#'   \item{mf}{Molecular function annotation.}
#'   \item{final_score}{Score derived from experimental data in Step 1.}
#'   \item{method}{Method used for protein inference.}
#'   \item{discordant}{Logical flag indicating inconsistency between
#'   inferred activity and final score.}
#' }
#'
#'
add_output_carnival_nodes_attributes <- function(carnival_result,
                                                 proteins_df,
                                                 organism,
                                                 with_atlas = TRUE,
                                                 direct = FALSE){

  # Get correct database
  if(organism == 'mouse'){
    if(with_atlas == TRUE){
      stop('If organism is \'mouse\' with_atlas parameter must be FALSE')
    }else{
      PKN_proteins <- SignalingProfiler::PKN_proteins_mouse
      if(direct == TRUE){
        db <- SignalingProfiler::db_mouse_dir
      }else{
        db <- SignalingProfiler::db_mouse_ind
      }
    }
  }else if(organism == 'human'){
    if(with_atlas == TRUE){
      PKN_proteins <- SignalingProfiler::PKN_proteins_human_atlas

      if(direct == TRUE){
        db <- SignalingProfiler::db_human_atlas_dir
      }else{
        db <- SignalingProfiler::db_human_atlas_ind
      }

    }else{
      PKN_proteins <- SignalingProfiler::PKN_proteins_human

      if(direct == TRUE){
        db <- SignalingProfiler::db_human_dir
      }else{
        db <- SignalingProfiler::db_human_ind
      }
    }
  }else{
    error('Please provide a valide organism')
  }

  nodes <- tibble::as_tibble(carnival_result$nodesAttributes)
  optimal_edges <- tibble::as_tibble(carnival_result$weightedSIF)

  nodes_igraph_ids <- unique(c(optimal_edges$Node1, optimal_edges$Node2))

  optimal_nodes <- nodes %>%
    dplyr::filter(Node %in% nodes_igraph_ids) %>%
    dplyr::select(Node, 'carnival_activity' = AvgAct) %>%
    dplyr::mutate_at('carnival_activity', as.numeric)

  # Take UNIPROT ID from Prior Knowledge Network
  nodes_df <- dplyr::left_join(optimal_nodes,
                               PKN_proteins,
                               by = c('Node' = 'ENTITY')) %>%
    dplyr::rename('gene_name' = 'Node',
                  'UNIPROT' = 'ID')

  nodes_df <- dplyr::left_join(
    nodes_df,
    proteins_df %>%
      dplyr::select(-UNIPROT) %>%
      dplyr::mutate(gene_name = stringr::str_to_upper(stringr::str_replace_all(gene_name, "[^[:alnum:]]", '_'))),
    by = c('gene_name')) %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score, method)

  # Add New Attributes
  nodes_df$discordant <- FALSE
  nodes_df$discordant[as.numeric(nodes_df$carnival_activity) * nodes_df$final_score < 0] <- TRUE

  nodes_df <- nodes_df %>%
    dplyr::relocate(gene_name) %>%
    dplyr::distinct()

  return(nodes_df)
}
