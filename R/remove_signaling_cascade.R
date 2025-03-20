#' Remove Signaling Cascade Regulators
#'
#' This function filters out proteins involved in signaling cascades
#' that regulate the same phenotype multiple times, ensuring only
#' direct regulatory relationships remain.
#'
#' @param input_prot_paths A dataframe containing protein-phenotypes paths,
#'   including regulatory effects (`Effect`), query proteins (`QueryNode`),
#'   and target phenotypes (`EndPathways`).
#' @param sp_graph An `igraph` object representing the signaling network
#'   from SignalingProfiler.
#'
#' @return A filtered dataframe containing only paths of
#'  proteins that are not part of redundant signaling cascades.
#'
#' @details
#' The function identifies proteins that regulate the same phenotype multiple times
#' and removes those that are downstream in the signaling cascade.
#' It evaluates the shortest path distance between proteins using the provided
#' `sp_graph`, considering both forward (A → B) and reverse (B → A) connections.
#'
#' If no redundant cascade regulators are found, the input dataframe
#' is returned unchanged.
#'
#' @export
#'
#' @examples
#' library(igraph)
#' library(dplyr)
#'
#' # Example input data
#' input_data <- data.frame(
#'   QueryNode = c("TP53", "AKT1", "MAPK1"),
#'   EndPathways = c("APOPTOSIS", "APOPTOSIS", "APOPTOSIS"),
#'   Effect = c("up-regulates", "up-regulates", "down-regulates")
#' )
#'
#' # Example graph
#' example_graph <- graph_from_edgelist(
#'   matrix(c("TP53", "AKT1", "AKT1", "MAPK1"), ncol = 2, byrow = TRUE),
#'   directed = TRUE
#' )
#'
#' # Run cascade removal
#' filtered_paths <- remove_signaling_cascade(input_data, example_graph)
#' print(filtered_paths)
#'
remove_signaling_cascade <- function(input_prot_paths, sp_graph){

  # filter phenotypes with effect multiply regulated
  INPUT.Phen.Paths_clean <- input_prot_paths %>%
    dplyr::filter(Effect != '-')

  INPUT.Phen.Paths_clean <- INPUT.Phen.Paths_clean %>%
    dplyr::mutate(key = ifelse(INPUT.Phen.Paths_clean$Effect == 'up-regulates',
                               paste0(INPUT.Phen.Paths_clean$EndPathways, '-up'),
                               paste0(INPUT.Phen.Paths_clean$EndPathways, '-down')))

  # Keep only true multiple regulators
  prot_to_phenotypes <- INPUT.Phen.Paths_clean %>%
    dplyr::select(key, QueryNode) %>%
    dplyr::distinct()

  multiple_regulated_phen <- prot_to_phenotypes %>%
    dplyr::count(key) %>%
    dplyr::filter(n>1)

  # Remove cascade only if multiple_regulated_phen exist
  if(nrow(multiple_regulated_phen) != 0){

    # Create a table with proteins on phenotypes multiply regulated
    prot_to_phenotypes_multiple <- prot_to_phenotypes %>%
      dplyr::filter(key %in% multiple_regulated_phen$key)

    # Select phenotypes
    phenotypes <- unique(prot_to_phenotypes_multiple$key)

    # Initialize empty list
    to_remove_list <- list()

    for(i_phen in c(1:length(phenotypes))){

      phenotype = phenotypes[i_phen]

      prot_to_phenotype <- prot_to_phenotypes_multiple %>%
        dplyr::filter(key == phenotype)

      combinations <- combn(unique(prot_to_phenotype$QueryNode), 2)
      # create a matrix for the distance and the flag
      dist_count <- matrix(nrow = 3, ncol = ncol(combinations))

      # forward run: check direction A --> B
      for(i in c(1:ncol(combinations))){

        # Handling NA in gene_name column
        true_table_from <- igraph::V(sp_graph)$name == combinations[1,i]
        true_table_to <- igraph::V(sp_graph)$name == combinations[2,i]

        true_table_from[is.na(true_table_from)] <- FALSE
        true_table_to[is.na(true_table_to)] <- FALSE

        if(sum(true_table_from)==0 |  sum(true_table_to) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network.')
          dist_count[1,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from]
          to = igraph::V(sp_graph)[true_table_to]
          # Compute the distance between two nodes
          # assign to the first row of the matrix the distance
          dist_count[1,i] <- unlist(igraph::distances(sp_graph, v = from,
                                                      to = to,
                                                      mode = "out"))

          # If the distance is infinite, NA node to remove
          if(dist_count[1,i] == Inf){
            dist_count[2,i] <- NA
          }else{ # If the nodes are connected, the second is more downstream and should be removed
            dist_count[2,i] <- combinations[2,i] # Opposite
          }
        }
      }
      # Reverse run: check the direction B --> A
      for(i in c(1:ncol(combinations))){
        true_table_from_rev <- igraph::V(sp_graph)$name == combinations[2,i]
        true_table_from_rev[is.na(true_table_from_rev)] <- FALSE

        true_table_to_rev <- igraph::V(sp_graph)$name == combinations[1,i]
        true_table_to_rev[is.na(true_table_to_rev)] <- FALSE

        if(sum(true_table_from_rev) == 0 | sum(true_table_to_rev) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network')
          dist_count[3,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from_rev]
          to = igraph::V(sp_graph)[true_table_to_rev]

          # Override the previous distance
          distance_rev <- unlist(igraph::distances(sp_graph, from, to = to, mode = "out"))

          # If it is infinite, not connected no removal
          if(distance_rev == Inf){
            dist_count[3,i] <- NA
          }else{ # If it is finite, remove the first node because it is more downstream
            if(distance_rev < dist_count[1,i]){
              # But if the new distance is shorter than the forward one, remove this
              dist_count[3,i] <- combinations[1,i] #opp
            }else{
              dist_count[3,i] <- combinations[2,i] #opp
            }
          }
        }
      }
      output <- dist_count[!(is.na(as.vector(dist_count)) | as.vector(dist_count) == Inf)]

      # if i binary remove all the proteins
      to_remove_list[[i_phen]] <- output
    }

    names(to_remove_list) <- phenotypes
    to_remove_list_clean <- to_remove_list[unlist(lapply(to_remove_list, function(x){!is.null(x)}))]

    to_remove_df <- tibble::tibble(phenotype = names(to_remove_list_clean),
                                   proteins = unlist(lapply(to_remove_list_clean, function(x){paste0(x, collapse = ',')})))
    to_remove_df <- to_remove_df %>% tidyr::separate_rows(proteins)

    prot_to_phenotypes_clean <- dplyr::anti_join(prot_to_phenotypes,
                                                 to_remove_df,
                                                 by = c('key' = 'phenotype',
                                                        'QueryNode' = 'proteins'))

    INPUT.Phen.Paths <- dplyr::inner_join(INPUT.Phen.Paths_clean,
                                          prot_to_phenotypes_clean,
                                          by = c('QueryNode', 'key'))

    return(INPUT.Phen.Paths)
  }
  return(input_prot_paths)
}
