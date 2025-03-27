#' Find All Paths Between Two Nodes in a Prior Knowledge Network
#'
#' Identifies all paths with a length up to `max_length` between
#' a starting node (`v_start`) and a target node (`v_end`) in a tabular
#' Prior Knowledge Network (`PKN_table`).
#'
#' @param v_start Character, Gene Symbol of the starting node.
#' @param v_end Character, Gene Symbol of the target node.
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param max_length Integer, maximum path length between start and end nodes (1 to 4).
#'
#'
#' @return A data frame containing interactions from all the shortest paths
#'   between `v_start` and `v_end`.
#'   - `ENTITYA`: Starting node.
#'   - `INTERACTION`: Interaction type.
#'   - `ENTITYB`: Target node.
#'
#' @export
#'
#' @examples
#' data('PKN_human_atlas_dir')
#' find_all_paths(v_start = 'FLT3',
#'                v_end = 'TFEB',
#'                max_length = 4,
#'                PKN_table = PKN_human_atlas_dir)
#'
find_all_paths <- function(v_start, v_end, PKN_table, max_length) {
  if (max_length > 4) {
    stop('Max length is 4!')
  }

  PKN_table <- PKN_table %>%
    dplyr::select('ENTITYA', 'INTERACTION', 'ENTITYB') %>%
    dplyr::distinct()

  ## transform the PKN in a set of vectors
  entitya = PKN_table$ENTITYA
  entityb = PKN_table$ENTITYB
  interactions = PKN_table$INTERACTION

  ## step 1
  i = 1

  type_1 <- interactions[entitya == v_start]
  interactors_1 <- entityb[entitya == v_start]

  if (v_end %in% interactors_1) {
    ## define paths
    type_ve <- interactions[entitya == v_start & entityb == v_end]
    paths <-
      tibble::tibble(ENTITYA = v_start,
                     INTERACTION = type_ve,
                     ENTITYB = v_end)
    paths <- paths %>% dplyr::distinct()

    return(paths)

  } else{
    if (i == max_length) {
      return(NULL)
    } else{
      i = i + 1

      ## position of interactors in ENTITYA
      entitya_2 <- entitya[entitya %in% unique(interactors_1)]
      type_2 <- interactions[entitya %in% unique(interactors_1)]
      interactors_2 <- entityb[entitya %in% unique(interactors_1)]

      if (v_end %in% interactors_2) {
        ## take the subset of entitya2, type2 and interactors representing v_end
        entitya_2 <- entitya_2[interactors_2 == v_end]
        type_2 <- type_2[interactors_2 == v_end]
        interactors_2 <-  interactors_2[interactors_2 == v_end]

        ## take the subset of interactors_1 representing  proteins regulating interactors2
        type_1 <- type_1[interactors_1  %in% unique(entitya_2)]
        interactors_1 <-
          interactors_1[interactors_1 %in% unique(entitya_2)]

        paths <-
          tibble::tibble(
            ENTITYA = c(rep(v_start, length(
              interactors_1
            )), entitya_2),
            INTERACTION = c(type_1, type_2),
            ENTITYB = c(interactors_1, interactors_2)
          )
        paths <- paths %>% dplyr::distinct()

        return(paths)
      } else {
        if (i == max_length) {
          return(NULL)

        } else{
          i = i + 1

          entitya_3 <- entitya[entitya %in% unique(interactors_2)]
          type_3 <- interactions[entitya %in% unique(interactors_2)]
          interactors_3 <-
            entityb[entitya %in% unique(interactors_2)]

          if (v_end %in% interactors_3) {
            ## take the subset of entitya3, type3
            ## and interactors3 representing v_end
            entitya_3 <- entitya_3[interactors_3 == v_end]
            type_3 <- type_3[interactors_3 == v_end]
            interactors_3 <-  interactors_3[interactors_3 == v_end]

            ## take the subset of interactors_2 representing
            ## proteins that are regulators of interactors3
            type_2 <- type_2[interactors_2 %in% unique(entitya_3)]
            entitya_2 <-
              entitya_2[interactors_2 %in% unique(entitya_3)]
            interactors_2 <-
              interactors_2[interactors_2 %in% unique(entitya_3)]

            ## take the subset of interactors_1 representing
            ## proteins that are regulators of interactors2
            type_1 <- type_1[interactors_1 %in% unique(entitya_2)]
            interactors_1 <-
              interactors_1[interactors_1 %in% unique(entitya_2)]

            paths <-
              tibble::tibble(
                ENTITYA = c(rep(
                  v_start, length(interactors_1)
                ), entitya_2, entitya_3),
                INTERACTION = c(type_1, type_2, type_3),
                ENTITYB = c(interactors_1, interactors_2, interactors_3)
              )

            paths <- paths %>% dplyr::distinct()
            return(paths)
          } else {
            if (i == max_length) {
              return(NULL)
            } else{
              i <- i + 1

              entitya_4 <-
                entitya[entitya %in% unique(interactors_3)]
              type_4 <-
                interactions[entitya %in% unique(interactors_3)]
              interactors_4 <-
                entityb[entitya %in% unique(interactors_3)]

              if (v_end %in% interactors_4) {
                ## define paths
                entitya_4 <- entitya_4[interactors_4 == v_end]
                type_4 <- type_4[interactors_4 == v_end]
                interactors_4 <-
                  interactors_4[interactors_4 == v_end]

                ## take the subset of entitya3, type3 and interactors3
                ## representing v_end
                type_3 <-
                  type_3[interactors_3 %in% unique(entitya_4)]
                entitya_3 <-
                  entitya_3[interactors_3 %in% unique(entitya_4)]
                interactors_3 <-
                  interactors_3[interactors_3 %in% unique(entitya_4)]

                ## take the subset of interactors_2 representing
                ## the proteins that are regulators of interactors3

                type_2 <-
                  type_2[interactors_2 %in% unique(entitya_3)]
                entitya_2 <-
                  entitya_2[interactors_2 %in% unique(entitya_3)]
                interactors_2 <-
                  interactors_2[interactors_2 %in% unique(entitya_3)]

                ## take the subset of interactors_1 representing the proteins
                ## regulating interactors2
                type_1 <-
                  type_1[interactors_1 %in% unique(entitya_2)]
                interactors_1 <-
                  interactors_1[interactors_1 %in% unique(entitya_2)]

                paths <-
                  tibble::tibble(
                    ENTITYA = c(
                      rep(v_start, length(interactors_1)),
                      entitya_2,
                      entitya_3,
                      entitya_4
                    ),
                    INTERACTION = c(type_1, type_2, type_3, type_4),
                    ENTITYB = c(
                      interactors_1,
                      interactors_2,
                      interactors_3,
                      interactors_4
                    )
                  )

                paths <- paths %>% dplyr::distinct()
                return(paths)
              } else {
                return(NULL)
              }
            }
          }
        }
      }
    }
  }
}

#' Find All Shortest Paths Between Sets of Nodes
#'
#' Applies `find_all_paths` iteratively to find all paths between
#' multiple starting (`start_nodes_gn`) and target (`target_nodes_gn`) nodes.
#'
#' @param start_nodes_gn Character vector, Gene Symbols of starting nodes.
#' @param target_nodes_gn Character vector, Gene Symbols of target nodes.
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param max_length Integer, maximum path length between start and end nodes (1 to 4).
#'
#' @return A data frame containing interactions in all shortest paths
#'   between the given sets of starting and target nodes.
#'
#' @export
#'
#' @examples
#' data('PKN_human_atlas_dir')
#' get_all_shortest_path_custom(
#'   start_nodes_gn = c("FLT3", "KRAS"),
#'   target_nodes_gn = c("STAT5A", "TFEB"),
#'   PKN_table = PKN_human_atlas_dir,
#'   max_length = 4)
#'
get_all_shortest_path_custom <-
  function(start_nodes_gn,
           target_nodes_gn,
           PKN_table,
           max_length) {
    for (i in c(1:length(start_nodes_gn))) {
      for (j in c(1:length(target_nodes_gn))) {
        paths_df <-
          find_all_paths(start_nodes_gn[i],
                         target_nodes_gn[j],
                         PKN_table = PKN_table,
                         max_length)

        if (i == 1 & j == 1) {
          all_paths_df <- paths_df
        } else{
          all_paths_df <-
            dplyr::bind_rows(all_paths_df, paths_df) %>%
            dplyr::distinct()
        }
      }
    }
    return(all_paths_df)
  }

#' Convert Paths Data to an igraph Network
#'
#' Converts a data frame of interactions (`all_paths_df`) into an `igraph` object,
#' optionally adding direct one-step connections between intermediate nodes.
#'
#' @param all_paths_df Data frame, interactions in all shortest paths.
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param connect_all Logical, whether to connect intermediate nodes of shortest paths. Default: `TRUE`.
#'
#' @return An `igraph` object representing the network with all paths.
#'
#'
create_graph_from_paths <-
  function(all_paths_df, PKN_table, connect_all = TRUE) {
    if (nrow(all_paths_df) == 0) {
      stop(
        'SignalingProfiler ERROR: No paths found for you analytes.
        Try to not preprocessing the PKN'
      )
    }

    if (connect_all == TRUE) {
      ## transform PKN table in a graph
      pkn_graph <- pkn_table_as_graph(PKN_table)

      ## get all entities in all_paths_df
      entities <-
        unique(c(all_paths_df$ENTITYA, all_paths_df$ENTITYB))

      ## create the subgraph
      subgraph <- igraph::induced.subgraph(pkn_graph,
                                           igraph::V(pkn_graph)$name[igraph::V(pkn_graph)$name %in% entities])

      return(subgraph)

    } else {
      edges_paths_df <- dplyr::left_join(all_paths_df,
                                         PKN_table,
                                         by = c('ENTITYA',
                                                'INTERACTION',
                                                'ENTITYB')) %>%
        dplyr::relocate(ENTITYA, ENTITYB, INTERACTION)

      nodes_paths_df <-
        dplyr::bind_rows(
          edges_paths_df %>% dplyr::select(
            ENTITY = ENTITYA,
            ID = IDA,
            TYPE = TYPEA
          ),
          edges_paths_df %>% dplyr::select(
            ENTITY = ENTITYB,
            ID = IDB,
            TYPE = TYPEB
          )
        ) %>%
        dplyr::distinct() %>%
        dplyr::group_by(ENTITY) %>%
        dplyr::reframe(ID = paste0(ID, collapse = ';'),
                       TYPE = paste0(unique(TYPE), collapse = ';')) %>%
        dplyr::relocate(ENTITY)

      graph <-
        igraph::graph_from_data_frame(d = edges_paths_df,
                                      vertices = nodes_paths_df)

      return(graph)
    }
  }

#' Convert a Prior Knowledge Network Table to an igraph Object
#'
#' Transforms a tabular Prior Knowledge Network (`PKN_table`) into an
#' `igraph` object, creating nodes and edges with SIGNOR-compliant attributes
#'
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#'
#' @return An `igraph` object representing the Prior Knowledge Network.
#'
pkn_table_as_graph <- function(PKN_table) {
  nodes_df <-
    tidyr::tibble(
      ENTITY = c(PKN_table$ENTITYA, PKN_table$ENTITYB),
      ID = c(PKN_table$IDA, PKN_table$IDB),
      TYPE = c(PKN_table$TYPEA, PKN_table$TYPEB),
      DATABASE = c(PKN_table$DATABASEA, PKN_table$DATABASEB)
    ) %>% dplyr::distinct()

  nodes_df <- nodes_df %>%
    dplyr::group_by(ENTITY) %>%
    dplyr::reframe(
      ID = paste0(ID, collapse = ';'),
      TYPE = paste0(TYPE, collapse = ';'),
      DATABASE = paste0(DATABASE, collapse = ';')
    )

  PKN_graph <-
    igraph::graph_from_data_frame(d = PKN_table,
                                  vertices = nodes_df)

  return(PKN_graph)
}
