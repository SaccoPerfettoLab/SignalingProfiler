
#' find_all_paths
#'
#' @param v_start gene_name of starting node
#' @param v_end gene_name of target node
#' @param PKN_table tibble of all causal interactions
#' @param max_length max_length of the desired path, max = 4
#'
#' @return tibble of all shortest paths among start and end nodes
#' @export
#'
#' @examples
find_all_paths <- function(v_start, v_end, PKN_table, max_length){

  # v_start = 'FLT3'
  # v_end = 'TFEB'
  # max_length = 4
  # PKN_table = PKN_human_atlas_dir
  #

  if(max_length > 4){
    stop('Max length is 4!')
  }

  PKN_table <- PKN_table %>% dplyr::select('ENTITYA', 'INTERACTION', 'ENTITYB') %>% dplyr::distinct()

  # Transform the PKN in a set of vectors
  entitya = PKN_table$ENTITYA
  entityb = PKN_table$ENTITYB
  interactions = PKN_table$INTERACTION

  # STEP 1
  i = 1

  type_1 <- interactions[entitya == v_start]
  interactors_1 <- entityb[entitya == v_start]

  if(v_end %in% interactors_1){
    # define in some way paths

    type_ve <- interactions[entitya == v_start & entityb == v_end]
    paths <- tibble(ENTITYA = v_start, INTERACTION = type_ve, ENTITYB = v_end)
    paths <- paths %>% distinct()

    return(paths)

  }else{
    if(i == max_length){
      return(NULL)
    }else{
      i = i+1

      # Position of interactors in ENTITYA
      entitya_2 <- entitya[entitya %in% unique(interactors_1)]
      type_2 <- interactions[entitya %in% unique(interactors_1)]
      interactors_2 <- entityb[entitya %in% unique(interactors_1)]

      if(v_end %in% interactors_2){

        # Take the subset of entitya2, type2 and interactors that represents v_end
        entitya_2 <- entitya_2[interactors_2 == v_end]
        type_2 <- type_2[interactors_2 == v_end]
        interactors_2 <-  interactors_2[interactors_2 == v_end]

        # Take the subset of interactors_1 representing the proteins that are regulators of interactors2
        type_1 <- type_1[interactors_1  %in% unique(entitya_2)]
        interactors_1 <- interactors_1[interactors_1 %in% unique(entitya_2)]

        paths <- tibble(ENTITYA = c(rep(v_start, length(interactors_1)), entitya_2),
                        INTERACTION = c(type_1, type_2),
                        ENTITYB = c(interactors_1, interactors_2))
        paths <- paths %>% distinct()

        return(paths)
      } else {

        if( i == max_length){
          return(NULL)

        }else{
          i = i + 1

          entitya_3 <- entitya[entitya %in% unique(interactors_2)]
          type_3 <- interactions[entitya %in% unique(interactors_2)]
          interactors_3 <- entityb[entitya %in% unique(interactors_2)]

          if(v_end %in% interactors_3){

            # Take the subset of entitya3, type3 and interactors3 that represents v_end
            entitya_3 <- entitya_3[interactors_3 == v_end]
            type_3 <- type_3[interactors_3 == v_end]
            interactors_3 <-  interactors_3[interactors_3 == v_end]

            # Take the subset of interactors_2 representing the proteins that are regulators of interactors3

            type_2 <- type_2[interactors_2 %in% unique(entitya_3)]
            entitya_2 <- entitya_2[interactors_2 %in% unique(entitya_3)]
            interactors_2 <- interactors_2[interactors_2 %in% unique(entitya_3)]


            # Take the subset of interactors_1 representing the proteins that are regulators of interactors2
            type_1 <- type_1[interactors_1 %in% unique(entitya_2)]
            interactors_1 <- interactors_1[interactors_1 %in% unique(entitya_2)]

            paths <- tibble(ENTITYA = c(rep(v_start, length(interactors_1)), entitya_2, entitya_3),
                            INTERACTION = c(type_1, type_2, type_3),
                            ENTITYB = c(interactors_1, interactors_2, interactors_3))

            paths <- paths %>% distinct()
            return(paths)
          } else {

            if ( i == max_length){
              return(NULL)
            }else{

              i = i + 1

              entitya_4 <- entitya[entitya %in% unique(interactors_3)]
              type_4 <- interactions[entitya %in% unique(interactors_3)]
              interactors_4 <- entityb[entitya %in% unique(interactors_3)]

              if(v_end %in% interactors_4){

                # define in some way paths
                entitya_4 <- entitya_4[interactors_4 == v_end]
                type_4 <- type_4[interactors_4 == v_end]
                interactors_4 <-  interactors_4[interactors_4 == v_end]

                # Take the subset of entitya3, type3 and interactors3 that represents v_end
                type_3 <- type_3[interactors_3 %in% unique(entitya_4)]
                entitya_3 <- entitya_3[interactors_3 %in% unique(entitya_4)]
                interactors_3 <- interactors_3[interactors_3 %in% unique(entitya_4)]

                # Take the subset of interactors_2 representing the proteins that are regulators of interactors3

                type_2 <- type_2[interactors_2 %in% unique(entitya_3)]
                entitya_2 <- entitya_2[interactors_2 %in% unique(entitya_3)]
                interactors_2 <- interactors_2[interactors_2 %in% unique(entitya_3)]


                # Take the subset of interactors_1 representing the proteins that are regulators of interactors2
                type_1 <- type_1[interactors_1 %in% unique(entitya_2)]
                interactors_1 <- interactors_1[interactors_1 %in% unique(entitya_2)]

                paths <- tibble(ENTITYA = c(rep(v_start, length(interactors_1)), entitya_2, entitya_3, entitya_4),
                                INTERACTION = c(type_1, type_2, type_3, type_4),
                                ENTITYB = c(interactors_1, interactors_2, interactors_3, interactors_4))

                paths <- paths %>% distinct()
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


  #' get_all_shortest_path_custom
  #'
  #' @param start_nodes_gn vector of gene_names of starting nodes
  #' @param target_nodes_gn vector of gene_names of target nodes
  #' @param PKN_table tibble of all causal interactions
  #' @param max_length integer, 1 to 4
  #'
  #' @return
  #' @export
  #'
  #' @examples
  get_all_shortest_path_custom <- function(start_nodes_gn, target_nodes_gn, PKN_table, max_length){
    # start_nodes_gn <- kinphos_gn
    # target_nodes_gn <- tfs_gn
    # path_length <- 'shortest'


    for(i in c(1:length(start_nodes_gn))){
      for(j in c(1:length(target_nodes_gn))){
        paths_df <- find_all_paths(start_nodes_gn[i], target_nodes_gn[j], PKN_table = PKN_table, max_length)

        if( i == 1 & j == 1){
          all_paths_df <- paths_df
        }else{
          all_paths_df <- dplyr::bind_rows(all_paths_df, paths_df) %>% dplyr::distinct()
        }
      }
    }
    return(all_paths_df)
  }
  #' create_graph_from_paths
  #'
  #' @param all_paths_df tibble with edges of shortest paths
  #' @param PKN_table tibble of all causal interactions
  #'
  #' @return igraph object of naive network
  #'
  #' @examples
  create_graph_from_paths <- function(all_paths_df, PKN_table){


    if(nrow(all_paths_df) == 0){
      stop('SignalingProfiler ERROR: No paths found for you analytes. Try to not preprocessing the PKN')
    }

    edges_paths_df <- dplyr::left_join(all_paths_df, PKN_table,
                                       by = c('ENTITYA', 'INTERACTION', 'ENTITYB')) %>%
      dplyr::relocate(ENTITYA, ENTITYB, INTERACTION)

    #edges_paths_df$PHOSPHO_KEY_GN_SEQ <- ''
    #edges_paths_df$PHOSPHO_KEY_GN_SEQ[!is.na(edges_paths_df$SEQUENCE)] <- paste0(edges_paths_df$ENTITYB[!is.na(edges_paths_df$SEQUENCE)], '-', edges_paths_df$SEQUENCE[!is.na(edges_paths_df$SEQUENCE)])

    nodes_paths_df <- dplyr::bind_rows(edges_paths_df %>% dplyr::select(ENTITY = ENTITYA, ID = IDA, TYPE = TYPEA),
                                       edges_paths_df %>% dplyr::select(ENTITY = ENTITYB, ID = IDB, TYPE = TYPEB)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(ENTITY) %>%
      dplyr::reframe(ID = paste0(ID, collapse = ';'), TYPE = paste0(unique(TYPE), collapse = ';')) %>%
      dplyr::relocate(ENTITY)

    graph <- igraph::graph_from_data_frame(d = edges_paths_df, vertices = nodes_paths_df)

    return(graph)
  }


