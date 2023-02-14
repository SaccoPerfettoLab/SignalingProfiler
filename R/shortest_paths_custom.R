
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

  # v_start = 'Cdk1'
  # v_end = 'Mapk3'
  # organism = 'mouse'
  # max_length = 5

  if(max_length > 4){
    stop('Max length is 4!')
  }

  PKN_table <- PKN_table %>% dplyr::select('ENTITYA', 'INTERACTION', 'ENTITYB')

  i = 1
  step1 <- PKN_table %>% dplyr::filter(ENTITYA == v_start)

  if(v_end %in% step1$ENTITYB){
    paths <- step1 %>% dplyr::filter(ENTITYB == v_end) %>% dplyr::distinct()
    return(paths)
  }else{
    if(i == max_length){ # stop if max_length = 1
      #message(paste0('No paths found of ', max_length, ' maximum length'))
      return(NULL)
    }else{
      i = i + 1
      step2 <- dplyr::left_join(step1,
                                PKN_table,
                                by = c('ENTITYB' = 'ENTITYA'),
                                suffix = c('step1', 'step2')) %>%
        dplyr::rename(ENTITYBstep1 = ENTITYB)

      if(v_end %in% step2$ENTITYBstep2){
        paths <- step2 %>% dplyr::filter(ENTITYBstep2 == v_end) %>% dplyr::distinct()

        step1_paths <- paths %>% dplyr::select("ENTITYA","INTERACTIONstep1","ENTITYBstep1")
        colnames(step1_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

        step2_paths <- paths %>% dplyr::select("ENTITYBstep1","INTERACTIONstep2","ENTITYBstep2")
        colnames(step2_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

        paths_df <- dplyr::bind_rows(step1_paths, step2_paths) %>%
          dplyr::distinct()

        return(paths_df)
      }else{
        if(i == max_length){ # stop if max_length == 2
          #message(paste0('No paths found of ', max_length, ' maximum length'))
          return(NULL)
        }else{
          i = i + 1
          step3 <- dplyr::left_join(step2, PKN_table,
                                    by = c('ENTITYBstep2' = 'ENTITYA'),
                                    suffix = c('step2', 'step3')) %>%
            dplyr::rename(ENTITYBstep3 = ENTITYB,
                          INTERACTIONstep3 = INTERACTION)

          if(v_end %in% step3$ENTITYBstep3){
            paths <- step3 %>% dplyr::filter(ENTITYBstep3 == v_end) %>% dplyr::distinct()

            step1_paths <- paths %>% dplyr::select("ENTITYA","INTERACTIONstep1","ENTITYBstep1")
            colnames(step1_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

            step2_paths <- paths %>% dplyr::select("ENTITYBstep1","INTERACTIONstep2","ENTITYBstep2")
            colnames(step2_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

            step3_paths <- paths %>% dplyr::select("ENTITYBstep2","INTERACTIONstep3","ENTITYBstep3")
            colnames(step3_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

            paths_df <- dplyr::bind_rows(step1_paths, step2_paths, step3_paths) %>%
              dplyr::distinct()

            return(paths_df)
          }else{
            if(i == max_length){ # stop if max_length == 3
              #message(paste0('No paths found of ', max_length, ' maximum length'))
              return(NULL)
            }else{
              i = i + 1
              step4 <- dplyr::left_join(step3, PKN_table,
                                        by = c('ENTITYBstep3' = 'ENTITYA'),
                                        suffix = c('step3', 'step4')) %>%
                dplyr::rename(ENTITYBstep4 = ENTITYB,
                              INTERACTIONstep4 = INTERACTION)

              if(v_end %in% step4$ENTITYBstep4){
                paths <- step4 %>% dplyr::filter(ENTITYBstep4 == v_end) %>% dplyr::distinct()
                step1_paths <- paths %>% dplyr::select("ENTITYA","INTERACTIONstep1","ENTITYBstep1")
                colnames(step1_paths) <-  c('ENTITYA', 'INTERACTION', 'ENTITYB')

                step2_paths <- paths %>% dplyr::select("ENTITYBstep1","INTERACTIONstep2", "ENTITYBstep2")
                colnames(step2_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

                step3_paths <- paths %>% dplyr::select("ENTITYBstep2","INTERACTIONstep3","ENTITYBstep3")
                colnames(step3_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

                step4_paths <- paths %>% dplyr::select("ENTITYBstep3","INTERACTIONstep4","ENTITYBstep4")
                colnames(step4_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')

                paths_df <- dplyr::bind_rows(step1_paths, step2_paths, step3_paths, step4_paths) %>%
                  dplyr::distinct()
                return(paths_df)
              }else{
                #message(paste0('No paths found from ', v_start, ' to ', v_end))
                return(NULL)
              }

                # }else{
                #   i = i + 1
                #   step5 <- dplyr::left_join(step4, PKN_table,
                #                             by = c('ENTITYBstep3' = 'ENTITYA'),
                #                             suffix = c('step4', 'step5')) %>%
                #     dplyr::rename(ENTITYBstep5 = ENTITYB,
                #                   INTERACTIONstep5 = INTERACTION)
                #   if(v_end %in% step5$ENTITYBstep5){
                #     paths <- step5 %>% dplyr::filter(ENTITYBstep5 == v_end) %>% dplyr::distinct()
                #     step1_paths <- paths %>% dplyr::select("ENTITYA","INTERACTIONstep1","ENTITYBstep1")
                #
                #     colnames(step1_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')
                #
                #     step2_paths <- paths %>% dplyr::select("ENTITYBstep1","INTERACTIONstep2", "ENTITYBstep2")
                #     colnames(step2_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')
                #
                #     step3_paths <- paths %>% dplyr::select("ENTITYBstep2","INTERACTIONstep3","ENTITYBstep3")
                #     colnames(step3_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')
                #
                #     step4_paths <- paths %>% dplyr::select("ENTITYBstep3","INTERACTIONstep4","ENTITYBstep4")
                #     colnames(step4_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')
                #
                #     step5_paths <- paths %>% dplyr::select("ENTITYBstep4","INTERACTIONstep5","ENTITYBstep5")
                #     colnames(step5_paths) <- c('ENTITYA', 'INTERACTION', 'ENTITYB')
                #
                #     paths_df <- dplyr::bind_rows(step1_paths, step2_paths, step3_paths, step4_paths, step5_paths) %>%
                #       dplyr::distinct()
                #     return(paths_df)
                #   }else{
                #     message(paste0('No paths found of ', max_length, ' maximum length'))
                #   }
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
        all_paths_df <- dplyr::bind_rows(all_paths_df, paths_df)
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

  edges_paths_df <- dplyr::left_join(all_paths_df, PKN_table, by = c('ENTITYA', 'INTERACTION', 'ENTITYB')) %>%
    dplyr::relocate(IDA, IDB)

  edges_paths_df$PHOSPHO_KEY_GN_SEQ <- ''
  edges_paths_df$PHOSPHO_KEY_GN_SEQ[!is.na(edges_paths_df$SEQUENCE)] <- paste0(edges_paths_df$ENTITYB[!is.na(edges_paths_df$SEQUENCE)], '-', edges_paths_df$SEQUENCE[!is.na(edges_paths_df$SEQUENCE)])

  nodes_paths_df <- dplyr::bind_rows(edges_paths_df %>% dplyr::select(ENTITY = ENTITYA, ID = IDA),
                                     edges_paths_df %>% dplyr::select(ENTITY = ENTITYB, ID = IDB)) %>%
    dplyr::distinct(ID,.keep_all = TRUE) %>% #
    dplyr::relocate(ID)

  graph <- igraph::graph_from_data_frame(d = edges_paths_df, vertices = nodes_paths_df)

  return(graph)
}


