# scripts for carnival analysis

#' generateTFList
#'
#' @param df df of inferred proteins
#' @param top integer
#' @param access_idx integer
#'
#' @return list of proteins
#'
#' @examples
generateTFList <- function (df = df, top = 50, access_idx = 1){
  if (top == "all") {
    top <- nrow(df)
  }
  if (top > nrow(df)) {
    warning("Number of to TF's inserted exceeds the number of actual TF's in the\n            data frame. All the TF's will be considered.")
    top <- nrow(df)
  }
  ctrl <- intersect(x = access_idx, y = 1:ncol(df))
  if (length(ctrl) == 0) {
    stop("The indeces you inserted do not correspond to \n              the number of columns/samples")
  }
  returnList <- list()
  for (ii in 1:length(ctrl)) {
    tfThresh <- sort(x = abs(df[, ctrl[ii]]), decreasing = TRUE)[top]
    temp <- which(abs(df[, ctrl[ii]]) >= tfThresh)
    currDF <- matrix(data = , nrow = 1, ncol = top)
    colnames(currDF) <- rownames(df)[temp[1:top]]
    currDF[1, ] <- df[temp[1:top], ctrl[ii]]
    currDF <- as.data.frame(currDF)
    returnList[[length(returnList) + 1]] <- currDF
  }
  names(returnList) <- colnames(df)[ctrl]
  return(returnList)
}

#' formatting_proteins_for_carnival
#'
#' @param proteins_df dataframe with inferred proteins with specific mf
#'
#' @return a list of proteins as desired by CARNIVAL
#'
#' @examples
formatting_proteins_for_carnival <- function(proteins_df){

  proteins_df <- proteins_df %>%
    dplyr::select(UNIPROT, t = final_score) %>%
    dplyr::distinct()
  proteins_activities <- data.frame(unique(proteins_df))
  rownames(proteins_activities) <- unique(unique(proteins_activities$UNIPROT))
  proteins_activities$UNIPROT <- NULL
  proteinList <- generateTFList(proteins_activities, top = 'all', access_idx = 1)

  return(proteinList)
}

#' add_output_carnival_nodes_attributes
#'
#' @param carnival_result runCARNIVAL output
#' @param proteins_df inferred proteins df
#' @param organism string, 'mouse' or 'human'
#'
#' @return dataframe with nodes and attributes
#' @export
#'
#' @examples
add_output_carnival_nodes_attributes <- function(carnival_result,
                                                 proteins_df,
                                                 organism){

  # get db and protein db for adding attributes
  if(organism == 'mouse'){
    PKN_proteins <- PKN_proteins_mouse
    db <- db_mouse
  }else if(organism == 'human'){
    PKN_proteins <- PKN_proteins_human
    db <- db_human
  }else{
    error('Please provide a valide organism')
  }


  nodes <- tibble::as_tibble(carnival_result$nodesAttributes)
  optimal_edges <- tibble::as_tibble(carnival_result$weightedSIF)

  nodes_igraph_ids <- unique(c(gsub('_', '-',optimal_edges$Node1),
                               gsub('_', '-',optimal_edges$Node2)))

  optimal_nodes <- nodes %>%
    dplyr::filter(gsub('_', '-', Node) %in% nodes_igraph_ids) %>%
    dplyr::select(Node, 'carnival_activity' = AvgAct) %>%
    dplyr::mutate_at('carnival_activity', as.numeric) %>%
    dplyr::mutate(Node = gsub('_', '-', Node))

  nodes_df <- dplyr::left_join(optimal_nodes, PKN_proteins, by = c('Node' = 'ID')) %>%
    dplyr::rename('UNIPROT' = 'Node','gene_name' = 'ENTITY')

  nodes_df <- dplyr::left_join(nodes_df, proteins_df, by = c('gene_name', 'UNIPROT')) %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score)

  # annotate missing genes in the network
  message('GO Molecular Function annotation of optimized nodes')
  nodes_df <- molecular_function_annotation(nodes_df)

  # add new attributes
  nodes_df$discordant <- FALSE

  nodes_df$discordant[as.numeric(nodes_df$carnival_activity) * nodes_df$final_score < 0] <- TRUE

  nodes_df <- nodes_df %>%
    dplyr::relocate(UNIPROT) %>%
    dplyr::distinct()

  nodes_df$UNIPROT <- gsub('_', '-',nodes_df$UNIPROT)

  return(nodes_df)
}

#' add_output_carnival_edges_attributes
#'
#' @param carnival_result runCARNIVAL output
#'
#' @return dataframe with edges attributes
#' @export
#'
#' @examples
add_output_carnival_edges_attributes <- function(carnival_result){

  optimal_edges <- tibble::as_tibble(carnival_result$weightedSIF)

  edges_df <- tibble::tibble(source = gsub('_', '-',optimal_edges$Node1),
                     target = gsub('_', '-',optimal_edges$Node2),
                     sign = optimal_edges$Sign,
                     carnival_weight = as.numeric(optimal_edges$Weight))
  return(edges_df)
}

#' create_discretized_initiators_for_carnival
#'
#' @param proteins_df tibble of initiator proteins
#'
#' @return initiators' tibble with discrtized activities
#' @export
#'
#' @examples
create_discretized_initiators_for_carnival <- function(proteins_df){

  proteins_df$discretized <- 0
  proteins_df$discretized[proteins_df$final_score < 0] <- -1
  proteins_df$discretized[proteins_df$final_score > 0] <- 1

  initiators_df <- proteins_df %>%
    dplyr::filter(discretized == 1 | discretized == -1) %>%
    dplyr::select(UNIPROT, gene_name, final_score = discretized, mf)

  return(initiators_df)
}


#' union_of_graphs
#'
#' @param graph_1 run1 igraph object of run_carnival_and_create_graph
#' @param graph_2 run2 igraph object of run_carnival_and_create_graph
#' @param proteins_df tibble of inferred proteins in the naive network
#' @param files boolean value, TRUE if you want output files
#' @param path_sif
#' @param path_rds
#'
#' @return
#' @export
#'
#' @examples
union_of_graphs <- function(graph_1, graph_2, proteins_df, files,
                            path_sif = './union_graph.sif',
                            path_rds = './union_graph.rds'){

  # graph_1 = output1$igraph_network
  # graph_2 = output2$igraph_network
  # proteins_df <- carnival_input_toy

  #nodes df
  nodes_rec_kin <- tibble::as_tibble(igraph::as_data_frame(graph_1, what = c('vertices'))) %>%
    dplyr::rename(UNIPROT = 'name') %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT)

  nodes_kin_tf <- tibble::as_tibble(igraph::as_data_frame(graph_2, what = c('vertices'))) %>%
    dplyr::rename(UNIPROT = 'name') %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT)

  nodes <- dplyr::bind_rows(nodes_rec_kin, nodes_kin_tf) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene_name)

  nodes <- dplyr::left_join(nodes, proteins_df, by = c('gene_name', 'UNIPROT')) %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score) %>%
    dplyr::relocate(UNIPROT)

  # se sono duplicati nei due network prendo la media
  nodes1 <- nodes %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(carnival_activity = mean(carnival_activity))

  nodes <- nodes %>%
    dplyr::select(UNIPROT, gene_name, mf, final_score) %>%
    dplyr::distinct()

  nodes <- dplyr::left_join(nodes1, nodes, by = c('gene_name')) %>%
    dplyr::relocate(UNIPROT)

  # adding attributes
  nodes$discordant <- FALSE
  nodes$discordant[nodes$carnival_activity * nodes$final_score < 0] <- TRUE

  #edges df
  edges_rec_kin <- tibble::as_tibble(igraph::as_data_frame(graph_1, what = c('edges')))
  edges_kin_tf <- tibble::as_tibble(igraph::as_data_frame(graph_2, what = c('edges')))

  edges <- dplyr::bind_rows(edges_rec_kin, edges_kin_tf)

  union_graph <- igraph::graph_from_data_frame(d = data.frame(edges), vertices = nodes, directed = T)

  if(files == TRUE){
    igraphToSif(union_graph, path_sif, 'sign')
    saveRDS(union_graph, path_rds)
  }

  return(list(network = union_graph,
              nodes_df = nodes,
              edges_df = edges))
}


#' Title
#'
#' @param source_df tibble with source nodes discretized among 1 and -1
#' @param target_df tibble with target nodes in a continuous range of activity
#' @param naive_network tibble with naive network in SIF format
#' @param solver_path path of the solver you use in CARNIVAL
#' @param solver name of solver
#' @param proteins_df all inferred proteins df
#' @param organism string, organism you are using
#' @param phospho_df tibble containing phosphoproteomics data
#' @param files boolean value, TRUE if you want output files
#' @param path_sif path of the sif output file of network
#' @param path_rds path of the rds output file of network
#'
#' @return list with igraph object, nodes df and edges df
#' @export
#'
#' @examples
run_carnival_and_create_graph <- function(source_df,
                                          target_df,
                                          naive_network,
                                          solver_path = NULL, solver = NULL,
                                          proteins_df,
                                          organism,
                                          phospho_df,
                                          files = TRUE,
                                          path_sif = './optimized_network.sif',
                                          path_rds = './optimized_network.RDS'){

  # discretize initiators
  source_df_disc <- create_discretized_initiators_for_carnival(source_df)

  carnival_result <- CARNIVAL::runCARNIVAL( inputObj = formatting_proteins_for_carnival(source_df_disc)$t,
                                            measObj = formatting_proteins_for_carnival(target_df)$t,
                                            netObj = unique(naive_network), #removing duplicated edges
                                            solverPath = solver_path,
                                            solver = solver,
                                            timelimit=7200,
                                            mipGAP=0.05,
                                            poolrelGAP=0.0001,
                                            betaWeight = 0.2)

  # organism <- 'mouse'
  # format result
  nodes_df <- add_output_carnival_nodes_attributes(carnival_result,
                                                   proteins_df,
                                                   organism) %>%
    dplyr::filter(carnival_activity != 0) %>%
    dplyr::distinct()


  if(nrow(nodes_df) == 0){
    message('No network found for your experiment')
    return(NULL)
  }else{

    edges_df <- add_output_carnival_edges_attributes(carnival_result) %>%
      dplyr::filter(carnival_weight != 0)

    CARNIVAL_igraph_network <- igraph::graph_from_data_frame(edges_df,
                                                             nodes_df,
                                                             directed = TRUE)

    if(files == TRUE){
      igraphToSif(CARNIVAL_igraph_network, path_sif, 'sign')
      saveRDS(CARNIVAL_igraph_network, path_rds)
    }

    # create carnival output as graph
    return(list(igraph_network = CARNIVAL_igraph_network,
                nodes_df = nodes_df,
                edges_df = edges_df))

  }
}

#' convert_output_nodes_in_next_input
#'
#' @param carnival_result output list of run_carnival_and_create_graph function
#'
#' @return formatted_nodes
#' @export
#'
#' @examples
convert_output_nodes_in_next_input <- function(carnival_result){
  nodes <- carnival_result$nodes_df
  formatted_nodes <- nodes %>% dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf)
  return(formatted_nodes)
}

#' expand_and_map_edges
#'
#' @param optimized_graph_rds # igraph object output of the run_carnival_and_create_graph
#' @param optimized_graph_sif # table output of the run_carnival_and_create_graph
#' @param organism #string, human or mouse
#' @param phospho_df # tibble of phosphoproteomics data
#' @param files # boolean value, TRUE if you want output files
#' @param path_sif # string of the path  of output sif file
#' @param path_rds # string of the path of RDS file
#'
#' @return # list with igraph object, optimized nodes df and optimized edges df
#' @export
#'
#' @examples
expand_and_map_edges <- function(optimized_graph_rds, optimized_graph_sif,
                                 organism, phospho_df, files, path_sif, path_rds){

  # organism = 'mouse'
  # phospho_df <- phospho_toy_df
  db <- choose_database_for_building(organism, format = 'table')

  optimized_graph_sif <- optimized_graph_sif %>% dplyr::mutate_at('interaction', as.character)
  db <- db %>% dplyr::mutate_at('INTERACTION', as.character)

  # create substring if necessary
  if(nchar(phospho_df$sequence_window[1]) != 15){
    center <- (nchar(phospho_df$sequence_window[1])+1)/2
    phospho_df <- phospho_df %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
  }else{ #if it is rename the column just to have a homogeneous variable
    phospho_df <- phospho_df %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }

  #phospho_df <- phospho_toy_df
  phospho_df <- phospho_df %>%
    dplyr::mutate(phosphositeID = paste0(gene_name, '-', sequence_window_sub)) %>%
    dplyr::select(gene_name, aminoacid, position, phosphositeID, difference, significant)

  # join db on network
  #taking from DBs the optimized interaction in order to get duplicated interaction
  merged_df <- dplyr::inner_join(db, optimized_graph_sif, by = c('IDA' = 'source',
                                                                 'INTERACTION' = 'interaction',
                                                                 'IDB' = 'target'))

  expanded_optimized_model <- unique(merged_df)

  # *** CREATION OF A GRAPH OBJECT ** #
  # nodes are the one in optimized network, while edges become the expanded ones

  # ** NODES DATAFRAME **
  nodes_df <- igraph::as_data_frame(optimized_graph_rds, what = 'vertice') %>%
    tibble::as_tibble()# nodes' dataframe already present
  edges_df <- igraph::as_data_frame(optimized_graph_rds, what = 'edges') %>%
    tibble::as_tibble()

  # ** EDGES DATAFRAME **
  edges_df_new <- expanded_optimized_model # should be created from expanded optimized model

  # adding some attributes: carnival_weight, is_quantified, is_significant, phosphosite_value
  edges_df_new$carnival_weight <- 0
  edges_df_new$phosphosite_value <- 0
  edges_df_new$is_significant <- FALSE
  edges_df_new$is_quantified <- FALSE

  # create a GN - SEQUENCE key
  edges_df_new$PHOSPHO_KEY_GN_SEQ <- paste0(edges_df_new$ENTITYB,'-', edges_df_new$SEQUENCE)


  #  ** setting $Is_Quantified **
  # setting Is_Quantified to TRUE if the combination GeneName - Seq is in experimental phosphosites (before with sequence)

  edges_df_new$is_quantified[edges_df_new$PHOSPHO_KEY_GN_SEQ %in% phospho_df$phosphositeID] <- TRUE

  # ** setting phosphosite_value
  for (i in c(1:length(edges_df_new$ENTITYA))){
    for( j in c(1:length(phospho_df$aminoacid))){
      if(is.na(edges_df_new$PHOSPHO_KEY_GN_SEQ[i]) == FALSE){
        if(edges_df_new$PHOSPHO_KEY_GN_SEQ[i] == phospho_df$phosphositeID[j]){

          edges_df_new$phosphosite_value[i] <- phospho_df$difference[j]

          if(is.na(phospho_df$significant[j]) == FALSE){
            edges_df_new$is_significant[i] <- TRUE
          }
        }
      }
    }
  }

  # ** setting $Weight_CARNIVAL **
  for (i in c(1:length(edges_df_new$IDA))){
    for (j in c(1:length(edges_df$from))){
      if(edges_df_new$IDA[i] == edges_df$from[j] &
         edges_df_new$IDB[i] == edges_df$to[j] &
         edges_df_new$INTERACTION[i] == edges_df$sign[j]){
        edges_df_new$carnival_weight[i] <- edges_df$carnival_weight[j]
      }
    }
  }

  edges_df_new <- edges_df_new %>%
    dplyr::select(source = IDA, target = IDB, sign = INTERACTION, carnival_weight, mechanism = MECHANISM,
                  residue = RESIDUE,sequence = SEQUENCE, is_significant,phosphosite_value) %>%
    dplyr::distinct()

  edges_df_new$residue[edges_df_new$is_significant == TRUE] <- paste0(edges_df_new$residue[edges_df_new$is_significant == TRUE], '*')
  edges_df_new$phosphosite_value[edges_df_new$is_significant == TRUE] <- paste0(edges_df_new$phosphosite_value[edges_df_new$is_significant == TRUE], '*')

  edges_df_new$phosphosite_value[edges_df_new$phosphosite_value == '0'] <- NA

  edges_df_new_with_amino <- edges_df_new %>%
    dplyr::filter(!is.na(phosphosite_value)) %>%
    dplyr::group_by(source, target, sign) %>%
    dplyr::summarise(aminoacid = paste0(residue,collapse = ';'),
              FC = paste0(as.character(phosphosite_value), collapse = ';'))

  edges_df_new_final <- dplyr::left_join(edges_df_new %>%
                                           dplyr::select(source, target, sign, carnival_weight) %>%
                                           dplyr::distinct(),
                                  edges_df_new_with_amino)

  edges_df_new_final$is_quantified <- FALSE
  edges_df_new_final$is_quantified[is.na(edges_df_new_final$FC) == FALSE] <- TRUE

  CARNIVAL_graph <- igraph::graph_from_data_frame(edges_df_new_final, nodes_df, directed = TRUE)

  unique(c(edges_df_new_final$source, edges_df_new_final$target)) %in% nodes_df$name

  if(files == TRUE){
    igraphToSif(CARNIVAL_graph, path_sif, 'sign')
    saveRDS(CARNIVAL_graph, path_rds)
  }

  # create carnival output as graph
  return(list(igraph_network = CARNIVAL_graph,
              nodes_df = nodes_df,
              edges_df = edges_df_new_final))
}

