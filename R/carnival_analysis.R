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
#' @export
#'
#' @examples
formatting_proteins_for_carnival <- function(proteins_df){

  proteins_df <- proteins_df %>%
    dplyr::select(gene_name, t = final_score) %>%
    dplyr::distinct()
  proteins_activities <- data.frame(unique(proteins_df))
  rownames(proteins_activities) <- unique(unique(proteins_activities$gene_name))
  proteins_activities$gene_name <- NULL
  proteinList <- generateTFList(proteins_activities, top = 'all', access_idx = 1)

  return(proteinList)
}

#' add_output_carnival_nodes_attributes
#'
#' @param carnival_result runCARNIVAL output
#' @param proteins_df inferred proteins df
#' @param organism string, 'mouse' or 'human'
#' @param with_atlas Boolean value, default TRUE, if FALSE exclude kinome atlas interactions
#' @param direct Boolean value, default FALSE, if TRUE uses direct interactions
#'
#' @return dataframe with nodes and attributes
#' @export
#'
#' @examples
add_output_carnival_nodes_attributes <- function(carnival_result,
                                                 proteins_df,
                                                 organism,
                                                 with_atlas = TRUE,
                                                 direct = FALSE){

  # get db and protein db for adding attributes
  # !!!!!!!!!!!!!!!!!! to add direct attribute
  if(organism == 'mouse'){
    if(with_atlas == TRUE){
      stop('If organism is \'mouse\' with_atlas parameter must be FALSE')
    }else{
      PKN_proteins <- PKN_proteins_mouse
      if(direct == TRUE){
        db <- db_mouse_dir
      }else{
        db <- db_mouse_ind
      }
    }
  }else if(organism == 'human'){
    if(with_atlas == TRUE){
      PKN_proteins <- PKN_proteins_human_atlas

      if(direct == TRUE){
        db <- db_human_atlas_dir
      }else{
        db <- db_human_atlas_ind
      }

    }else{
      PKN_proteins <- PKN_proteins_human

      if(direct == TRUE){
        db <- db_human_dir
      }else{
        db <- db_human_ind
      }
    }
  }else{
    error('Please provide a valide organism')
  }

  nodes <- tibble::as_tibble(carnival_result$nodesAttributes)
  optimal_edges <- tibble::as_tibble(carnival_result$weightedSIF)

  nodes_igraph_ids <- unique(c(optimal_edges$Node1,
                              optimal_edges$Node2))

  optimal_nodes <- nodes %>%
    dplyr::filter(Node %in% nodes_igraph_ids) %>%
    dplyr::select(Node, 'carnival_activity' = AvgAct) %>%
    dplyr::mutate_at('carnival_activity', as.numeric) #%>%
    #dplyr::mutate(Node = gsub('_', '-', Node))

  # Taking UNIPROTs from PKN
  nodes_df <- dplyr::left_join(optimal_nodes, PKN_proteins, by = c('Node' = 'ENTITY')) %>%
    dplyr::rename('gene_name' = 'Node','UNIPROT' = 'ID')

  nodes_df <- dplyr::left_join(nodes_df, proteins_df %>% dplyr::select(-UNIPROT), by = c('gene_name')) %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score, method)

  # annotate missing genes in the network
  #message('GO Molecular Function annotation of optimized nodes')
  #nodes_df <- molecular_function_annotation(nodes_df)

  # add new attributes
  nodes_df$discordant <- FALSE

  nodes_df$discordant[as.numeric(nodes_df$carnival_activity) * nodes_df$final_score < 0] <- TRUE

  nodes_df <- nodes_df %>%
    dplyr::relocate(gene_name) %>%
    dplyr::distinct()

  #nodes_df$UNIPROT <- gsub('_', '-',nodes_df$UNIPROT)

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

  edges_df <- tibble::tibble(source = optimal_edges$Node1,
                             target = optimal_edges$Node2,
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
    dplyr::select(gene_name, final_score = discretized, mf)

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

  # graph_1 = transc_graph
  # graph_2 = sign_graph
  # proteins_df = proteins_df
  # files = TRUE
  # path_sif = paste0(carnival_dir, patient, '_union_model.sif')
  # path_rds = paste0(carnival_dir, patient, '_union_object.RDS')

  # graph_1 = output1$igraph_network
  # graph_2 = output2$igraph_network
  # proteins_df <- carnival_input_toy

  if(sum(unlist(grepl('UNIPROT', colnames(proteins_df)))) == 1){
    proteins_df$UNIPROT <- NULL #remove old UNIPROT and take the one of database
  }

  #nodes df
  nodes_rec_kin <- tibble::as_tibble(igraph::as_data_frame(graph_1, what = c('vertices'))) %>%
    dplyr::rename(gene_name = 'name') %>%
    dplyr::select(gene_name, UNIPROT, carnival_activity)

  nodes_kin_tf <- tibble::as_tibble(igraph::as_data_frame(graph_2, what = c('vertices'))) %>%
    dplyr::rename(gene_name = 'name') %>%
    dplyr::select(gene_name, UNIPROT, carnival_activity)

  nodes <- dplyr::bind_rows(nodes_rec_kin, nodes_kin_tf) %>%
    dplyr::distinct() %>%
    dplyr::arrange(gene_name)

  nodes <- dplyr::left_join(nodes, proteins_df, by = c('gene_name')) %>%
    dplyr::select(gene_name, carnival_activity, UNIPROT, mf, final_score, method) %>%
    dplyr::relocate(gene_name)

  # se sono duplicati nei due network prendo la media
  nodes1 <- nodes %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(carnival_activity = mean(carnival_activity))

  nodes <- nodes %>%
    dplyr::select(gene_name, UNIPROT, mf, final_score, method) %>%
    dplyr::distinct()

  nodes <- dplyr::left_join(nodes1, nodes, by = c('gene_name')) %>%
    dplyr::relocate(gene_name)

  # adding attributes
  nodes$discordant <- FALSE
  nodes$discordant[nodes$carnival_activity * nodes$final_score < 0] <- TRUE

  #edges df
  edges_rec_kin <- tibble::as_tibble(igraph::as_data_frame(graph_1, what = c('edges')))
  edges_kin_tf <- tibble::as_tibble(igraph::as_data_frame(graph_2, what = c('edges')))

  edges <- dplyr::bind_rows(edges_rec_kin, edges_kin_tf)

  union_graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = T)

  sp_object <- list(igraph_network = union_graph,
                    nodes_df = nodes,
                    edges_df = edges)

  if(files == TRUE){
    igraphToSif(union_graph, path_sif, 'sign')
    saveRDS(sp_object, path_rds)
  }

  return(sp_object)
}

#' default_CARNIVAL_options
#'
#' @param solver string specifying solver
#'
#' @return list of CARNIVAL options
#'
#' @examples
#' @export
default_CARNIVAL_options = function(solver = NULL){

  if(is.null(solver)){
    stop("please call default_CARNIVAL_options(solver = 'cplex') with a
        specific solver argument. Valid solvers are: 'lpSolve','cplex' or 'cbc'.")
  }
  solver_options = c("lpSolve","cplex","cbc")
  solver <- match.arg(solver,choices = solver_options)

  if(solver == "lpSolve"){
    opts = CARNIVAL::defaultLpSolveCarnivalOptions()
  }else if(solver == "cplex"){
    opts = CARNIVAL::defaultCplexCarnivalOptions()
  }else if(solver == "cbc"){
    opts = CARNIVAL::defaultCbcSolveCarnivalOptions()
  }
  opts$keepLPFiles = FALSE

  return(opts)
}

#' check_CARNIVAL_inputs
#'
#' @param source_df tibble with source nodes discretized among 1 and -1
#' @param target_df tibble with target nodes in a continuous range of activity
#' @param naive_network tibble with naive network in SIF format
#' @param proteins_df all inferred proteins df
#' @param organism string, organism you are usinggit
#'
#' @return TRUE if everything is correct
#'
#' @examples
check_CARNIVAL_inputs <- function(source_df, target_df,
                                  naive_network, proteins_df,
                                  organism){

  # check the input proteins
  if(!is.null(source_df)){
    stopifnot(is.data.frame(source_df))
    stopifnot(all(c("gene_name","mf","final_score") %in% names(source_df)))
    stopifnot(ncol(source_df)>4)
  }

  # check the input proteins
  stopifnot(is.data.frame(target_df))
  stopifnot(all(c("gene_name","mf","final_score") %in% names(target_df)))
  stopifnot(ncol(target_df)>4)

  # check naive network
  stopifnot(is.data.frame(naive_network))
  stopifnot(all(c("source","interaction","target" ) %in% names(naive_network)))
  stopifnot(ncol(naive_network)==3)

  # check meta informations for mapping attributes
  stopifnot(is.data.frame(proteins_df))
  stopifnot(all(c("gene_name","mf","final_score", "method") %in% names(proteins_df)))
  stopifnot(ncol(proteins_df)==5)

  # check organism
  stopifnot(typeof(organism) == 'character')
  stopifnot(organism %in% c('mouse', 'human'))

  return(TRUE)
}

#' run_carnival_and_create_graph
#'
#' @param source_df if NULL, inverse CARNIVAL,
#' otherwise, tibble with source nodes discretized among 1 and -1
#' @param target_df tibble with target nodes in a continuous range of activity
#' @param naive_network tibble with naive network in SIF format
#' @param carnival_options list of options returned by default_CARNIVAL_options
#' @param proteins_df all inferred proteins df
#' @param organism string, organism you are using
#' @param files boolean value, TRUE if you want output files
#' @param path_sif path of the sif output file of network
#' @param path_rds path of the rds output file of network
#' @param with_atlas Boolean value, default TRUE, if FALSE excludes Kinome Altas derived regulons
#' @param direct Boolean value, default FALSE, if TRUE uses only direct interactions
#' @param topbottom Boolean value, default FALSE, if TRUE optimization priorities sources
#'
#' @return list with igraph object, nodes df and edges df
#' @export
#'
#' @examples
run_carnival_and_create_graph <- function(source_df,
                                          target_df,
                                          naive_network,
                                          carnival_options,
                                          proteins_df,
                                          organism,
                                          topbottom = FALSE,
                                          with_atlas = TRUE,
                                          direct = FALSE,
                                          files = TRUE,
                                          path_sif = './optimized_network.sif',
                                          path_rds = './optimized_SP_oject.RDS'){

  if(is.null(source_df)){
    message(' ** Running inverse CARNIVAL (No perturbations) ** ')
    message('Credits to Prof. Julio Saez-Rodriguez. For more information visit: https://saezlab.github.io/CARNIVAL/ ')
    check_CARNIVAL_inputs(source_df = source_df,
                          target_df = target_df,
                          naive_network = naive_network,
                          proteins_df = proteins_df,
                          organism = organism)

    carnival_result <- CARNIVAL::runInverseCarnival(measurements = unlist(formatting_proteins_for_carnival(target_df)$t),
                                                    priorKnowledgeNetwork = unique(naive_network),
                                                    carnivalOptions = carnival_options)
  }else{
    if( topbottom == TRUE ){

      # Invert naive network
      naive_network <- naive_network %>%
        dplyr::rename('target' = 'source', 'source' = 'target') %>%
        dplyr::relocate(source, interaction, target)

      # Inverte source proteins with target proteins
      source_appo <- target_df
      target_df <- source_df
      source_df <- source_appo

    }

    message(' ** Running vanilla CARNIVAL (With perturbations) ** ')
    message('Credits to Prof. Julio Saez-Rodriguez. For more information visit: https://saezlab.github.io/CARNIVAL/ ')

    # check inputs
    check_CARNIVAL_inputs(source_df = source_df,
                          target_df = target_df,
                          naive_network = naive_network,
                          proteins_df = proteins_df,
                          organism = organism)

    # keep only present perturbation
    source_df <- keep_only_present_perturbation(source_df, naive_network)
    target_df <- keep_only_present_perturbation(target_df, naive_network)

    # naive_graph <- igraph::graph_from_data_frame(naive_network %>% dplyr::relocate(source, target))
    #
    # source_df_present <- target_df %>%
    #   dplyr::filter(gene_name %in% igraph::V(naive_graph)$name)

    # discretize initiators
    source_df_disc <- create_discretized_initiators_for_carnival(source_df)

    carnival_result <- CARNIVAL::runVanillaCarnival(perturbations = unlist(formatting_proteins_for_carnival(source_df_disc)$t),
                                                    measurements = unlist(formatting_proteins_for_carnival(target_df)$t),
                                                    priorKnowledgeNetwork = unique(naive_network),
                                                    carnivalOptions = carnival_options)
  }

  if(nrow(carnival_result$sifAll[[1]]) == 0){
    message('No network found for your experiment')
    return(NULL)
  }

  if( topbottom == TRUE ){

    carnival_result$weightedSIF <- carnival_result$weightedSIF %>% dplyr::rename('Node2' = 'Node1', 'Node1' = 'Node2') %>%
      dplyr::relocate(Node1, Sign, Node2)

  }

  # add attributes to the edges
  edges_df <- add_output_carnival_edges_attributes(carnival_result) %>%
    dplyr::filter(carnival_weight != 0)

  # keep only nodes not 0 or 0 involved in relations in the network
  nodes_df <- add_output_carnival_nodes_attributes(carnival_result,
                                                   proteins_df,
                                                   organism,
                                                   with_atlas = with_atlas,
                                                   direct = direct) %>%
    #keep nodes that have an activity OR that are 0 but are involved in some interactions
    # dplyr::filter(carnival_activity != 0 |
    #                 (carnival_activity == 0 & (gene_name %in% edges_df$source | gene_name %in% edges_df$target))) %>%
    dplyr::filter(gene_name %in% edges_df$source | gene_name %in% edges_df$target) %>%
    dplyr::distinct()

  CARNIVAL_igraph_network <- igraph::graph_from_data_frame(edges_df,
                                                           nodes_df,
                                                           directed = TRUE)

  SP_object <- list(igraph_network = CARNIVAL_igraph_network,
                    nodes_df = nodes_df,
                    edges_df = edges_df)

  if(files == TRUE){
    message(paste0('Creating optimized network file in tabular form in ', path_sif))
    igraphToSif(CARNIVAL_igraph_network, path_sif, 'sign')

    message(paste0('Creating an RDS object with SP optimization in ', path_rds))
    saveRDS(SP_object, path_rds)
  }

  # create carnival output as graph
  return(SP_object)
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
  formatted_nodes <- nodes %>% dplyr::select(gene_name, final_score = carnival_activity, mf, method)
  return(formatted_nodes)
}

#' expand_and_map_edges
#'
#' @param optimized_object SP object containing optimized graphs, nodes and edges of run_carnival_and_create_graph
#' @param organism #string, human or mouse
#' @param phospho_df # tibble of phosphoproteomics data
#' @param files # boolean value, TRUE if you want output files
#' @param path_sif # string of the path of output sif file (network)
#' @param path_rds # string of the path of RDS file (SP object)
#' @param with_atlas # Boolean value, FALSE default, if TRUE uses integrated regulons
#' @param direct  # Boolean value, FALSE default, if TRUE uses indirect interactions for mapping
#'
#' @return # list with igraph object, optimized nodes df and optimized edges df
#' @export
#'
#' @examples
expand_and_map_edges <- function(optimized_object,
                                 organism,
                                 phospho_df,
                                 files,
                                 with_atlas = FALSE,
                                 direct = FALSE,
                                 path_sif,
                                 path_rds){

  # Input of the function

  # Takes as input the CARNIVAL object
  # optimized_object = readRDS(paste0(carnival_dir, patient, '_union_carnival_object.RDS'))
  #
  # phos_df <- read_tsv(paste0(input_dir, 'phos/', patient, '_for_sp.tsv'))
  # phospho_df <- phos_df
  #
  # organism <- 'human'
  # with_atlas = TRUE
  #
  # files = TRUE
  # path_sif = paste0(carnival_dir, patient, '_union_carnival_validated.sif')
  # path_rds = paste0(carnival_dir, patient, '_union_carnival_validated.RDS')

  # Nodes df
  nodes_df <- optimized_object$nodes_df

  # Edges df
  edges_df <- optimized_object$edges_df %>%
    dplyr::rename('source' = 'from', 'target' = 'to', 'interaction' = 'sign') %>%
    dplyr::mutate_at('interaction', as.character)

  # output <- readRDS(paste0(carnival_dir, patient, '_union_carnival_object.RDS'))
  # optimized_graph_rds <- output$igraph_network
  # optimized_graph_sif <- output$edges_df
  # optimized_graph_sif$carnival_weight <- NULL
  # optimized_graph_sif <- optimized_graph_sif %>% rename( )

  # Derive the PKN to expand the edges_df
  # direct = TRUE because we want to select only direct phosphorylation events

  db <- choose_database_for_building(organism = organism,
                                     format = 'table',
                                     with_atlas = with_atlas,
                                     direct = FALSE) %>%
    dplyr::mutate_at('INTERACTION', as.character)


  # ========================================================================== #
  # Parse phosphoproteomics data
  # ========================================================================== #

  # Create substring of phosphopeptide sequence in phospho_df
  # for mapping to PKN, if necessary
  if(nchar(phospho_df$sequence_window[1]) != 15){
    center <- (nchar(phospho_df$sequence_window[1])+1)/2
    phospho_df <- phospho_df %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
  }else{ #if it is rename the column just to have a homogeneous variable
    phospho_df <- phospho_df %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }

  phospho_df <- phospho_df %>%
    dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(gene_name, '-', sequence_window_sub)) %>%
    dplyr::select(gene_name, aminoacid, position, PHOSPHO_KEY_GN_SEQ, difference, significant)

  # ========================================================================== #
  # Expand edges_df
  # ========================================================================== #
  # Join DB on network
  # taking from DBs the optimized interaction in order to expand interactions
  # on the phosphorylation events
  # and retrieve information about the mechanism, DIRECT or INDIRECT

  edges_df_new <- dplyr::left_join(edges_df, db, by = c('source' = 'ENTITYA',
                                                        'interaction' = 'INTERACTION',
                                                        'target' = 'ENTITYB')) %>%
    dplyr::distinct()

  # Map the expanded edges to experimental data
  edges_df_new <- dplyr::left_join(edges_df_new,
                                   phospho_df %>%
                                     dplyr::select(-c('gene_name')),
                                   by = c('PHOSPHO_KEY_GN_SEQ'))

  # Check if the mapped is coherent with the proteins' activity status
  dplyr::inner_join(dplyr::inner_join(edges_df_new,
                                      nodes_df %>% dplyr::select(gene_name, act_source = carnival_activity),
                                      by = c('source' = 'gene_name')), nodes_df %>% dplyr::select(gene_name, act_target = carnival_activity),
                    by = c('target' = 'gene_name')) -> edges_df_new_activity

  edges_df_new_activity %>%
    dplyr::mutate(keep = ifelse(as.numeric(interaction) * difference * act_target > 0, '+', NA)) -> edges_df_new_activity

  # Rename and add some columns as flags like
  # is_quantified, is_significant (derived from significant in phospho_df),
  # phosphosite_value is the abundance in experimental data
  edges_df_new <- edges_df_new %>%
    dplyr::mutate(is_quantified = ifelse(!is.na(edges_df_new$difference) & !is.na(edges_df_new_activity$keep), TRUE, FALSE),
                  is_significant = ifelse(!is.na(significant) & !is.na(edges_df_new_activity$keep), TRUE, FALSE)) %>%
    dplyr::rename('phosphosite_value' = 'difference')

  # Remove phosphosite values mapped but not coherent with protein activities
  edges_df_new <- edges_df_new %>%
    dplyr::mutate(phosphosite_value = ifelse(is_quantified == TRUE, phosphosite_value, NA),
                  aminoacid = ifelse(is_quantified == TRUE, aminoacid, NA),
                  position = ifelse(is_quantified == TRUE, position, NA))

  # Reorganize the edges_df
  edges_df_new <- edges_df_new %>%
    dplyr::select(source, target, sign = interaction, carnival_weight, mechanism = MECHANISM,
                  residue = RESIDUE,sequence = SEQUENCE, is_significant, phosphosite_value, DIRECT) %>%
    dplyr::distinct() %>%
    dplyr::mutate(residue = ifelse(is_significant == TRUE, paste0(residue, '*'), residue),
                  phosphosite_value = ifelse(is_significant == TRUE, paste0(round(phosphosite_value, 2), '*'), round(phosphosite_value,2)))

  # For rows mapped in the phospho_df
  # collapse the same edge mapping to different sites
  edges_df_new_with_amino <- edges_df_new %>%
    dplyr::filter(!is.na(phosphosite_value)) %>%
    dplyr::group_by(source, target, sign) %>%
    dplyr::reframe(aminoacid = paste0(residue,collapse = ';'),
                   FC = paste0(as.character(phosphosite_value), collapse = ';'),
                   direct = paste0(unique(DIRECT), collapse = ';'),
                   mechanism = paste0(unique(mechanism), collapse = ';')) %>%
    dplyr::mutate_at('direct', as.logical)

  edges_df_new_final <- dplyr::left_join(edges_df_new %>%
                                           dplyr::select(source, target, sign, carnival_weight, direct = DIRECT, mechanism) %>%
                                           dplyr::distinct(),
                                         edges_df_new_with_amino,
                                         by = c('source', 'target', 'sign', 'direct', 'mechanism'))

  edges_df_new_final <- edges_df_new_final %>%
    dplyr::mutate(is_quantified = ifelse(is.na(FC), FALSE, TRUE),
                  is_significant = ifelse(grepl('\\*', FC), TRUE, FALSE))

  # When there are duplicated interactions because of different mechanism
  # paste them all together

  edges_df_new_final <- edges_df_new_final %>% dplyr::group_by(source, target, sign, carnival_weight) %>%
    dplyr::reframe(direct = paste0(unique(direct), collapse = ';'),
                   mechanism = paste0(na.omit(mechanism), collapse = ';'),
                   aminoacid = paste0(na.omit(aminoacid), collapse = ';'),
                   is_quantified = paste0(unique(is_quantified), collapse = ';'),
                   is_significant = paste0(unique(is_significant), collapse = ';'),
                   FC = paste0(na.omit(FC), collapse = ';'))

  # Create igraph object
  CARNIVAL_graph <- igraph::graph_from_data_frame(edges_df_new_final, nodes_df, directed = TRUE)

  SP_object <- list(igraph_network = CARNIVAL_graph,
                    nodes_df = nodes_df,
                    edges_df = edges_df_new_final)

  if(files == TRUE){
    message(paste0('Creating validated optimized network file in tabular form in ', path_sif))
    igraphToSif(CARNIVAL_graph, path_sif, 'sign')

    message(paste0('Creating an RDS object with SP validated optimization in ', path_rds))
    saveRDS(SP_object, path_rds)
  }

  # create carnival output as graph
  return(SP_object)
}

#' keep_only_present_perturbation
#'
#' @param source_df dataframe with receptors
#' @param naive_network dataframe with naive network interactions
#'
#' @return dataframe with receptors present in the naive network
#' @export
#'
#' @examples
keep_only_present_perturbation <- function(source_df, naive_network){

  naive_graph <- igraph::graph_from_data_frame(naive_network %>%
                                                 dplyr::relocate(source, target))

  source_df_present <- source_df %>%
    dplyr::filter(gene_name %in% igraph::V(naive_graph)$name)

  return(source_df_present)
}



