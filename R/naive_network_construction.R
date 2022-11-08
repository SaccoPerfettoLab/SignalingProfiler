# scripts for building a naive network containing all inferred proteins

# function which converts an igraph object in a sif file
#' igraphToSif
#'
#' @param inGraph igraph object
#' @param outfile output file in sif format
#' @param edgeLabel label to use as edge
#'
#' @return igraph object in sif format
#'
#' @examples
igraphToSif <- function(inGraph, outfile="output.sif", edgeLabel="label") {

  if(file.exists(outfile)){
    command1 <- paste0('rm ', outfile)
    system(command1)
  }

  singletons <- as.list(igraph::get.vertex.attribute(inGraph, "name"))
  edgeList <- igraph::get.edgelist(inGraph, names=FALSE)
  nodeNames <- igraph::get.vertex.attribute(inGraph, "name")
  edgeAttribute <- igraph::get.edge.attribute(inGraph, edgeLabel)
  numE <- igraph::ecount(inGraph)

  for (i in 1:numE) {
    node1 <- edgeList[i,1]
    node2 <- edgeList[i,2]
    singletons <- singletons[which(singletons != nodeNames[node1])]
    singletons <- singletons[which(singletons != nodeNames[node2])]
    write(paste0(nodeNames[node1], "\t", edgeAttribute[i], "\t", nodeNames[node2], "\n"),
          outfile, append = TRUE)
  }

  for (single in singletons) {
    write(paste0(single, "\n"),
          outfile, append = TRUE)
  }
}


##### --------- #####


#' choose_database_for_building
#'
#' @param organism 'human' or 'mouse''
#' @param format 'igraph', 'table'
#'
#' @return built-in database as igraph object or tibble
#' @export
#'
#' @examples
choose_database_for_building <- function(organism ,format){
  # db choosing
  if(organism == 'mouse'){
    if(format == 'igraph'){
      pkn <- db_mouse
    }else if(format == 'table'){
      pkn <- PKN_mouse
    }else{error('Please provide a valid format among igraph or table')}
  }else if(organism == 'human'){
    if(format == 'igraph'){
      pkn <- db_human
    }else if(format == 'table'){
      pkn <- PKN_human
    }else{error('Please provide a valid format among igraph or table')}
  }else{error('Please provide a valid organism')}
  return(pkn)
}

##### --------- #####

#' coverage_of_inferred_proteins_in_db
#'
#' @param prediction_output concatenation of inferred proteins from omics data
#' @param organism 'human' or 'mouse'
#' @param report boolean value, if TRUE returns a file containing
#' the report of the proteins' coverage analysis, if FALSE print it
#'
#' @return print information about coverage or annotates them in a file
#' @export
#'
#' @examples
coverage_of_inferred_proteins_in_db <- function(prediction_output, organism, report){

  # organism <- 'mouse'
  # report <- TRUE
  pkn <- choose_database_for_building(organism, 'igraph')

  if(report == TRUE){sink('report.txt')}

  for(m in c('tf', 'kin', 'phos', 'other')){

    #m <- 'tf'
    prot_subset <- prediction_output %>%
      dplyr::filter(mf == m)

    found <- length(igraph::V(pkn)[ENTITY %in% prot_subset$gene_name])
    total <- nrow(prot_subset)

    sentence <- paste0('For ', m, ' molecular function have been found ', found, ' proteins out of ', total)
    message(sentence)
    if(report == TRUE){cat(paste0(sentence, '\n'))}
  }
  if(report == TRUE){sink()}
}

##### --------- #####

#' one_layer_naive_network
#'
#' @param starts_gn gene names of starting nodes (e.g. receptors)
#' @param targets_gn gene names of target nodes (e.g. transcription factors)
#' @param organism string, 'human' or 'mouse'
#' @param max_length integer, 1 to 4 for max_length connecting start to end
#' @param rds_path path of network rds file
#' @param sif_path path of network sif file

#'
#' @return naive network
#' @export
#'
#' @examples
one_layer_naive_network <- function(starts_gn, targets_gn, organism, max_length,
                                    rds_path = 'one_layer_naive.RDS',
                                    sif_path = 'one_layer_naive.sif'){
  message('One layer: shortest paths from receptor(s) to all proteins')
  all_paths_df <- get_all_shortest_path_custom(starts_gn, targets_gn, organism, max_length)
  network <- create_graph_from_paths(all_paths_df, organism)

  # set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[ENTITY %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[ENTITY %in% targets_gn]$mf_naive <- 'target'

  # save files
  message(paste0('Writing in ', getwd(), ' sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')

  return(network)
}

#' two_layer_naive_network
#'
#' @param starts_gn gene names of starting nodes
#' @param intermediate_gn gene names of intermediate nodes (e.g. kins and phos)
#' @param targets_gn gene names of target nodes (e.g. transcription factors)
#' @param organism string, 'human' or 'mouse'
#' @param max_length_1 max_length of shortest path from start to intermediates
#' @param max_length_2 max_length of shortest path from intermediates to targets
#' @param rds_path path of network rds file
#' @param sif_path path of network sif file
#'
#' @return naive network
#' @export
#'
#' @examples
two_layer_naive_network <- function(starts_gn, intermediate_gn, targets_gn,
                                    organism, max_length_1, max_length_2,
                                    rds_path = 'two_layer_naive.RDS',
                                    sif_path = 'two_layer_naive.sif'){


  message(paste0('First layer: ', max_length_1, ' length paths from starting nodes to intermediates'))
  all_paths_layer1_df <- get_all_shortest_path_custom(starts_gn, intermediate_gn, organism, max_length_1)

  message(paste0('Second layer: ', max_length_2, ' length paths from intermediates nodes to targets'))
  all_paths_layer2_df <- get_all_shortest_path_custom(intermediate_gn, targets_gn, organism, max_length_2)

  all_paths_df <- dplyr::bind_rows(all_paths_layer1_df, all_paths_layer2_df)
  network <- create_graph_from_paths(all_paths_df, organism)

  # set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[ENTITY %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[ENTITY %in% intermediate_gn]$mf_naive <- 'intermediate'
  igraph::V(network)[ENTITY %in% targets_gn]$mf_naive <- 'target'

  # save files
  message(paste0('Writing in ', getwd(), ' sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')

  return(network)
}

#' three_layer_naive_network
#'
#' @param starts_gn gene names of starting nodes
#' @param intermediate1_gn gene names of intermediates1 nodes
#' @param intermediate2_gn gene names of intermediates2 nodes
#' @param targets_gn gene names of targets nodes
#' @param organism string, 'human' or 'mouse'
#' @param max_length_1 max_length of shortest path from start to intermediates1
#' @param max_length_2 max_length of shortest path from intermediates1 to intermediates2
#' @param max_length_3  max_length of shortest path from intermediates2 to targets
#' @param rds_path path of network rds file
#' @param sif_path path of network sif file
#' @param both_intermediates boolean, TRUE if you want both intermediates as the starting point of the third layer,
#' FALSE if you want to use only intermediates2 (not suggested)
#'
#' @return naive network
#' @export
#'
#' @examples
#'
three_layer_naive_network <- function(starts_gn, intermediate1_gn, intermediate2_gn,
                                      targets_gn, organism,
                                      max_length_1, max_length_2, max_length_3,
                                      both_intermediates = TRUE,
                                      rds_path = 'three_layer_naive.RDS',
                                      sif_path = 'three_layer_naive.sif'){

  # create network
  message(paste0('First layer: ', max_length_1, ' length paths from starting nodes to intermediates1'))
  all_paths_layer1_df <- get_all_shortest_path_custom(starts_gn, intermediate1_gn, organism, max_length_1)

  message(paste0('Second layer: ', max_length_2, ' length paths from intermediates1 nodes to intermediates'))
  all_paths_layer2_df <- get_all_shortest_path_custom(intermediate1_gn, targets_gn, organism, max_length_2)

  # third layer
  if(both_intermediates == TRUE){
    message(paste0('Third layer: ', max_length_3, ' length paths from intermediates2 nodes to targets'))
    all_paths_layer3_df <- get_all_shortest_path_custom(c(intermediate1_gn, intermediate2_gn),
                                                        targets_gn, organism, max_length_3)

  }else if(both_intermediates == FALSE){
    message(paste0('Third layer: ', max_length_3, ' length paths from intermediates2 nodes to targets'))
    all_paths_layer3_df <- get_all_shortest_path_custom(intermediate2_gn, targets_gn, organism, max_length_3)
  }else{
    error('Please provide a boolean value for both_intermediates parameter')
  }


  all_paths_df <- dplyr::bind_rows(all_paths_layer1_df, all_paths_layer2_df, all_paths_layer3_df)
  network <- create_graph_from_paths(all_paths_df, organism)

  # set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[ENTITY %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[ENTITY %in% intermediate1_gn]$mf_naive <- 'intermediate1'
  igraph::V(network)[ENTITY %in% intermediate2_gn]$mf_naive <- 'intermediate2'
  igraph::V(network)[ENTITY %in% targets_gn]$mf_naive <- 'target'

  # save files
  message(paste0('Writing in ', getwd(), ' sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')

  return(network)
}

#' prepare_carnival_input
#'
#' @param naive_network naive network connecting inferred proteins
#' @param prediction_output inferred proteins from experimental data
#' @param recept_list list of receptors with their desired activity or NULL if you don't have any
#'
#' @return inferred protein filtered for presence in naive network
#' @export
#'
#' @examples
prepare_carnival_input <- function(naive_network, prediction_output, recept_list){

  # filter prediction
  prediction_output_filt <- prediction_output %>%
    dplyr::filter(gene_name %in% igraph::V(naive_network)$ENTITY) %>%
    dplyr::arrange(gene_name)

  if(!is.null(recept_list)){
    # create receptor list
    recept_df <- tibble::tibble(gene_name = names(recept_list),
                                mf = 'rec',
                                method = 'user',
                                final_score = unlist(recept_list))

    recept_df <- convert_gene_name_in_uniprotid(recept_df, 'mouse') %>%
      dplyr::relocate('UNIPROT')

    # unify elements
    prediction_output_filt_rec <- dplyr::bind_rows(recept_df, prediction_output_filt)

    return(prediction_output_filt_rec)

  }
  # create receptors tibble
  # recept_list <- list('Flt3' = 1)
  # prediction_output <- toy_prot_activity_df
  # naive_network <- one_layer_toy

  return(prediction_output_filt)
}


