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

  sink(outfile)

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
    cat(nodeNames[node1], "\t", edgeAttribute[i], "\t", nodeNames[node2], "\n")
  }

  for (single in singletons) {
    cat(single, "\n")
  }

  sink()
}

##### --------- #####

#' get all edge ids
#'
#' @param all_shortest_paths takes in input all shortest paths of igraph command
#' @param db fb used for the analysis
#'
#' @return all edges id along the shortest path
#'
#' @examples
get_all_edge_ids <- function(all_shortest_paths, db){
  shortest_paths <- lapply(all_shortest_paths$res,
                           FUN = function(x){names(unlist(x))})

  edges_ids_list <- c()

  for(i_path in c(1:length(shortest_paths))){
    # unlisting each path
    # e.g ['A', 'B', 'C', 'D]
    path = unlist(shortest_paths[i_path])

    # getting for each couple (AB, BC, CD) all edges ids
    all_couple_edges_ids <- c()

    if(length(path) > 1){
      for( node_idx in c( 1:(length(path)-1)) ){
        # selecting node A
        node1 = path[ node_idx ]
        # selecting node B
        node2 = path[ node_idx + 1 ]

        out_edges <- igraph::incident_edges( db, node1, mode = 'out' )
        in_edges <- igraph::incident_edges( db, node2, mode = 'in' )
        couple_edge_ids <- intersect( out_edges[[node1]],
                                      in_edges[[node2]] )

        all_couple_edges_ids <- c( all_couple_edges_ids, couple_edge_ids )
      }

      edges_ids_list <- c(edges_ids_list, all_couple_edges_ids)
    }

  }
  return(unique(edges_ids_list))
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

#' get_all_shortest_path_custom
#'
#' @param start_nodes_gn vector of gene names of starting nodes
#' @param target_nodes_gn vector of gene names of target nodes
#' @param db igraph object representing chosen built in databas
#' @param path_length 'one' or 'shortest' to choose the length of the path connecting
#' proteins: 1 step or the shortest
#'
#' @return a vector with edges ids along the shortest paths
#' @export
#'
#' @examples
get_all_shortest_path_custom <- function(start_nodes_gn, target_nodes_gn, db, path_length){
  message('Warning message from developer: this analysis can take even hours according to the number of nodes')
  # start_nodes_gn <- kinphos_gn
  # target_nodes_gn <- tfs_gn
  # path_length <- 'shortest'

  all_edges_ids <- c()

  for(protein in start_nodes_gn){
    #protein <- 'FLT3'

    source2target <- igraph::get.all.shortest.paths(graph = db,
                                                    from = igraph::V(db)[ENTITY %in% c(protein)],
                                                    to = igraph::V(db)[ENTITY %in% target_nodes_gn],
                                                    mode = 'out')
    if(path_length == 'one'){

      if(length(source2target$res) == 0){
        next
      }else{
        for (i in c(1:length(source2target$res))){
          if (length(source2target$res[[i]]) == 2){
            node1 = names(unlist(source2target$res[[i]]))[1]
            # selecting node B
            node2 = names(unlist(source2target$res[[i]]))[2]

            out_edges <- igraph::incident_edges( db, node1, mode = 'out' )
            in_edges <- igraph::incident_edges( db, node2, mode = 'in' )
            edge_ids <- intersect( out_edges[[node1]],
                                   in_edges[[node2]] )
            all_edges_ids <- c(all_edges_ids, edge_ids)
          }
        }
      }

    }else if(path_length == 'shortest'){
      if(length(source2target$res) == 0){
        next
      }else{
        all_edges_ids <- c(all_edges_ids, get_all_edge_ids(source2target, db))
      }

    }else{
      error('Please provide a valid type of analysis')
    }
  }
  return(all_edges_ids)
}

#' one_layer_naive_network
#'
#' @param receptors_gn gene names of receptors
#' @param targets_gn gene names of target nodes (transcription factors)
#' @param db built in causal relation database
#'
#' @return naive network
#' @export
#'
#' @examples
one_layer_naive_network <- function(receptors_gn, targets_gn, db,
                                    rds_path = 'one_layer_naive.RDS',
                                    sif_path = 'one_layer_naive.sif'){
  message('One layer: shortest paths from receptor(s) to all proteins')
  edges_ids <- get_all_shortest_path_custom(receptors_gn, targets_gn, db, 'shortest')
  network <- igraph::subgraph.edges(db, edges_ids)

  # set node attributes
  igraph::V(network)$mf <- 'unknown'
  igraph::V(network)[ENTITY %in% receptors_gn]$mf <- 'recept'
  igraph::V(network)[ENTITY %in% targets_gn]$mf <- 'target'

  # save files
  message(paste0('Writing in ', getwd(), 'sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')

  return(network)
}

#' two_layer_naive_network
#'
#' @param receptors_gn gene names of receptors
#' @param kinphos_gn gene names of kinases and phosphatases
#' @param tfs_gn gene names of transcription factors
#' @param db  built in causal relation database
#'
#' @return naive network
#' @export
#'
#' @examples
two_layer_naive_network <- function(receptors_gn, kinphos_gn, tfs_gn, db,
                                    rds_path = 'two_layer_naive.RDS',
                                    sif_path = 'two_layer_naive.sif'){
  # receptors_gn <- c('FLT3')
  # kinphos_gn <- toy_kin$gene_name
  # tfs_gn <- toy_tf$gene_name

  message('First layer: shortest paths from receptor(s) to kinases, phosphatases and others')
  rec_kins_edges_ids <- get_all_shortest_path_custom(receptors_gn, kinphos_gn, db, 'shortest')

  message('Second layer: shortest paths from kinases, phosphtases and others to transcription factors')
  kins_tfs_edges_ids <-  get_all_shortest_path_custom(kinphos_gn, tfs_gn, db, 'shortest')
  network <- igraph::subgraph.edges(db, c(rec_kins_edges_ids, kins_tfs_edges_ids))

  # set node attributes
  igraph::V(network)$mf <- 'unknown'
  igraph::V(network)[ENTITY %in% receptors_gn]$mf <- 'recept'
  igraph::V(network)[ENTITY %in% kinphos_gn]$mf <- 'kinphos'
  igraph::V(network)[ENTITY %in% tfs_gn]$mf <- 'tf'

  # save files
  message(paste0('Writing in ', getwd(), 'sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')


  return(network)
}

#' three_layer_naive_network
#'
#' @param receptors_gn gene names of receptors
#' @param kinphos_gn gene names of kinases and phosphatases
#' @param subs_gn gene names of substrates
#' @param tfs_gn gene names of transcription factors
#' @param db  built in causal relation database
#'
#' @return naive network
#' @export
#'
#' @examples
three_layer_naive_network <- function(receptors_gn, kinphos_gn, subs_gn, tfs_gn, db,
                                      rds_path = 'three_layer_naive.RDS',
                                      sif_path = 'three_layer_naive.sif'){

  # create network
  message('First layer: shortest paths from receptor(s) to kinases and phosphatases')
  rec_kins_edges_ids <- get_all_shortest_path_custom(receptors_gn, kinphos_gn, db, 'shortest')

  message('Second layer: one step paths from kinases and phosphatases to others')
  kins_subs_edges_ids <- get_all_shortest_path_custom(kinphos_gn, subs_gn, db, 'one')

  message('Third layer: shortest paths from kinases, phosphatases and others to transcription factors')
  kins_tfs_edges_ids <-  get_all_shortest_path_custom(c(kinphos_gn, subs_gn), tfs_gn, db, 'shortest')

  network <- igraph::subgraph.edges(db, c(rec_kins_edges_ids, kins_subs_edges_ids, kins_tfs_edges_ids))

  # set node attributes
  igraph::V(network)$mf <- 'unknown'
  igraph::V(network)[ENTITY %in% receptors_gn]$mf <- 'recept'
  igraph::V(network)[ENTITY %in% kinphos_gn]$mf <- 'kinphos'
  igraph::V(network)[ENTITY %in% subs_gn]$mf <- 'other'
  igraph::V(network)[ENTITY %in% tfs_gn]$mf <- 'tf'

  # save files
  message(paste0('Writing in ', getwd(), 'sif and RDS file of the naive network'))
  saveRDS(network, rds_path)
  igraphToSif(network, outfile = sif_path, edgeLabel = 'SIGN')

  return(network)
}

#' filter_inferred_protein_for_presence_in_naive_network
#'
#' @param naive_network naive network connecting inferred proteins
#' @param prediction_output inferred proteins from experimental data
#'
#' @return inferred protein filtered for presence in naive network
#' @export
#'
#' @examples
filter_inferred_protein_for_presence_in_naive_network <- function(naive_network, prediction_output){
  prediction_output_filt <- prediction_output %>%
    dplyr::filter(gene_name %in% igraph::V(naive_network)$ENTITY) %>%
    dplyr::arrange(gene_name)
  return(prediction_output_filt)
}


