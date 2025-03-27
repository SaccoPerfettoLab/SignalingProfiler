#' Get Prior Knowledge Network in multiple formats
#'
#' Returns the built-in Prior Knowledge Network (PKN) from SIGNOR and PhosphoSitePlus
#' in different formats, with optional kinase-substrate predicted interactions
#' from Kinome Atlases (PMIDs: 36631611, 38720073).

#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param format Character string, format of the output; choose between `"igraph"` or `"table"`.
#' @param with_atlas Logical, whether to integrate inferred kinase-substrate
#'   interactions from the Ser/Thr and Tyr Kinome Atlas (PMIDs: 36631611, 38720073); default: `FALSE`.
#' @param direct Logical, whether to retain only direct interactions among biological entities; default `TRUE`.
#' @param custom Logical, whether the user wants to provide a custom PKN.
#'   If `TRUE`, the function will stop with an error message; default `FALSE`.
#'
#' @return The prior knowledge network in the specified format:
#'   - `"igraph"` → An `igraph` object.
#'   - `"table"` → A `data.frame`.

#' @details
#' The function selects the appropriate dataset based on user-defined parameters.
#' If `custom = TRUE`, the function stops since custom PKNs must be handled separately.
#' For `"mouse"`, kinase-substrate interactions (`with_atlas = TRUE`) are not supported.
#'
#' @export
#'
#' @examples
#' # Retrieve human PKN in table format with direct interactions
#' choose_database_for_building(organism='human',with_atlas=FALSE,direct=TRUE,custom=FALSE,format='table')
#'
#' # Retrieve mouse PKN as an igraph object (without kinase-substrate interactions)
#' PKN_mouse <- choose_database_for_building(organism = "mouse", format = "igraph")
choose_database_for_building <- function(organism,
                                         with_atlas = FALSE,
                                         direct = FALSE,
                                         custom = FALSE,
                                         format) {

  if (custom) {
    stop("This function is not available for custom PKNs.")
  }

  if (!organism %in% c("human", "mouse")) {
    stop("Please provide a valid organism: 'human' or 'mouse'.")
  }

  if (!format %in% c("igraph", "table")) {
    stop("Please provide a valid format: 'igraph' or 'table'.")
  }

  if (organism == "mouse" && with_atlas) {
    stop("If organism is 'mouse', 'with_atlas' must be FALSE.")
  }

  ## built-in object selection based on name
  dataset_prefix <- ifelse(organism == "mouse", "mouse", "human")
  atlas_suffix <- ifelse(with_atlas, "_atlas", "")
  direct_suffix <- ifelse(direct, "_dir", "_ind")

  dataset_name <- paste0(ifelse(format == "igraph", "db_", "PKN_"),
                         dataset_prefix, atlas_suffix, direct_suffix)

  e <- new.env()
  data(list = dataset_name, package = "SignalingProfiler", envir = e)

  return(e[[dataset_name]])
}

#' Generate a One-Layer Naive Network
#'
#' Constructs a one-layer naive network by identifying paths up to `max_length`
#' between starting (`starts_gn`) and target (`targets_gn`) nodes in a
#' Prior Knowledge Network (`PKN_table`).
#' The `connect_all = TRUE` parameter ensures that intermediate nodes in independent
#' shortest paths are connected by one-step interactions.
#'
#' @param starts_gn Character vector, Gene Symbols of starting nodes (e.g., receptors).
#' @param targets_gn Character vector, Gene Symbols of target nodes (e.g., transcription factors).
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param max_length Integer, maximum path length between start and end nodes (1 to 4).
#' @param rds_path Character, path to save the naive network as an RDS file. Default: `"one_layer_naive.RDS"`.
#' @param sif_path Character, path to save the naive network in SIF format. Default: `"one_layer_naive.sif"`.
#' @param connect_all Logical, whether to connect intermediate nodes of independent shortest paths.
#'   Default is `TRUE`.
#' @param files Logical, whether to save the naive network as RDS and SIF files. Default is `TRUE`.
#'
#' @return An `igraph` object representing the naive network.
#'
#' @details
#' The function:
#' - Identifies all shortest paths (up to `max_length`) between `starts_gn` and `targets_gn`.
#' - Creates an `igraph` object representing the naive network.
#' - Optionally saves the network as an RDS file (`rds_path`) and in SIF format (`sif_path`).
#'
#' @export
#'
#' @examples
#' PKN = choose_database_for_building(organism='mouse',with_atlas=FALSE,direct=TRUE,custom=FALSE,format='table')
#' one_layer_naive_network(starts_gn=c('Flt3'), targets_gn=c('Mapk1', 'Mapk3','Stat5a'),PKN_table=PKN,max_length=4)
one_layer_naive_network <- function(starts_gn,
                                    targets_gn,
                                    PKN_table,
                                    max_length,
                                    files = TRUE,
                                    rds_path = 'one_layer_naive.RDS',
                                    sif_path = 'one_layer_naive.sif',
                                    connect_all = TRUE) {
  message(paste0("One-layer network: Finding shortest paths (maximum = ", max_length,
                 "), from receptors to all proteins."))

  all_paths_df <-
    get_all_shortest_path_custom(starts_gn, targets_gn, PKN_table, max_length)

  network <-
    create_graph_from_paths(all_paths_df, PKN_table, connect_all = connect_all)

  ## set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[name %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[name %in% targets_gn]$mf_naive <- 'target'

  ## save files
  if (files) {
    saveRDS(network, rds_path)
    igraphToSif(network, outfile = sif_path, edgeLabel = 'INTERACTION')
  }

  return(network)
}

#' Generate a Two-Layers Naive Network
#'
#' Constructs a two-layer naive network by identifying paths up to `max_length_1`
#' between `starts_gn` and `intermediate_gn`, and up to `max_length_2` between
#' `intermediate_gn` and `targets_gn` in a Prior Knowledge Network (`PKN_table`).
#' The `connect_all = TRUE` parameter ensures that intermediate nodes in independent
#' shortest paths are connected by one-step interactions.
#'
#' @param starts_gn Character vector, Gene Symbols of starting nodes (e.g., receptors).
#' @param intermediate_gn Character vector, Gene Symbols of intermediate nodes (e.g., kinases, phosphatases).
#' @param targets_gn Character vector, Gene Symbols of target nodes (e.g., transcription factors).
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param max_length_1 Integer, maximum path length (1 to 4) from start to intermediate nodes. Default: `3`.
#' @param max_length_2 Integer, maximum path length (1 to 4) from intermediate to target nodes. Default: `4`.
#' @param rds_path Character, file path to save the naive network as an RDS file. Default: `"two_layer_naive.RDS"`.
#' @param sif_path Character, file path to save the naive network in SIF format. Default: `"two_layer_naive.sif"`.
#' @param connect_all Logical, whether to connect intermediate nodes from independent shortest paths. Default: `TRUE`.
#' @param files Logical, whether to save the naive network in RDS and SIF formats. Default: `TRUE`.
#'
#' @return An `igraph` object representing the two-layer naive network.
#'
#' The function:
#' - Identifies all shortest paths (up to `max_length_1`) between `starts_gn` and `intermediate_gn`.
#' - Identifies all shortest paths (up to `max_length_2`) between `intermediate_gn` and `targets_gn`.
#' - Constructs an `igraph` object representing the naive network.
#' - Optionally saves the network as an RDS file (`rds_path`) and in SIF format (`sif_path`).
#'
#' @export
#'
#' @examples
#' PKN <- choose_database_for_building(organism='mouse',format='table')
#' two_layer_naive_network(
#'   starts_gn=c('Flt3'),
#'   intermediate_gn=c('Mapk1', 'Mapk3'),
#'   targets_gn=c('Stat5a'),
#'   PKN_table=PKN,
#'   max_length_1=3,
#'   max_length_2=4)
#'
two_layer_naive_network <- function(starts_gn,
                                    intermediate_gn,
                                    targets_gn,
                                    PKN_table,
                                    max_length_1 = 3,
                                    max_length_2 = 4,
                                    files = TRUE,
                                    rds_path = 'two_layer_naive.RDS',
                                    sif_path = 'two_layer_naive.sif',
                                    connect_all = TRUE) {

  message(paste0("First layer: Identifying paths up to ", max_length_1,
                 " steps from start to intermediate nodes."))

  all_paths_layer1_df <- get_all_shortest_path_custom(starts_gn,
                                                      intermediate_gn,
                                                      PKN_table,
                                                      max_length_1)

  message(paste0("Second layer: Identifying paths up to ", max_length_2,
                 " steps from intermediate to target nodes."))
  all_paths_layer2_df <-
    get_all_shortest_path_custom(intermediate_gn,
                                 targets_gn,
                                 PKN_table,
                                 max_length_2)

  ## combine and ensure uniqueness
  all_paths_df <- dplyr::bind_rows(all_paths_layer1_df,
                                   all_paths_layer2_df) %>%
    dplyr::distinct()

  ## create the network graph
  network <- create_graph_from_paths(all_paths_df,
                                     PKN_table,
                                     connect_all = connect_all)

  ## set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[name %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[name %in% intermediate_gn]$mf_naive <-
    'intermediate'
  igraph::V(network)[name %in% targets_gn]$mf_naive <- 'target'

  ## save files
  if (files) {
    saveRDS(network, rds_path)
    igraphToSif(network, outfile = sif_path, edgeLabel = 'INTERACTION')
  }

  return(network)
}

#' Generate a Three-Layer Naive Network
#'
#' Generate a three-layered naive network by looking for all paths
#' equal or shorter than `max_length_1`, `max_length_2`, and `max_length_3` in
#' the Prior Knowledge Network (`PKN_table`).

#' - The first layer links starting nodes (`starts_gn`) to intermediate nodes (`intermediate1_gn`).
#' - The second layer links `intermediate1_gn` to themselves and additional intermediate nodes (`intermediate2_gn`).
#'   If `keep_only_connected = TRUE`, only `intermediate1_gn` nodes that were connected in the first layer are considered.
#' - The third layer links intermediate nodes to target nodes (`targets_gn`).
#'
#' If `both_intermediates = TRUE`, both `intermediate1_gn` and `intermediate2_gn` are used as
#' starting points in the third layer; otherwise, only `intermediate2_gn` is considered.
#'
#' The `connect_all = TRUE` parameter ensures that intermediate nodes from independent shortest paths
#' are connected by one-step interactions.
#'
#' @param starts_gn Character vector, Gene Symbols of starting nodes (e.g., receptors).
#' @param intermediate1_gn Character vector, Gene Symbols of first set of intermediate nodes (e.g., kinases, phosphatases).
#' @param intermediate2_gn Character vector, Gene Symbols of second set of intermediate nodes (e.g., other phosphorylated proteins).
#' @param targets_gn Character vector, Gene Symbols of target nodes (e.g., transcription factors).
#' @param PKN_table Data frame, Prior Knowledge Network of causal (signed and oriented) interactions.
#' @param max_length_1 Integer, maximum path length (1 to 4) from starting nodes to `intermediate1_gn`. Default: `3`.
#' @param max_length_2 Integer, maximum path length (1 to 4) from `intermediate1_gn` to `intermediate2_gn`. Default: `1`.
#' @param max_length_3 Integer, maximum path length (1 to 4) from intermediate nodes to `targets_gn`. Default: `4`.
#' @param rds_path Character, file path to save the naive network as an RDS file. Default: `"three_layer_naive.RDS"`.
#' @param sif_path Character, file path to save the naive network in SIF format. Default: `"three_layer_naive.sif"`.
#' @param both_intermediates Logical, if `TRUE`, both `intermediate1_gn` and `intermediate2_gn` are used
#'   as the starting point of the third layer; otherwise, only `intermediate2_gn` is used. Default: `TRUE`.
#' @param keep_only_connected Logical, if `TRUE`, only `intermediate1_gn` nodes connected in the first layer are
#'   considered in the second layer. Default: `FALSE`.
#' @param connect_all Logical, whether to connect intermediate nodes of independent shortest paths. Default: `TRUE`.
#' @param files Logical, whether to save the naive network as RDS and SIF files. Default: `TRUE`.
#'
#' @return An `igraph` object representing the three-layer naive network.
#'
#' @details
#' The function:
#' - Identifies shortest paths up to `max_length_1` between `starts_gn` and `intermediate1_gn`.
#' - Identifies shortest paths up to `max_length_2` between `intermediate1_gn` and `intermediate2_gn`.
#' - Identifies shortest paths up to `max_length_3` between intermediate nodes and `targets_gn`.
#' - Constructs an `igraph` object representing the naive network.
#' - Optionally saves the network as an RDS file (`rds_path`) and in SIF format (`sif_path`).
#'
#' @export
#'
#' @examples
#' PKN = choose_database_for_building(organism='mouse',format='table')
#' three_layer_naive_network(starts_gn=c('Flt3'),
#'                           intermediate1_gn=c('Mapk1', 'Mapk3'),
#'                           intermediate2_gn=c('Hox9'),
#'                           targets_gn=c('Stat5a'),
#'                           PKN_table=PKN,
#'                           max_length_1=3,
#'                           max_length_2=1,
#'                           max_length_3=4)
#'
three_layer_naive_network <- function(starts_gn,
                                      intermediate1_gn,
                                      intermediate2_gn,
                                      targets_gn,
                                      PKN_table,
                                      max_length_1,
                                      max_length_2,
                                      max_length_3,
                                      both_intermediates = TRUE,
                                      keep_only_connected = FALSE,
                                      files = TRUE,
                                      rds_path = 'three_layer_naive.RDS',
                                      sif_path = 'three_layer_naive.sif',
                                      connect_all = TRUE) {

  message(paste0("First layer: Finding paths up to ", max_length_1,
                 " steps from starting nodes to intermediate1 nodes."))

  ## create first layer
  all_paths_layer1_df <- get_all_shortest_path_custom(starts_gn,
                                                      intermediate1_gn,
                                                      PKN_table,
                                                      max_length_1)

  ## filter intermediate1_gn if keep_only_connected is TRUE
  if (keep_only_connected) {
    connected_intermediate1 <- unique(c(all_paths_layer1_df$ENTITYA, all_paths_layer1_df$ENTITYB))
    removed <- setdiff(intermediate1_gn, connected_intermediate1)

    if (length(removed) > 0) {
      warning("These intermediate1 nodes were removed as they were not connected in the first layer: ",
              paste(removed, collapse = " | "))
    } else {
      warning("No intermediate1 nodes were connected in the first layer; all are retained.")
    }
    intermediate1_gn <- intersect(intermediate1_gn, connected_intermediate1)
  }

  ## create second layer
  message(paste0("Second layer: Finding paths up to ", max_length_2,
                 " steps from intermediate1 to intermediate2."))

  intermediates <- c(intermediate1_gn, intermediate2_gn)
  all_paths_layer2_df <-
    get_all_shortest_path_custom(intermediate1_gn,
                                 intermediates,
                                 PKN_table,
                                 max_length_2)


  if (keep_only_connected) {
    connected_intermediate2 <- unique(c(all_paths_layer2_df$ENTITYA, all_paths_layer2_df$ENTITYB))
    removed2 <- setdiff(intermediates, connected_intermediate2)

    if (length(removed2) > 0) {
      warning("These intermediate2 nodes were removed as they were not connected in the second layer: ",
              paste(removed2, collapse = " | "))
    } else {
      warning("No intermediate2 nodes were connected in the second layer; all are retained.")
    }
    intermediates <- intersect(intermediates, connected_intermediate2)
  }

  ## create third layer
  if (both_intermediates) {
    message(paste0("Third layer: Finding paths up to ", max_length_3,
                   " steps from all intermediate nodes to targets."))

    all_paths_layer3_df <- get_all_shortest_path_custom(
      intermediates, targets_gn, PKN_table, max_length_3
    )
  } else {
    message(paste0("Third layer: Finding paths up to ", max_length_3,
                   " steps from intermediate2 nodes to targets."))

    intermediate2_gn <- intersect(intermediate2_gn, intermediates)
    all_paths_layer3_df <- get_all_shortest_path_custom(
      intermediate2_gn, targets_gn, PKN_table, max_length_3
    )
  }

  ## combine all paths and ensure uniqueness
  all_paths_df <- dplyr::bind_rows(all_paths_layer1_df,
                                   all_paths_layer2_df,
                                   all_paths_layer3_df) %>%
    dplyr::distinct()

  network <- create_graph_from_paths(all_paths_df,
                                     PKN_table,
                                     connect_all = connect_all)

  ## set node attributes
  igraph::V(network)$mf_naive <- 'unknown'
  igraph::V(network)[name %in% starts_gn]$mf_naive <- 'start'
  igraph::V(network)[name %in% intermediate1_gn]$mf_naive <-
    'intermediate1'
  igraph::V(network)[name %in% intermediate2_gn]$mf_naive <-
    'intermediate2'
  igraph::V(network)[name %in% targets_gn]$mf_naive <- 'target'

  ## save files if required
  if (files) {
    saveRDS(network, rds_path)
    igraphToSif(network, outfile = sif_path, edgeLabel = 'INTERACTION')
  }

  return(network)
}

#' Prepare Carnival Input for Signaling Network Optimization
#'
#' Creates a protein activity table to be used as input for signaling network optimization.
#' This function:
#' \enumerate{
#'   \item Optionally adds a set of user-defined perturbed proteins (receptors) with specified activity.
#'   \item Removes predicted proteins that are not connected in the naive network.
#' }
#'
#' @param naive_network An \code{igraph} object representing the naive network connecting inferred proteins.
#' @param prediction_output A data frame of inferred proteins obtained from experimental data.
#' @param recept_list A named list of receptors with their desired activity values. If \code{NULL}, no receptors are added.
#'   The names of the list elements should be the gene symbols.
#' @param organism Character string specifying the organism; either \code{"human"} or \code{"mouse"}.
#'
#' @return A data frame of inferred proteins containing the following columns:
#'   \item{UNIPROT}{Protein Uniprot ID.}
#'   \item{gene_name}{Gene Symbol of the protein.}
#'   \item{mf}{Molecular function (e.g., transcription factor (\code{tf}), kinase (\code{kin}),
#'     phosphatase (\code{phos}), or other phosphorylated proteins (\code{other}).}
#'   \item{final_score}{Activity score for the protein.}
#'   \item{method}{Method(s) used for the computation of \code{final_score}.}
#'
#' @details
#' The function first filters the \code{prediction_output} to retain only those proteins that
#' are present in the naive network. If a receptor list is provided, it converts the list into a
#' data frame, adds the corresponding Uniprot IDs using \code{convert_gene_name_in_uniprotid()},
#' and removes any predicted proteins whose gene symbols conflict with the user-specified receptors.
#' Finally, the receptor data and the filtered predictions are combined and returned.
#'
#' @export
#'
#' @examples
#' # Load example data
#' data('toy_prot_activity_df')
#'
#' # Define a receptor list with a desired activity value
#' receptor_list <- list('MTOR' = -1, 'AMPK' = 1)
#' 
#' # Create a toy igraph object as naive network
#' naive_network <- data.frame(source = c('MTOR', 'AMPK'), target = c('ATF2', 'MYCN'), sign = c(1, -1))
#' two_layers_toy <- graph_from_data_frame(naive_network)
#' 
#' carnival_input_toy <- prepare_carnival_input(two_layers_toy, 
#'                                             toy_prot_activity_df, 
#'                                              receptor_list, 
#'                                              organism = 'human')
#'
prepare_carnival_input <- function(naive_network,
                                   prediction_output,
                                   recept_list = NULL,
                                   organism) {

  ## filter prediction output: keep only proteins present in the naive network
  prediction_output_filt <- prediction_output %>%
    dplyr::filter(gene_name %in% igraph::V(naive_network)$name) %>%
    dplyr::arrange(gene_name)

  if (!is.null(recept_list)) {
    ## transform receptor list into a data frame
    recept_df <- tibble::tibble(
      gene_name = names(recept_list),
      mf = 'rec',
      method = 'user',
      final_score = unlist(recept_list)
    )

    ## add Uniprot IDs using the conversion function and relocate the 'UNIPROT' column
    recept_df <-
      convert_gene_name_in_uniprotid(recept_df, organism) %>%
      dplyr::relocate('UNIPROT')

    ## remove predicted nodes that are specified in the receptor list
    prediction_output_filt <- prediction_output_filt %>%
      dplyr::filter(!gene_name %in% recept_df$gene_name)

    ## combine receptor data with filtered prediction output
    prediction_output_filt_rec <- dplyr::bind_rows(recept_df,
                                                   prediction_output_filt)
    return(prediction_output_filt_rec)
  }
  return(prediction_output_filt)
}
