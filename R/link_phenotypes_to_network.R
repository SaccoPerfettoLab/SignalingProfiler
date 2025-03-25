#' Link Phenotypes to the Signaling Network
#'
#' This function integrates phenotype regulatory relationships into the
#' existing signaling network from SignalingProfiler, creating an extended
#' network with phenotype nodes and phenotype-regulator edges.
#'
#' @param phenotype_regulators A dataframe containing phenotype-regulator
#'   interactions with columns `EndPathways`, `Effect`, and `regulators` (Gene Symbol).
#' @param phenoscore_df A dataframe containing computed phenoscores for
#'   each phenotype, including `EndPathways` and `phenoscore`.
#' @param sp_graph An `igraph` object representing the signaling network
#'   from SignalingProfiler network generation step.
#'
#' @return A list containing:
#'   \item{igraph_network}{An `igraph` object representing the extended network, including phenotype nodes.}
#'   \item{nodes_df}{A dataframe containing node attributes.}
#'   \item{edges_df}{A dataframe containing edge relationships, including phenotype-regulator interactions.}
#'
#' @details
#' This function expands the existing signaling network by:
#' - Extracting regulators for each phenotype.
#' - Creating phenotype nodes with computed phenoscores.
#' - Connecting phenotype nodes to their associated regulators.
#' - Merging the new nodes and edges with the existing `sp_graph`.
#'
#' The resulting network allows for the analysis of how protein regulators
#' influence phenotypic pathways.
#'
#' @examples
#' library(igraph)
#' library(dplyr)
#'
#' # Example phenotype regulator data
#' phenotype_regulators <- data.frame(
#'   EndPathways = c("Apoptosis", "Cell Growth"),
#'   Effect = c("up-regulates", "down-regulates"),
#'   regulators = c("P53;AKT1", "MAPK1")
#' )
#'
#' # Example phenoscore dataframe
#' phenoscore_df <- data.frame(
#'   EndPathways = c("Apoptosis", "Cell Growth"),
#'   phenoscore = c(2.5, -1.8)
#' )
#'
#' # Example signaling network
#' example_graph <- graph_from_edgelist(
#'   matrix(c("P53", "AKT1", "AKT1", "MAPK1"), ncol = 2, byrow = TRUE),
#'   directed = TRUE
#' )
#'
#' # Generate phenotype-augmented network
#' sp_object <- link_phenotypes_to_network(phenotype_regulators, phenoscore_df, example_graph)
link_phenotypes_to_network <- function(phenotype_regulators, phenoscore_df, sp_graph){

  # Extract phenotypes regulators
  phenotype_regulators <- phenotype_regulators %>%
    dplyr::select(EndPathways, Effect, regulators) %>%
    tidyr::separate_rows(regulators, sep = ';') %>%
    dplyr::mutate(Effect = ifelse(Effect == 'down-regulates', -1, 1)) %>%
    dplyr::rename('source' = 'regulators' , 'sign' = 'Effect', 'target' = 'EndPathways')

  # Create a phenotype table in SignalingProfiler compliant format

  pheno_nodes <- tidyr::tibble(gene_name = stringr::str_replace_all(str_to_upper(phenoscore_df$EndPathways), "[^[:alnum:]]", '_'),
                               carnival_activity = NA,
                               final_score = phenoscore_df$phenoscore,
                               method = 'phenoscore',
                               discordant = FALSE
  )

  # Add UNIPROT ID and molecular function to phenotypes
  pheno_nodes <- convert_gene_name_in_uniprotid(bio_dataset = pheno_nodes, organism = 'human')
  pheno_nodes <- molecular_function_annotation(inferred_proteins_dataset = pheno_nodes, organism = 'human')

  pheno_nodes <- pheno_nodes %>%
    dplyr::mutate(carnival_activity = ifelse(final_score < 0, -100, 100))

  # Get SignalingProfiler network proteins
  nodes_df <- igraph::as_data_frame(sp_graph, what = 'vertices')
  colnames(nodes_df)[1] <- 'gene_name'

  # Create a nodes table
  node_df_pheno <- dplyr::bind_rows(nodes_df, pheno_nodes) %>%
    dplyr::distinct()

  pheno_edges <- phenotype_regulators %>% dplyr::mutate(target = stringr::str_replace_all(str_to_upper(target), "[^[:alnum:]]", '_'))
  pheno_edges <- tidyr::tibble(target = pheno_edges$target,
                               sign = as.character(pheno_edges$sign),
                               source = pheno_edges$source,
                               carnival_weight = 100,
                               direct = 'FALSE',
                               aminoacid = '',
                               is_quantified = 'FALSE',
                               is_significant = 'FALSE',
                               FC = '' )

  # Get SignalingProfiler network edges
  edges_df <- tidyr::as_tibble(igraph::as_data_frame(sp_graph, what = 'edges'))
  colnames(edges_df)[1:2] <- c('source', 'target')
  edges_df$sign <- as.character(edges_df$sign)
  edges_df_pheno <- dplyr::bind_rows(edges_df, pheno_edges)

  setdiff(unique(c(pheno_edges$source, pheno_edges$target)), node_df_pheno$gene_name)

  # Create igraph object
  pheno_graph <- igraph::graph_from_data_frame(edges_df_pheno,
                                               vertices = node_df_pheno %>% dplyr::distinct())

  sp_object <- list(igraph_network = pheno_graph,
                    nodes_df = node_df_pheno,
                    edges_df = edges_df_pheno)

  return(sp_object)
}
