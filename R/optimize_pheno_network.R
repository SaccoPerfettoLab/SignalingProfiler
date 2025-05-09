#' Optimize Phenotype Network with CARNIVAL
#'
#' This function optimizes the protein-to-phenotype network by retrieving only
#' sign-coherent protein-to-phenotype edges with phenotypes activity.
#' Optionally it calls [expand_and_map_edges()] to map edges using phosphoproteomics data.
#'
#' @param sp_object A list containing the SignalingProfiler network, nodes, and edges.
#' @param organism A character string specifying the organism (`"human"` or `"mouse"`).
#' @param phospho_df An optional tibble containing phosphoproteomics data. Default is `NULL`.
#' @param carnival_options A list of options returned by `default_CARNIVAL_options()`.
#' @param files Logical, whether to generate output files (`TRUE` or `FALSE`).
#' @param direct Logical, whether to use only direct interactions (`TRUE` or `FALSE`, default: `FALSE`).
#' @param with_atlas Logical, whether to include Kinome Atlas regulons (`TRUE` or `FALSE`, default: `TRUE`).
#' @param path_sif A character string specifying the path for the SIF output file. Default is `'./pheno_opt_graph.sif'`.
#' @param path_rds A character string specifying the path for the RDS output file. Default is `'./pheno_opt_graph.rds'`.
#'
#' @return A modified SP object containing the optimized phenotype network.
#'
#' @details
#' This function:
#' - Extracts phenotype and non-phenotype nodes from the input network.
#' - Prepares the edge list in SIF format.
#' - Runs CARNIVAL-based network optimization.
#' - If phosphoproteomics data is provided, maps validated edges.
#' - Returns the updated SP object with the optimized phenotype network.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data('toy_sp_object')
#' carnival_options <- default_CARNIVAL_options(solver = "cplex")
#' optimize_pheno_network(sp_object = toy_sp_object, organism = 'human', carnival_options)
#' }
optimize_pheno_network <- function(sp_object,
                                   organism,
                                   phospho_df = NULL,
                                   carnival_options,
                                   files,
                                   direct = FALSE,
                                   with_atlas = FALSE,
                                   path_sif = './pheno_opt_graph.sif',
                                   path_rds = './pheno_opt_graph.rds'){

  # Extract phenotype and non-phenotype nodes
  start_df <- sp_object$sp_object_phenotypes$nodes_df %>%
    dplyr::filter(mf != 'phenotype') %>%
    dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf, method)

  pheno_df <- sp_object$sp_object_phenotypes$nodes_df %>%
    dplyr::filter(mf == 'phenotype') %>%
    dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf, method)

  # Transform edges into SIF format
  pheno_naive_df <- sp_object$sp_object_phenotypes$edges_df %>%
    dplyr::select(source, target, interaction = sign)

  # Run CARNIVAL optimization
  pheno_out <- run_carnival_and_create_graph(source_df = start_df,
                                             target_df = pheno_df,
                                             naive_network = unique(pheno_naive_df),
                                             proteins_df = sp_object$sp_object_phenotypes$nodes_df %>%
                                               dplyr::select(UNIPROT, gene_name, final_score, mf, method),
                                             organism = organism,
                                             carnival_options = carnival_options,
                                             files = files,
                                             direct = direct,
                                             with_atlas = with_atlas,
                                             path_sif = path_sif,
                                             path_rds = path_rds)

  # Validate network with phosphoproteomics if provided
  sp_object$sp_object_phenotypes <- if(!is.null(phospho_df)){
    expand_and_map_edges(optimized_object = pheno_out,
                         organism = organism,
                         phospho_df = phospho_df,
                         files = files,
                         direct = direct,
                         with_atlas = with_atlas,
                         path_sif = path_sif,
                         path_rds = path_rds)
  }else{
    pheno_out
  }

  return(sp_object)
}
