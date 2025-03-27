#' Generate Transcription Factor (TF) List
#'
#' Transform dataframe of proteins with their inferred activity in a list
#' compliant with R CARNIVAL analysis. Optionally, it extracts
#' the top-ranked proteins (TFs).
#' Function re-used from Saez-lab (https://saezlab.github.io/).
#'
#' @param df A dataframe of inferred proteins with their activity scores.
#' @param top An integer specifying the number of TFs to return (default: 50).
#'   If `"all"`, returns all proteins in `df`.
#' @param access_idx An integer index referring to the column containing activity values.
#'
#' @return A list of ranked proteins with their scores.
generateTFList <- function (df = df, top = 'all', access_idx = 1){
  if (top == "all")  top <- nrow(df)

  if (top > nrow(df)) {
    warning("Number of to TF's inserted exceeds the number of actual TF's in the\n            data frame. All the TF's will be considered.")
    top <- nrow(df)
  }
  ctrl <- intersect(x = access_idx, y = seq_len(ncol(df)))
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

#' Format Inferred Proteins for CARNIVAL Analysis
#'
#' Converts inferred protein scores into a format compatible with CARNIVAL.
#'
#' @param proteins_df A dataframe with inferred protein activities and scores.
#'
#' @return A formatted list of protein activities for CARNIVAL.
#' @export
#'
#' @examples
#' data <- data.frame(gene_name = c("TP53", "MYC"), final_score = c(2.5, -1.8))
#' formatted_proteins <- formatting_proteins_for_carnival(data)
#'
formatting_proteins_for_carnival <- function(proteins_df){

  proteins_df <- proteins_df %>%
    dplyr::select(gene_name, t = final_score) %>%
    dplyr::distinct()

  proteins_activities <- data.frame(unique(proteins_df))
  rownames(proteins_activities) <- unique(proteins_activities$gene_name)
  proteins_activities$gene_name <- NULL
  proteinList <- generateTFList(proteins_activities,
                                top = 'all',
                                access_idx = 1)
  return(proteinList)
}

#' Discretize Initiator Proteins for CARNIVAL Analysis
#'
#' This function discretizes protein activity scores into categorical values
#' (-1, or 1) for use in `CARNIVAL optimization`. It filters proteins
#' that have either positive or negative inferred activity while excluding neutral ones.
#'
#' @param proteins_df A tibble containing initiator proteins with their
#'        inferred activity scores. The dataframe must include:
#'        \describe{
#'          \item{gene_name}{Character. Gene symbol of the protein.}
#'          \item{final_score}{Numeric. The continuous activity score of the protein.}
#'          \item{mf}{Character. Molecular function category (e.g., "kin", "tf").}
#'        }
#'
#' @return A tibble containing discretized initiators with the following columns:
#' \describe{
#'   \item{gene_name}{Character. Gene symbol of the protein.}
#'   \item{final_score}{Integer. Discretized activity score (-1 for inhibition, 1 for activation).}
#'   \item{mf}{Character. Molecular function of the protein.}
#' }
#'
#' @export
#'
#' @examples
#' # Example dataframe
#' proteins_df <- tibble::tibble(
#'   gene_name = c("ProteinA", "ProteinB", "ProteinC"),
#'   final_score = c(0.8, -0.5, 0),
#'   mf = c("kin", "tf", "phos")
#' )
#'
#' # Convert to discretized format
#' initiators <- create_discretized_initiators_for_carnival(proteins_df)
#' print(initiators)
#'
create_discretized_initiators_for_carnival <- function(proteins_df) {
  # Initialize discretized column with zero (neutral state)
  proteins_df$discretized <- 0

  # Assign -1 for negative activity and 1 for positive activity
  proteins_df$discretized[proteins_df$final_score < 0] <- -1
  proteins_df$discretized[proteins_df$final_score > 0] <- 1

  # Filter for only significantly modulated proteins
  initiators_df <- proteins_df %>%
    dplyr::filter(discretized %in% c(-1, 1)) %>%
    dplyr::select(gene_name, final_score = discretized, mf)

  return(initiators_df)
}

#' Generate Default Parameters to Run CARNIVAL
#'
#' This function provides the default configuration settings for the
#' CARNIVAL optimization pipeline, selecting the appropriate Integer
#' Linear Programming solver and returning the corresponding default parameters.
#'
#' @param solver A string specifying the solver to be used for optimization.
#'        Valid options are:
#'        \itemize{
#'          \item `"lpSolve"` - Uses `lpSolve` as the solver.
#'          \item `"cplex"` - Uses IBM CPLEX as the solver.
#'          \item `"cbc"` - Uses the CBC (Coin-or branch and cut) solver.
#'           \item `"gurobi"` - Uses Gurobi as the solver.
#'        }
#' @param gurobi_path A string specifying the 'gurobi' solver path. To use only if gurobi is chosen.
#' 
#' @return A named list containing the default CARNIVAL parameters for the
#'         specified solver.
#'
#' @export
#'
#' @examples
#' # Generate default options for CPLEX solver
#' cplex_options <- default_CARNIVAL_options(solver = "cplex")
#'
#' # Generate default options for CBC solver
#' cbc_options <- default_CARNIVAL_options(solver = "cbc")
#'
#'  # Generate default options for Gurobi solver
#' gurobi_options <- default_CARNIVAL_options(solver = "gurobi",
#'                                            gurobi_path = './gurobi.exe')
default_CARNIVAL_options <- function(solver = NULL, gurobi_path = NULL) {
  # Validate solver input
  if (is.null(solver)) {
    stop("Please call default_CARNIVAL_options(solver) with a
          specific solver argument. Valid solvers are: 'lpSolve', 'cplex', 'gurobi', or 'cbc'.
         If solver='gurobi' provide gurobi_path")
  }

  solver_options <- c("lpSolve", "cplex", "cbc", "gurobi")
  solver <- match.arg(solver, choices = solver_options)

  # Retrieve default options based on selected solver
  opts <- switch(solver,
                 "lpSolve" = CARNIVAL::defaultLpSolveCarnivalOptions(),
                 "cplex" = CARNIVAL::defaultCplexCarnivalOptions(),
                 "cbc" = CARNIVAL::defaultCbcSolveCarnivalOptions(),
                 "gurobi" = CARNIVAL::defaultCplexCarnivalOptions())

  # Trick to use gurobi with cplex settings
  if(solver=='gurobi'){
    if(is.null(gurobi_path)){
      stop('Please provide the gurobi.exe file path!')
    }
    opts$solverPath <- gurobi_path
    opts$solver <- 'gurobi'
  }

  # Ensure LP files are not retained
  opts$keepLPFiles <- FALSE

  return(opts)
}

#' Validate CARNIVAL Input Data
#'
#' This function verifies the integrity of input data provided to CARNIVAL,
#' ensuring that all required elements are correctly formatted and present.
#' It checks the source nodes, target nodes, naive network structure,
#' inferred protein data, and organism type.
#'
#' @param source_df A tibble containing source nodes (perturbations) with discretized values of 1 and -1.
#'        Can be `NULL` if running an inverse CARNIVAL analysis.
#' @param target_df A tibble containing target nodes (e.g., transcription factors, kinases)
#'        with continuous activity scores.
#' @param naive_network A tibble representing the naive network in SIF format, containing
#'        three columns: `source`, `interaction`, and `target`.
#' @param proteins_df A tibble of inferred proteins with associated metadata, including
#'        molecular function (*mf*) and method of protein inference (*method*).
#' @param organism A string specifying the organism being analyzed.
#'        Valid options are `"mouse"` or `"human"`.
#'
#' @return Returns `TRUE` if all input data passes validation.
#' Otherwise, stops execution with an error message.
#'
#' @examples
#' # Example dataset of perturbations (source nodes)
#' toy_carnival_input <- data.frame(UNIPROT = c('P42345', 'SIGNOR-C15', 'P31749', 'P15336'),
#'           gene_name = c('MTOR', 'AMPK', 'AKT1', 'ATF2'),
#'           mf = c('rec', 'rec', 'kin', 'tf'),
#'           method = c('user', 'user', 'VIPER', 'VIPER'),
#'           final_score = c(-1, 1, -2.64, 4.48))
#'           
#' source_df <- toy_carnival_input %>% dplyr::filter(mf == 'rec')
#'
#' # Example dataset of targets (continuous values)
#' target_df <- toy_carnival_input %>% dplyr::filter(mf %in% c('kin', 'phos', 'other'))
#'
#' # Example naive network in SIF format
#' naive_network <- data.frame(source = c('MTOR', 'AMPK'), interaction = c(1, -1), target = c('ATM', 'CDK2'))
#'
#' # Example inferred protein dataset
#' data(toy_prot_activity_df)
#' proteins_df <- toy_prot_activity_df
#'
#' # Run validation
#' check_CARNIVAL_inputs(source_df, target_df, naive_network, proteins_df, organism = "human")
#' @export
check_CARNIVAL_inputs <- function(source_df, target_df,
                                  naive_network, proteins_df,
                                  organism){

  # Validate source perturbations
  if (!is.null(source_df)) {
    stopifnot(is.data.frame(source_df))
    stopifnot(all(c("gene_name", "mf", "final_score") %in% names(source_df)))
    stopifnot(ncol(source_df) >= 3)
  }

  # Validate target nodes
  stopifnot(is.data.frame(target_df))
  stopifnot(all(c("gene_name", "mf", "final_score") %in% names(target_df)))
  stopifnot(ncol(target_df) >= 3)

  # Validate naive network structure
  stopifnot(is.data.frame(naive_network))
  stopifnot(all(c("source", "interaction", "target") %in% names(naive_network)))
  stopifnot(ncol(naive_network) == 3)

  # Validate inferred protein dataset
  stopifnot(is.data.frame(proteins_df))
  stopifnot(all(c("gene_name", "mf", "final_score", "method") %in% names(proteins_df)))
  stopifnot(ncol(proteins_df) >= 4)

  # Validate organism selection
  stopifnot(is.character(organism))
  stopifnot(organism %in% c("mouse", "human"))

  return(TRUE)
}


#' Filter Source Nodes Present in the Naive Network
#'
#' This function filters the source dataframe (`source_df`) to keep only
#' the receptors (source nodes) that are present in the naive network.
#'
#' @param source_df A dataframe containing receptor nodes with a column `gene_name`.
#' @param naive_network A dataframe representing the naive network, containing `source` and `target` columns.
#'
#' @return A dataframe containing only the receptors that are present in the naive network.
#' @export
#'
#' @seealso [one_layer_naive_network()] [two_layer_naive_network()] [three_layer_naive_network()]
#'
#' @examples
#' # Example naive network
#' naive_network <- data.frame(source = c('MTOR', 'MTOR'), interaction = c(1, -1), target = c('ATM', 'CDK2'))
#'
#' # Example source dataframe
#' toy_carnival_input <- data.frame(UNIPROT = c('P42345', 'SIGNOR-C15', 'P31749', 'P15336'),
#'           gene_name = c('MTOR', 'AMPK', 'AKT1', 'ATF2'),
#'           mf = c('rec', 'rec', 'kin', 'tf'),
#'           method = c('user', 'user', 'VIPER', 'VIPER'),
#'           final_score = c(-1, 1, -2.64, 4.48))
#' source_df <- toy_carnival_input %>% dplyr::filter(mf == 'rec')
#'
#' # Filter source nodes present in naive network
#' filtered_source_df <- keep_only_present_perturbation(source_df, naive_network)
#' print(filtered_source_df)
#'
keep_only_present_perturbation <- function(source_df, naive_network) {

  # Ensure the naive network has the required structure
  if (!all(c("source", "target") %in% colnames(naive_network))) {
    stop("The naive_network dataframe must contain 'source' and 'target' columns.")
  }

  # Convert naive network to igraph object
  naive_graph <- igraph::graph_from_data_frame(naive_network, directed = TRUE)

  # Filter source_df to keep only nodes present in the naive network
  source_df_present <- dplyr::filter(source_df, gene_name %in% igraph::V(naive_graph)$name)

  return(source_df_present)
}

#' Run CARNIVAL and Create an Optimized Network Graph
#'
#' This function runs CARNIVAL, either in **inverse mode** (without perturbations)
#' or **vanilla mode** (with perturbations), and generates an optimized network
#' using **igraph**. The output network is based on the given naive regulatory
#' interactions and inferred protein activity scores.
#'
#' @param source_df A tibble containing source nodes (perturbations) with discretized values (1 or -1).
#'        If `NULL`, the function runs **Inverse CARNIVAL** (i.e., without perturbations).
#' @param target_df A tibble containing target nodes (e.g., transcription factors, kinases)
#'        with continuous activity scores.
#' @param naive_network A tibble representing the naive network in SIF format, containing
#'        three columns: `source`, `interaction`, and `target`.
#' @param carnival_options A list of optimization options generated by `default_CARNIVAL_options()`.
#' @param proteins_df A tibble of inferred proteins with associated metadata, including
#'        molecular functions and the method of inference.
#' @param organism A string specifying the organism being analyzed. Valid options: `"mouse"` or `"human"`.
#' @param with_atlas Logical, `TRUE` (default) includes a Prior Knowledge Network integrated with Kinome Atlas-derived regulons, `FALSE` excludes them.
#' @param direct Logical, `FALSE` (default) allows both direct and indirect interactions; `TRUE` uses only direct interactions.
#' @param files Logical, `TRUE` (default) saves output files, `FALSE` runs analysis without saving files.
#' @param path_sif A string specifying the output file path for the optimized network in **SIF format**.
#' @param path_rds A string specifying the output file path for the optimized **RDS object**.

#' @return *SP_object*, a list  containing:
#' \itemize{
#'   \item `igraph_network` - An igraph object representing the optimized network.
#'   \item `nodes_df` - A dataframe with node attributes, including inferred activity scores.
#'   \item `edges_df` - A dataframe with edge attributes, including regulatory interactions and confidence scores.
#' }
#'
#' @seealso [add_output_carnival_nodes_attributes()] [add_output_carnival_edges_attributes()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' toy_carnival_input <- data.frame(UNIPROT = c('P42345', 'SIGNOR-C15', 'P31749', 'P15336'),
#'           gene_name = c('MTOR', 'AMPK', 'AKT1', 'ATF2'),
#'           mf = c('rec', 'rec', 'kin', 'tf'),
#'           method = c('user', 'user', 'VIPER', 'VIPER'),
#'           final_score = c(-1, 1, -2.64, 4.48))

#' # Example dataset of perturbations (source nodes)
#' receptors_df <- toy_carnival_input %>% dplyr::filter(mf == 'rec')
#'
#' # Example dataset of targets (continuous values)
#' target_df <- toy_carnival_input %>%
#' dplyr::filter(mf %in% c('kin', 'phos', 'other'))
#'
#' # Example naive network in SIF format
#' toy_naive_network_sif <- data.frame(from = c("MTOR", "AKT1", "AMPK"), 
#'                                     interaction = c('-1', '1', '1'), 
#'                                     to = c("AKT1", "ATF2", "MTOR"))
#'
#' # Example inferred protein dataset
#' data(toy_prot_activity_df)
#'
#' # Load default CARNIVAL options
#' carnival_options <- default_CARNIVAL_options(solver = "cplex")
#'
#' # Run CARNIVAL and generate optimized network
#' optimized_network <- run_carnival_and_create_graph(
#'   source_df, target_df, toy_naive_network, carnival_options,
#'   proteins_df, organism = "human",
#'   files = FALSE)
#' }
run_carnival_and_create_graph <- function(source_df,
                                          target_df,
                                          naive_network,
                                          carnival_options,
                                          proteins_df,
                                          organism,
                                          with_atlas = TRUE,
                                          direct = FALSE,
                                          files = TRUE,
                                          path_sif = './optimized_network.sif',
                                          path_rds = './optimized_SP_object.RDS'){

  message('Running CARNIVAL Optimization...')

  if (is.null(source_df)) {
    message(' ** Running Inverse CARNIVAL (No perturbations) ** ')
    message('Credits to Prof. Julio Saez-Rodriguez. More info: https://saezlab.github.io/CARNIVAL/')

    check_CARNIVAL_inputs(source_df, target_df, naive_network, proteins_df, organism)

    carnival_result <- CARNIVAL::runInverseCarnival(
      measurements = unlist(formatting_proteins_for_carnival(target_df)$t),
      priorKnowledgeNetwork = unique(naive_network),
      carnivalOptions = carnival_options
    )
  } else {
    message(' ** Running Vanilla CARNIVAL (With Perturbations) ** ')

    check_CARNIVAL_inputs(source_df, target_df, naive_network, proteins_df, organism)

    # Keep only present perturbations
    source_df <- keep_only_present_perturbation(source_df, naive_network)
    target_df <- keep_only_present_perturbation(target_df, naive_network)

    # Discretize initiators
    source_df_disc <- create_discretized_initiators_for_carnival(source_df)

    carnival_result <- CARNIVAL::runVanillaCarnival(
      perturbations = unlist(formatting_proteins_for_carnival(source_df_disc)$t),
      measurements = unlist(formatting_proteins_for_carnival(target_df)$t),
      priorKnowledgeNetwork = unique(naive_network),
      carnivalOptions = carnival_options
    )
  }

  if(nrow(carnival_result$sifAll[[1]]) == 0){
    message('No network found for your experiment')
    return(NULL)
  }

  # Add attributes to the edges
  edges_df <- add_output_carnival_edges_attributes(carnival_result) %>%
    dplyr::filter(carnival_weight != 0)

  # Keep only nodes that are active or involved in relations
  nodes_df <- add_output_carnival_nodes_attributes(carnival_result, proteins_df, organism, with_atlas, direct) %>%
    dplyr::filter(gene_name %in% c(edges_df$source, edges_df$target)) %>%
    dplyr::distinct()

  # Create igraph object
  CARNIVAL_igraph_network <- igraph::graph_from_data_frame(edges_df,
                                                           nodes_df,
                                                           directed = TRUE)

  # Create a list representing the SignalingProfiler object
  SP_object <- list(igraph_network = CARNIVAL_igraph_network,
                    nodes_df = nodes_df,
                    edges_df = edges_df)

  if (files) {
    message(paste0('Saving optimized network to: ', path_sif))
    igraphToSif(CARNIVAL_igraph_network, path_sif, 'sign')

    message(paste0('Saving optimized RDS object to: ', path_rds))
    saveRDS(SP_object, path_rds)
  }

  return(SP_object)
}


