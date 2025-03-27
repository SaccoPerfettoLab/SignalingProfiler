#' Expand and Map Edges with Experimental Phosphoproteomics Data
#'
#' This function integrates CARNIVAL-optimized edges with experimental phosphoproteomics data.
#' It maps regulatory interactions to phosphorylation events incorporating quantitative phosphoproteomics information.
#'
#' @param optimized_object A list containing the optimized CARNIVAL network, including:
#'   \itemize{
#'     \item `nodes_df`: A dataframe with node attributes.
#'     \item `edges_df`: A dataframe with regulatory interactions.
#'     \item `igraph_network`: An igraph object representing the network.
#'   }
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism used.
#' @param phospho_df A tibble containing phosphoproteomics data with the following required columns:
#'   \itemize{
#'     \item `gene_name` - Gene symbol of the phosphoprotein.
#'     \item `sequence_window` - The phosphopeptide sequence window surrounding the phosphorylated site.
#'     \item `aminoacid` - The phosphorylated amino acid (e.g., `"S"`, `"T"`, `"Y"`).
#'     \item `position` - The numeric position of the phosphosite in the protein sequence.
#'     \item `difference` - The fold-change in phosphorylation between experimental conditions.
#'     \item `significant` - Logical, indicating whether the phosphosite is significantly modulated.
#'   }
#' @param with_atlas Logical (`TRUE` by default). If `TRUE`, includes the Kinome Atlas-derived regulons
#'        in the Prior Knowledge Network (PKN); if `FALSE`, excludes them.
#' @param direct Logical (`FALSE` by default). If `TRUE`, uses only direct interactions in the Prior
#'        Knowledge Network; if `FALSE`, includes both direct and indirect interactions.
#' @param files Logical (`TRUE` by default). If `TRUE`, saves the resulting expanded network to disk;
#'        if `FALSE`, the function runs without writing files.
#' @param path_sif Character string specifying the output file path for the optimized network in
#'        **SIF format**.
#' @param path_rds Character string specifying the output file path for the optimized **RDS object**.
#'
#' @return A list (*SP_object*) containing:
#' \item{igraph_network}{An igraph object representing the expanded network, incorporating phosphoproteomics data.}
#' \item{nodes_df}{A dataframe with node attributes, including regulatory activity scores and mapped metadata.}
#' \item{edges_df}{A dataframe with edge attributes, including:}
#'   \itemize{
#'     \item `source` - The source node of the regulatory interaction.
#'     \item `target` - The target node of the regulatory interaction.
#'     \item `sign` - The interaction sign (`+1` or `-1`).
#'     \item `carnival_weight` - The weight assigned to the interaction by CARNIVAL.
#'     \item `mechanism` - The molecular mechanism of the interaction annotated in SIGNOR.
#'     \item `direct` - Logical, indicating whether the interaction is direct.
#'     \item `residue` - The putative phosphosite modulated in the interaction.
#'     \item `is_quantified` - Logical, indicating whether the phosphosite was quantified in the dataset.
#'     \item `is_significant` - Logical, indicating whether the phosphosite is significantly modulated.
#'     \item `FC` - Character, fold-change in phosphorylation (starred values are significant).
#'   }
#'
#' @seealso [add_output_carnival_nodes_attributes()], [add_output_carnival_edges_attributes()], [molecular_function_annotation()]
#'
#' @export
#'
#' @examples
#' # Load an optimized CARNIVAL network
#' data('toy_opt_network')
#' toy_opt_network$edges_df <- toy_opt_network$edges_df[, c('source', 'target', 'sign', 'carnival_weight')]
#'
#' # Load phosphoproteomics dataset
#' data('phospho_toy_df')
#'
#' expand_and_map_edges(
#'   optimized_object = toy_opt_network,
#'   organism = "human",
#'   phospho_df = phospho_toy_df,
#'   files = TRUE,
#'   with_atlas = TRUE,
#'   direct = FALSE,
#'   path_sif = "expanded_network.sif",
#'   path_rds = "expanded_network.rds"
#' )
#'
expand_and_map_edges <- function(optimized_object,
                                 organism,
                                 phospho_df,
                                 files,
                                 with_atlas = FALSE,
                                 direct = FALSE,
                                 path_sif,
                                 path_rds) {

  # Validate input organism
  if (!organism %in% c("human", "mouse")) {
    stop("Invalid organism. Please specify 'human' or 'mouse'.")
  }

  # Extract nodes and annotate molecular functions
  nodes_df <- optimized_object$nodes_df
  nodes_df <- molecular_function_annotation(inferred_proteins_dataset = nodes_df, organism = organism)
  nodes_df$method[is.na(nodes_df$final_score)] <- "CARNIVAL"

  # Extract edges and ensure correct data types
  edges_df <- optimized_object$edges_df
  colnames(edges_df) <- c('source', 'target', 'interaction','carnival_weight')
  edges_df <- edges_df %>%
    dplyr::mutate_at('interaction', as.character)

  # Load the prior knowledge network (PKN) for interactions mapping
  db <- choose_database_for_building(
    organism = organism,
    format = "table",
    with_atlas = with_atlas,
    direct = direct
  ) %>%
    dplyr::mutate(INTERACTION = as.character(INTERACTION))

  # ============================ #
  # Process phosphoproteomics data
  # ============================ #

  # Create mapping key in phosphoproteomics for joining with PKN
  # by transforming every special character in the 'gene_name' column as '_'

  if(organism == 'human'){
    phospho_df <- phospho_df %>%
      tidyr::separate_rows(gene_name, sep = ';') %>%
      dplyr::mutate(gene_name = stringr::str_to_upper(stringr::str_replace_all(gene_name, "[^[:alnum:]]", '_')))
  }else{
    phospho_df <- phospho_df %>%
      tidyr::separate_rows(gene_name, sep = ';') %>%
      dplyr::mutate(gene_name = ifelse(grepl("[^[:alnum:]]", gene_name),
                                       stringr::str_to_upper(stringr::str_replace_all(gene_name, "[^[:alnum:]]", '_')),
                                       stringr::str_replace_all(gene_name, "[^[:alnum:]]", '_')))
  }

  if('sequence_window' %in% colnames(phospho_df)){
    # Adjust phosphopeptide sequences for mapping consistency
    if (nchar(phospho_df$sequence_window[1]) > 15) {
      # Transform longer peptides in 15mers
      center <- (nchar(phospho_df$sequence_window[1]) + 1) / 2
      phospho_df <- phospho_df %>%
        dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
    } else if (nchar(phospho_df$sequence_window[1]) < 15) {
      # Transform PKN in peptide short as phosphoproteomics
      offset <- (nchar(phospho_df$sequence_window[1]) - 1) / 2
      db <- db %>%
        tidyr::separate(PHOSPHO_KEY_GN_SEQ, into = c("gene_name", "seq"), sep = "-") %>%
        dplyr::filter(nchar(seq) == 15) %>%
        dplyr::mutate(seq = stringr::str_sub(seq, 7 - offset + 1, 7 + offset + 1)) %>%
        tidyr::unite("PHOSPHO_KEY_GN_SEQ", gene_name, seq, sep = "-")
      phospho_df <- phospho_df %>%
        dplyr::rename(sequence_window_sub = sequence_window)
    } else {
      # Rename sequence_window for consistency
      phospho_df <- phospho_df %>%
        dplyr::rename(sequence_window_sub = sequence_window)
    }

    # Create a mapping key: GeneSymbol-SequenceWindow
    phospho_df <- phospho_df %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(gene_name, "-", sequence_window_sub)) %>%
      dplyr::select(gene_name, aminoacid, position, PHOSPHO_KEY_GN_SEQ, difference, significant)
  }else{
    # Create a mapping key: GeneSymbol-SequenceWindow
    phospho_df <- phospho_df %>%
      dplyr::mutate(aminoacid = dplyr::case_when(
        aminoacid == 'Y' ~ "Tyr",
        aminoacid == 'T' ~ "Thr",
        aminoacid == 'S' ~ "Ser",
        TRUE ~ NA_character_
      )) %>%
      dplyr::mutate(RESIDUE = paste0(aminoacid, position)) %>%
      dplyr::select(gene_name, aminoacid, position, RESIDUE, difference, significant)
  }

  # ============================ #
  # Expand network edges
  # ============================ #

  if('PHOSPHO_KEY_GN_SEQ' %in% colnames(phospho_df)){
    # Map edges to known phosphorylation interactions
    edges_df_expanded <- edges_df %>%
      dplyr::left_join(db, by = c("source" = "ENTITYA", "interaction" = "INTERACTION", "target" = "ENTITYB"),
                       relationship = 'many-to-many') %>%
      dplyr::distinct() %>%
      dplyr::left_join(phospho_df %>% dplyr::select(-gene_name), by = "PHOSPHO_KEY_GN_SEQ",
                       relationship = 'many-to-many')
  }else{
    edges_df_expanded <- edges_df %>%
      dplyr::left_join(db, by = c("source" = "ENTITYA", "interaction" = "INTERACTION", "target" = "ENTITYB"),
                       relationship = 'many-to-many') %>%
      dplyr::distinct() %>%
      dplyr::left_join(phospho_df, by = c("target" = 'gene_name' , "RESIDUE"),
                       relationship = 'many-to-many')
  }

  # Validate interaction coherence with protein activities
  edges_df_expanded <- edges_df_expanded %>%
    dplyr::left_join(nodes_df %>% dplyr::select(gene_name, act_source = carnival_activity), by = c("source" = "gene_name")) %>%
    dplyr::left_join(nodes_df %>% dplyr::select(gene_name, act_target = carnival_activity), by = c("target" = "gene_name")) %>%
    dplyr::mutate(keep = ifelse(as.numeric(interaction) * difference * act_target > 0, "+", NA))

  # Assign metadata flags
  edges_df_expanded <- edges_df_expanded %>%
    dplyr::mutate(
      is_quantified = !is.na(difference) & !is.na(keep),
      is_significant = !is.na(significant) & !is.na(keep),
      phosphosite_value = ifelse(is_quantified, difference, NA)
    )

  # Reorganize and consolidate edge attributes
  edges_df_final <- edges_df_expanded %>%
    dplyr::select(
      source, target, sign = interaction, carnival_weight, mechanism = MECHANISM,
      residue = RESIDUE, sequence = SEQUENCE, is_significant, phosphosite_value, DIRECT
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      residue = ifelse(is_significant, paste0(residue, "*"), residue),
      phosphosite_value = ifelse(is_significant, paste0(round(phosphosite_value, 2), "*"), round(phosphosite_value, 2))
    )

  # Merge multiple phosphosites mapping to the same interaction
  edges_df_final_amino <- edges_df_final %>%
    dplyr::filter(!is.na(phosphosite_value)) %>%
    dplyr::group_by(source, target, sign) %>%
    dplyr::reframe(aminoacid = paste0(residue,collapse = ';'),
                   FC = paste0(as.character(phosphosite_value), collapse = ';'),
                   direct = paste0(unique(DIRECT), collapse = ';'),
                   mechanism = paste0(unique(mechanism), collapse = ';')) %>%
    dplyr::mutate_at('direct', as.logical)

  edges_df_final <- dplyr::left_join(edges_df_final %>%
                                           dplyr::select(source,
                                                         target,
                                                         sign,
                                                         carnival_weight,
                                                         direct = DIRECT,
                                                         mechanism) %>%
                                           dplyr::distinct(),
                                         edges_df_final_amino,
                                         by = c('source', 'target', 'sign', 'direct', 'mechanism'))

  edges_df_final <- edges_df_final %>%
    dplyr::mutate(is_quantified = ifelse(is.na(FC), FALSE, TRUE),
                  is_significant = ifelse(grepl('\\*', FC), TRUE, FALSE))

  # When there are duplicated interactions because of different mechanisms
  # annotated in the PKN, paste them all together
  edges_df_final <- edges_df_final %>%
    dplyr::group_by(source, target, sign, carnival_weight) %>%
    dplyr::reframe(direct = paste0(unique(direct), collapse = ';'),
                   mechanism = paste0(na.omit(mechanism), collapse = ';'),
                   aminoacid = paste0(na.omit(aminoacid), collapse = ';'),
                   is_quantified = paste0(unique(is_quantified), collapse = ';'),
                   is_significant = paste0(unique(is_significant), collapse = ';'),
                   FC = paste0(na.omit(FC), collapse = ';'))

  # Create igraph object
  CARNIVAL_graph <- igraph::graph_from_data_frame(edges_df_final, nodes_df, directed = TRUE)

  # Prepare output
  SP_object <- list(igraph_network = CARNIVAL_graph, nodes_df = nodes_df, edges_df = edges_df_final)

  # Save output files if required
  if (files) {
    message("Saving annotated network with phosphoproteomics to ", path_sif)
    igraphToSif(CARNIVAL_graph, path_sif, "sign")
    message("Saving annotated network with phosphoproteomics as RDS object to ", path_rds)
    saveRDS(SP_object, path_rds)
  }

  return(SP_object)
}
