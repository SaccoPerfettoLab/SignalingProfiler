#' Compute Phenoscore for a Given Protein Dataset
#'
#' This function computes the activity modulation of cellular phenotypes (*Phenoscore*)
#' by analyzing the interaction between proteins and phenotypes in SIGNOR and
#' the protein activity modulation in a SignalingProfiler generated network.
#' For each phenotype, it calculates enrichment scores, statistical significance,
#' and optionally visualizes the results.
#'
#' @param proteins_df A dataframe containing differentially expressed proteins
#'   (upregulated or downregulated) with associated activity scores.
#' @param desired_phenotypes A vector of phenotype names to investigate.
#'   If `NULL`, all available phenotypes are used.
#' @param sp_graph An `igraph` object representing the SignalingProfiler
#'   signaling network.
#' @param pheno_distances_table A dataframe containing protein-phenotype distances.
#'   If `NULL`, the full distance table is used.
#' @param path_length An integer (2 to 4) specifying the maximum signaling path
#'   length considered for phenotype association. Default is `4`.
#' @param stat A character string, either `"mean"` or `"median"`, indicating
#'   how the phenoscore should be calculated. Default is `"mean"`.
#' @param zscore_threshold A numeric threshold for filtering significant
#'   pathways based on the computed z-score. Default is `-1.96`.
#' @param n_random Number of randomization iterations for statistical
#'   testing. Default is `1000`.
#' @param pvalue_threshold A numeric value specifying the significance threshold
#'   for filtering pathways. Default is `0.05`.
#' @param remove_cascade Logical, whether to consider only model proteins
#'   independently regulating the phenotype. Default is `TRUE`.
#' @param node_idx Logical, whether to weight the contribution of
#'   each protein on the phenotype based on the number of regulatory paths
#'   (redundancy) to the phenotype. Default is `FALSE`.
#' @param use_carnival_activity Logical, whether to use CARNIVAL activity
#'   scores instead of final scores. Default is `FALSE`.
#' @param create_pheno_network Logical, whether to integrate phenotypes
#'   in the SignalingProfiler model. Default is `TRUE`.
#'
#' @return A list containing:
#'   \item{barplot}{A ggplot2 object representing phenotype modulation scores.}
#'   \item{table_regulators}{A dataframe summarizing the proteic regulators of each phenotypes.}
#'   \item{table_phenotypes}{A dataframe containing the computed Phenoscore values.}
#'   \item{sp_object_phenotypes}{(Optional) A list with phenotype network objects if `create_pheno_network = TRUE`.}
#'
#' @details
#' This function performs the following steps:
#' - Checks input data validity and processes gene names.
#' - Identifies significant paths from proteins to phenotypes through `script.py` in SIGNOR database.
#' - Identifies significantly close protein to user-selected phenotypes via `z-test`.
#' - Computes Phenotypes enrichment scores based on randomization and statistical tests via `t-test`.
#' - Filters results based on user-defined thresholds.
#' - (Optionally) It integrates phenotypes within the network.
#'
#' @export
#'
#' @examples
#' # Example dataset
#' data('prot_toy_df')
#' data('phospho_toy_df')
#' data('toy_opt_network')
#' desired_phenotypes <- c("APOPTOSIS", "PROLIFERATION")
#'
#' results <- phenoscore_computation(proteins_df = toy_opt_network$nodes_df,
#'                                   desired_phenotypes = c('APOPTOSIS', 'PROLIFERATION'),
#'                                   sp_graph = toy_opt_network,
#'                                   create_pheno_network = TRUE)
phenoscore_computation <- function(proteins_df,
                                   desired_phenotypes = NULL,
                                   sp_graph,
                                   pheno_distances_table = NULL,
                                   path_length = 4,
                                   stat = 'mean',
                                   zscore_threshold = -1.96,
                                   n_random = 1000,
                                   pvalue_threshold = 0.05,
                                   remove_cascade = TRUE,
                                   node_idx = FALSE,
                                   use_carnival_activity = FALSE,
                                   create_pheno_network = TRUE){

  # Determine organism type based on gene naming convention
  if (sum(proteins_df$gene_name == stringr::str_to_title(proteins_df$gene_name)) == nrow(proteins_df)) {
    organism_type <- 'mouse'
    proteins_df <- proteins_df %>%
      dplyr::mutate(gene_name = stringr::str_to_upper(gene_name))
    igraph::V(sp_graph)$name <- stringr::str_to_upper(igraph::V(sp_graph)$name)
  } else {
    organism_type <- 'human'
  }

  message('** PHENOSCORE ANALYSIS **')

  ##############################################################################
  # PARAMETERS INPUT CHECK #
  ##############################################################################

  # Validate path_length
  if (!(path_length %in% 2:4)) {
    stop("path_length should be an integer between 2 and 4!")
  }

  # Validate z-score threshold
  if (!is.numeric(zscore_threshold)) {
    stop("zscore_threshold should be numeric, e.g., -1.96 or -2.58")
  }

  if (!stat %in% c('mean', 'median')) {
    stop("Please insert a valid stat: mean or median")
  }

  # Validate p-value threshold
  if (!is.numeric(pvalue_threshold)) {
    stop("pvalue_threshold should be numeric")
  }

  # Filter proteins with missing final_score if use_carnival_activity is FALSE
  if (!use_carnival_activity) {
    message('Removing proteins without a final_score value')
    proteins_df <- proteins_df %>%
      dplyr::filter(!is.na(final_score))
  }

  # If desired_phenotypes is NULL, use all available phenotypes
  if (is.null(desired_phenotypes)) {
    message('No desired_phenotypes specified. Using all available phenotypes.')
    desired_phenotypes <- unique(pheno_distances_table$EndPathways)
  }

  ##############################################################################
  # BUILD PATHS TABLE OF SIGNIFICANTLY CLOSE PROTEINS TO PHENOTYPES #
  ##############################################################################
  message('Building significant paths to phenotypes table')

  # Use full phenotype distances table if none is provided
  if (is.null(pheno_distances_table)) {
    message('No custom pheno_distance_table provided. Using full dataset without path filtering.')
    pheno_distances_table <- get(data('phenoscore_distances_table'))
  }

  # Clean phenotype distance table by replacing non-alphanumeric characters
  pheno_distances_table <- pheno_distances_table %>%
    dplyr::mutate(QueryNode = stringr::str_replace_all(QueryNode, "[^[:alnum:]]", '_'))

  # Filter paths based on user-defined path_length
  filtered_paths <- pheno_distances_table %>%
    dplyr::filter(Path_Length <= path_length & EndPathways %in% desired_phenotypes) %>%
    dplyr::select(1:8) %>%
    dplyr::distinct(Path_String, QueryNode, EndNode, EndPathways, .keep_all = TRUE)

  filtered_paths$Effect <- '-'
  filtered_paths$Effect[ filtered_paths$Final_Effect == 1 ] <- 'up-regulates'
  filtered_paths$Effect[ filtered_paths$Final_Effect == 0 ] <- '-'
  filtered_paths$Effect[ filtered_paths$Final_Effect == -1 ] <- 'down-regulates'

  # Compute path statistics
  path_summary <- filtered_paths %>%
    dplyr::filter(Effect != "-") %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(n = n(),
                     mean = mean(Path_Score),
                     median = median(Path_Score),
                     sd = sd(Path_Score),
                     .groups = "drop")

  # Merge path data with computed statistics
  path_data <- dplyr::left_join(filtered_paths %>%
                                  dplyr::filter(Effect != "-"),
                                path_summary, by = c("EndPathways", "Effect"))

  # Compute z-score with user-defined stat
  path_data <- path_data %>%
    dplyr::mutate(zscore = (Path_Score - dplyr::case_when(
      stat == "mean" ~ mean,
      stat == "median" ~ median
    )) / sd) %>%
    dplyr::filter(zscore <= zscore_threshold) %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways, .keep_all = TRUE)

  # If no significant pathways, return NULL early
  if (nrow(path_data) == 0) {
    message("No significantly close proteins found! Try adjusting max_length parameters.")
    return(NULL)
  }

  ##############################################################################
  # COMPUTE PHENOTYPES ENRICHMENT
  ##############################################################################

  # Count the number of activating/inhibiting paths on each phenotype
  enrichment_results <- path_data %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(total_paths = n(), .groups = "drop")

  # Perform randomization for statistical significance
  random_distributions <- tibble::tibble()
  for (i in seq_len(n_random)) {
    random_sample <- sample(get(data('background_phenoscore')), size = nrow(proteins_df))
    randomized <- path_data %>%
      dplyr::filter(QueryNode %in% random_sample) %>%
      dplyr::count(EndPathways, Effect, name = "hits") %>%
      dplyr::left_join(enrichment_results, by = c("EndPathways", "Effect")) %>%
      dplyr::mutate(Frac_rand = dplyr::case_when(total_paths > 0 ~ hits / total_paths,
                                          total_paths <= 0 ~ 0))

    random_distributions <- dplyr::bind_rows(random_distributions, randomized)
  }

  # Compute mean/sd of randomized fractions
  random_stats <- random_distributions %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(mean_rand = mean(Frac_rand), sd_rand = sd(Frac_rand), .groups = "drop")

  # Compare paths from input proteins against random distributions
  observed_results <- path_data %>%
    # Count #Paths from input proteins
    dplyr::filter(QueryNode %in% proteins_df$gene_name) %>%
    dplyr::count(EndPathways, Effect, name = "hits_INPUT") %>%
    dplyr::left_join(enrichment_results, by = c("EndPathways", "Effect")) %>%
    dplyr::left_join(random_stats, by = c("EndPathways", "Effect")) %>%
    dplyr::mutate(Fraction = hits_INPUT / total_paths,
           t_score = (Fraction - mean_rand) / (sd_rand + 1e-15),
           pvalue = pnorm(t_score, lower.tail = FALSE),
           pvalue_adj = p.adjust(pvalue, method = "BH"))

  # Filter by p-value threshold
  significant_results <- observed_results %>%
    dplyr::filter(pvalue_adj < pvalue_threshold, Effect != "-") %>%
    dplyr::mutate(Significance = dplyr::case_when(
      pvalue_adj < 0.00001 ~ "****",
      pvalue_adj < 0.0001 ~ "***",
      pvalue_adj < 0.001 ~ "**",
      pvalue_adj < 0.05 ~ "*",
      TRUE ~ "!"
    ))

  # Return only desired phenotypes
  if(!is.null(desired_phenotypes)){
    significant_results <- significant_results %>%
      dplyr::filter(EndPathways %in% desired_phenotypes)
  }

  if (nrow(significant_results) == 0) {
    message("No significantly modulated phenotypes found! Try adjusting parameters.")
    return(NULL)
  }

  # Ensure data types
  path_data <- path_data %>%
    dplyr::mutate_at(c('EndPathways', 'QueryNode'), as.character)

  # Map phenotypes' enrichment on significantly close paths
  path_data_enriched <- dplyr::inner_join(path_data,
                                          significant_results, by = c('EndPathways', 'Effect'))

  # Retrieve paths invol
  input_prot_paths <- path_data_enriched %>%
    dplyr::filter(QueryNode %in% proteins_df$gene_name)

  # (Optional) Remove cascade
  if(remove_cascade == TRUE & length(unique(input_prot_paths$QueryNode)) != 1 ){
    message('Removing signaling cascade regulators')
    input_prot_paths <- remove_signaling_cascade(input_prot_paths = input_prot_paths,
                                                 sp_graph = sp_graph)
  }

  # Count the # Paths with a specific effect per phenotype regulator
  count_query <- input_prot_paths %>%
    dplyr::count(QueryNode, EndPathways, Effect, name = "n")

  # Aggregate regulators and compute node_idx
  # the weight of that regulator over the phenotype
  ProteinsPaths <- count_query %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(
      regulators = paste(QueryNode, collapse = ";"),
      node_idx = paste(round(n / sum(n), 2), collapse = ";"),
      .groups = "drop") %>%
    dplyr::mutate(
      EndPathways = stringr::str_replace_all(tolower(EndPathways), "_", " ")
    )

  # Compute log10 p-values and format the results table
  significant_results <- significant_results %>%
    dplyr::mutate(
      Log10_p_value_plot = -log10(pvalue),
      EndPathways = stringr::str_replace_all(tolower(EndPathways), "_", " ")
    ) %>%
    dplyr::arrange(Effect, Log10_p_value_plot) %>%
    tibble::as_tibble()


  # Compute node_idx
  results.table_reg <- dplyr::left_join(significant_results,
                                        ProteinsPaths,
                                        by = c("EndPathways", "Effect")) %>%
    tidyr::separate_rows(regulators, node_idx, sep = ";") %>%
    dplyr::mutate(node_idx = as.numeric(node_idx))

  ##############################################################################
  # Compute PhenoScore
  ##############################################################################

  # Merge with the enriched phenotypes regulators in the SignalingProfiler
  # network with the protein activity in the mode,
  prot_pheno_act <- dplyr::inner_join(results.table_reg,
                                      proteins_df,
                                      by = c("regulators" = "gene_name")) %>%
    dplyr::mutate(Sign = ifelse(Effect == "up-regulates", 1, -1))

  # Compute PhenoScore according to user-specified parameters
  if(node_idx == TRUE){
    if(use_carnival_activity == TRUE){
      phenoscore_df <- prot_pheno_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Sign*Log10_p_value_plot*node_idx))

    }else{
      phenoscore_df <- prot_pheno_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(final_score*Sign*Log10_p_value_plot*node_idx))
    }
  }else{
    if(use_carnival_activity == TRUE){
      phenoscore_df <- prot_pheno_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Sign*Log10_p_value_plot))

    }else{
      phenoscore_df <- prot_pheno_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(final_score*Sign*Log10_p_value_plot))
    }
  }

  phenoscore_df_filt <- phenoscore_df %>%
    dplyr::filter(phenoscore != 0)

  if(nrow(phenoscore_df_filt) != nrow(phenoscore_df)){
    message("Some phenotypes are missing because had just one regulator annotated
            both as positive and negative actor")
  }

  ##############################################################################
  # BarPlot generation
  ##############################################################################

  message('PhenoScore Barplot generation')

  phenoscore_df_filt <- phenoscore_df_filt %>%
    dplyr::mutate(reg = ifelse(phenoscore < 0, 'down', 'up'))

  color_list <- list('down' = '#407F7F', 'up' = '#D46A6A')

  barplot_phenotypes <- ggplot2::ggplot(phenoscore_df_filt,
                                        ggplot2::aes(x = forcats::fct_reorder(EndPathways, phenoscore),
                                                     y = phenoscore, fill = reg))+
    ggplot2::geom_bar(stat = 'identity', alpha = 0.8)+
    ggplot2::scale_fill_manual(values = color_list, labels = names(color_list)) +
    ggplot2::ggtitle(paste0('Phenoscore')) +
    ggplot2::ylab("phenotype modulation") +
    ggplot2::xlab("") +
    ggplot2::theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          plot.title=element_text(hjust = 0.5),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text( size = 14, face = 'bold'),
          axis.text.x=element_text(size=10, angle=0, vjust=0.5, hjust=0.5))+
    ggplot2::coord_flip()

  if(organism_type == 'mouse'){
    prot_pheno_act$regulators <- stringr::str_to_title(prot_pheno_act$regulators)
    igraph::V(sp_graph)$name <- stringr::str_to_title(igraph::V(sp_graph)$name)
  }

  if(create_pheno_network){
    # Link phenotypes in the SignalingProfiler network
    pheno_graph_object <- link_phenotypes_to_network(phenotype_regulators = prot_pheno_act,
                                                     phenoscore_df = phenoscore_df_filt,
                                                     sp_graph = sp_graph)

    # Return with network object
    return(list(barplot = barplot_phenotypes,
                table_regulators = results.table_reg,
                table_phenotypes = phenoscore_df,
                sp_object_phenotypes = pheno_graph_object))
  }

  # Return with network object
  return(list(barplot = barplot_phenotypes,
              table_regulators = results.table_reg,
              table_phenotypes = phenoscore_df))
}
