
#' compute_phenoscore
#'
#' @param sp_output Signaling Profiler output list
#' @param desired_phenotypes NULL if you want all phenotypes, vector of phenotypes name in uppercase
#' @param compact TRUE if you want a compact result, FALSE otherwise
#'
#' @return phenoSCORE table
#' @export
#'
#' @examples
compute_phenoscore <- function(sp_output, desired_phenotypes, compact = TRUE){

  # filter out discordant nodes and the ones not coming from experimental data
  nodes_to_consider <- sp_output$nodes_df %>%
    dplyr::filter(!is.na(final_score) & !discordant) %>%
    dplyr::mutate_at('gene_name', toupper)


  if(!is.null(desired_phenotypes)){
    paths_to_pheno <- paths_to_pheno %>%
      dplyr::filter(EndPathways %in% desired_phenotypes)
  }

  prot_to_pheno <- dplyr::inner_join(nodes_to_consider,
                                     paths_to_pheno, by = c('gene_name' = 'QueryNode')) %>%
    dplyr::mutate(signed_exp_score = final_score * Final_Effect)

  # compute phenotype score
  phenotype_score <- prot_to_pheno %>%
    dplyr::group_by(EndPathways) %>%
    dplyr::summarise(phenoSCORE = sum(signed_exp_score))

  extended_result <- dplyr::inner_join(prot_to_pheno, phenotype_score, by = 'EndPathways') %>%
    dplyr::select(name, gene_name, mf, final_score, EndPathways, phenoSCORE, Path_String, Path_Length)

  compact_result <- extended_result %>%
    dplyr::group_by(EndPathways, phenoSCORE) %>% dplyr::summarise(Regulators = paste0(unique(gene_name), collapse = ';'),
                                                                  Regulators_Paths = paste0(Path_String, collapse = ';'),
                                                                  N_Paths = dplyr::n()) %>%
    dplyr::ungroup()

  if(compact == TRUE){
    return(compact_result)
  }else{
    return(extended_result)
  }
}
