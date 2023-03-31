
#' compute_phenoscore
#'
#' @param sp_output Signaling Profiler output list
#' @param desired_phenotypes NULL if you want all phenotypes, vector of phenotypes name in uppercase
#' @param compact TRUE if you want a compact result, FALSE otherwise
#' @param max_path_length integer, from 1 to 4, path length of
#' @param remove_cascade default TRUE
#'
#' @return phenoSCORE table
#' @export
#'
#' @examples
compute_phenoscore <- function(sp_output,
                               desired_phenotypes,
                               compact = TRUE,
                               max_path_length = 4,
                               remove_cascade = TRUE){

  # filter out discordant nodes and the ones not coming from experimental data
  nodes_to_consider <- sp_output$nodes_df %>%
    dplyr::filter(!is.na(final_score) & !discordant) %>%
    dplyr::mutate_at('gene_name', toupper)


  if(!is.null(desired_phenotypes)){
    paths_to_pheno <- paths_to_pheno %>%
      dplyr::filter(EndPathways %in% desired_phenotypes)
  }

  prot_to_pheno <- dplyr::inner_join(nodes_to_consider,
                                     paths_to_pheno %>% dplyr::filter(Path_Length <= max_path_length),
                                     by = c('gene_name' = 'QueryNode')) %>%
    dplyr::mutate(signed_exp_score = final_score * Final_Effect)


  if(remove_cascade == TRUE){
    g <- sp_output$igraph_network
    multiple_regulated_phen <- prot_to_pheno %>%
      dplyr::count(EndPathways) %>%
      dplyr::filter(n>1)

    # select phenotypes
    phenotypes <- multiple_regulated_phen$EndPathways

    # initialize empty list
    to_remove_list <- list()

    for(i_phen in c(1:length(phenotypes))){
      phenotype = phenotypes[i_phen]

      prot_to_phenotype <- prot_to_pheno %>%
        dplyr::filter(EndPathways == phenotype)

      combinatios <- combn(prot_to_phenotype$name, 2)
      dist_count <- matrix(nrow = 3, ncol = ncol(combinatios))

      # forward run
      for(i in ncol(combinatios)){
        #i = 1
        dist_count[1,i] <- unlist(igraph::distances(g, combinatios[1,i], to = combinatios[2,i], mode = "out"))
        if(dist_count[1,i] == Inf){
          dist_count[2,i] <- NA
        }else{
          dist_count[2,i] <- combinatios[2,i]
        }
      }
      # reverse run
      for(i in ncol(combinatios)){
        #i = 1
        dist_count[1,i] <- unlist(igraph::distances(g, combinatios[2,i], to = combinatios[1,i], mode = "out"))
        if(dist_count[1,i] == Inf){
          dist_count[3,i] <- NA
        }else{
          dist_count[3,i] <- combinatios[1,i]
        }
      }
      output <- dist_count[!(is.na(as.vector(dist_count)) | as.vector(dist_count) == Inf)]

      to_remove_list[[i_phen]] <- output
    }

    names(to_remove_list) <- phenotypes

    to_remove_list_clean <- to_remove_list[unlist(lapply(to_remove_list, function(x){!is.null(x)}))]

    to_remove_df <- tibble::tibble(phenotype = names(to_remove_list_clean),
                                   proteins = unlist(lapply(to_remove_list_clean, function(x){paste0(x, collapse = ',')})))
    to_remove_df <- to_remove_df %>% tidyr:::separate_rows(proteins)

    prot_to_pheno <- dplyr::anti_join(prot_to_pheno,
                                            to_remove_df,
                                            by = c('EndPathways' = 'phenotype',
                                                   'name' = 'proteins'))
  }

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
