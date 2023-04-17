# this set of function manages go terms

#' Title
#'
#' @param inferred_proteins_dataset dataset of inferred proteins with UNIPROT
#' @param organism string, 'mouse', or 'human'
#' @return dataset of inferred proteins with molecular function column
#' @export
#'
#' @examples
molecular_function_annotation <- function(inferred_proteins_dataset, organism){

  # ADD A CHECK FOR UNIPROT COLUMN PRESENCE
  inferred_proteins_dataset <- inferred_proteins_dataset %>%
    dplyr::filter(!is.na(UNIPROT)) %>%
    tidyr::separate_rows('UNIPROT', sep = ';') %>%
    dplyr::distinct()

  #get uniprot annotation
  inferred_proteins_dataset_f <- inferred_proteins_dataset %>%
    dplyr::filter(!grepl('*SIGNOR*', UNIPROT))


  if(organism == 'human'){
    inferred_proteins_dataset_f <- dplyr::left_join(inferred_proteins_dataset_f,
                                                    gomf_human,
                                                    by = c('UNIPROT'))
  }else if(organism == 'mouse' | organism == 'hybrid'){
    inferred_proteins_dataset_f <- dplyr::left_join(inferred_proteins_dataset_f,
                                                    gomf_mouse,
                                                    by = c('UNIPROT'))
  }else{
    stop('please provide a valid organism')
  }

  inferred_proteins_dataset_f$mf[is.na(inferred_proteins_dataset_f$mf)] <- 'other'

  missing_genes <- inferred_proteins_dataset %>%
    dplyr::filter(!gene_name %in% inferred_proteins_dataset_f$gene_name)

  inferred_proteins_dataset <- dplyr::bind_rows(inferred_proteins_dataset_f, missing_genes)

  inferred_proteins_dataset <- inferred_proteins_dataset %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(UNIPROT = paste0(UNIPROT, collapse = ';')) %>%
    dplyr::ungroup()

  return(inferred_proteins_dataset)
}

