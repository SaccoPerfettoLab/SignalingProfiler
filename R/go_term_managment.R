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

  if('mf' %in% colnames(inferred_proteins_dataset)){
    inferred_proteins_dataset$mf <- NULL
  }

  if(organism == 'human'){
    inferred_proteins_dataset <- dplyr::left_join(inferred_proteins_dataset,
                                                  gomf_human %>% dplyr::select(gene_name, mf),
                                                  by = c('gene_name'))

  }else if(organism == 'mouse' | organism == 'hybrid'){
    inferred_proteins_dataset <- dplyr::left_join(inferred_proteins_dataset,
                                                  gomf_mouse %>% dplyr::select(gene_name, mf),
                                                  by = c('gene_name'))

  }else{
    stop('please provide a valid organism')
  }

  inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf) & grepl('SIGNOR-C', inferred_proteins_dataset$UNIPROT)] <- 'complex'
  inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf) & grepl('SIGNOR-F', inferred_proteins_dataset$UNIPROT)] <- 'fusion protein'
  inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf)] <- 'other'
  if ("method" %in% colnames(inferred_proteins_dataset)){
    inferred_proteins_dataset$mf[inferred_proteins_dataset$method == 'user'] <- 'rec'
  }

  return(inferred_proteins_dataset)
}

