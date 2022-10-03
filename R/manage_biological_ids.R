# these functions aim at managing annoying biological IDs
# passing from gene_name to UNIPROT ID and viceversa
# these functions exploit SIGNOR database
# since all the package aim at using proteins annotated in SIGNOR

#' add uniprot id in dataset from gene name column
#'
#' @param bio_dataset biological dataset without uniprot id
#' @param organism organism human or mouse
#'
#' @return biological dataset with uniprot id from SIGNOR added
#' @export
#'
#' @examples
convert_gene_name_in_uniprotid <- function(bio_dataset, organism){

  # bio_dataset <- f
  # organism <- 'mouse'
  if(organism == 'human'){
    db <- PKN_proteins_human
  }else if(organism == 'mouse'){
    db <- PKN_proteins_mouse
  }else{
    stop('please provide a valid organism')
  }


  bio_dataset_with_id <- dplyr::left_join(bio_dataset,
                         db %>% dplyr::select(ID, ENTITY),
                         by = c('gene_name' = 'ENTITY')) %>%
    dplyr::rename(UNIPROT=ID) %>%
    dplyr::filter(!is.na(UNIPROT)) %>%
    dplyr::mutate(UNIPROT = sub('-\\d', '', UNIPROT)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(UNIPROT = paste0(UNIPROT, collapse = ';')) %>%
    dplyr::distinct()

  return(bio_dataset_with_id)
}
