# this set of function manages go terms

#' Title
#'
#' @param inferred_proteins_dataset dataset of inferred proteins with UNIPROT
#'
#' @return dataset of inferred proteins with molecular function column
#' @export
#'
#' @examples
molecular_function_annotation <- function(inferred_proteins_dataset){

  #inferred_proteins_dataset <- output_uniprot
  # ADD A CHECK FOR UNIPROT COLUMN PRESENCE
  inferred_proteins_dataset <- inferred_proteins_dataset %>%
    dplyr::filter(!is.na(UNIPROT)) %>%
    tidyr::separate_rows('UNIPROT', sep = ';') %>%
    dplyr::distinct()

  #get uniprot annotation
  annotation <- UniprotR::GetProteinGOInfo(inferred_proteins_dataset$UNIPROT)

  # extract GO term for each protein
  ovl <- annotation$Gene.Ontology..molecular.function.
  pattern <- "(\\[.*?\\])"
  matches <- gregexpr(pattern, ovl)
  overlap <- regmatches(ovl, matches)
  GOterm <- lapply(overlap, function(x){stringr::str_sub(x, 2,-2)})

  # set ontology as MF
  GOSim::setOntology(ont = "MF", loadIC=TRUE, DIR=NULL)

  # get ancestors for each GO term
  GOanc <- GOSim::getAncestors()

  inferred_proteins_dataset$MF <- NA
  for(i in c(1:length(inferred_proteins_dataset$gene_name))){

    # transcription factor
    if(sum(unlist(lapply(GOterm[[i]], function(x){'GO:0140110' %in% GOanc[[x]] |
        'GO:0140110' == x})))>=1){
      inferred_proteins_dataset$MF[i] <- 'tf'
      # kinase
    }else if(sum(unlist(lapply(GOterm[[i]], function(x){'GO:0016301' %in% GOanc[[x]] |
        'GO:0016301' == x})))>=1){
      inferred_proteins_dataset$MF[i] <- 'kin'
      # phosphatase
    }else if(sum(unlist(lapply(GOterm[[i]], function(x){'GO:0004725' %in% GOanc[[x]] |
        'GO:0004725' == x})))>=1){
      inferred_proteins_dataset$MF[i] <- 'phos'
      # other types
    }else{
      inferred_proteins_dataset$MF[i] <- 'other'
    }
  }
  inferred_proteins_dataset <- inferred_proteins_dataset %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(UNIPROT = paste0(UNIPROT, collapse = ';')) %>%
    dplyr::ungroup()

  return(inferred_proteins_dataset)
}
