


#' preprocess_PKN
#'
#' @param omics_data list of all the omics data you have
#' @param organism string, 'mouse' or 'human'
#'
#' @return table of PKN with interactions involving expressed proteins in your system
#' @export
#'
#' @examples
preprocess_PKN <- function(omics_data, organism){

  # omics_data <- list(tr_toy_df, prot_toy_df, phospho_toy_df)
  # organism = 'mouse'

  if(organism == 'mouse'){
    PKN_table <- PKN_mouse
  }else if(organism == 'human'){
    PKN_table <- PKN_human
  }else{
    stop('Please provide a valid organism')
  }

  genes <- unique(unlist(lapply(omics_data, function(x){x$gene_name})))

  PKN_expressed <- PKN_table %>%
    dplyr::filter(ENTITYA %in% genes & ENTITYB %in% genes)

  return(PKN_expressed)
}
