


#' preprocess_PKN
#'
#' @param omics_data list of all the omics data you have
#' @param PKN_table table derived by choose_PKN
#'
#' @return table of PKN with interactions between proteins expressed in your system
#' @export
#'
#' @examples
preprocess_PKN <- function(omics_data, PKN_table){

  # omics_data <- list(tr_toy_df, prot_toy_df, phospho_toy_df)
  # organism = 'mouse'

  genes <- unique(unlist(lapply(omics_data, function(x){x$gene_name})))

  PKN_expressed <- PKN_table %>%
    dplyr::filter(ENTITYA %in% genes & ENTITYB %in% genes)

  return(PKN_expressed)
}


#' choose_PKN
#'
#' @param organism string, 'mouse' or 'human'
#' @param with_atlas Boolean value, default FALSE, if TRUE you can use PKN integrated with kinome atlas data
#' @param custom Boolean value, default FALSE, if TRUE you can provide the path for your custom PKN
#' @param custom_path string, specify the path of the custom PKN in tab separated format
#'
#' @return PKN_table, table of PKN
#' @export
#'
#' @examples
choose_PKN <- function(organism,
                       with_atlas = FALSE,
                       custom = FALSE,
                       custom_path = NULL){

  if(organism == 'mouse'){
    if(with_atlas == TRUE){
      stop('If organism is \'mouse\' with_atlas parameter must be FALSE')
    }else{
      if(custom == TRUE){
        message(paste0('Reading your custom ', organism, ' PKN'))
        PKN_table <- readr::read_tsv(custom_path)

      }else{
        PKN_table <- PKN_mouse
      }
    }
  }else if(organism == 'human'){
    if(with_atlas == TRUE){
      PKN_table <- PKN_human_atlas
    }else{
      if(custom == TRUE){
        message(paste0('Reading your custom ', organism, ' PKN'))
        PKN_table <- readr::read_tsv(custom_path)
      }else{
        PKN_table <- PKN_human
      }
    }
  }else{
    stop('Please provide a valid organism')
  }

  return(PKN_table)
}
