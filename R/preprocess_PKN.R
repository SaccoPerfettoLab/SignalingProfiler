


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
  genes <- stringr::str_replace_all(genes, "[^[:alnum:]]", "_")

  PKN_expressed <- PKN_table %>%
    dplyr::filter((TYPEA == 'protein' & ENTITYA %in% genes) |
                    (TYPEB == 'protein' & ENTITYB %in% genes) |
                    (TYPEA != 'protein' & TYPEB != 'protein'))

  return(PKN_expressed)
}


#' choose_PKN
#'
#' @param organism string, 'mouse' or 'human'
#' @param with_atlas Boolean value, default TRUE, if FALSE you exclude the relations in the PKN integrated with kinome atlas data
#' @param custom Boolean value, default FALSE, if TRUE you can provide the path for your custom PKN
#' @param custom_path string, specify the path of the custom PKN in tab separated format
#' @param direct Boolean value, default FALSE, if TRUE you choose only direct interactions between entities
#'
#' @return PKN_table, table of PKN
#' @export
#'
#' @examples
choose_PKN <- function(organism,
                       with_atlas = TRUE,
                       direct = FALSE,
                       custom = FALSE,
                       custom_path = NULL){

  if(custom == TRUE){
    message(paste0('Reading your custom ', organism, ' PKN'))
    if(is.null(custom_path)){
      stop('Please provide a path for your custom PKN')
    }else{
      PKN_table <- readr::read_tsv(custom_path)
    }
  }else{
    if(organism == 'mouse'){
      message('Ignoring with atlas parameter since it is only for human')
      if(direct == TRUE){
        PKN_table <- PKN_mouse_dir
      }else if(direct == FALSE){
        PKN_table <- PKN_mouse_ind
      }else{
        stop('direct parameter must be TRUE or FALSE')
      }
    }else if(organism == 'human'){
      if(with_atlas == TRUE){
        if(direct == TRUE){
          PKN_table <- PKN_human_atlas_dir
        }else if(direct == FALSE){
          PKN_table <- PKN_human_atlas_ind
        }else{
          stop('direct parameter must be TRUE or FALSE')
        }
      }else{
        if(direct == TRUE){
          PKN_table <- PKN_human_dir
        }else if(direct == FALSE){
          PKN_table <- PKN_human_ind
        }else{
          stop('direct parameter must be TRUE or FALSE')
        }
      }
    }else{
      stop('Please provide a valid organism')
    }
  }

  return(PKN_table)
}
