#' Filter Prior Knowledge Network Based on Expressed Proteins
#'
#' Filters a Prior Knowledge Network (PKN) to retain only interactions
#' involving proteins expressed in the provided omics datasets.
#'
#' @param omics_data A list of data frames, each containing at least a `gene_name` column
#'   representing expressed proteins.
#' @param PKN_table A data frame representing the Prior Knowledge Network,
#'   derived from `choose_PKN()`.
#'
#' @return A filtered PKN table containing only interactions involving
#'   proteins detected in the provided omics datasets.
#'
#' @details
#' - The function extracts **unique gene names** from all omics datasets.
#' - It filters the **PKN_table** to retain:
#'   - Interactions where `ENTITYA` or `ENTITYB` is a detected protein.
#'   - Non-protein interactions (e.g., transcription factors, complexes) remain unchanged.
#'
#' @export
#'
#' @examples
#' data('tr_toy_df')
#' data('prot_toy_df')
#' data('phospho_toy_df')
#' omics_list <- list(tr_toy_df, prot_toy_df, phospho_toy_df)
#' filtered_PKN <- preprocess_PKN(omics_data = omics_list, PKN_table = PKN_human_atlas_dir)
#'
preprocess_PKN <- function(omics_data, PKN_table){

  # Extract and clean unique gene names from omics data
  genes <- unique(unlist(lapply(omics_data, `[[`, "gene_name"))) %>%
    stringr::str_split(";") %>%
    unlist() %>%
    stringr::str_replace_all("[^[:alnum:]]", "_") %>%
    unique()

  # Filter PKN to retain only relevant interactions
  PKN_expressed <- PKN_table %>%
    dplyr::filter((TYPEA == 'protein' & ENTITYA %in% genes) |
                    (TYPEB == 'protein' & ENTITYB %in% genes) |
                    (TYPEA != 'protein' & TYPEB != 'protein'))

  return(PKN_expressed)
}

#' Select a Prior Knowledge Network (PKN)
#'
#' Retrieves a predefined or custom Prior Knowledge Network (PKN)
#' based on user-defined parameters. The PKN is a structured dataset
#' representing known interactions among signaling molecules.
#'
#' @param organism  Character. The organism for which to retrieve the PKN.
#'. Must be either `"human"` or `"mouse"`.
#' @param with_atlas Logical. If `TRUE`, integrates inferred kinase-substrate
#'   interactions from the Ser/Thr and Tyr Kinome Atlas (PMIDs: 36631611, 38720073).
#'   Applicable only to human PKNs. Default: `FALSE`.
#' @param direct Logical. If `TRUE`, retains only direct interactions among biological entities;
#'  otherwise, includes indirect interactions. Default `TRUE`.
#' @param custom Logical. If `TRUE`, the user must provide a path to a custom
#'   PKN file via the `custom_path` parameter. Default: `FALSE`.
#' @param custom_path Character. The file path of the custom PKN, expected
#'   in tab-separated format. Must be provided if `custom = TRUE`. Default: `NULL`.
#'
#' @return A `data.frame` containing the selected Prior Knowledge Network (PKN).
#'   The columns may vary based on the chosen dataset but typically include
#'   source and target nodes, interaction types, and supporting evidence.
#'
#' @seealso [SignalingProfiler::PKN_human_dir],
#'   [SignalingProfiler::PKN_mouse_dir]
#'
#' @export
#'
#' @examples
#' # Retrieve the human PKN without kinase-substrate interactions and with direct interactions
#' human_pkn <- choose_PKN(organism = "human", with_atlas = FALSE, direct = TRUE)
#'
#' # Retrieve the mouse PKN with indirect interactions
#' mouse_pkn <- choose_PKN(organism = "mouse", direct = FALSE)
#'
#' \dontrun{
#' # Load a custom PKN from a file
#' custom_pkn <- choose_PKN(organism = "human", custom = TRUE, custom_path = "your_pkn.tsv")
#' }
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
        PKN_table <- SignalingProfiler::PKN_mouse_dir
      }else if(direct == FALSE){
        PKN_table <- SignalingProfiler::PKN_mouse_ind
      }else{
        stop('direct parameter must be TRUE or FALSE')
      }
    }else if(organism == 'human'){
      if(with_atlas == TRUE){
        if(direct == TRUE){
          PKN_table <- SignalingProfiler::PKN_human_atlas_dir
        }else if(direct == FALSE){
          PKN_table <- SignalingProfiler::PKN_human_atlas_ind
        }else{
          stop('direct parameter must be TRUE or FALSE')
        }
      }else{
        if(direct == TRUE){
          PKN_table <- SignalingProfiler::PKN_human_dir
        }else if(direct == FALSE){
          PKN_table <- SignalingProfiler::PKN_human_ind
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
