
#' Compute ProteoScore from Proteomic Data
#'
#' This function processes proteomic data to compute a **ProteoScore**, which represents
#' protein activity measures based on differential expression. The function first filters
#' significant proteins, converts gene names to UniProt IDs, and annotates molecular function.
#'
#' @param prot_df A dataframe containing proteomic data, including gene names and expression differences.
#' @param organism A character string specifying the organism. Supported values are `"human"` and `"mouse"`.
#'
#' @details
#' - The function filters proteins marked as significant (`significant == '+'`).
#' - Converts gene symbols to UniProt IDs using `convert_gene_name_in_uniprotid()`.
#' - Annotates molecular function using `molecular_function_annotation()`.
#'
#' @seealso [molecular_function_annotation]
#' @return A tibble containing proteomic activity measures with UniProt identifiers and molecular function annotations.
#'
#' @export
#'
#' @examples
#' data('prot_toy_df')
#' # Compute ProteoScore for human proteomic data
#' proteo_score <- activity_from_proteomics(prot_df = prot_toy_df,
#'                                          organism = "human")
#'
activity_from_proteomics <- function(prot_df, organism){

  message('Computing ProteoScore from proteomics data')

  prot_df <- prot_df %>%
    dplyr::filter(significant == '+') %>%
    dplyr::select(gene_name, difference)

  prot_df <- convert_gene_name_in_uniprotid(prot_df, organism)

  prot_df <- molecular_function_annotation(prot_df, organism)

  return(prot_df)
}




