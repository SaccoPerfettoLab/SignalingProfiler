
#' Filters PKN for analytes present in experimental data
#'
#' @param PKN Igraph object containing the built-in Prior Knowledge Network
#' @param transcriptomics A tibble containing transcriptomics data
#' @param proteomics A tibble containing proteomics data
#' @param phosphoproteomics A tibble containing phosphoproteomics data
#'
#' @return igraph object containing only analytes present in experimental data
#' @export
#'
#' @examples
pkn_filtering_for_expression <- function(PKN, transcriptomics, proteomics, phosphoproteomics){
  genes <- c(transcriptomics$gene_name, proteomics$gene_name, phosphoproteomics$gene_name)
  nodes_PKN <- igraph::as_data_frame(PKN, what = c('vertices')) %>% dplyr::filter(ENTITY %in% toupper(genes))
  vertices <- nodes_PKN$name
  PKN_expressed <- igraph::subgraph(PKN, igraph::V(PKN)[name %in% vertices])
  return(PKN_expressed)
}
