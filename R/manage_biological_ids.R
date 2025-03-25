#' Convert Gene Symbol to UniProt IDs
#'
#' Maps Gene Symbols in a biological dataset to their corresponding
#' UniProt IDs based on SIGNOR and PhosphoSitePlus databases.
#'
#' @param bio_dataset  Data frame, biological dataset containing a `gene_name` column.
#' @param organism Character string specifying the organism; either \code{"human"} or \code{"mouse"}.
#'
#' @return A data frame identical to `bio_dataset`, but with an additional `UNIPROT` column
#'   containing the corresponding UniProt ID(s) retrieved from SignalingProfiler databases.
#'   If multiple UniProt IDs are associated with a single gene, they are concatenated with `";"`.
#'
#' @details
#' - The function uses predefined SIGNOR-based protein tables:
#'   - `"human"` → `PKN_proteins_human_atlas`
#'   - `"mouse"` → `PKN_proteins_mouse`
#' - The function removes isoform-specific UniProt annotations (e.g., `P12345-1` → `P12345`).
#' - Any genes not found in the SIGNOR or PhosphoSitePlus database are excluded from the output.
#'
#' @export
#'
#' @examples
#' bio_data <- data.frame(gene_name = c("TP53", "KRAS", "MAPK1"))
#'
#' # Convert gene names to UniProt IDs for human proteins
#' converted_data <- convert_gene_name_in_uniprotid(bio_data, organism = "human")
#'
convert_gene_name_in_uniprotid <- function(bio_dataset, organism) {
  if (organism == 'human') {
    db <- get(data('PKN_proteins_human'))
  } else if (organism == 'mouse' | organism == 'hybrid') {
    db <- get(data('PKN_proteins_mouse'))
  } else{
    stop("Please provide a valid organism: 'human' or 'mouse'.")
  }

  bio_dataset_with_id <- dplyr::left_join(bio_dataset,
                                          db %>% dplyr::select(ID, ENTITY),
                                          by = c('gene_name' = 'ENTITY')) %>%
    dplyr::rename(UNIPROT = ID) %>%
    dplyr::filter(!is.na(UNIPROT)) %>%
    dplyr::mutate(UNIPROT = sub('-\\d', '', UNIPROT)) %>%
    tidyr::separate_rows(UNIPROT, sep = ';') %>%
    dplyr::distinct() %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(UNIPROT = paste0(unique(UNIPROT), collapse = ';')) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  return(bio_dataset_with_id)
}
