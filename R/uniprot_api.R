#' Retrieve UniProt Protein Information
#'
#' Queries the UniProt REST API to retrieve protein metadata for a given set of UniProt IDs.
#'
#' @param id_input Character vector. List of UniProt accession IDs.
#' @param batch_size Integer. Number of IDs per batch query to UniProt. Default: `400`.
#'
#' @return A `data.frame` containing UniProt metadata, including accession ID, protein name,
#'   Gene Name (Primary), organism, sequence length, and modification date.
#'
#' @export
#'
#' @examples
#' result <- query_uniprot_proteins(ids=c("P12345", "Q9Y6X1", "A0A0B4J2D5"))
#'
query_uniprot_proteins <- function(id_input, batch_size = 400) {
  ## Keep only ID input that are UNIPROT IDs
  id_input <- unique(id_input)
  id_input <-
    id_input[grepl(
      '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',
      id_input
    )]
  id_input <-
    id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id_input)]

  ## Define column names for output data frame
  column_names <-
    c(
      "Entry",
      "Reviewed",
      "Entry Name",
      "Protein names",
      "Gene Names (primary)",
      "Organism",
      "Length",
      "Sequence",
      "Date_modified"
    )

  ## Batch query processing
  for (i in seq(1, length(id_input), by = batch_size)) {
    batch_ids <- id_input[i:min(i + batch_size - 1, length(id_input))]

    query <-
      paste0("accession%3A",
             paste(batch_ids, collapse = "%20OR%20accession%3A"))
    url <-
      paste0(
        "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=",
        query
      )

    response <- httr::GET(url)
    content <- httr::content(response, "text", encoding = "UTF-8")

    batch_df <-
      readr::read_delim(
        content,
        delim = "\t",
        skip_empty_rows = TRUE,
        show_col_types = FALSE,
        locale = readr::locale(encoding = "UTF-8")
      )
    if (i == 1) {
      result_df <- batch_df
      colnames(result_df) <- column_names
    } else{
      result_df <- dplyr::bind_rows(result_df, batch_df)
    }
  }
  return(dplyr::distinct(result_df))
}


#' Harmonize UniProt IDs Including Isoforms
#'
#' Queries UniProt for metadata of given UniProt IDs, processing isoforms separately.
#'
#' @param id_input Character vector. List of UniProt accession IDs.
#' @param batch_size Integer. Number of IDs per batch query to UniProt. Default: 400.
#'
#' @return  A `data.frame` containing UniProt metadata for both standard proteins and isoforms.
#'
#' @seealso [query_uniprot_proteins]
#'
#' @export
#'
#' @examples
#' # Query metadata including isoforms
#' ids <- c("P12345", "Q9Y6X1-2", "A0A0B4J2D5")
#' result <- query_uniprot(ids)
#' print(result)
#'
#' # Query only proteins (without isoforms)
#' result_proteins <- query_uniprot(c("P12345", "Q9Y6X1"), batch_size = 200)
#' print(result_proteins)
#'
query_uniprot <- function(id_input, batch_size = 400) {
  ## Filter valid UniProt IDs
  id_input <-
    unique(unlist(stringr::str_split(id_input, pattern = ';')))
  id_input <-
    id_input[grepl(
      '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',
      id_input
    )]
  id_input <-
    id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606|VAR|REV',
                    id_input)]

  ## Separate isoforms for independent querying
  id_isoforms <- unique(id_input[grepl('-', id_input)])
  result_iso <- if (length(id_isoforms) > 0) {
    query_uniprot_proteins(id_isoforms, min(batch_size, length(id_isoforms)))
  } else {
    NULL
  }

  ## Query for proteins
  result_proteins <-
    query_uniprot_proteins(id_input, batch_size) # Query for not isoforms

  ## Merge Results
  if (exists('result_iso')) {
    return(dplyr::bind_rows(result_iso, result_proteins))
  }

  return(result_proteins)
}
