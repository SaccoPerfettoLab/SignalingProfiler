#' Retrieve UniProt Protein Information
#'
#' Queries the UniProt REST API to retrieve protein metadata for a given set of UniProt IDs.
#'
#' @param id_input Character vector. List of UniProt accession IDs.
#' @param batch_size Integer. Number of IDs per batch query to UniProt. Default: `90`.
#'
#' @return A `data.frame` containing UniProt metadata, including accession ID, protein name,
#'   Gene Name (Primary), organism, sequence length, and modification date.
#'
#' @export
#'
#' @examples
#' result <- query_uniprot_proteins(id_input=c("P12345", "Q9Y6X1", "A0A0B4J2D5"))
#' print(result)
#' 
chunk_vec <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

query_uniprot_proteins <- function(id_input, batch_size = 90){
  
  # Keep only ID input that are UNIPROT IDs
  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id_input)]
  if (length(id_input) == 0) return(tibble::tibble())
  
  
  chunks <- chunk_vec(id_input, min(batch_size, 90))
  
  header <- c("Entry","Reviewed","Entry Name","Protein names",
              "Gene Names (primary)","Organism","Length","Sequence","Date_modified")
  out <- vector("list", length(chunks))
  k <- 1
  
  for (ids in chunks) {
    q <- paste0('accession%3A', paste0(ids, collapse = '%20OR%20accession%3A'))
    url <- paste0(
      'https://rest.uniprot.org/uniprotkb/stream?',
      'fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified',
      '&format=tsv&query=', q
    )
    
    resp <- httr::GET(url)
    txt  <- as.character(httr::content(resp, "text"))
    
    #
    if (!grepl("^Entry\\t", txt)) {
      warning("UniProt error on chunk: ", substr(txt, 1, 120))
      next
    }
    
    
    df <- readr::read_delim(
      txt,
      delim = '\t',
      skip_empty_rows = TRUE,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      locale = readr::locale(encoding = "UTF-8")
    )
    
    
    if ("Date sequence modified" %in% names(df) && !("Date_modified" %in% names(df))) {
      names(df)[names(df) == "Date sequence modified"] <- "Date_modified"
    }
    
    
    header <- c("Entry","Reviewed","Entry Name","Protein names",
                "Gene Names (primary)","Organism","Length","Sequence","Date_modified")
    missing_cols <- setdiff(header, names(df))
    if (length(missing_cols)) {
      for (mc in missing_cols) df[[mc]] <- NA_character_
    }
    df <- df[, header]
    
    
    suppressWarnings({
      df$Length <- as.integer(df$Length)
    })
    
    
    out[[k]] <- df
    k <- k + 1
  }
  
  dplyr::bind_rows(out) %>% dplyr::distinct()
}

#' Harmonize UniProt IDs Including Isoforms
#'
#' Queries UniProt for metadata of given UniProt IDs, processing isoforms separately.
#'
#' @param id_input Character vector. List of UniProt accession IDs.
#' @param batch_size Integer. Number of IDs per batch query to UniProt. Default: 90.
#'
#' @return  A `data.frame` containing UniProt metadata for both standard proteins and isoforms.
#'
#' @seealso [query_uniprot_proteins]
#'
#' @export
#'
#' @examples
#' # Toy Query
#' ids <- c("P12345", "Q9Y6X1", "A0A0B4J2D5")
#' result <- query_uniprot(ids)
#' print(result)
query_uniprot <- function(id_input, batch_size = 90) {
  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606|VAR|REV', id_input)]
  if (length(id_input) == 0) return(tibble::tibble())
  
  id_iso  <- unique(id_input[grepl('-', id_input)]) # Query for isoforms
  id_core <- setdiff(id_input, id_iso) # Query for not isoforms
  
  res_iso  <- if (length(id_iso))  query_uniprot_proteins(id_iso,  batch_size = min(batch_size, 90)) else NULL
  res_core <- if (length(id_core)) query_uniprot_proteins(id_core, batch_size = min(batch_size, 90)) else NULL
  
  result_total <- dplyr::bind_rows(res_iso, res_core) %>% dplyr::distinct()
  
  return(result_total)
}