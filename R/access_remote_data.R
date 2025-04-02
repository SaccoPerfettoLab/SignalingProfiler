
#' Access Remote Data from PerfettoLab Server
#'
#' Downloads and reads a `.tsv` file hosted on the PerfettoLab server 
#' (https://perfettolab.bio.uniroma1.it/PerfettoLabData/SignalingProfiler).
#' This function is useful for accessing publicly available signaling datasets 
#' for downstream analyses in the context of SignalingProfiler.
#'
#' @usage access_remote_file(file, dir)
#'
#' @param file Character. The name of the file to download (e.g. "example.tsv").
#' @param dir Character. The subdirectory within the SignalingProfiler directory where the file is located (e.g. "PKN").
#'
#' @return A tibble containing the parsed contents of the downloaded TSV file.
#'
#' @examples
#' file <- "SerThrAtlas_PsP_regulons_filtered_99.tsv"
#' access_remote_file(file = file, dir = "PKN")
#'
#' @export
access_remote_file <- function(file, dir){
  
  base_url <- "https://perfettolab.bio.uniroma1.it/PerfettoLabData/SignalingProfiler"
  
  file_dir <- file.path(base_url, dir, file)
  
  file_path <- file.path(file_dir)
  
  response <- httr::GET(file_path)
  
  if (httr::status_code(response) == 200) {
    return(readr::read_tsv(httr::content(response, "text"),show_col_types = FALSE))
  } else {
    stop(paste("Unable to download from PerfettoLabServer:", file))
  }
}
