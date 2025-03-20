#' Preprocess Network Data for Phenoscore Analysis
#'
#' This function processes proteomics and phosphoproteomics data to generate
#' a filtered interaction network from proteins to phenotypes based on expressed proteins.
#' It writes the input data to temporary files, runs a Python script
#' for processing, and returns the cleaned SIGNOR network.
#'
#' @param proteomics A dataframe containing processed proteomics data for SignalingProfiler.
#' @param phospho A dataframe containing processed phosphoproteomics data for SignalingProfiler.
#' @param local Logical. If `TRUE`, for *development purposes*, it uses a local path (`./inst/`) for the Python script;
#'   otherwise, the function attempts to find the script in the package library path.
#'
#' @return A dataframe containing the cleaned and filtered SIGNOR interaction network.
#'
#' @details
#' The function writes the provided proteomics and phosphoproteomics data to TSV files,
#' executes a Python script for processing (available at `python/script.py`),
#' and retrieves the output. If the output file (`Global_result_final_table_minimized.txt`)
#' is successfully generated, it is read and returned
#' as a dataframe. If the process fails, an error is raised.
#'
#' @export
#'
#' @examples
#' data('toy_prot_df')
#' data('toy_phospho_df')
#'
#' signor_network <- phenoscore_network_preprocessing(proteomics_data = toy_prot_df,
#'                                                    phospho_data = toy_phos_df)
phenoscore_network_preprocessing <- function(proteomics, phospho,
                                             local = FALSE){

  if(local == TRUE){
    path_package <- './inst/'
  }else{
    path_package <- paste0(.libPaths(), '/SignalingProfiler/')
  }

  # Set HOME dir according to OS
  if (.Platform$OS.type == "windows") {
    home_dir <- normalizePath(Sys.getenv("USERPROFILE"), winslash = "/")
  } else {
    home_dir <- normalizePath(Sys.getenv("HOME"), winslash = "/")
  }

  readr::write_tsv(proteomics, paste0(home_dir, '/proteomics.tsv'))
  readr::write_tsv(phospho, paste0(home_dir, '/phosphoproteomics.tsv'))

  for(path in path_package){
    result <<- tryCatch({
      reticulate::py_run_file(paste0(path, "/python/script.py"))
    }, error = function(e) {
      NA # Return NA if an error occurs
    })
  }

  if(file.exists(paste0(home_dir, '/Global_result_final_table_minimized.txt'))){
    signor_filtered <- readr::read_tsv(paste0(home_dir, '/Global_result_final_table_minimized.txt'),
                                       show_col_types = FALSE)
    file.remove(paste0(home_dir, '/Global_result_final_table_minimized.txt'))
    file.remove(paste0(home_dir, '/proteomics.tsv'))
    file.remove(paste0(home_dir, '/phosphoproteomics.tsv'))
    return(signor_filtered)
  }else{
    stop('Something went wrong with Python calling from R!')
  }
}
