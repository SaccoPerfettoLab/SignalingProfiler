#' Coverage of Inferred Proteins in Prior Knowledge Network
#'
#' This function evaluates the presence of inferred proteins from omics data
#' in a selected Prior Knowledge Network (PKN) based on the chosen organism.
#'
#' @param prediction_output `data.frame` containing inferred proteins with a column named `gene_name` and `mf` (molecular function).
#' @param organism Character. Either `'human'` or `'mouse'`.
#' @param with_atlas Logical. Default `TRUE`. If `TRUE`, considers PKN with Kinase-Substrate Atlas.
#' @param direct Logical. Default `FALSE`. If `TRUE`, considers PKN with only direct interactions.
#' @param custom Logical. Default `FALSE`. If `TRUE`, allows custom database selection.
#' @param custom_path Logical or character. Default `FALSE`. If a path is provided, it uses a custom PKN file.
#' @param report Logical. Default `FALSE`. If `TRUE`, writes results to a file (`report.txt`); otherwise, prints them.
#'
#' @return Prints coverage statistics or writes them to a file.
#'
#' @export
#'
#' @examples
#'
#' inferred_proteins <- data.frame(
#'   gene_name = c("TP53", "EGFR", "MYC"),
#'   mf = c("tf", "kin", "rec")
#' )
#' # Control coverage in human PKN with direct interactions
#' coverage_of_inferred_proteins_in_db(inferred_proteins, organism = "human", with_atlas = F, direct = T)
#'
#' # Control coverage in human PKN with indirect interactions
#' coverage_of_inferred_proteins_in_db(inferred_proteins, organism = "human", with_atlas = F, direct = F)
#'
coverage_of_inferred_proteins_in_db <- function(prediction_output,
                                                organism,
                                                with_atlas = TRUE,
                                                direct = FALSE,
                                                custom = FALSE,
                                                custom_path = FALSE,
                                                report = FALSE) {

  # Validate organism argument
  organism <- match.arg(c('human', 'mouse'))

  # Select the appropriate database
  pkn <- choose_database_for_building(organism, with_atlas, direct, custom, format = "igraph")

  # Initialize a list to store report content if needed
  report_lines <- c()

  # Iterate through molecular function categories
  molecular_functions <- c("tf", "kin", "phos", "other", "rec")

  for (m in molecular_functions) {

    # Filter proteins for the current molecular function
    prot_subset <- dplyr::filter(prediction_output, mf == m)

    # Count found proteins in the network
    found <- sum(igraph::V(pkn)$name %in% prot_subset$gene_name)
    total <- nrow(prot_subset)

    # Generate output message
    sentence <- sprintf("For %s molecular function: %d proteins found out of %d", m, found, total)
    message(sentence)

    # If report is enabled, store the result
    if (report) {
      report_lines <- c(report_lines, sentence)
    }
  }
  if (report) {
    writeLines(report_lines, "report.txt")
  }
}
