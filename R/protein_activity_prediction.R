
#' Combine Footprint-Based and PhosphoScore Analyses
#'
#' This function integrates footprint-based analysis (`VIPER`) with
#' regulatory phosphosite-based (`PhosphoScore`) activity inference.
#' It computes a **final score** for inferred proteins by merging the results
#' of both approaches.
#'
#' @param footprint_output A dataframe containing inferred protein activity from the **VIPER footprint-based** method.
#' @param phosphoscore_df A dataframe containing inferred protein activity from the **PhosphoScore phosphosite-based** method.
#' @param analysis A string specifying the analysis type (`"tfea"` for transcription factor enrichment or `"ksea"` for kinase enrichment).
#'
#' @return A dataframe containing:
#' \describe{
#'   \item{`gene_name`}{The name of the inferred protein.}
#'   \item{`UNIPROT`}{The UniProt ID of the inferred protein.}
#'   \item{`mf`}{Molecular function (`"tf"` for transcription factors, `"kin"` or `"phos"` for kinases and phosphatases).}
#'   \item{`final_score`}{The combined activity score from both methods.}
#'   \item{`method`}{Indicates whether the score is derived from `"VIPER"`, `"PhosphoScore"`, or `"VIPER+PhosphoScore"`.}
#' }
#'
#' @examples
#' # Example footprint-based output
#' footprint_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   UNIPROT = c("P04637", "P01106"),
#'   mf = c("tf", "tf"),
#'   weightedNES = c(2.5, -1.8)
#' )
#'
#' # Example phosphoscore output
#' phosphoscore_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   UNIPROT = c("P04637", "P01106"),
#'   mf = c("tf", "tf"),
#'   phosphoscore = c(1.9, -1.5)
#' )
#'
#' # Combine footprint and phosphoscore analyses
#' combined_results <- combine_footprint_and_phosphoscore(footprint_data, phosphoscore_data, "tfea")
#' head(combined_results)
#'
#' @export
combine_footprint_and_phosphoscore <- function(footprint_output, phosphoscore_df, analysis) {

  message("** Combining Footprint and PhosphoScore Analyses **")

  # Validate analysis type
  analysis <- match.arg(analysis, c("tfea", "ksea"))

  # Filter PhosphoScore dataset based on analysis type
  if ("mf" %in% colnames(phosphoscore_df)) {
    phosphoscore_df <- dplyr::filter(phosphoscore_df,
                                     mf == ifelse(analysis == "tfea", "tf", c("kin", "phos")))
  }

  # Ensure footprint_output contains `weightedNES` (rename if necessary)
  if (!"weightedNES" %in% colnames(footprint_output)) {
    footprint_output <- dplyr::rename(footprint_output, weightedNES = NES)
    nes_was_renamed <- TRUE
  } else {
    nes_was_renamed <- FALSE
  }

  # Merge datasets
  merge_cols <- if ("mf" %in% colnames(footprint_output) & "mf" %in% colnames(phosphoscore_df)) {
    c("gene_name", "UNIPROT", "mf")
  } else {
    c("gene_name", "UNIPROT")
  }

  combined_df <- dplyr::full_join(footprint_output, phosphoscore_df, by = merge_cols)

  # Compute the final score and method
  combined_df <- combined_df %>%
    dplyr::mutate(
      final_score = dplyr::case_when(
        !is.na(weightedNES) & !is.na(phosphoscore) ~ rowMeans(dplyr::select(., weightedNES, phosphoscore), na.rm = TRUE),
        !is.na(weightedNES) & is.na(phosphoscore) ~ weightedNES,
        is.na(weightedNES) & !is.na(phosphoscore) ~ phosphoscore,
        TRUE ~ NA_real_
      ),
      method = dplyr::case_when(
        !is.na(weightedNES) & !is.na(phosphoscore) ~ "VIPER+PhosphoScore",
        !is.na(weightedNES) & is.na(phosphoscore) ~ "VIPER",
        is.na(weightedNES) & !is.na(phosphoscore) ~ "PhosphoScore",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(gene_name, UNIPROT, mf, final_score, method, everything()) %>%
    dplyr::arrange(dplyr::desc(final_score))

  # Restore original `NES` column name if necessary
  if (nes_was_renamed) {
    combined_df <- dplyr::rename(combined_df, NES = weightedNES)
  }

  return(combined_df)
}

#' Combine Activity Score and ProteoScore
#'
#' This function merges activity measures derived from **footprint-based** and
#' **PhosphoScore** techniques (`activity_score`) with significant
#' **proteomics** modulations (`proteo_score`), computing a **final score**
#' that integrates both data sources.
#'
#' @param activity_score A dataframe containing activity measures derived from phosphoproteomics.
#' @param proteo_score A dataframe containing activity measures derived from proteomics.
#'
#' @return A dataframe containing:
#' \describe{
#'   \item{`gene_name`}{The name of the inferred protein.}
#'   \item{`UNIPROT`}{The UniProt ID of the inferred protein.}
#'   \item{`mf`}{Molecular function (`"tf"` for transcription factors, `"kin"` or `"phos"` for kinases and phosphatases).}
#'   \item{`final_score`}{The combined activity score from phosphoproteomics and proteomics.}
#'   \item{`method`}{Indicates whether the score is derived from `"ActivityScore"`, `"ProteoScore"`, or both.}
#' }
#'
#' @examples
#' # Example phosphoproteomics activity scores
#' activity_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   UNIPROT = c("P04637", "P01106"),
#'   mf = c("tf", "tf"),
#'   final_score = c(2.5, -1.8)
#' )
#'
#' # Example proteomics-based activity scores
#' proteo_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   UNIPROT = c("P04637", "P01106"),
#'   mf = c("tf", "tf"),
#'   difference = c(1.9, -1.5)
#' )
#'
#' # Combine activity scores
#' combined_results <- combine_activityscore_proteoscore(activity_data, proteo_data)
#' head(combined_results)
#'
#' @export
combine_activityscore_proteoscore <- function(activity_score, proteo_score) {

  message("** Combining Activity Score and ProteoScore **")

  # Merge datasets based on common identifiers
  combined_score <- dplyr::full_join(
    activity_score, proteo_score,
    by = c("gene_name", "UNIPROT", "mf")
  ) %>%
    dplyr::rename(
      activity_score = final_score,
      proteo_score = difference
    )

  # Compute the final score based on available data
  combined_score <- combined_score %>%
    dplyr::mutate(
      final_score = dplyr::coalesce(activity_score, proteo_score),  # Use activity_score if available, otherwise proteo_score
      method = dplyr::case_when(
        !is.na(activity_score) & !is.na(proteo_score) ~ "ActivityScore+ProteoScore",
        !is.na(activity_score) ~ "ActivityScore",
        !is.na(proteo_score) ~ "ProteoScore",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-c(activity_score, proteo_score)) %>%
    dplyr::arrange(gene_name)

  return(combined_score)
}


