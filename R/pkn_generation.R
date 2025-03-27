#' Create a Prior Knowledge Network (PKN)
#'
#' This function constructs a Prior Knowledge Network (PKN) by integrating
#' interactions from user-specified databases, including SIGNOR, OmniPath,
#' PhosphoSitePlus (PsP), and the Ser/Thr Kinome Atlas.
#'
#' @param database Character vector. Specifies which databases to integrate. Options include:
#'   - `"SIGNOR"`: Regulatory interactions from SIGNOR.
#'   - `"Omnipath"`: Interactions from the OmniPath database.
#'   - `"PsP"`: Phosphorylation-based interactions from PhosphoSitePlus.
#'   - `"SerThr_Atlas"`: Kinome Atlas interactions (requires `"PsP"`).
#'   Default: `c("SIGNOR", "Omnipath", "PsP", "SerThr_Atlas")`
#'
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param direct Logical. If `TRUE`, retains only direct interactions between proteins
#'   or indirect interactions involving stimuli and phenotypes. Default: `FALSE`.
#' @param file_path Character. Optional. Path to save the output as an RDS file. Default: `NULL`.
#' @param omnipath_resources Character vector. List of OmniPath resources to retrieve interactions from.
#'   Default: `c("SignaLink3", "PhosphoSite", "SIGNOR")`.
#' @param psp_reg_site_path Character. Path to the PhosphoSitePlus "Regulatory_sites" file.
#'   Required if `"PsP"` is included in `database`.
#' @param psp_kin_sub_path Character. Path to the PhosphoSitePlus  "Kinases_Substrate_Dataset" file.
#'   Required if `"PsP"` is included in `database`.
#'
#' @return A `list` containing:
#'   - `igraph_PKN`: An `igraph` object representing the PKN.
#'   - `interactions`: A `data.frame` of parsed interactions.
#'   - `entities`: A `data.frame` of molecular entities included in the PKN.
#'
#' @details
#' The function:
#' - **Retrieves interactions** from the selected databases.
#' - **Maps interactions to a unified format**, ensuring consistent entity naming.
#' - **Resolves molecular complexes** by mapping OmniPath interactions to SIGNOR.
#' - **Updates UniProt IDs and gene names** by querying the UniProt API.
#' - **Builds an `igraph` object**, allowing network analysis.
#'
#' If `direct = TRUE`, only direct interactions (or indirect stimulus/phenotype relations) are retained.
#'
#' If `file_path` is provided, the function saves the output as an `.rds` file.
#'
#' @seealso [signor_parsing], [omnipath_parsing], [psp_parsing]
#'
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a human PKN integrating SIGNOR and PsP
#' pkn_data <- create_PKN(database = c("SIGNOR", "PsP", "SerThr_Atlas"))
#'
#' # Create a mouse PKN including PsP interactions
#' pkn_mouse <- create_PKN(
#'   database = c("SIGNOR", "PsP"),
#'   organism = "mouse",
#'   psp_reg_site_path = "./Kinases_Substrates",
#'   psp_kin_sub_path = "./Regulatory_sites")
#'
#' # Save the PKN to a file
#' create_PKN(
#'   database = c("SIGNOR", "Omnipath", "PsP"),
#'   file_path = "output/pkn_network.rds")
#' }

create_PKN <- function(database = c("SIGNOR", "Omnipath", "PsP", "SerThr_Atlas"),
                       organism = "human",
                       direct = FALSE,
                       file_path = NULL,
                       omnipath_resources = c("SignaLink3", "PhosphoSite", "SIGNOR"),
                       psp_reg_site_path = NULL,
                       psp_kin_sub_path = NULL) {

  PKN_list <- list()

  if (("SIGNOR" %in% database | "PsP" %in% database) & is.null(organism)) {
    stop("Please specify an organism between 'human' and 'mouse'")
  }

  if ("SIGNOR" %in% database) {
    message("Querying SIGNOR database")
    SIGNOR <- signor_parsing(organism, direct)$interactions
    SIGNOR$SOURCE <- "SIGNOR"
    PKN_list[[length(PKN_list) + 1]] <- SIGNOR
  }

  if ("PsP" %in% database) {
    if (is.null(psp_reg_site_path) | is.null(psp_kin_sub_path)) {
      stop("Please provide paths to PhosphoSitePlus Regulatory_sites and Kinases_Substrates files.")
    } else {
      message("Parsing PhosphoSitePlus data")
      PsP <- psp_parsing(reg_site_path = psp_reg_site_path,
                         kin_sub_path = psp_kin_sub_path,
                         organism = organism,
                         with_atlas = "SerThr_Atlas" %in% database)
      PsP$SOURCE <- ifelse("SerThr_Atlas" %in% database, "PsP+Atlas", "PsP")
      PsP$DIRECT <- TRUE
      PKN_list[[length(PKN_list) + 1]] <- PsP
    }
  }

  if ("Omnipath" %in% database) {
    if (is.null(omnipath_resources)) {
      stop('Plese provide a vector of resources to download from Omnipath')
    } else {
      message('Querying OmniPath database using OmniPathR package')
      omnipath <- omnipath_parsing(resources = omnipath_resources)
      omnipath$SOURCE <- "OmniPath"
      omnipath$DIRECT <- FALSE
      PKN_list[[length(PKN_list) + 1]] <- omnipath
    }
  }

  # Combine all interactions
  PKN <- do.call("bind_rows", PKN_list) %>%
    dplyr::filter(!is.na(ENTITYA) & !is.na(ENTITYB))

  # Harmonize mouse gene names
  if (organism == "mouse") {
    PKN <- PKN %>%
      dplyr::mutate(ENTITYA = ifelse(TYPEA == "protein", stringr::str_to_title(ENTITYA), ENTITYA),
                    ENTITYB = ifelse(TYPEB == "protein", stringr::str_to_title(ENTITYB), ENTITYB))
  }

  # Query UniProt for updated gene names
  message("Querying UniProt API for gene name updates")
  res <- query_uniprot(id_input = unique(c(PKN$IDA, PKN$IDB)), batch_size = 400)
  res_sub <- res %>% dplyr::select(ID = Entry, ENTITY = `Gene Names (primary)`) %>% dplyr::distinct()

  # ** Create a unified node table **
  all_proteins <- tibble::tibble(ID = c(PKN$IDA, PKN$IDB),
                                 ENTITY = c(PKN$ENTITYA, PKN$ENTITYB),
                                 TYPE = c(PKN$TYPEA, PKN$TYPEB),
                                 DATABASE = c(PKN$DATABASEA, PKN$DATABASEB)) %>%
    dplyr::distinct(ID, TYPE, DATABASE, .keep_all = TRUE)

  # If a gene name has multiple UNIPROT, TYPES or DATABASE are collapsed in one string
  all_proteins <- all_proteins %>% dplyr::group_by(ENTITY) %>%
    dplyr::mutate(ID = paste0(unique(ID), collapse = ';'),
                  TYPE = paste0(unique(TYPE), collapse = ';'),
                  DATABASE = paste0(unique(DATABASE), collapse = ';')) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  # Since we want to use the UNIPROT for mapping, we reseparate the rows on the ID
  # obtaining a table which separates only UNIPROTs like this:
  #   P33800 ACAB1 protein;complex UNIPROT;SIGNOR
  #   SIGNOR-3 ACAB1 protein;complex UNIPROT;SIGNOR
  all_proteins <- all_proteins %>% tidyr::separate_rows(ID, sep = ';')

  analytes <- dplyr::inner_join(res_sub, all_proteins %>% dplyr::select(-ENTITY), by = c("ID"))

  # Take the proteins not found in UNIPROT because are chemicals, complex, fusion proteins
  not_uniprot <- dplyr::anti_join(all_proteins, res_sub, by = "ID")

  # Merge them
  analytes <- dplyr::bind_rows(analytes, not_uniprot)

  # Check for duplicates and remove Ubiquitin proteins
  # because they have the same UNIPROT ID
  # keeping only the entry with all gene names
  analytes %>% dplyr::distinct(ID, TYPE, DATABASE, .keep_all = TRUE) -> analytes

  # reconcatenate the UNIPROTs and recreate analytes table like before
  # P33800;SIGNOR-3 PIPPO protein;complex UNIPROT;SIGNOR
  analytes_UNI <- analytes %>%
    dplyr::group_by(ENTITY) %>%
    dplyr::reframe(ID = paste0(ID, collapse =';'))

  analytes <- dplyr::inner_join(analytes %>%
                                  dplyr::select(-ID),
                                analytes_UNI,
                                by = 'ENTITY') %>%
    dplyr::distinct()

  # ** Create a unified edges table **

  # Re-separate UNIPROT IDs
  analytes_for_edges <- analytes %>%
    tidyr::separate_rows(ID, sep = ';') %>%
    dplyr::distinct()

  # Update ENTITYA, TYPEA, DATABASEA with Primary Gene Name
  # and concatenation of types and databases
  ENTITYA_update <- dplyr::inner_join(PKN %>%
                                        dplyr::select(-c('ENTITYA', 'TYPEA', 'DATABASEA')),
                                      analytes_for_edges,
                                      by = c('IDA' = 'ID')) %>%
    dplyr::rename('ENTITYA' = 'ENTITY', 'TYPEA' = 'TYPE', 'DATABASEA' = 'DATABASE')

  # update ENTITYB, TYPEB, DATABASEB with Primary Gene Name and concatenation of types and dbs
  ENTITYB_update <- dplyr::inner_join(ENTITYA_update %>%
                                        dplyr::select(-c('ENTITYB', 'TYPEB', 'DATABASEB')),
                                      analytes_for_edges,
                                      by = c('IDB' = 'ID')) %>%
    dplyr::rename('ENTITYB' = 'ENTITY', 'TYPEB' = 'TYPE', 'DATABASEB' = 'DATABASE')

  # Parse interactions file
  PKN_final <- ENTITYB_update %>%
    dplyr::relocate(ENTITYA, ENTITYB, INTERACTION)

  PKN_final$SCORE <- NULL
  PKN_final$GROUP <- NULL
  PKN_final$MECHANISM[is.na(PKN_final$MECHANISM)] <- ''
  PKN_final$INTERACTION <- PKN_final$INTERACTION %>% as.character #setting 1 and -1 as char

  # Add key Phospho-Key-GeneName-Sequence necessary in SignalingProfiler
  PKN_final$PHOSPHO_KEY_GN_SEQ <- ''
  PKN_final$PHOSPHO_KEY_GN_SEQ[is.na(PKN_final$SEQUENCE) == FALSE] <- paste0(PKN_final$ENTITYB[is.na(PKN_final$SEQUENCE) == FALSE], '-',
                                                                             PKN_final$SEQUENCE[is.na(PKN_final$SEQUENCE) == FALSE])
  # ** Final adjustments ** #
  #  Modify BCR-ABL format of PsP in BCR-ABL format of SIGNOR
  # remove NA as gene names (immunoglobulins) because interfere with any filtering on the key
  PKN_final <- PKN_final %>% dplyr::filter(!is.na(ENTITYB) & !is.na(ENTITYA))

  PKN_final  <-  PKN_final %>% dplyr::mutate(ENTITYA = ifelse(ENTITYA == 'BCR/ABL fusion', 'BCR-ABL', ENTITYA),
                                             IDA = ifelse(ENTITYA == 'BCR-ABL', 'SIGNOR-FP6', IDA),
                                             TYPEA = ifelse(ENTITYA == 'BCR-ABL', 'fusion protein', TYPEA),
                                             DATABASEA = ifelse(ENTITYA == 'BCR-ABL', 'SIGNOR', DATABASEA),

                                             ENTITYB = ifelse(ENTITYB == 'BCR/ABL fusion', 'BCR-ABL', ENTITYB),
                                             IDB = ifelse(ENTITYB == 'BCR-ABL', 'SIGNOR-FP6', IDB),
                                             TYPEB = ifelse(ENTITYB == 'BCR-ABL', 'fusion protein', TYPEB),
                                             DATABASEB = ifelse(ENTITYB == 'BCR-ABL', 'SIGNOR', DATABASEB))

  # Remove NA as gene names (Immunoglobulins) also in the nodes table
  analytes <- analytes %>% dplyr::filter(!is.na(ENTITY))
  analytes <- analytes %>% dplyr::filter(ENTITY != 'BCR/ABL fusion')

  # Group according to the source
  PKN_final <- PKN_final %>% dplyr::group_by(ENTITYA,
                                             ENTITYB,
                                             INTERACTION,
                                             IDA,
                                             IDB,
                                             MECHANISM,
                                             RESIDUE,
                                             SEQUENCE,
                                             DIRECT,
                                             TYPEA, DATABASEA,
                                             TYPEB, DATABASEB,
                                             PHOSPHO_KEY_GN_SEQ) %>%
    dplyr::reframe(SOURCE = paste0(unique(SOURCE), collapse = ';'))

  # Remove duplicated interactions because of DIRECT attributes
  PKN_final <- PKN_final %>%
    dplyr::mutate(DIRECT = ifelse(DIRECT == TRUE, 1, 0)) %>%
    dplyr::group_by(ENTITYA, ENTITYB, INTERACTION, IDA, IDB,
                    MECHANISM, RESIDUE, SEQUENCE,
                    TYPEA, DATABASEA, TYPEB, DATABASEB, PHOSPHO_KEY_GN_SEQ, SOURCE) %>%
    dplyr::reframe(DIRECT = sum(DIRECT)) %>%
    dplyr::mutate(DIRECT = ifelse(DIRECT == 1, TRUE, FALSE))

  # This modality includes only direct interactions between proteins
  # or indirect interactions between stimuli and proteins or proteins and
  # phenotypes

  if(direct){

    # Edges table
    PKN_prot <- PKN_final %>% dplyr::filter(DIRECT == TRUE &
                                              TYPEA == 'protein' &
                                              TYPEB == 'protein' &
                                              !grepl('transcr|transl', MECHANISM))

    PKN_pheno <- PKN_final %>%
      dplyr::filter(DIRECT == FALSE & (TYPEB == 'phenotype' | TYPEA == 'stimulus'))

    # Nodes table
    analytes <- analytes %>%
      dplyr::filter(ENTITY %in% PKN_final$ENTITYA | ENTITY %in% PKN_final$ENTITYB)

    # If we want the INDIRECT relation, keep only 'transcriptional regulations'
  }else{

    PKN_final <- PKN_final %>%
      dplyr::filter(DIRECT == TRUE | (DIRECT == FALSE & grepl('^transc', MECHANISM)))

    # Nodes table
    analytes <- analytes %>%
      dplyr::filter(ENTITY %in% PKN_final$ENTITYA | ENTITY %in% PKN_final$ENTITYB)
  }

  # ** Harmonize elements ** #
  # This chunck transforms '-' and '/' in underscores
  # and if ENTITIES start with a number concatenate 'X' to the name
  analytes1 <- analytes %>% dplyr::mutate(ENTITY = stringr::str_replace_all(ENTITY, "[^[:alnum:]]", '_'))
  analytes1$ENTITY[grepl('^\\d', analytes1$ENTITY)] <- paste0('X', analytes1$ENTITY[grepl('^\\d', analytes1$ENTITY)])

  analytes1 <- analytes1 %>% dplyr::group_by(ENTITY) %>%
    dplyr::reframe(ID = paste0(unique(ID), collapse =';'),
                   TYPE = paste0(unique(TYPE), collapse = ';'),
                   DATABASE = paste0(unique(DATABASE), collapse = ';')) %>%
    dplyr::distinct()

  PKN_final1 <- PKN_final %>% dplyr::mutate(ENTITYA = stringr::str_replace_all(ENTITYA, "[^[:alnum:]]", '_'),
                                            ENTITYB = stringr::str_replace_all(ENTITYB, "[^[:alnum:]]", '_'))
  PKN_final1$ENTITYA[grepl('^\\d', PKN_final1$ENTITYA)] <- paste0('X', PKN_final1$ENTITYA[grepl('^\\d', PKN_final1$ENTITYA)])
  PKN_final1$ENTITYB[grepl('^\\d', PKN_final1$ENTITYB)] <- paste0('X', PKN_final1$ENTITYB[grepl('^\\d', PKN_final1$ENTITYB)])

  # ** Create an igraph object **
  PKN_graph <- igraph::graph_from_data_frame(d = PKN_final1, vertices = analytes1)

  final_list <- list(igraph_PKN = PKN_graph,
                     interactions = PKN_final1,
                     entities = analytes1)

  if(!is.null(file_path)){
    write_rds(final_list, file_path)
  }
  return(final_list)
}
