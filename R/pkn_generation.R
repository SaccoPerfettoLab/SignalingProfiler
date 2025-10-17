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
create_PKN <- function(
    database = c("SIGNOR","Omnipath","PsP","SerThr_Atlas"),
    organism = "human",
    direct   = FALSE,
    file_path = NULL,
    omnipath_resources = c("SignaLink3","PhosphoSite","SIGNOR"),
    psp_reg_site_path = NULL,
    psp_kin_sub_path  = NULL
){
  PKN_list <- list()
  
  # -- Internal helper: enforce consistent types before binding -----------------
  # Ensures: INTERACTION is character; DIRECT is logical; MECHANISM is character.
  .coerce_edge_types <- function(df){
    if (!"INTERACTION" %in% names(df)) df$INTERACTION <- NA_character_
    df$INTERACTION <- as.character(df$INTERACTION)
    
    if (!"DIRECT" %in% names(df)) df$DIRECT <- NA
    df$DIRECT <- as.logical(df$DIRECT)
    
    if ("MECHANISM" %in% names(df)) {
      df$MECHANISM <- as.character(df$MECHANISM)
    } else {
      df$MECHANISM <- NA_character_
    }
    df
  }
  # ---------------------------------------------------------------------------
  
  # Basic input checks
  if (("SIGNOR" %in% database | "PsP" %in% database) && is.null(organism)) {
    stop("Please specify an organism between 'human' and 'mouse'")
  }
  
  # SIGNOR
  if ("SIGNOR" %in% database) {
    message("Querying SIGNOR database")
    SIGNOR <- signor_parsing(organism, direct)$interactions
    SIGNOR$SOURCE <- "SIGNOR"
    SIGNOR <- .coerce_edge_types(SIGNOR)
    PKN_list[[length(PKN_list) + 1]] <- SIGNOR
  }
  
  # PhosphoSitePlus (optionally with Ser/Thr Atlas merge inside psp_parsing)
  if ("PsP" %in% database) {
    if (is.null(psp_reg_site_path) || is.null(psp_kin_sub_path)) {
      stop("Please provide paths to PhosphoSitePlus Regulatory_sites and Kinases_Substrates files.")
    } else {
      message("Parsing PhosphoSitePlus data")
      PsP <- psp_parsing(
        reg_site_path = psp_reg_site_path,
        kin_sub_path  = psp_kin_sub_path,
        organism      = organism,
        with_atlas    = "SerThr_Atlas" %in% database
      )
      PsP$SOURCE <- ifelse("SerThr_Atlas" %in% database, "PsP+Atlas", "PsP")
      PsP$DIRECT <- TRUE
      PsP <- .coerce_edge_types(PsP)
      PKN_list[[length(PKN_list) + 1]] <- PsP
    }
  }
  
  # OmniPath
  if ("Omnipath" %in% database) {
    if (is.null(omnipath_resources)) {
      stop("Please provide a vector of resources to download from Omnipath")
    } else {
      message("Querying OmniPath database using OmniPathR package")
      # If your omnipath_parsing supports 'organism', prefer passing organism = 9606.
      # omnipath <- omnipath_parsing(resources = omnipath_resources, organism = 9606)
      omnipath <- omnipath_parsing(resources = omnipath_resources)
      omnipath$SOURCE <- "OmniPath"
      omnipath$DIRECT <- FALSE
      omnipath <- .coerce_edge_types(omnipath)
      PKN_list[[length(PKN_list) + 1]] <- omnipath
    }
  }
  
  # Combine all interactions (safe after type coercion)
  PKN <- do.call(dplyr::bind_rows, PKN_list) |>
    dplyr::filter(!is.na(ENTITYA) & !is.na(ENTITYB))
  
  # Harmonize mouse gene names
  if (organism == "mouse") {
    PKN <- PKN |>
      dplyr::mutate(
        ENTITYA = ifelse(TYPEA == "protein", stringr::str_to_title(ENTITYA), ENTITYA),
        ENTITYB = ifelse(TYPEB == "protein", stringr::str_to_title(ENTITYB), ENTITYB)
      )
  }
  
  # Query UniProt for updated primary gene names
  message("Querying UniProt API for gene name updates")
  res <- query_uniprot(id_input = unique(c(PKN$IDA, PKN$IDB)), batch_size = 400)
  res_sub <- res |>
    dplyr::select(ID = Entry, ENTITY = `Gene Names (primary)`) |>
    dplyr::distinct()
  
  # Build a unified node table (collapse multi-mapped IDs per entity, then re-split)
  all_proteins <- tibble::tibble(
    ID = c(PKN$IDA, PKN$IDB),
    ENTITY = c(PKN$ENTITYA, PKN$ENTITYB),
    TYPE = c(PKN$TYPEA, PKN$TYPEB),
    DATABASE = c(PKN$DATABASEA, PKN$DATABASEB)
  ) |>
    dplyr::distinct(ID, TYPE, DATABASE, .keep_all = TRUE) |>
    dplyr::group_by(ENTITY) |>
    dplyr::mutate(
      ID = paste0(unique(ID), collapse = ";"),
      TYPE = paste0(unique(TYPE), collapse = ";"),
      DATABASE = paste0(unique(DATABASE), collapse = ";")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct() |>
    tidyr::separate_rows(ID, sep = ";")
  
  # Map UniProt IDs to entities; keep non-UniProt analytes (chemicals, complexes, fusions)
  analytes <- dplyr::inner_join(res_sub, all_proteins |> dplyr::select(-ENTITY), by = c("ID"))
  not_uniprot <- dplyr::anti_join(all_proteins, res_sub, by = "ID")
  analytes <- dplyr::bind_rows(analytes, not_uniprot) |>
    dplyr::distinct(ID, TYPE, DATABASE, .keep_all = TRUE)
  
  # Re-concatenate IDs per entity for the final analyte table
  analytes_UNI <- analytes |>
    dplyr::group_by(ENTITY) |>
    dplyr::reframe(ID = paste0(ID, collapse = ";"))
  
  analytes <- dplyr::inner_join(analytes |> dplyr::select(-ID), analytes_UNI, by = "ENTITY") |>
    dplyr::distinct()
  
  # -- Unified edges: re-split IDs and update ENTITYA/ENTITYB metadata ----------
  analytes_for_edges <- analytes |>
    tidyr::separate_rows(ID, sep = ";") |>
    dplyr::distinct()
  
  ENTITYA_update <- PKN |>
    dplyr::select(-c("ENTITYA","TYPEA","DATABASEA")) |>
    dplyr::inner_join(analytes_for_edges, by = c(IDA = "ID")) |>
    dplyr::rename(ENTITYA = "ENTITY", TYPEA = "TYPE", DATABASEA = "DATABASE")
  
  ENTITYB_update <- ENTITYA_update |>
    dplyr::select(-c("ENTITYB","TYPEB","DATABASEB")) |>
    dplyr::inner_join(analytes_for_edges, by = c(IDB = "ID")) |>
    dplyr::rename(ENTITYB = "ENTITY", TYPEB = "TYPE", DATABASEB = "DATABASE")
  
  PKN_final <- ENTITYB_update |>
    dplyr::relocate(ENTITYA, ENTITYB, INTERACTION)
  
  # Cleanup and keys
  PKN_final$SCORE <- NULL
  PKN_final$GROUP <- NULL
  PKN_final$MECHANISM[is.na(PKN_final$MECHANISM)] <- ""
  PKN_final$INTERACTION <- as.character(PKN_final$INTERACTION) # ensure char "1"/"-1"
  
  # Phospho key used downstream
  PKN_final$PHOSPHO_KEY_GN_SEQ <- ""
  PKN_final$PHOSPHO_KEY_GN_SEQ[!is.na(PKN_final$SEQUENCE)] <-
    paste0(PKN_final$ENTITYB[!is.na(PKN_final$SEQUENCE)], "-", PKN_final$SEQUENCE[!is.na(PKN_final$SEQUENCE)])
  
  # Canonicalize BCR-ABL; drop NA entities
  PKN_final <- PKN_final |>
    dplyr::filter(!is.na(ENTITYB) & !is.na(ENTITYA)) |>
    dplyr::mutate(
      ENTITYA   = ifelse(ENTITYA == "BCR/ABL fusion", "BCR-ABL", ENTITYA),
      IDA       = ifelse(ENTITYA == "BCR-ABL", "SIGNOR-FP6", IDA),
      TYPEA     = ifelse(ENTITYA == "BCR-ABL", "fusion protein", TYPEA),
      DATABASEA = ifelse(ENTITYA == "BCR-ABL", "SIGNOR", DATABASEA),
      ENTITYB   = ifelse(ENTITYB == "BCR/ABL fusion", "BCR-ABL", ENTITYB),
      IDB       = ifelse(ENTITYB == "BCR-ABL", "SIGNOR-FP6", IDB),
      TYPEB     = ifelse(ENTITYB == "BCR-ABL", "fusion protein", TYPEB),
      DATABASEB = ifelse(ENTITYB == "BCR-ABL", "SIGNOR", DATABASEB)
    ) |>
    dplyr::group_by(
      ENTITYA, ENTITYB, INTERACTION, IDA, IDB, MECHANISM,
      RESIDUE, SEQUENCE, DIRECT, TYPEA, DATABASEA, TYPEB, DATABASEB, PHOSPHO_KEY_GN_SEQ
    ) |>
    dplyr::reframe(SOURCE = paste0(unique(SOURCE), collapse = ";"))
  
  # Remove duplicates inflated by DIRECT flags; make DIRECT logical again
  PKN_final <- PKN_final |>
    dplyr::mutate(DIRECT = ifelse(DIRECT == TRUE, 1, 0)) |>
    dplyr::group_by(
      ENTITYA, ENTITYB, INTERACTION, IDA, IDB, MECHANISM, RESIDUE, SEQUENCE,
      TYPEA, DATABASEA, TYPEB, DATABASEB, PHOSPHO_KEY_GN_SEQ, SOURCE
    ) |>
    dplyr::reframe(DIRECT = sum(DIRECT)) |>
    dplyr::mutate(DIRECT = DIRECT == 1)
  
  # Edge/node selection by 'direct' mode
  if (direct) {
    PKN_prot <- PKN_final |>
      dplyr::filter(DIRECT == TRUE & TYPEA == "protein" & TYPEB == "protein" & !grepl("transcr|transl", MECHANISM))
    PKN_pheno <- PKN_final |>
      dplyr::filter(DIRECT == FALSE & (TYPEB == "phenotype" | TYPEA == "stimulus"))
    analytes <- analytes |>
      dplyr::filter(ENTITY %in% PKN_final$ENTITYA | ENTITY %in% PKN_final$ENTITYB)
  } else {
    PKN_final <- PKN_final |>
      dplyr::filter(DIRECT == TRUE | (DIRECT == FALSE & grepl("^transc", MECHANISM)))
    analytes <- analytes |>
      dplyr::filter(ENTITY %in% PKN_final$ENTITYA | ENTITY %in% PKN_final$ENTITYB)
  }
  
  # Harmonize entity strings (replace non-alnum, prefix numerics with 'X')
  analytes1 <- analytes |>
    dplyr::mutate(ENTITY = stringr::str_replace_all(ENTITY, "[^[:alnum:]]", "_"))
  analytes1$ENTITY[grepl("^\\d", analytes1$ENTITY)] <- paste0("X", analytes1$ENTITY[grepl("^\\d", analytes1$ENTITY)])
  analytes1 <- analytes1 |>
    dplyr::group_by(ENTITY) |>
    dplyr::reframe(
      ID = paste0(unique(ID), collapse = ";"),
      TYPE = paste0(unique(TYPE), collapse = ";"),
      DATABASE = paste0(unique(DATABASE), collapse = ";")
    ) |>
    dplyr::distinct()
  
  PKN_final1 <- PKN_final |>
    dplyr::mutate(
      ENTITYA = stringr::str_replace_all(ENTITYA, "[^[:alnum:]]", "_"),
      ENTITYB = stringr::str_replace_all(ENTITYB, "[^[:alnum:]]", "_")
    )
  PKN_final1$ENTITYA[grepl("^\\d", PKN_final1$ENTITYA)] <- paste0("X", PKN_final1$ENTITYA[grepl("^\\d", PKN_final1$ENTITYA)])
  PKN_final1$ENTITYB[grepl("^\\d", PKN_final1$ENTITYB)] <- paste0("X", PKN_final1$ENTITYB[grepl("^\\d", PKN_final1$ENTITYB)])
  
  # Build igraph
  PKN_graph <- igraph::graph_from_data_frame(d = PKN_final1, vertices = analytes1)
  
  final_list <- list(
    igraph_PKN = PKN_graph,
    interactions = PKN_final1,
    entities = analytes1
  )
  
  # Optional save
  if (!is.null(file_path)) {
    readr::write_rds(final_list, file_path)
  }
  return(final_list)
}
