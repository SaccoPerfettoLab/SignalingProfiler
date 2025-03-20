#' Retrieve PhosphoScore Information
#'
#' This function compiles phosphorylation site regulatory data from SIGNOR and PhosphoSitePlus (PsP).
#' The resulting dataset identifies sites with consistent activation or inhibition effects.
#'
#' @param resources Character vector specifying the regulatory phosphosite databases to query. Allowed values are `"SIGNOR"` and `"PsP"`.
#' @param organism Character string specifying the organism. Supported values are `"human"` and `"mouse"`.
#' @param psp_reg_site_path Character string specifying the file path to the PhosphoSitePlus 'Regulatory Sites' data.
#'                          Required if `"PsP"` is included in `resources`.
#' @param psp_kin_sub_path Character. Path to the PhosphoSitePlus "Kinases_Substrate_Dataset" file manually downloaded from it.
#' @param only_activatory Logical, indicating whether to filter for only phosphosites that regulate protein activity. Default is `TRUE`.
#' @param aa_pos Logical, whether using in combination with Gene Symbol the amino acid residue position (`RESIDUE`)
#'               instead of 15-mer peptide (`SEQUENCE`) for phosphosite identification. Default is `FALSE`.
#'
#' @param local Logical, indicating whether to use a local or package-installed version of the SIGNOR dictionary. Default is `FALSE`.
#'
#' @details
#' - If `"SIGNOR"` is included in `resources`, the function extracts phosphorylation site interactions from SIGNOR.
#' - If `"PsP"` is included, the function extracts phosphorylation sites from PhosphoSitePlus.
#' - The function removes phosphosites with conflicting regulatory effects (e.g., activating in one context but inhibiting in another).
#' - If `aa_pos = TRUE`, phosphosites are identified using their residue position (e.g., `"MAPK1-T185"`),
#'   otherwise, their phosphopeptide context (e.g., `ABL1-RLMTGDTYTAHAGAK `) is used.
#'
#' @return A tibble containing regulatory phosphosites with the following columns:
#' - `PHOSPHO_KEY_GN_SEQ`: Unique phosphosite identifier (`Gene-Residue` or `Gene-Sequence`).
#' - `UNIPROT`: UniProt accession ID.
#' - `ACTIVATION`: Regulatory effect (`1` for activation, `-1` for inhibition).
#'
#' @details
#' - If `"SIGNOR"` is included in `resources`, the function extracts phosphorylation site interactions from SIGNOR.
#' - If `"PsP"` is included, the function extracts phosphorylation sites from PhosphoSitePlus.
#' - The function removes phosphosites with conflicting regulatory effects (e.g., activating in one context but inhibiting in another).
#'
#' @examples
#'
#' # Retrieve phosphosite information from SIGNOR
#' phosphoscore_info <- get_phosphoscore_info(resources = c("SIGNOR"),
#'                                            organism = "human")
#'
#' \dontrun{
#' # Retrieve phosphosite information from SIGNOR and PsP for human using sequence-based identifiers
#' phosphoscore_info <- get_phosphoscore_info(resources = c("SIGNOR", "PsP"),
#'                                            organism = "human",
#'                                            psp_reg_site_path = "./Regulatory_sites")
#'
#' # Retrieve phosphosite information using amino acid residue position instead of sequence
#' phosphoscore_info_residue <- get_phosphoscore_info(resources = c("SIGNOR", "PsP"),
#'                                                    organism = "human",
#'                                                    psp_reg_site_path = "./Regulatory_sites",
#'                                                    aa_pos = TRUE)
#' }
#'
#' @export
#'
get_phosphoscore_info <- function(resources = c('SIGNOR', 'PsP'),
                                  organism,
                                  psp_reg_site_path = NULL,
                                  psp_kin_sub_path = NULL,
                                  only_activatory = TRUE,
                                  aa_pos = FALSE,
                                  local = FALSE){

  phosphoscore_list <- list()

  if('SIGNOR' %in% resources){
    message('Querying SIGNOR database')
    SIGNOR <- signor_parsing(organism, direct = TRUE,
                             only_activatory = only_activatory)$interactions
    DBs <- SIGNOR %>% dplyr::filter(INTERACTION != '0')

    # Create the key for phosphosite unique identification
    if(aa_pos){
      DBs <- DBs %>%
        dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITYB, '-', RESIDUE)) %>%
        dplyr::rename('ENTITY' = 'ENTITYB') %>%
        dplyr::mutate(ENTITY = stringr::str_to_upper(stringr::str_replace_all(ENTITY, "[^[:alnum:]]", '_')))
    }else{
      DBs <- DBs %>%
        dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITYB, '-', SEQUENCE)) %>%
        dplyr::rename('ENTITY' = 'ENTITYB') %>%
        dplyr::mutate(ENTITY = stringr::str_to_upper(stringr::str_replace_all(ENTITY, "[^[:alnum:]]", '_')))
    }

    phos_mech <- DBs %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM) %>%
      dplyr::filter(grepl('*phos*', DBs$MECHANISM)) %>%
      dplyr::arrange(PHOSPHO_KEY_GN_SEQ) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ != '') %>%
      dplyr::distinct()

    phosphoscore_list[[length(phosphoscore_list) + 1]] <- phos_mech
  }

  if('PsP' %in% resources){
    if(is.null(psp_reg_site_path) | is.null(psp_kin_sub_path)){
      stop('Please provide a valid path to Regulatory_Sites and
           Kinases_Substrate_Dataset information of PhosphoSitePlus database')
    }
    message('Reading PsP Regulatory_Site file')

    psp_db <- psp_parsing(reg_site_path = psp_reg_site_path,
                          kin_sub_path = psp_kin_sub_path,
                          organism = organism,
                          only_activatory = only_activatory,
                          with_atlas = FALSE,
                          local = local)

    psp_db_phosphoscore <- psp_db %>%
      dplyr::select(ENTITY = ENTITYB, IDB, INTERACTION, RESIDUE, SEQUENCE) %>%
      dplyr::mutate(ENTITY = stringr::str_replace_all(ENTITY, "[^[:alnum:]]", '_')) %>%
      dplyr::mutate(ENTITY = ifelse(ENTITY == 'BCR_ABL', 'BCR_ABL1', ENTITY))

    if(organism == 'mouse'){
      psp_db_phosphoscore <- psp_db_phosphoscore %>%
        dplyr::mutate(ENTITY = ifelse(grepl("[^[:alnum:]]", psp_db_phosphoscore$ENTITY),
                                      ENTITY,
                                      stringr::str_to_title(psp_db_phosphoscore)))

    }else{
      psp_db_phosphoscore <- psp_db_phosphoscore %>%
        dplyr::mutate(ENTITY = stringr::str_to_upper(ENTITY))
    }

    if(aa_pos){
      phosphosite_mech <- psp_db_phosphoscore %>%
        dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITY, '-', RESIDUE),
                      IDB,
                      INTERACTION,
                      MECHANISM = 'phosphorylation') %>%
        dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM)
    }else{
      phosphosite_mech <- psp_db_phosphoscore %>%
        dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITY, '-', SEQUENCE),
                      IDB,
                      INTERACTION,
                      MECHANISM = 'phosphorylation') %>%
        dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM)
    }
    phosphoscore_list[[length(phosphoscore_list) + 1]] <- phosphosite_mech
  }

  # Dynamically bind the rows of all the entities
  phos_mech <- do.call("bind_rows", phosphoscore_list) %>% dplyr::distinct()

  # Add activation row
  phos_mech$ACTIVATION <- ''
  phos_mech$ACTIVATION[(phos_mech$MECHANISM == 'dephosphorylation' & phos_mech$INTERACTION == '-1') |
                         (phos_mech$MECHANISM == 'phosphorylation' & phos_mech$INTERACTION == '1')] <- '1'
  phos_mech$ACTIVATION[(phos_mech$MECHANISM == 'phosphorylation' & phos_mech$INTERACTION == '-1') |
                         (phos_mech$MECHANISM == 'dephosphorylation' & phos_mech$INTERACTION == '1')] <- '-1'

  # Select only phosphosites without controversial regulatory role in different contexts
  good_phos_df <- phos_mech %>%
    dplyr::distinct(PHOSPHO_KEY_GN_SEQ, ACTIVATION, .keep_all = TRUE) %>%
    dplyr::count(PHOSPHO_KEY_GN_SEQ) %>%
    dplyr::filter(n == 1) %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ)

  good_phos <- as.vector(unlist(good_phos_df))

  phosphoscore_table_final <- phos_mech %>%
    dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% good_phos) %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, ACTIVATION) %>%
    dplyr::rename(UNIPROT = IDB) %>%
    dplyr::distinct()

  return(phosphoscore_table_final)
}
