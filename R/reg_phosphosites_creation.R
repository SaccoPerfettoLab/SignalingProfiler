#' Retrieve regulatory phosphosites without upstream kinases
#'
#' This function keeps only phosphosites without an upstream kinase, hence
#' absent in PhosphoSitePlus prior knowledge network.
#'
#' @param organism Character string specifying the organism. Supported values are `"human"` and `"mouse"`.
#' @param reg_site_path Character string specifying the file path to the PhosphoSitePlus 'Regulatory Sites' data.
#'                          Required if `"PsP"` is included in `resources`.
#' @param only_activatory Logical, indicating whether to filter for only phosphosites that regulate protein activity. Default is `TRUE`.
#' @param local Logical, indicating whether to use a local or package-installed version of the SIGNOR dictionary. Default is `FALSE`.
#' @param psp_db PhosphoSitePlus database returned from the `psp_parsing` function.
#'
#' @return A tibble containing regulatory phosphosites with the following columns:
#' - `ENTITY`: Protein gene symbol.
#' - `IDB`: UniProt accession ID.
#' - `RESIDUE`: Modified residue (aminoacid and position).
#' - `INTERACTION`: Regulatory effect (`1` for activation, `-1` for inhibition).
#'
#'
integrate_psp_without_kin <- function(organism,
                                      reg_site_path,
                                      psp_db,
                                      only_activatory = TRUE,
                                      local = FALSE){

  regulatory_sites <- tibble::tibble(read.delim2(reg_site_path, sep ='\t', skip = 2))
  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == organism & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  # Read SIGNOR dictionary
  path_package <- ifelse(local == TRUE, './inst/',
                         paste0(.libPaths()[1], '/SignalingProfiler/'))


  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))


  # PSP
  phosphoscore_table <- regulatory_org %>%
    dplyr::select(c('ENTITY' = 'PROTEIN',
                    'IDB' = 'ACC_ID',
                    'RESIDUE' = 'MOD_RSD',
                    'SEQUENCE' = 'SITE_...7_AA',
                    'ON_FUNCTION')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'S', 'Ser')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'T', 'Thr')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'Y', 'Tyr')) %>%
    dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE)) %>%
    tidyr::separate_rows('ON_FUNCTION', sep = '; ')

  phosphoscore_table <- dplyr::left_join(phosphoscore_table,
                                         phospho2SIGNOR,
                                         by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  # Retrieve regulated phosphosites without an upstream kinases
  phosphoscore_table <- dplyr::anti_join(phosphoscore_table,
                                         psp_db,
                                         by = c('ENTITY' = 'ENTITYB',
                                                'RESIDUE',
                                                'SEQUENCE'))

  phosphoscore_table <- phosphoscore_table %>%
    dplyr::mutate(ENTITY = stringr::str_replace_all(ENTITY, "[^[:alnum:]]", '_')) %>%
    dplyr::mutate(ENTITY = ifelse(ENTITY == 'BCR_ABL', 'BCR_ABL1', ENTITY))

  if(only_activatory){
    phosphoscore_table <- phosphoscore_table %>%
      dplyr::filter(grepl('activity', SIGNOR)) %>%
      dplyr::select(ENTITY, IDB, RESIDUE, SEQUENCE, INTERACTION = sign)
  }else{
    phosphoscore_table <- phosphoscore_table %>%
      dplyr::select(ENTITY, IDB, RESIDUE, SEQUENCE, INTERACTION = sign)
  }

  return(phosphoscore_table)
}

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
    SIGNOR <- signor_parsing(organism,
                             direct = TRUE,
                             only_activatory = only_activatory)$interactions
    DBs <- SIGNOR %>%
      dplyr::filter(!is.na(SEQUENCE)) %>%
      dplyr::filter(INTERACTION != '0')

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
      dplyr::distinct() %>%
      dplyr::mutate_at('INTERACTION', as.character)

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


    psp_db_phospho_int <- integrate_psp_without_kin(organism = organism,
                              reg_site_path = psp_reg_site_path,
                              psp_db = psp_db,
                              only_activatory = only_activatory,
                              local = local)

    psp_db_phosphoscore <- dplyr::bind_rows(psp_db_phospho_int,
                                            psp_db_phosphoscore)

    # Bring to upper complexes
    psp_db_phosphoscore <- psp_db_phosphoscore %>%
      dplyr::mutate(ENTITY = ifelse(grepl("[^[:alnum:]]", ENTITY),
                                   stringr::str_to_upper(ENTITY),
                                   ENTITY))

    if(organism == 'mouse'){
      psp_db_phosphoscore <- psp_db_phosphoscore %>%
        dplyr::mutate(ENTITY = ifelse(grepl("[^[:alnum:]]", ENTITY),
                                      ENTITY,
                                      stringr::str_to_title(ENTITY)))

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
        dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM) %>%
        dplyr::distinct() %>%
        dplyr::mutate_at('INTERACTION', as.character)
    }else{
      phosphosite_mech <- psp_db_phosphoscore %>%
        dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITY, '-', SEQUENCE),
                      IDB,
                      INTERACTION,
                      MECHANISM = 'phosphorylation') %>%
        dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM) %>%
        dplyr::distinct() %>%
        dplyr::mutate_at('INTERACTION', as.character)
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
