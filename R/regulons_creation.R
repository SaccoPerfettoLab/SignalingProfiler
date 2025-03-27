#' Retrieve SIGNOR Regulons
#'
#' This function extracts transcription factor regulons or kinase-substrate interactions from the SIGNOR database.
#' The function supports two types of analysis: Transcription Factor Enrichment Analysis (TFEA) and Kinase-Substrate Enrichment Analysis (KSEA).
#'
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param analysis A character string specifying the type of analysis. Can be `"tfea"` for transcriptional regulation or `"ksea"` for kinase-substrate interactions.
#'
#' @details
#' - If `analysis = "tfea"`, the function retrieves transcriptional regulation interactions, returning a data frame of transcription factors and their target genes.
#' - If `analysis = "ksea"`, the function retrieves kinase-substrate phosphorylation interactions, returning a data frame of kinases and their phosphorylation targets.
#'
#' @return A tibble containing the regulon data. The structure of the output depends on the analysis type:
#' - For `"tfea"`: Columns include `tf` (transcription factor), `target` (regulated gene), `mor` (mode of regulation), and `confidence` (set to `"A"`).
#' - For `"ksea"`: Columns include `tf` (kinase), `target` (substrate in `Gene-Residue-Position` format), and `mor` (mode of regulation).
#'
#' @export
#'
#' @examples
#' # Retrieve transcription factor regulons for human
#' tfea_regulons <- get_signor_regulons(organism = "human", analysis = "tfea")
#'
#' # Retrieve kinase-substrate interactions for mouse
#' ksea_regulons <- get_signor_regulons(organism = "mouse", analysis = "ksea")
#'
get_signor_regulons <- function(organism, analysis){

  if(analysis == 'tfea'){
    signor_db <- signor_parsing(organism)
    SIGNOR_tfs <- signor_db$interactions %>%
      dplyr::filter(MECHANISM == 'transcriptional regulation')
    SIGNOR_regulons <- SIGNOR_tfs %>%
      dplyr::select(tf = ENTITYA, target = ENTITYB, mor = INTERACTION) %>%
      dplyr::distinct() %>%
      dplyr::mutate(confidence = 'A')

  }else if(analysis == 'ksea'){
    signor_db <- signor_parsing(organism, direct = TRUE)
    SIGNOR_kin_phos <- signor_db$interactions %>%
      dplyr::filter(grepl('phos', MECHANISM))

    SIGNOR_ksea <- SIGNOR_kin_phos %>%
      dplyr::select(ENTITYA, IDB, RESIDUE, SEQUENCE, mor = INTERACTION) %>%
      dplyr::filter(!is.na(RESIDUE)) %>%
      dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'Thr', 'T'),
                    RESIDUE = stringr::str_replace(RESIDUE, 'Tyr', 'Y'),
                    RESIDUE = stringr::str_replace(RESIDUE, 'Ser', 'S'))

    separated_tibble <- tidyr::separate(SIGNOR_ksea, RESIDUE,
                                        into = c("AMINO", "POSITION"),
                                        sep = "(?<=[A-Za-z])(?=[0-9])",
                                        remove = FALSE)

    SIGNOR_regulons <- separated_tibble %>%
      dplyr::mutate(target = paste0(IDB, '-', AMINO, '-', POSITION)) %>%
      dplyr::select(tf = ENTITYA, target, mor)
  }else{
    stop('Please provide a valid analysis type between "tfea" and "ksea"')
  }

  if(organism == 'mouse'){
    # Put in start case only proteins (tf without strange characters)
    SIGNOR_regulons <- SIGNOR_regulons %>%
      dplyr::mutate(tf = ifelse(grepl("[^[:alnum:]]", SIGNOR_regulons$tf), tf, stringr::str_to_title(tf)))

    if(analysis == 'tfea'){
      SIGNOR_regulons <- SIGNOR_regulons %>%
        dplyr::mutate(target = ifelse(grepl("[^[:alnum:]]", SIGNOR_regulons$target), target, stringr::str_to_title(target)))
    }
  }

  # Correct complexes removing special characters
  SIGNOR_regulons <- SIGNOR_regulons %>%
    dplyr::mutate(tf = stringr::str_replace_all(tf, "[^[:alnum:]]", '_'))

  if(analysis == 'tfea'){
    SIGNOR_regulons <- SIGNOR_regulons %>%
      dplyr::mutate(target = stringr::str_replace_all(target, "[^[:alnum:]]", '_'))
  }

  return(SIGNOR_regulons)
}


#' Retrieve PhosphoSitePlus Regulons
#'
#' This function extracts kinase-substrate interactions from the PsP user-provided files:
#' Kinase_Substrate_Dataset and Regulatory_sites.
#'
#' @param reg_site_path Character. Path to the PhosphoSitePlus "Regulatory_sites" file.
#' @param kin_sub_path Character. Path to the PhosphoSitePlus  "Kinases_Substrate_Dataset" file.
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param with_atlas Logical. If `TRUE`, integrates Ser/Thr/Tyr Kinome Atlas interactions (human only). Default: `FALSE`.
#'
#' @return A tibble containing the kinase-substrate relationships data for `"ksea"`
#' with three columns:`tf` (kinase), `target` (substrate in `Gene-Residue-Position` format), and `mor` (mode of regulation).
#' `mor` is -1, 1, and is continous from 0 to 1 for Ser/Thr Kinome Atlas relationships.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' psp_db <- psp_parsing(reg_site_path = './Regulatory_sites',
#'                       kin_sub_path = './Kinase_Substrate_Dataset',
#'                       organism = 'human',
#'                       with_atlas = FALSE,
#'                       local = FALSE)
#' }
#'
get_psp_regulons <- function(reg_site_path,
                             kin_sub_path,
                             organism,
                             with_atlas = FALSE){

  psp_db <- psp_parsing(reg_site_path = reg_site_path,
                        kin_sub_path = kin_sub_path,
                        organism = organism,
                        with_atlas = with_atlas,
                        local = FALSE)

  psp_ksea <- psp_db %>%
    dplyr::select(ENTITYA, IDB, RESIDUE, SEQUENCE, mor = INTERACTION) %>%
    dplyr::filter(!is.na(RESIDUE)) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'Thr', 'T'),
                  RESIDUE = stringr::str_replace(RESIDUE, 'Tyr', 'Y'),
                  RESIDUE = stringr::str_replace(RESIDUE, 'Ser', 'S'))

  separated_tibble <- tidyr::separate(psp_ksea, RESIDUE,
                                      into = c("AMINO", "POSITION"),
                                      sep = "(?<=[A-Za-z])(?=[0-9])",
                                      remove = FALSE)

  psp_regulons <- separated_tibble %>%
    dplyr::mutate(target = paste0(IDB, '-', AMINO, '-', POSITION)) %>%
    dplyr::select(tf = ENTITYA, target, mor)

  # Correct complexes removing special characters
  psp_regulons <- psp_regulons %>%
    dplyr::mutate(tf = stringr::str_replace_all(tf, "[^[:alnum:]]", '_')) %>%
    dplyr::mutate(tf = ifelse(tf == 'BCR_ABL', 'BCR_ABL1', tf))

  if(organism == 'mouse'){
    psp_regulons <- psp_regulons %>%
      dplyr::mutate(tf = ifelse(grepl("[^[:alnum:]]", psp_regulons$tf), tf, stringr::str_to_title(tf)))
  }

  return(psp_regulons)
}

#' Retrieve Kinase-Substrate Interaction Data from OmniPath
#'
#' This function queries the OmniPath database for kinase/phosphatase-substrate interactions,
#' retrieving phosphorylation and dephosphorylation sites for a specified organism.
#'
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param resources A character vector specifying the OmniPath resources to query.
#'                  Default sources include `"PhosphoSite"`, `"SIGNOR"`, `"Reactome_ProtMapper"`, `"HPRD"`, `"DEPOD"`, and others.
#'                  For a complete list, refer to the OmniPath documentation.
#'
#' @details
#' - If no `resources` are specified, the function retrieves interactions from all available OmniPath resources.
#' - The function combines phosphorylation (`phosphorylation`) and dephosphorylation (`dephosphorylation`) interactions.
#' - The mode of regulation (`mor`) is assigned as `1` for phosphorylation and `-1` for dephosphorylation.
#'
#' @return A tibble containing kinase/phosphatase interactions with the following columns:
#' - `tf`: Kinase or phosphatase.
#' - `target`: Substrate, formatted as `Gene-Residue-Position` (e.g., `"MAPK1-T-185"`).
#' - `mor`: Mode of regulation (`1` for phosphorylation, `-1` for dephosphorylation).
#'
#' @export
#'
#' @examples
#' # Retrieve OmniPath kinase-substrate interactions for human
#' omnipath_regulons <- get_omnipath(organism = "human")
#'
#' # Retrieve interactions for mouse using specific resources
#' omnipath_mouse <- get_omnipath(organism = "mouse", resources = c("SIGNOR", "DEPOD"))
#
get_omnipath <- function(organism,
                         resources = c('PhosphoSite', 'PhosphoSite_ProtMapper',
                                       'SIGNOR', 'SIGNOR_ProtMapper',
                                       'Reactome_ProtMapper', 'phosphoELM',
                                       'HPRD', 'DEPOD')){

  if(organism == 'human'){
    ID = 9606
  }else if (organism == 'mouse'){
    ID = 10090
  }else{
    stop('Please provide a valid organism name between "human" and "mouse"')
  }

  phos_df <- read.table(paste0('https://omnipathdb.org/ptms?genesymbols=1&types=phosphorylation&organisms=', ID, '&fields=sources,references'),header = TRUE,fill = TRUE)
  dephos_df <- read.table(paste0('https://omnipathdb.org/ptms?genesymbols=1&types=dephosphorylation&organisms=', ID, '&fields=sources,references'),header = TRUE, fill = TRUE)

  mod_df <- rbind(phos_df, dephos_df)

  if(is.null(resources)){
    message('Using all OmniPath resources')
    mod_df_filt <- mod_df
  }else{
    mod_df_filt <- mod_df %>%
      tidyr::separate_rows(sources, sep = ';') %>%
      dplyr::filter(sources %in% resources)
  }

  mod_df_filt <- mod_df_filt %>% dplyr::select(-c(sources, references)) %>% dplyr::distinct()
  df_regulons <- tibble::tibble(
    tf = mod_df_filt$enzyme_genesymbol,
    target = paste0(mod_df_filt$substrate,
                    '-',mod_df_filt$residue_type,
                    '-',mod_df_filt$residue_offset), # targets
    mor = mod_df_filt$modification)

  df_regulons$mor[df_regulons$mor=='phosphorylation'] <- 1
  df_regulons$mor[df_regulons$mor=='dephosphorylation'] <- -1

  df_regulons<- df_regulons %>%
    dplyr::mutate(mor = as.numeric(mor), tf = tf)

  return(df_regulons)
}

#' Generate Transcription Factor Enrichment Analysis (TFEA) Regulons
#'
#' This function retrieves transcription factor (TF) regulons
#' from user-specified resources, including SIGNOR, Dorothea, and Collectri.
#' The function supports both human and mouse organisms.
#'
#' @param resources Character vector specifying the regulatory databases to query.
#' Allowed values are `"SIGNOR"`, `"Dorothea"`, and `"Collectri"`.
#' Multiple resources can be provided.
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#'
#' @details
#' - If `"SIGNOR"` is included in `resources`, the function retrieves TF-target interactions from the SIGNOR database.
#' - If `"Dorothea"` is included (without `"Collectri"`), the function retrieves curated TF regulons from the Dorothea database.
#' - If `"Collectri"` is included, the function retrieves regulatory interactions from Collectri, replacing hyphens with underscores in TF and target names for consistency.
#' The resulting table merges all retrieved TF-target interactions from the selected resources.
#'
#' @return A tibble containing the compiled TF-target gene regulatory interactions. The output includes:
#' - `tf`: Transcription factor.
#' - `target`: Regulated gene.
#' - `mor`: Mode of regulation.
#' - `confidence`: Confidence level (applies to Collectri).
#'
#' @export
#'
#' @examples
#' # Retrieve TF regulons from SIGNOR and Dorothea for human
#' regulons <- create_tfea_regulons(resources = c("SIGNOR", "Dorothea"), organism = "human")
#'
#' # Retrieve TF regulons from Collectri for mouse
#' regulons_mouse <- create_tfea_regulons(resources = "Collectri", organism = "mouse")
#'
create_tfea_regulons <- function(resources, organism){

  regulons_list <- list()

  if('SIGNOR' %in% resources){
    message('Querying SIGNOR database')
    regulons_list[[length(regulons_list)+1]] <- get_signor_regulons(organism, analysis = 'tfea')
  }

  if('Dorothea' %in% resources & !('Collectri' %in% resources)){
    message('Querying Dorothea database')
    dorothea_regulons <- decoupleR::get_dorothea(organism = organism, levels = c('A'))
    dorothea_regulons <- dorothea_regulons %>% dplyr::rename('tf' = 'source')
    regulons_list[[length(regulons_list)+1]] <- dorothea_regulons
  } else if ( 'Collectri' %in% resources ){
    message('Querying Collectri database')
    collectri <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
    colnames(collectri) <- c('tf', 'target', 'mor')
    collectri <-  collectri %>%
      dplyr::mutate(tf = stringr::str_replace_all(tf, '-', '_'),
                    target = stringr::str_replace_all(target, '-', '_'),
                    confidence = 'A')
    regulons_list[[length(regulons_list)+1]] <- collectri
  }

  # Dynamically bind the rows of all the entities
  regulons <- do.call("bind_rows", regulons_list) %>% dplyr::distinct()

  return(regulons)
}

#' Generate Kinase-Substrate Enrichment Analysis (KSEA) Regulons
#'
#' This function retrieves kinase-substrate interaction data from the SIGNOR and OmniPath databases,
#' creating a unified dataset of regulatory interactions for Kinase-Substrate Enrichment Analysis (KSEA).
#'
#' @param resources A character vector specifying the regulatory databases to query.
#'                  Allowed values are `"SIGNOR"`, `"PsP"`, `"Omnipath"` and `"Atlas"` (to integrate Kinome-Atlas).
#'                  Default is `c("SIGNOR", "Omnipath")`.
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param reg_site_path Character. Path to the PhosphoSitePlus "Regulatory_sites" file.
#' @param kin_sub_path Character. Path to the PhosphoSitePlus  "Kinases_Substrate_Dataset" file.
#' @param omni_resources A character vector specifying the OmniPath sources to use.
#'                       Default sources include `"PhosphoSite"`, `"SIGNOR"`, `"HPRD"`, `"DEPOD"`, and others.
#'                       For a complete list, refer to the OmniPath documentation.
#' @param local Logical. If `TRUE`, uses local package paths; otherwise, uses system-installed paths. Default: `FALSE`.
#'
#' @details
#' - If `"SIGNOR"` is included in `resources`, the function retrieves kinase-substrate interactions from the SIGNOR database.
#' - If `"Omnipath"` is included, the function retrieves kinase-substrate interactions from the OmniPath database.
#' - The resulting table merges interactions from the selected resources into a single dataset.
#'
#' @return A tibble containing the compiled kinase-substrate regulatory interactions with the following columns:
#' - `tf`: Kinase or phosphatase.
#' - `target`: Substrate, formatted as `Gene-Residue-Position` (e.g., `"AKT1-S-473"`).
#' - `mor`: Mode of regulation (`1` for phosphorylation, `-1` for dephosphorylation).
#'
#' @examples
#' # Retrieve KSEA regulons from SIGNOR and OmniPath for human
#' ksea_regulons <- create_ksea_regulons(resources = c("SIGNOR", "Omnipath"), organism = "human")
#'
#' # Retrieve KSEA regulons from OmniPath for mouse using specific resources
#' ksea_mouse <- create_ksea_regulons(resources = "Omnipath", organism = "mouse", omni_resources = c("SIGNOR", "DEPOD"))
#'
#' @export
#'
create_ksea_regulons <- function(resources = c('SIGNOR', 'PsP', 'Omnipath'),
                                 organism,
                                 reg_site_path = NULL,
                                 kin_sub_path = NULL,
                                 omni_resources = c('PhosphoSite',
                                                    'PhosphoSite_ProtMapper',
                                                    'SIGNOR',
                                                    'SIGNOR_ProtMapper',
                                                    'Reactome_ProtMapper',
                                                    'phosphoELM',
                                                    'phosphoELM',
                                                    'HPRD',
                                                    'DEPOD'),
                                 local = FALSE){
  regulons_list <- list()

  if('SIGNOR' %in% resources){
    message('Querying SIGNOR database')
    regulons_list[[length(regulons_list)+1]] <- get_signor_regulons(organism, analysis = 'ksea')
  }

  if('Omnipath' %in% resources){
    message('Querying OmniPath database')
    regulons_list[[length(regulons_list)+1]] <- get_omnipath(organism, resources = omni_resources)
  }

  if('PsP' %in% resources){
    message('Reading PsP database files')
    regulons_list[[length(regulons_list)+1]] <- get_psp_regulons(organism,
                                                                 reg_site_path = reg_site_path,
                                                                 kin_sub_path = kin_sub_path,
                                                                 with_atlas = 'Atlas' %in% resources)
  }

  if('Atlas' %in% resources){
    if(organism == 'mouse'){
      warning('Atlas was not integrated because organism is mouse!')
    }else {
      message('Integrating Serine/Threonine and Tyrosine Kinome Atlas')

      # Define path for package resources, in this case take it from server
      path_file <- ifelse(local == TRUE,
                          './data-raw/yaffe_integration/',
                          paste0(.libPaths()[1], '/SignalingProfiler/'))

      # Load Kinome Atlas dataset
      regulons_list[[length(regulons_list)+1]] <- readr::read_tsv(paste0(path_file, 'yaffe_regulons_ksea_STY.tsv'),
                                                                  show_col_types = FALSE)
    }
  }

  # Dynamically bind the rows of all the entities
  regulons <- do.call("bind_rows", regulons_list) %>%
    dplyr::distinct(tf, target, mor, .keep_all = TRUE)

  return(regulons)
}

