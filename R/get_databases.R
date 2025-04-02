#' Parse SIGNOR Network Data
#'
#' Retrieves and processes interaction data from the SIGNOR database for the specified organism.
#'
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param direct Logical, whether to retain only direct interactions among biological entities; default `TRUE`.
#' @param file_path Character. Path to save the SIGNOR processed files. If `NULL`, results are not saved. Default: `NULL`.
#' @param only_activatory Logical. If `TRUE`, only interactions regulating protein activity are included. Default: `FALSE`.
#'
#' @return A list containing:
#'   - `proteins`: A data frame of SIGNOR entities (protein IDs and names).
#'   - `interactions`: A data frame of parsed SIGNOR interactions with regulatory effects.
#'
#' @details
#' `only_activatory==TRUE`is used to create information for the PhosphoScore algorithm
#'
#' @export
#'
#' @examples
#' # Parse human SIGNOR interactions, including indirect interactions
#' signor_data <- signor_parsing("human")
#'
#' # Parse mouse SIGNOR interactions, including only direct interactions
#' signor_data_direct <- signor_parsing("mouse", direct = TRUE)
#'
#' # Save parsed results to a specific file path
#' signor_data_saved <- signor_parsing("human", file_path = "./", only_activatory = TRUE)
signor_parsing <- function(organism, direct = FALSE, file_path = NULL, only_activatory = FALSE) {

  # Validate organism input
  if(organism == 'human'){
    taxID <- 9606
  } else if (organism == 'mouse') {
    taxID <- 10090
  } else {
    stop ('Please provide a valid organism between human and mouse')
  }

  # Load SIGNOR data
  url <-
    paste0("https://signor.uniroma2.it/getData.php?organism=", taxID)
  query <- readr::read_tsv(
    url,
    col_names = c(
      "ENTITYA",
      "TYPEA",
      "IDA",
      "DATABASEA",
      "ENTITYB",
      "TYPEB",
      "IDB",
      "DATABASEB",
      "EFFECT",
      "MECHANISM",
      "RESIDUE",
      "SEQUENCE",
      "TAX_ID",
      "CELL_DATA",
      "TISSUE_DATA",
      "MODULATOR_COMPLEX",
      "TARGET_COMPLEX",
      "MODIFICATIONA",
      "MODASEQ",
      "MODIFICATIONB",
      "MODBSEQ",
      "PMID",
      "DIRECT",
      "NOTES",
      "ANNOTATOR",
      "SENTENCE",
      "SIGNOR_ID",
      "SCORE"
    ),
    show_col_types = FALSE
  )

  ## If the used selected DIRECT interactions:
  if (direct) {

    ## keep only phenotype-related indirect interactions
    query_sub_pheno <- query %>%
      dplyr::filter(DIRECT == FALSE & (TYPEB == 'phenotype' | TYPEA == 'stimulus')) %>%
      dplyr::select(c('ENTITYA','TYPEA','IDA','DATABASEA', 'EFFECT',
                      'ENTITYB','TYPEB','IDB','DATABASEB','MECHANISM',
                      'RESIDUE', 'SEQUENCE', 'SCORE')) %>%
      dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE),
                    ENTITYA = stringr::str_to_upper(ENTITYA),
                    ENTITYB = stringr::str_to_upper(ENTITYB))

    query_sub_pheno$INTERACTION <- 0
    query_sub_pheno$INTERACTION[grepl('down*', query_sub_pheno$EFFECT)] <- -1
    query_sub_pheno$INTERACTION[grepl('up*', query_sub_pheno$EFFECT)] <- 1
    query_sub_pheno$EFFECT <- NULL

    # Filter the direct interactions
    query <- query %>% dplyr::filter(DIRECT == TRUE)
  }

  # Remove interactions involving protein families
  query_sub_prot <- query %>%
    dplyr::select(c('ENTITYA','TYPEA','IDA','DATABASEA', 'EFFECT',
                    'ENTITYB','TYPEB','IDB','DATABASEB','MECHANISM',
                    'RESIDUE', 'SEQUENCE', 'SCORE', 'DIRECT')) %>%
    dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE),
                  ENTITYA = stringr::str_to_upper(ENTITYA),
                  ENTITYB = stringr::str_to_upper(ENTITYB)) %>%
    dplyr::filter(TYPEA != 'proteinfamily' & TYPEB != 'proteinfamily')

  # Remove residues annotations not regarding (de)phosphorylation events
  query_sub_prot[!grepl('*phos*', query_sub_prot$MECHANISM) &
                   (grepl('down*', query_sub_prot$EFFECT) | grepl('up*', query_sub_prot$EFFECT)),c('RESIDUE', 'SEQUENCE')] <-NA

  query_sub_prot$INTERACTION <- 0

  # Prepare SIGNOR entities table (PKN nodes)
  entities <- dplyr::bind_rows(query_sub_prot %>%
                                 dplyr::select(IDA, ENTITYA),
                               query_sub_prot %>%
                                 dplyr::select(IDA = IDB, ENTITYA = ENTITYB)) %>%
    dplyr::distinct()

  # Save entities file if path is provided
  if(!is.null(file_path)){
    readr::write_tsv(entities, file.path(file_path, paste0('SIGNOR_entities_', Sys.Date(), '.tsv')))
  }

  # Adjust interactions classification
  if(only_activatory){
    query_sub_prot$INTERACTION[grepl('*up-regulates activity*', query_sub_prot$EFFECT)] <- 1
    query_sub_prot$INTERACTION[grepl('*down-regulates activity*', query_sub_prot$EFFECT)] <- -1
  }else{
    # Match to 1, -1 notation
    query_sub_prot$INTERACTION[grepl('*up*', query_sub_prot$EFFECT)] <- 1
    query_sub_prot$INTERACTION[grepl('*down*', query_sub_prot$EFFECT)] <- -1
    query_sub_prot$INTERACTION[grepl('*form*', query_sub_prot$EFFECT)] <- 1
  }

  query_sub_prot$EFFECT <- NULL

  # Merge indirect phenotype interactions if undirect interactions were removed
  if(direct & !only_activatory){
    query_sub <- dplyr::bind_rows(query_sub_prot, query_sub_pheno)
  }else{
    query_sub <- query_sub_prot
  }

  # Save interactions file if path is provided
  if(!is.null(file_path)){
    message(paste0('SIGNOR PKN created at ', file_path))
    readr::write_tsv(query_sub, file.path(file_path, paste0('SIGNOR_PKN_', Sys.Date(), '.tsv')))
  }

  return(list(proteins = entities, interactions = query_sub))
}

#' Integrate Ser/Thr Kinome Atlas and PhosphoSitePlus Data into PKN
#'
#' This function integrate kinase-substrate interactions from PhosphoSitePlus
#' and inferred in the Serine/Threonine and Tyrosine Kinome Atlases (PMIDs: 36631611, 38720073) into
#' a Prior Knowledge Network (PKN).
#'
#' @param reg_site_path Character. Path to the PhosphoSitePlus "Regulatory_sites" file manually downloaded from it.
#' @param kin_sub_path Character. Path to the PhosphoSitePlus "Kinases_Substrate_Dataset" file manually downloaded from it.
#' @param file_path Character. Optional. Path to save the output file. If `NULL`, results are not saved. Default: `NULL`.
#' @param local Logical. If `TRUE`, uses local package paths; otherwise, uses system-installed paths. Default: `FALSE`.
#'
#' @return A `data.frame` containing integrated kinase-substrate interactions with
#'   regulatory effects mapped to SIGNOR conventions.
#'
#' @details
#' The function reads and processes:
#' - The **Kinome Atlas** dataset (filtered for 99% confidence that a kinase phosphorylate the substrate)
#' - The **PhosphoSitePlus Kinases_Substrates** file (manually downloaded)
#' - The **PhosphoSitePlus Regulatory_sites** file (manually downloaded)
#'
#' The output includes interactions where phosphorylation sites have known functional effects
#' (e.g., activation or inhibition of substrates), mapped to the **SIGNOR** regulatory format.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with local PhosphoSitePlus files
#' result <- integrate_atlas_in_pkn(
#'   reg_site_path = "./Regulatory_sites",
#'   kin_sub_path = "./Kinases_Substrate_Dataset",
#'
#' # Save results to a specific file path
#' integrate_atlas_in_pkn(
#'   reg_site_path = "./Regulatory_sites",
#'   kin_sub_path = "./Kinases_Substrate_Dataset",
#'   file_path = "output/")
#' }
integrate_atlas_in_pkn <- function(reg_site_path,
                                   kin_sub_path,
                                   file_path = NULL,
                                   local = FALSE){

  # Define path for package resources
  path_package <- ifelse(local == TRUE,
                         './inst/',
                         paste0(.libPaths()[1], '/SignalingProfiler/'))

  # Load Kinome Atlas dataset
  kinome_atlas_regulons <- access_remote_file(file = 'SerThrAtlas_PsP_regulons_filtered_99.tsv', dir = 'PKN')
                                          
  kinome_atlas_regulons$target <- NULL
  colnames(kinome_atlas_regulons) <- c('KIN_ACC_ID', 'KINASE', 'kin_percentile', 'SUB_ACC_ID',
                                       'SUB_GENE', 'sub_promiscuity_index',
                                       'SUB_MOD_RSD', 'SITE_...7_AA')

  # Load PhosphoSitePlus Kinase-Substrate dataset (filter human interactions)
  kinase_subs <- read.delim2(kin_sub_path, sep ='\t', skip = 2)
  kinase_org <- kinase_subs %>%
    dplyr::filter(KIN_ORGANISM == 'human' & SUB_ORGANISM == 'human')

  # Load SIGNOR dictionary mapping phosphosites to regulatory functions
  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))

  # Load PhosphoSitePlus Regulatory Sites dataset (filter human phosphosites)
  regulatory_sites <- read.delim2(reg_site_path, sep ='\t', skip = 2)
  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == 'human' & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  # Merge Kinome Atlas kinase-substrate interactions with regulatory data
  kin2res <- dplyr::inner_join(kinome_atlas_regulons,
                               regulatory_org,
                               by = c('SUB_ACC_ID' = 'ACC_ID',
                                      'SUB_MOD_RSD' = 'MOD_RSD')) %>%
    dplyr::select(c('KIN_ACC_ID',
                    'KINASE',
                    'SUB_ACC_ID', 'SUB_GENE',
                    'SUB_MOD_RSD',
                    'SITE_...7_AA.x',
                    'ON_FUNCTION')) %>%
    dplyr::filter(!(is.na(ON_FUNCTION) | ON_FUNCTION == ''))

  # Standardize phosphosite notation
  kin2res <- kin2res %>%
    dplyr::filter(!is.na(SUB_MOD_RSD)) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'S', 'Ser')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'T', 'Thr')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'Y', 'Tyr'))

  # Separate regulatory functions and map to SIGNOR interactions
  kin2res <- kin2res %>%
    tidyr::separate_rows('ON_FUNCTION', sep = '; ') %>%
    dplyr::left_join(phospho2SIGNOR,
                     by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  # Assign SIGNOR-like database fields
  kin2res <- kin2res %>%
    dplyr::mutate(TYPEA = 'protein',
                  TYPEB = 'protein',
                  DATABASEA = 'UNIPROT',
                  DATABASEB = 'UNIPROT',
                  MECHANISM = 'phosphorylation',
                  GROUP = '')

  kin2res_mapped_SIGNOR <- kin2res %>%
    dplyr::select(ENTITYA = KINASE,
                  TYPEA,
                  IDA = KIN_ACC_ID,
                  DATABASEA,
                  INTERACTION = sign,
                  ENTITYB = SUB_GENE,
                  TYPEB,
                  IDB = SUB_ACC_ID,
                  DATABASEB,
                  MECHANISM,
                  GROUP,
                  RESIDUE = SUB_MOD_RSD,
                  SEQUENCE = SITE_...7_AA.x
    ) %>%
    dplyr::mutate(ENTITYA = stringr::str_to_upper(ENTITYA),
                  SEQUENCE = stringr::str_to_upper(SEQUENCE),
                  ENTITYB = stringr::str_to_upper(ENTITYB)) %>%
    dplyr::distinct()

  # Save output if file_path is provided
  if(!is.null(file_path)){
    readr::write_tsv(kin2res_mapped_SIGNOR,
                     file.path(file_path, paste0('atlas_derived_edges_', Sys.Date(), '.tsv'))
    )
  }

  return(kin2res_mapped_SIGNOR)
}

#' Parse PhosphoSitePlus Data and (optionally) Integrate it with Kinome Atlas
#'
#' This function processes phosphorylation site interactions from PhosphoSitePlus
#' and optionally integrates the Serine/Threonine and Tyrosine Kinome
#' Atlases (PMIDs: 36631611, 38720073) data into the Prior Knowledge Network (PKN).
#'
#' @param reg_site_path Character. Path to the PhosphoSitePlus "Regulatory_sites" file manually downloaded from it.
#' @param kin_sub_path Character. Path to the PhosphoSitePlus "Kinases_Substrate_Dataset" file manually downloaded from it.
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param with_atlas Logical. If `TRUE`, integrates Ser/Thr/Tyr Kinome Atlas interactions (human only). Default: `FALSE`.
#' @param file_path Character. Optional. Path to save the output file. If `NULL`, results are not saved. Default: `NULL`.
#' @param only_activatory Logical. Indicates whether to filter for only phosphosites that regulate protein activity (not abundance). Default is `TRUE`.
#' @param local Logical. If `TRUE`, uses local package paths; otherwise, uses system-installed paths. Default: `FALSE`.
#'
#' @return A `data.frame` containing phosphorylation-based interactions, formatted for integration into a Prior Knowledge Network (PKN).
#'
#' @details
#' The function reads and processes:
#' - The **PhosphoSitePlus Kinases_Substrates** file (manually downloaded)
#' - The **PhosphoSitePlus Regulatory_sites** file (manually downloaded)
#' - If `with_atlas = TRUE`, integrates **Ser/Thr and Tyr Kinome Atlas interactions** (human only)
#' - If `only_activatory = TRUE`, only interactions regulating protein activity are retained. Used in PhosphoScore information creation.
#'
#' It filters and maps phosphosites to known regulatory effects using a SIGNOR-based dictionary.
#' The resulting table includes phosphorylation events that regulate protein activity.
#'
#' @seealso [integrate_atlas_in_pkn]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse PhosphoSitePlus data for human without Kinome Atlas integration
#' psp_data <- psp_parsing(
#'   reg_site_path = "./Regulatory_sites",
#'   kin_sub_path = "./Kinases_Substrate_Dataset",
#'   organism = "human")
#'
#' # Parse and integrate Kinome Atlas interactions (for human only)
#' psp_atlas_data <- psp_parsing(
#'   reg_site_path = "./Regulatory_sites",
#'   kin_sub_path = "./Kinases_Substrate_Dataset",
#'   organism = "human",
#'   with_atlas = TRUE)
#'
#' # Save parsed results to a specific file path
#' psp_parsing(
#'   reg_site_path = "./Regulatory_sites",
#'   kin_sub_path = "./Kinases_Substrate_Dataset",
#'   organism = "human",
#'   file_path = "output/")
#' }
#'
psp_parsing <- function(reg_site_path,
                        kin_sub_path,
                        organism,
                        only_activatory = FALSE,
                        with_atlas = FALSE,
                        file_path = NULL,
                        local = FALSE){

  # Load PhosphoSitePlus regulatory and kinase-substrate datasets
  regulatory_sites <- tibble::tibble(read.delim2(reg_site_path, sep ='\t', skip = 2))
  kinase_subs <- tibble::tibble(read.delim2(kin_sub_path, sep ='\t', skip = 2))

  # Define path for package resources
  path_package <- ifelse(local == TRUE, './inst/',
                         paste0(.libPaths()[1], '/SignalingProfiler/'))

  # Load phosphosite-to-SIGNOR mapping dictionary
  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))

  # Filter datasets for the selected organism and phosphosites
  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == organism & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  kinase_org <- kinase_subs %>%
    dplyr::filter(KIN_ORGANISM == organism & SUB_ORGANISM == organism)

  # Join kinase-substrate interactions with regulatory information
  kin2res <- dplyr::left_join(kinase_org,
                              regulatory_org,
                              by = c('SUB_ACC_ID' = 'ACC_ID',
                                     'SUB_MOD_RSD' = 'MOD_RSD'))

  # Select relevant columns
  kin2res <- kin2res %>% dplyr::select(c('GENE.x', 'KINASE', 'KIN_ACC_ID','SUBSTRATE',
                                         'SUB_GENE_ID', 'SUB_ACC_ID', 'SUB_GENE',
                                         'SUB_MOD_RSD', 'SITE_GRP_ID.x', 'SITE_...7_AA.x',
                                         'DOMAIN.x', 'IN_VIVO_RXN', 'IN_VITRO_RXN',
                                         'PROT_TYPE', 'HU_CHR_LOC', 'ON_FUNCTION'))

  # Standardize phosphosite notation (S, T, Y -> Ser, Thr, Tyr)
  kin2res <- kin2res %>% dplyr::filter(!is.na(SUB_MOD_RSD)) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'S', 'Ser')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'T', 'Thr')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'Y', 'Tyr'))

  # Separate regulatory functions and map them to SIGNOR interactions
  kin2res <- kin2res %>% tidyr::separate_rows('ON_FUNCTION', sep = '; ')
  kin2res_mapped <- dplyr::left_join(kin2res, phospho2SIGNOR, by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  if(only_activatory){
    kin2res_mapped <- kin2res_mapped %>% dplyr::filter(grepl('activity', kin2res_mapped$ON_FUNCTION))
  }
  # Assign SIGNOR-like database fields
  kin2res_mapped_SIGNOR <- kin2res_mapped %>%
    dplyr::mutate(TYPEA = 'protein',
                  TYPEB = 'protein',
                  DATABASEA = 'UNIPROT',
                  DATABASEB = 'UNIPROT',
                  MECHANISM = 'phosphorylation',
                  GROUP = '') %>%
    dplyr::select(ENTITYA = GENE.x,
                  TYPEA,
                  IDA = KIN_ACC_ID,
                  DATABASEA,
                  INTERACTION = sign,
                  ENTITYB = SUB_GENE,
                  TYPEB,
                  IDB = SUB_ACC_ID,
                  DATABASEB,
                  MECHANISM,
                  GROUP,
                  RESIDUE = SUB_MOD_RSD,
                  SEQUENCE = SITE_...7_AA.x) %>%
    dplyr::mutate(ENTITYA = stringr::str_to_upper(ENTITYA),
                  SEQUENCE = stringr::str_to_upper(SEQUENCE),
                  ENTITYB = stringr::str_to_upper(ENTITYB))

  # Optionally integrate Ser/Thr and Tyr Kinome Atlas interactions (human only)
  if(organism == 'human' & with_atlas){
    kinome_atlas_PKN <- integrate_atlas_in_pkn(reg_site_path, kin_sub_path)
    kin2res_mapped_SIGNOR <- dplyr::bind_rows(kin2res_mapped_SIGNOR, kinome_atlas_PKN) %>%
      dplyr::distinct()
  }

  # Save output if file_path is provided
  if(!is.null(file_path)){
    readr::write_tsv(kin2res_mapped_SIGNOR,
                     file.path(file_path, paste0('PhosphoSitePlus_PKN_', Sys.Date(), '.tsv'))
    )
  }

  return(kin2res_mapped_SIGNOR)
}

#' Parse OmniPath Interaction Data and Map to SIGNOR Format
#'
#' This function retrieves, processes, and maps protein-protein interactions from OmniPath
#' to the SIGNOR format, handling 'protein complexes' misaligned notation and integrating
#' them into a Prior Knowledge Network (PKN).
#'
#' @param resources Character vector. List of OmniPath resources to retrieve interactions from.
#' @param file_path Character. Optional. Path to save the output file. If `NULL`, results are not saved. Default: `NULL`.
#'
#' @return A `data.frame` containing OmniPath-derived protein interactions, mapped to the SIGNOR format.
#'
#' @details
#' The function:
#' - **Retrieves interactions** from OmniPath using `OmnipathR::omnipath_interactions()`.
#' - **Filters and formats interactions** into a SIGNOR-like structure.
#' - **Handles molecular complexes**, mapping them using SIGNOR internal data.
#' - **Resolves complex interactions** into individual protein-protein relationships.
#'
#' If `file_path` is provided, the resulting table is saved as a `.tsv` file.
#'
#' @seealso [OmnipathR::omnipath_interactions]
#'
#' @export
#'
#' @examples
#' # Retrieve and process OmniPath interactions
#' omni_data <- omnipath_parsing(resources = c("SIGNOR", "SignaLink3"))
#' 
omnipath_parsing <- function(resources, file_path = NULL){

  # Retrieve OmniPath interactions
  interactions <- OmnipathR::omnipath_interactions(resources=resources)

  # Map interactions to SIGNOR format
  omni_interactions <- interactions %>%
    dplyr::mutate(INTERACTION = ifelse(is_stimulation == 1 & is_inhibition == 0, 1,
                                       ifelse(is_stimulation == 0 & is_inhibition == 1, -1, NA))) %>%
    dplyr::filter(!is.na(INTERACTION))

  omni_interactions_clean <- omni_interactions %>%  dplyr::select(IDA = source,
                                                                  IDB = target,
                                                                  ENTITYA = source_genesymbol,
                                                                  ENTITYB = target_genesymbol,
                                                                  INTERACTION)


  # Map OmniPath interactions to SIGNOR internal database
  omni2signor <- dplyr::left_join(omni_interactions_clean %>% 
                                    dplyr::mutate_at('INTERACTION', as.character), 
                                  get(data('PKN_human_atlas_ind')),
                                  by = c('IDA', 'IDB', 'ENTITYA', 'ENTITYB', 'INTERACTION'))

  # Identify and process interactions involving Protein Complexes
  complex_interactions <- omni_interactions %>%
    dplyr::filter(grepl('COMPLEX', source ) | grepl('COMPLEX', target ))

  complex_table <- tibble::tibble(
    ID_omni = c(complex_interactions$source[grepl('COMPLEX', complex_interactions$source)],
                complex_interactions$target[grepl('COMPLEX', complex_interactions$target)]),
    gn_omni = c(complex_interactions$source_genesymbol[grepl('COMPLEX', complex_interactions$source)],
                complex_interactions$target_genesymbol[grepl('COMPLEX', complex_interactions$target)])) %>%
    dplyr::distinct() %>%
    dplyr::mutate(ID_omni1 = stringr::str_remove_all(ID_omni, 'COMPLEX:')) #complex_table1

  # Retrieve SIGNOR complexes for remapping
  query_signor <- readr::read_tsv("https://signor.uniroma2.it/getDataInternal.php?complexes=all%E2%80%99",
                                  show_col_types = FALSE)
  query_signor$...4 <- NULL
  query_signor <- query_signor %>%
    dplyr::mutate(COMPONENTS1 = stringr::str_replace_all(COMPONENTS, ',', '_'),
                  COMPONENTS2 = unlist(split_and_sort(COMPONENTS1)))

  complex_table1 <- complex_table %>%
    dplyr::mutate(ID_omni2 = unlist(split_and_sort(ID_omni1)))
  
  omni_to_signor <- dplyr::inner_join(complex_table1,
                                      query_signor,
                                      by = c('ID_omni2' = 'COMPONENTS2')) %>%
    dplyr::select(ID_omni, gn_omni, SIG_ID, COMPLEX_NAME)

  data('PKN_proteins_human')
  dictionary_sp_omni_signor_dic <- dplyr::inner_join(PKN_proteins_human,
                                                     omni_to_signor,
                                                     by = c('ID' = 'SIG_ID'))

  # Resolve complex interactions into individual protein-protein interactions
  complex_interactions_simple <- complex_interactions %>%
    dplyr::select(source, target, source_genesymbol, target_genesymbol, INTERACTION)

  entitya_fixed <- dplyr::left_join(complex_interactions_simple,
                                    dictionary_sp_omni_signor_dic,
                                    by = c('source' = 'ID_omni')) %>%
    dplyr::mutate(source = ifelse(is.na(ENTITY), source, NA),
                  source_genesymbol = ifelse(is.na(ENTITY),
                                             source_genesymbol, NA)) %>%
    dplyr::mutate(IDA = dplyr::coalesce(source, ID),
                  ENTITYA = dplyr::coalesce(source_genesymbol, ENTITY)) %>%
    dplyr::select(target,
                  target_genesymbol,
                  INTERACTION,
                  IDA,
                  ENTITYA)


  entityab_fixed <- dplyr::left_join(entitya_fixed,
                                     dictionary_sp_omni_signor_dic,
                                     by = c('target' = 'ID_omni'),
                                     relationship = 'many-to-many') %>%
    dplyr::mutate(target = ifelse(is.na(ENTITY), target, NA),
                  target_genesymbol = ifelse(is.na(ENTITY), target_genesymbol, NA)) %>%
    dplyr::mutate(IDB = dplyr::coalesce(target, ID),
                  ENTITYB = dplyr::coalesce(target_genesymbol, ENTITY)) %>%
    dplyr::select(IDA, ENTITYA, INTERACTION, IDB, ENTITYB) %>%
    dplyr::distinct()

  # Remove OmniPath complex notation
  no_complex_notation <- entityab_fixed %>%
    dplyr::filter(!(grepl('COMPLEX', IDA ) | grepl('COMPLEX', IDB ))) %>%
    dplyr::mutate_at('INTERACTION', as.character)

  # Unite interactions of complexes and proteins
  no_complex <- omni2signor %>% dplyr::filter(!(grepl('COMPLEX', IDA ) | grepl('COMPLEX', IDB )))

  omni2signor_final <- dplyr::bind_rows(no_complex, no_complex_notation) %>%
    dplyr::mutate(TYPEA = ifelse(grepl('SIGNOR-C', IDA), 'complex', 'protein'),
                  TYPEB = ifelse(grepl('SIGNOR-C', IDB), 'complex', 'protein'),
                  DATABASEA = ifelse(grepl('SIGNOR-C', IDA), 'SIGNOR', 'UNIPROT'),
                  DATABASEB = ifelse(grepl('SIGNOR-C', IDA), 'SIGNOR', 'UNIPROT'))

  # Save output if file_path is provided
  if(!is.null(file_path)){
    readr::write_tsv(omni2signor_final, file.path(file_path, 'PKN_Omnipath_interactions.tsv'))
  }

  return(omni2signor_final)
}

