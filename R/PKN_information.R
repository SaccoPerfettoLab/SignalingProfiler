
# ==================================================================== #
# Custom PKN and regulons generation
# ==================================================================== #


# ** Functions for the harmonization of entities IDs **
# updating UNIPROT IDs and Gene Names (Primary)
query_uniprot_proteins <- function(id_input, batch_size = 400){

  # Keep only ID input that are UNIPROT IDs
  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id_input)]

  header_df_uni2seq_fin <- c("Entry", "Reviewed", "Entry Name",
                             "Protein names", "Gene Names (primary)",
                             "Organism", "Length", "Sequence", "Date_modified")
  df_uni2seq_fin <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(df_uni2seq_fin) <- header_df_uni2seq_fin

  for (i in seq(from= 1, to= length(id_input)-(batch_size-1), by = batch_size)){
    print (i)
    id= unique(id_input[i:(i+(batch_size-1))])

    query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
    url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=',query_test)
    result <- httr::GET(url_test)
    as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
    df_uni2seq <- readr::read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE, show_col_types =  F)
    df_uni2seq_fin <- rbind(df_uni2seq_fin, df_uni2seq )
  }

  id= unique(id_input[i:length(id_input)])
  query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
  url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence%2Cdate_sequence_modified&format=tsv&query=',query_test)
  result <- httr::GET(url_test)
  as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
  df_uni2seq <- readr::read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE,show_col_types = F)
  df_uni2seq_fin <- rbind(df_uni2seq_fin,df_uni2seq ) %>% dplyr::distinct()

  return(df_uni2seq_fin)
}
query_uniprot <- function(id_input, batch_size = 400) {

  id_input <- unique(id_input)
  id_input <- id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]
  id_input <- id_input[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606|VAR|REV', id_input)]

  # Get isoforms
  id_input_iso <- unique(id_input[grepl('-', id_input)])

  if(length(id_input_iso) != 0){ # If there are isoforms do a separate query for them
    if(length(id_input_iso) > 400) {
      result_iso <- query_uniprot_proteins(id_input = id_input_iso, batch_size = batch_size)
    }else{
      result_iso <- query_uniprot_proteins(id_input = id_input_iso, batch_size = length(id_input_iso)-1)
    }
  }

  result_proteins <- query_uniprot_proteins(id_input, batch_size) # Query for not isoforms

  if(exists('result_iso')){
    result_total <- dplyr::bind_rows(result_iso, result_proteins)
  }else{
    result_total <- result_proteins
  }

  return(result_total)
}


# ======================================================================= #
# ** Functions for the Prior Knowledge Network Generation ** #
# ======================================================================= #

# =======>> using SIGNOR

#' signor_parsing
#'
#' @param organism string, 'human' or 'mouse'
#' @param direct Boolean, if TRUE only direct interactions are included
#' @param file_path string, path for file saving
#' @param only_activatory Boolean, if TRUE, only interactions regulating protein activity are included
#'
#' @return list of SIGNOR entities and their interactions
#' @export
#'
#' @examples
signor_parsing <- function(organism, direct = FALSE, file_path = NULL, only_activatory = FALSE){

  if(organism == 'human'){
    taxID = 9606
  } else if (organism == 'mouse') {
    taxID = 10090 }else{
      stop ('Please provide a valid organism between human and mouse')
    }


  query <- readr::read_tsv(paste0("https://signor.uniroma2.it/getData.php?organism=", taxID),
                           col_names = c("ENTITYA",
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
                                         "SCORE"))


  if(direct){
    # And create a tibble with indirect interactions only for phenotypes
    query_sub_pheno <- query %>%
      dplyr::filter(DIRECT == FALSE & (TYPEB == 'phenotype' | TYPEA == 'stimulus')) %>%
      dplyr::select(c('ENTITYA','TYPEA','IDA','DATABASEA', 'EFFECT',
                      'ENTITYB','TYPEB','IDB','DATABASEB','MECHANISM',
                      'RESIDUE', 'SEQUENCE', 'SCORE')) %>%
      dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE),
                    ENTITYA = stringr::str_to_uppertoupper(ENTITYA),
                    ENTITYB = stringr::str_to_uppertoupper(ENTITYB))

    query_sub_pheno$INTERACTION <- 0
    query_sub_pheno$INTERACTION[grepl('down*', query_sub_pheno$EFFECT)] <- -1
    query_sub_pheno$INTERACTION[grepl('up*', query_sub_pheno$EFFECT)] <- 1
    query_sub_pheno$EFFECT <- NULL

    # Then filter the query
    query <- query %>% dplyr::filter(DIRECT == TRUE)
  }

  query_sub_prot <- query %>%
    dplyr::select(c('ENTITYA','TYPEA','IDA','DATABASEA', 'EFFECT',
                    'ENTITYB','TYPEB','IDB','DATABASEB','MECHANISM',
                    'RESIDUE', 'SEQUENCE', 'SCORE', 'DIRECT')) %>%
    dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE),
                  ENTITYA = stringr::str_to_upper(ENTITYA),
                  ENTITYB = stringr::str_to_upper(ENTITYB)) %>%
    # removing protein families
    dplyr::filter(TYPEA != 'proteinfamily' & TYPEB != 'proteinfamily')

  # non annoto il residuo e la sequenza di eventi diversi dalla fosforilazione o defosforilazione
  query_sub_prot[!grepl('*phos*', query_sub_prot$MECHANISM) &
                   (grepl('down*', query_sub_prot$EFFECT) | grepl('up*', query_sub_prot$EFFECT)),c('RESIDUE', 'SEQUENCE')] <-NA

  query_sub_prot$INTERACTION <- 0

  entities <- dplyr::bind_rows(query_sub_prot %>%
                                 dplyr::select(IDA, ENTITYA),
                               query_sub_prot %>%
                                 dplyr::select(IDA = IDB, ENTITYA = ENTITYB)) %>%
    dplyr::distinct()

  if(!is.null(file_path)){
    readr::write_tsv(entities, paste0(file_path, 'SIGNOR_entities_', Sys.Date(), '.tsv'))
  }

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

  if(direct & !only_activatory){
    query_sub <- dplyr::bind_rows(query_sub_prot, query_sub_pheno)
  }else{
    query_sub <- query_sub_prot
  }

  if(!is.null(file_path)){
    message(paste0('SIGNOR PKN created at ', file_path))
    readr::write_tsv(query_sub, paste0(file_path, 'SIGNOR_PKN_', Sys.Date(), '.tsv'))
  }

  return(list(proteins = entities, interactions = query_sub))

}

# =======>> using Ser/Thr Kinome Atlas and PhosphoSitePlus

#' integrate_atlas_in_pkn
#'
#' @param reg_site_path path to the PhosphoSitePlus 'Regulatory Sites' file
#' @param kin_sub_path path to the PhosphoSitePlus 'Kinases-Substrates' file
#' @param file_path string, path for file saving
#'
#' @return table of interactions with Ser/Thr-derived interactions included
#' @export
#'
#' @examples
integrate_atlas_in_pkn <- function(reg_site_path,
                                   kin_sub_path, file_path = NULL){

  path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')

  kinome_atlas_regulons <- readr::read_tsv(paste0(path_package, 'extdata/SerThrAtlas_PsP_regulons_filtered_99.tsv'),
                                           show_col_types = FALSE)
  kinome_atlas_regulons$target <- NULL
  colnames(kinome_atlas_regulons) <- c('KIN_ACC_ID', 'KINASE', 'kin_percentile', 'SUB_ACC_ID',
                                       'SUB_GENE', 'sub_promiscuity_index',
                                       'SUB_MOD_RSD', 'SITE_...7_AA')

  # Read kin substrate dataset
  kinase_subs <- read.delim2(kin_sub_path, sep ='\t', skip = 2)
  kinase_org <- kinase_subs %>%
    dplyr::filter(KIN_ORGANISM == 'human' & SUB_ORGANISM == 'human')

  path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')

  # Read SIGNOR dictionary
  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))

  #Read PsP regulatory sites and filter human
  regulatory_sites <- read.delim2(reg_site_path, sep ='\t', skip = 2)
  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == 'human' & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  # Merge kinome_atlas KINASE to SUB_GENE with regulatory role of phosphosites
  kin2res <- dplyr::inner_join(kinome_atlas_regulons,
                               regulatory_org,
                               by = c('SUB_ACC_ID' = 'ACC_ID','SUB_MOD_RSD' = 'MOD_RSD'))


  # Create a smaller dataset removing percentiles and promiscuity index and
  # mergend rows without regulatory roles
  kin2res1 <- kin2res %>% dplyr::select(c('KIN_ACC_ID',
                                          'KINASE',
                                          'SUB_ACC_ID', 'SUB_GENE',
                                          'SUB_MOD_RSD',
                                          'SITE_...7_AA.x',
                                          'ON_FUNCTION')) %>%
    dplyr::filter(!(is.na(ON_FUNCTION) | ON_FUNCTION == ''))

  # Change S, T, Y notation
  kin2res1 <- kin2res1 %>%
    dplyr::filter(!is.na(SUB_MOD_RSD)) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'S', 'Ser')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'T', 'Thr')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'Y', 'Tyr'))

  # Split ON FUNCTION and select only activatiory or inhibitory relations
  # and transform according to SIGNOR_to_PsP dictionary
  kin2res1 <- kin2res1 %>%
    tidyr::separate_rows('ON_FUNCTION', sep = '; ')

  kin2res_mapped <- dplyr::left_join(kin2res1, phospho2SIGNOR,
                                     by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  # Assign the fields of SIGNOR database
  kin2res_mapped <- kin2res_mapped %>%
    dplyr::mutate(TYPEA = 'protein',
                  TYPEB = 'protein',
                  DATABASEA = 'UNIPROT',
                  DATABASEB = 'UNIPROT',
                  MECHANISM = 'phosphorylation',
                  GROUP = '')

  kin2res_mapped_SIGNOR <- kin2res_mapped %>%
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
                  SEQUENCE = SITE_...7_AA.x) %>%

    dplyr::mutate(ENTITYA = stringr::str_to_upper(ENTITYA),
                  SEQUENCE = stringr::str_to_upper(SEQUENCE),
                  ENTITYB = stringr::str_to_upper(ENTITYB)) %>%
    dplyr::distinct()

  if(!is.null(file_path)){

    readr::write_tsv(kin2res_mapped_SIGNOR,
                     paste0(file_path, 'atlas_derived_edges_', Sys.Date(), '.tsv'))
  }

  return(kin2res_mapped_SIGNOR)

}


#' psp_parsing
#'
#' @param reg_site_path path to the PhosphoSitePlus 'Regulatory Sites' file
#' @param kin_sub_path path to the PhosphoSitePlus 'Kinases-Substrates' file
#' @param organism string, 'human' or 'mouse'
#' @param with_atlas Boolean, if TRUE Ser/Thr Kinome Atlas interactions are included
#' @param file_path string, path for file saving
#'
#' @return table of interactions derived from PhosphoSitePlus
#' @export
#'
#' @examples
psp_parsing <- function(reg_site_path,
                        kin_sub_path,
                        organism,
                        with_atlas = FALSE,
                        file_path = NULL){

  regulatory_sites <- tibble::tibble(read.delim2(reg_site_path, sep ='\t', skip = 2))
  kinase_subs <- tibble::tibble(read.delim2(kin_sub_path, sep ='\t', skip = 2))

  path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')

  # dictionary for phosphosite-SIGNOR conversion internal in SignalingProfiler
  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))

  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == organism & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  kinase_org <- kinase_subs %>%
    dplyr::filter(KIN_ORGANISM == organism & SUB_ORGANISM == organism)

  # join of two tables
  kin2res <- dplyr::left_join(kinase_org,
                              regulatory_org,
                              by = c('SUB_ACC_ID' = 'ACC_ID',
                                     'SUB_MOD_RSD' = 'MOD_RSD'))

  # create a smaller df
  kin2res <- kin2res %>% dplyr::select(c('GENE.x', 'KINASE', 'KIN_ACC_ID','SUBSTRATE',
                                         'SUB_GENE_ID', 'SUB_ACC_ID', 'SUB_GENE',
                                         'SUB_MOD_RSD', 'SITE_GRP_ID.x', 'SITE_...7_AA.x',
                                         'DOMAIN.x', 'IN_VIVO_RXN', 'IN_VITRO_RXN',
                                         'PROT_TYPE', 'HU_CHR_LOC', 'ON_FUNCTION'))

  # change S, T, Y notation
  kin2res <- kin2res %>% dplyr::filter(!is.na(SUB_MOD_RSD)) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'S', 'Ser')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'T', 'Thr')) %>%
    dplyr::mutate(SUB_MOD_RSD = stringr::str_replace(SUB_MOD_RSD, 'Y', 'Tyr'))

  # split ON FUNCTION and select only activatiory or inhibitory relations
  kin2res <- kin2res %>% tidyr::separate_rows('ON_FUNCTION', sep = '; ')
  kin2res_mapped <- dplyr::left_join(kin2res, phospho2SIGNOR, by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  if(only_activatory){
    phosphoscore_table <- phosphoscore_table %>%
      dplyr::filter(grepl('activity', SIGNOR)) %>%
      dplyr::select(ENTITY, ID, RESIDUE, SEQUENCE, sign)
  }

  kin2res_mapped <- kin2res_mapped %>%
    dplyr::mutate(TYPEA = 'protein',
                  TYPEB = 'protein',
                  DATABASEA = 'UNIPROT',
                  DATABASEB = 'UNIPROT',
                  MECHANISM = 'phosphorylation',
                  GROUP = '')


  kin2res_mapped_SIGNOR <- kin2res_mapped %>%
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

  if(organism == 'human' & with_atlas){
    kinome_atlas_PKN <- integrate_atlas_in_pkn(reg_site_path, kin_sub_path)
    kin2res_mapped_SIGNOR <- dplyr::bind_rows(kin2res_mapped_SIGNOR, kinome_atlas_PKN) %>%
      dplyr::distinct()
  }

  if(!is.null(file_path)){
    readr::write_tsv(kin2res_mapped_SIGNOR, paste0(file_path, 'PhosphoSitePlus_PKN_', Sys.Date(), '.tsv'))
  }

  return(kin2res_mapped_SIGNOR)
}

# =======>> using OmniPath

#' omnipath_parsing
#'
#' @param resources vector of resources (see OmniPath tutorial for the complete list)
#' @param file_path string, path for file saving
#'
#' @return table of interactions derived from OmniPath
#' @export
#'
#' @examples
omnipath_parsing <- function(resources, file_path = NULL){

  split_and_sort <- function(x) {
    x <- strsplit(x, "_")
    lapply(x, function(y) paste(sort(unique(y)), collapse = "_"))
  }

  interactions <-
    OmnipathR::import_omnipath_interactions(resources=resources)

  # Edge format like SIGNOR
  omni_interactions <- interactions %>%
    dplyr::mutate(INTERACTION = ifelse(is_stimulation == 1 & is_inhibition == 0, 1,
                                       ifelse(is_stimulation == 0 & is_inhibition == 1, -1, NA))) %>%
    dplyr::filter(!is.na(INTERACTION))

  omni_interactions_clean <- omni_interactions %>%  dplyr::select(IDA = source,
                                                                  IDB = target,
                                                                  ENTITYA = source_genesymbol,
                                                                  ENTITYB = target_genesymbol,
                                                                  INTERACTION)


  # Map Omnipath interactions on SignalingProfiler internal database
  omni2signor <- dplyr::left_join(omni_interactions_clean, SignalingProfiler::PKN_human_atlas_ind,
                                  by = c('IDA', 'IDB', 'ENTITYA', 'ENTITYB', 'INTERACTION'))

  # Parsing of molecular entities like COMPLEXES
  complex_interactions <- omni_interactions %>% dplyr::filter(grepl('COMPLEX', source ) | grepl('COMPLEX', target ))

  tibble::tibble(ID_omni = c(complex_interactions$source[grepl('COMPLEX', complex_interactions$source)],
                             complex_interactions$target[grepl('COMPLEX', complex_interactions$target)]),
                 gn_omni = c(complex_interactions$source_genesymbol[grepl('COMPLEX', complex_interactions$source)],
                             complex_interactions$target_genesymbol[grepl('COMPLEX', complex_interactions$target)])) %>%
    dplyr::distinct() -> complex_table

  complex_table %>% dplyr::mutate(ID_omni1 = stringr::str_remove_all(ID_omni, 'COMPLEX:')) -> complex_table1

  # Get SIGNOR complexes for the remapping
  query_signor <- readr::read_tsv("https://signor.uniroma2.it/getDataInternal.php?complexes=all%E2%80%99", show_col_types = FALSE)
  query_signor$...4 <- NULL
  query_signor <- query_signor %>% dplyr::mutate(COMPONENTS1 = stringr::str_replace_all(COMPONENTS, ',', '_'))

  complex_table1 <- complex_table1 %>% dplyr::mutate(ID_omni2 = unlist(split_and_sort(ID_omni1)))
  query_signor <-  query_signor %>% dplyr::mutate(COMPONENTS2 = unlist(split_and_sort(COMPONENTS1)))

  omni_to_signor <- dplyr::inner_join(complex_table1, query_signor, by = c('ID_omni2' = 'COMPONENTS2'))
  omni_to_signor <- omni_to_signor %>% dplyr::select(ID_omni, gn_omni, SIG_ID, COMPLEX_NAME)

  dictionary_sp_omni_signor_dic <- dplyr::inner_join(PKN_proteins_human, omni_to_signor, by = c('ID' = 'SIG_ID'))

  complex_interactions_simple <- complex_interactions %>%
    dplyr::select(source, target, source_genesymbol, target_genesymbol, INTERACTION)

  dplyr::left_join(complex_interactions_simple,
                   dictionary_sp_omni_signor_dic,
                   by = c('source' = 'ID_omni')) %>%
    dplyr::mutate(source = ifelse(is.na(ENTITY), source, NA),
                  source_genesymbol = ifelse(is.na(ENTITY), source_genesymbol, NA)) %>%
    dplyr::mutate(IDA = dplyr::coalesce(source, ID),
                  ENTITYA = dplyr::coalesce(source_genesymbol, ENTITY)) -> entitya_fixed

  entitya_fixed$source <- NULL
  entitya_fixed$source_genesymbol <- NULL
  entitya_fixed <- entitya_fixed %>%
    dplyr::select(target, target_genesymbol,
                  INTERACTION, IDA, ENTITYA)#, TYPEA = TYPE, DATABASEA = DATABASE)

  dplyr::left_join(entitya_fixed,
                   dictionary_sp_omni_signor_dic,
                   by = c('target' = 'ID_omni'),
                   relationship = 'many-to-many') %>%
    dplyr::mutate(target = ifelse(is.na(ENTITY), target, NA),
                  target_genesymbol = ifelse(is.na(ENTITY), target_genesymbol, NA)) %>%
    dplyr::mutate(IDB = dplyr::coalesce(target, ID),
                  ENTITYB = dplyr::coalesce(target_genesymbol, ENTITY)) -> entityab_fixed

  complex_interactions_fixed <- entityab_fixed %>%
    dplyr::select(IDA, ENTITYA, INTERACTION, IDB, ENTITYB) %>% distinct()
  complex_interactions_fixed %>%
    dplyr::filter(!(grepl('COMPLEX', IDA ) | grepl('COMPLEX', IDB ))) -> no_complex_notation

  # Unite interactions of complexes and proteins
  omni2signor %>% dplyr::filter(!(grepl('COMPLEX', IDA ) | grepl('COMPLEX', IDB ))) -> no_complex
  dplyr::bind_rows(no_complex, no_complex_notation) -> omni2signor_final


  if(!is.null(file_path)){
    readr::write_tsv(omni2signor_final, paste0(file_path, 'PKN_Omnipath_interactions.tsv'))
  }

  return(omni2signor_final)
}

# =======>> function to integrate all the user-defined resources

#' create_PKN
#'
#' @param database vector, resources can be SIGNOR, OmniPath, PsP, SerThrAtlas
#' @param file_path string, path for file saving
#' @param organism string, 'human' or 'mouse'
#' @param direct Boolean, if TRUE only direct interactions are included
#' @param omnipath_resources vector, resources included in OmniPath database (see OmniPath tutorial for the complete list)
#' @param reg_site_path path to the PhosphoSitePlus 'Regulatory Sites' file
#' @param kin_sub_path path to the PhosphoSitePlus 'Kinases-Substrates' file
#'
#' @return table of interactions derived from the user-specified resources
#' @export
#'
#' @examples
create_PKN <- function(database = c('SIGNOR', 'OmniPath',
                                    'PsP', 'SerThr_Atlas'),
                       organism = 'human',
                       direct = FALSE,
                       file_path = NULL,
                       omnipath_resources = NULL,
                       reg_site_path = NULL,
                       kin_sub_path = NULL){

  # organism <- 'human'
  # direct <- TRUE
  # omnipath_resources <-  c("SignaLink3","PhosphoSite","SIGNOR")
  # database <- c('SIGNOR', 'OmniPath')
  # reg_site_path = './revisions/input/Regulatory_sites_2023-08-24'
  # kin_sub_path = './revisions/input/Kinase_Substrate_Dataset_2023-08-24'

  PKN_list <- list()

  if(('SIGNOR' %in% database | 'PsP' %in% database) & is.null(organism)){
    stop('Please specifify an organism between human and mouse')
  }

  if('SIGNOR' %in% database){
    SIGNOR <- signor_parsing(organism, direct)$interactions
    SIGNOR$SOURCE <- 'SIGNOR'

    PKN_list[[length(PKN_list) + 1]] <-  SIGNOR
  }


  if('PsP' %in% database){
    if((is.null(reg_site_path) | is.null(kin_sub_path))){
      stop('Please provide the paths to Kinase-Substrates and Regulatory Sites of
         PhosphoSitePlus database available at
         https://www.phosphosite.org/homeAction.action under registration')
    }else{
      if('SerThr_Atlas' %in% database){
        PsP <- psp_parsing(reg_site_path = reg_site_path,
                           kin_sub_path = kin_sub_path, organism = organism,
                           with_atlas = TRUE)
        PsP$SOURCE <- 'PsP+Atlas'
      }else{
        PsP <- psp_parsing(reg_site_path = reg_site_path,
                           kin_sub_path = kin_sub_path, organism = organism,
                           with_atlas = FALSE)
        PsP$SOURCE <- 'PsP'
        PsP$DIRECT <- TRUE
      }

      PKN_list[[length(PKN_list) + 1]] <- PsP
    }
  }else{
    if('SerThr_Atlas' %in% database){
      stop('You should provide PsP as a database to use Ser/Thr Atlas!')
    }
  }


  if('Omnipath' %in% database){
    if(!is.null(omnipath_resources)){
      stop('Plese provide a vector of resources to download from Omnipath')
    }else{
      omnipath <- omnipath_parsing(resources = omnipath_resources)
      omnipath$SOURCE <- 'OmniPath'
      omnipath$DIRECT <- FALSE
      PKN_list[[length(PKN_list) + 1]] <- omnipath
    }

  }

  # Dynamically bind the rows of all the entities
  PKN <- do.call("bind_rows", PKN_list)

  # Clean PKN_global
  PKN <- PKN %>% dplyr::filter(!(is.na(ENTITYA) | is.na(ENTITYB)))

  if(organism == 'mouse'){
    PKN <- PKN %>% dplyr::mutate(ENTITYA = ifelse(TYPEA == 'protein', stringr::str_to_title(ENTITYA), ENTITYA),
                                 ENTITYB = ifelse(TYPEB == 'protein', stringr::str_to_title(ENTITYB), ENTITYB))
  }

  # Call the function for the update of UNIPROT IDs
  message('** Query UniProt database for UNIPROT IDs and Gene Names (Primary) update ** ')
  res <- query_uniprot(id_input = unique(c(PKN$IDA, PKN$IDB)),
                       batch_size = 400)

  # Subset the result removing SEQUENCE column
  res_sub <- res %>%
    dplyr::select(ID = Entry, ENTITY = `Gene Names (primary)`) %>%
    dplyr::distinct()

  # ** Create a tibble with all the proteins of the PKN **
  all_proteins <- tibble::tibble( ID = c(PKN$IDA, PKN$IDB),
                                  ENTITY = c(PKN$ENTITYA, PKN$ENTITYB),
                                  TYPE = c(PKN$TYPEA, PKN$TYPEB),
                                  DATABASE = c(PKN$DATABASEA, PKN$DATABASEB)) %>%
    # if there are IDs associated to multiple gene names sample at random
    # it isn't a problem since the gene name will be replaced
    dplyr::distinct(ID, TYPE, DATABASE, .keep_all = TRUE)

  # if a gene name has multiple UNIPROT, TYPES or DATABASE are collapsed in one string
  all_proteins %>% dplyr::group_by(ENTITY) %>%
    dplyr::mutate(ID = paste0(unique(ID), collapse = ';'),
                  TYPE = paste0(unique(TYPE), collapse = ';'),
                  DATABASE = paste0(unique(DATABASE), collapse = ';')) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() -> all_proteins

  # since we want to use the UNIPROT for mapping, we reseparate the rows on the ID
  # obtaining a table which separates only UNIPROTs like this:
  #   P33800 PIPPO protein;complex UNIPROT;SIGNOR
  #   SIGNOR-3 PIPPO protein;complex UNIPROT;SIGNOR
  all_proteins <- all_proteins %>% tidyr::separate_rows(ID, sep = ';')

  # ** Analytes table creation (nodes_df) **

  # take from all proteins the one successfully found in UNIPROT
  analytes <- dplyr::inner_join(res_sub,
                                all_proteins %>% dplyr::select(-c('ENTITY')),
                                by = c('ID'))

  # take the proteins not found in UNIPROT because are chemicals, complex, fusion proteins
  not_uniprot <- dplyr::anti_join(all_proteins, res_sub, by = 'ID')

  # concatenate
  analytes <- dplyr::bind_rows(analytes, not_uniprot)

  # check for duplicates and remove Ubiquitin proteins
  # because they have the same UNIPROT ID
  # keeping only the entry with all gene names
  analytes %>% dplyr::count(ID, TYPE, DATABASE) %>% dplyr::filter(n>1) -> dup
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

  # ** Interaction table creation (edges_df) **

  # reseparate UNIPROT IDs
  analytes_for_edges <- analytes %>%
    tidyr::separate_rows(ID, sep = ';') %>%
    dplyr::distinct()

  # update ENTITYA, TYPEA, DATABASEA with Primary Gene Name and concatenation of types and dbs
  dplyr::inner_join(PKN %>%
                      dplyr::select(-c('ENTITYA', 'TYPEA', 'DATABASEA')),
                    analytes_for_edges,
                    by = c('IDA' = 'ID')) %>%
    dplyr::rename('ENTITYA' = 'ENTITY', 'TYPEA' = 'TYPE', 'DATABASEA' = 'DATABASE') -> ENTITYA_update

  # update ENTITYB, TYPEB, DATABASEB with Primary Gene Name and concatenation of types and dbs
  dplyr::inner_join(ENTITYA_update %>%
                      dplyr::select(-c('ENTITYB', 'TYPEB', 'DATABASEB')),
                    analytes_for_edges,
                    by = c('IDB' = 'ID')) %>%
    dplyr::rename('ENTITYB' = 'ENTITY', 'TYPEB' = 'TYPE', 'DATABASEB' = 'DATABASE') -> ENTITYB_update

  # Parse interactions file
  PKN_final <- ENTITYB_update %>% dplyr::relocate(ENTITYA, ENTITYB, INTERACTION)

  PKN_final$SCORE <- NULL
  PKN_final$GROUP <- NULL
  PKN_final$MECHANISM[is.na(PKN_final$MECHANISM)] <- ''
  PKN_final$INTERACTION <- PKN_final$INTERACTION %>% as.character #setting 1 and -1 as char

  # add key Phospho-Key-GeneName-Sequence
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

  # remove NA as gene names (immunoglobulins) also in the nodes table
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
                   TYPE = paste0(TYPE, collapse = ';'),
                   DATABASE = paste0(DATABASE, collapse = ';')) %>%
    dplyr::distinct()

  PKN_final1 <- PKN_final %>% dplyr::mutate(ENTITYA = stringr::str_replace_all(ENTITYA, "[^[:alnum:]]", '_'),
                                            ENTITYB = stringr::str_replace_all(ENTITYB, "[^[:alnum:]]", '_'))
  PKN_final1$ENTITYA[grepl('^\\d', PKN_final1$ENTITYA)] <- paste0('X', PKN_final1$ENTITYA[grepl('^\\d', PKN_final1$ENTITYA)])
  PKN_final1$ENTITYB[grepl('^\\d', PKN_final1$ENTITYB)] <- paste0('X', PKN_final1$ENTITYB[grepl('^\\d', PKN_final1$ENTITYB)])

  # ** Creath an igraph object **
  PKN_graph <- graph_from_data_frame(d = PKN_final1, vertices = analytes1)

  return(list(igraph_PKN = PKN_graph,
              interactions = PKN_final1,
              entities = analytes1))

}

# ======================================================================= #
# ** Functions for TFEA/KSEA regulons generation  ** #
# ======================================================================= #

#' get_signor_regulons
#'
#' @param organism string, 'human' and 'mouse'
#' @param analysis string, 'tfea' and 'ksea'
#'
#' @return table of regulons derived from SIGNOR
#' @export
#'
#' @examples
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

    SIGNOR_ksea <- SIGNOR_kin_phos %>% dplyr::select(ENTITYA, IDB, RESIDUE, SEQUENCE, mor = INTERACTION) %>%
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
    stop('Please provide a valid analysis type between tfea and ksea')
  }

  return(SIGNOR_regulons)
}

# =======>> TFEA regulons generation

#' create_tfea_regulons
#'
#' @param resources vector, allowed resources are SIGNOR, Dorothea or Collectri
#' @param organism string, 'human' and 'mouse'
#'
#' @return table of TF-target genes derived from user-defined resources
#' @export
#'
#' @examples
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
    collectri <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
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


# =======>> KSEA regulons generation

#' get_omnipath
#'
#' @param organism string, 'human' or 'mouse'
#' @param resources vector, see OmniPath tutorial for a complete list of kin-substrates resources
#'
#' @return a table of KIN/PHOS-target sites from user-defined OmniPath sources
#' @export
#'
#' @examples
get_omnipath <- function(organism,
                         resources = c('PhosphoSite', 'PhosphoSite_ProtMapper',
                                       'SIGNOR', 'SIGNOR_ProtMapper',
                                       'Reactome_ProtMapper', 'phosphoELM',
                                       'phosphoELM', 'HPRD', 'DEPOD')){

  if(organism == 'human'){
    ID = 9606
  }else if (organism == 'mouse'){
    ID = 10090
  }else{
    stop('Please provide a valide organism name between human and mouse')
  }

  phos_df <- read.table(paste0('https://omnipathdb.org/ptms?genesymbols=1&types=phosphorylation&organisms=', ID, '&fields=sources,references'),header = T,fill = T)
  dephos_df <- read.table(paste0('https://omnipathdb.org/ptms?genesymbols=1&types=dephosphorylation&organisms=', ID, '&fields=sources,references'),header = T, fill = T)

  mod_df <- rbind(phos_df, dephos_df)

  if(is.null(resources)){
    message('Using all OmniPath resources')
    mod_df_filt <- mod_df
  }else{
    mod_df_filt <- mod_df %>%
      tidyr::separate_rows(sources, sep = ';') %>%
      dplyr::filter(sources %in% accepted_sources$sources)
  }

  mod_df_filt <- mod_df_filt %>% dplyr::select(-c(sources, references)) %>% dplyr::distinct()
  df_regulons <- tibble::tibble(tf = mod_df_filt$enzyme_genesymbol, #kinase/phosphatase
                                target = paste0(mod_df_filt$substrate, #key: UNIPROT
                                                '-',mod_df_filt$residue_type,
                                                '-',mod_df_filt$residue_offset), #targets
                                mor = mod_df_filt$modification)

  df_regulons$mor[df_regulons$mor=='phosphorylation'] <- 1
  df_regulons$mor[df_regulons$mor=='dephosphorylation'] <- -1

  df_regulons<- df_regulons %>%
    dplyr::mutate(mor = as.numeric(mor), tf = tf)

  return(df_regulons)
}

#' create_ksea_regulons
#'
#' @param resources vector, sources such as SIGNOR or OmniPath
#' @param organism string, 'human' or 'mouse'
#' @param omni_resources vector, see OmniPath tutorial for a complete list of kin-substrates resources
#'
#' @return a table of KIN/PHOS-target sites from user-defined OmniPath sources
#' @export
#'
#' @examples
create_ksea_regulons <- function(resources = c('SIGNOR', 'Omnipath'), organism,
                                 omni_resources = c('PhosphoSite',
                                                    'PhosphoSite_ProtMapper',
                                                    'SIGNOR',
                                                    'SIGNOR_ProtMapper',
                                                    'Reactome_ProtMapper',
                                                    'phosphoELM',
                                                    'phosphoELM',
                                                    'HPRD',
                                                    'DEPOD')){


  regulons_list <- list()

  if('SIGNOR' %in% resources){
    message('Querying SIGNOR database')
    regulons_list[[length(regulons_list)+1]] <- get_signor_regulons(organism, analysis = 'ksea')

  }
  if('Omnipath' %in% resources){
    regulons_list[[length(regulons_list)+1]] <- get_omnipath(organism, resources = omni_resources)
  }

  # Dynamically bind the rows of all the entities
  regulons <- do.call("bind_rows", regulons_list) %>% dplyr::distinct()

  return(regulons)
}


# =======>> PhosphoScore uniform regulatory sites information

#' psp_phospho_parsing
#'
#' @param organism string, 'human' or 'mouse'
#' @param reg_site_path path to the PhosphoSitePlus 'Regulatory Sites' file
#' @param only_activatory Boolean, if TRUE, only interactions regulating protein activity are included
#'
#' @return table of regulatory phosphosites derived from PsP
#' @export
#'
#' @examples
psp_phospho_parsing <- function(organism, reg_site_path, only_activatory = TRUE){

  regulatory_sites <- tibble::tibble(read.delim2(reg_site_path, sep ='\t', skip = 2))

  regulatory_org <- regulatory_sites %>%
    dplyr::filter(ORGANISM == organism & #select organism of interest
                    grepl('p$', regulatory_sites$MOD_RSD)) %>% #select only phosphosites
    dplyr::mutate(MOD_RSD = stringr::str_remove(MOD_RSD, '-p')) # remove -p in the column

  # Read SIGNOR dictionary
  phospho2SIGNOR <- read.delim2(paste0(path_package, 'extdata/diz_FUNCTION_phosphosite.txt'),
                                sep ='\t',header=FALSE,
                                col.names = c('phospho', 'SIGNOR','sign'))
  # PSP
  phosphoscore_table <- regulatory_org %>%
    dplyr::select(c('ENTITY' = 'PROTEIN', 'ID' = 'ACC_ID',
                    'RESIDUE' = 'MOD_RSD', 'SEQUENCE' = 'SITE_...7_AA',
                    'ON_FUNCTION')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'S', 'Ser')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'T', 'Thr')) %>%
    dplyr::mutate(RESIDUE = stringr::str_replace(RESIDUE, 'Y', 'Tyr')) %>%
    dplyr::mutate(SEQUENCE = stringr::str_to_upper(SEQUENCE)) %>%
    tidyr::separate_rows('ON_FUNCTION', sep = '; ')

  phosphoscore_table <- dplyr::left_join(phosphoscore_table,
                                         phospho2SIGNOR, by = c('ON_FUNCTION' = 'phospho')) %>%
    dplyr::filter(!is.na(sign) & (sign == 1 | sign == -1))

  if(only_activatory){
    phosphoscore_table <- phosphoscore_table %>%
      dplyr::filter(grepl('activity', SIGNOR)) %>%
      dplyr::select(ENTITY, ID, RESIDUE, SEQUENCE, sign)

  }else{
    phosphoscore_table <- phosphoscore_table %>%
      dplyr::select(ENTITY, ID, RESIDUE, SEQUENCE, sign)
  }
  return(phosphoscore_table)
}

#' get_phosphoscore_info
#'
#' @param resources vector, source for regulatory phosphosites information
#' @param organism string, 'human' or 'mouse'
#' @param psp_reg_site_path path to the PhosphoSitePlus 'Regulatory Sites' file
#' @param only_activatory Boolean, if TRUE, only interactions regulating protein activity are included
#'
#' @returntable of regulatory phosphosites derived user-defined sources
#' @export
#'
#' @examples
get_phosphoscore_info <- function(resources = c('SIGNOR', 'PsP'), organism,
                                  psp_reg_site_path = NULL, only_activatory = TRUE){

  phosphoscore_list <- list()

  if('SIGNOR' %in% resources){
    SIGNOR <- signor_parsing(organism, direct = TRUE, only_activatory = only_activatory)$interactions
    DBs <- SIGNOR %>% dplyr::filter(INTERACTION != '0')
    # create a key GN-SEQ in PKN
    DBs <- DBs %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITYB, '-', SEQUENCE))
    phos_mech <- DBs %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM) %>%
      dplyr::filter(grepl('*phos*', DBs$MECHANISM)) %>%
      dplyr::arrange(PHOSPHO_KEY_GN_SEQ) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ != '') %>%
      dplyr::distinct()

    phosphoscore_list[[length(phosphoscore_list) + 1]] <- phos_mech
  }

  if('PsP' %in% resources){
    if(is.null(reg_site_path)){
      stop('Please provide a valid path to RegulatorySites information of PhosphoSitePlus database')
    }else{
      PsP <- psp_phospho_parsing(organism = organism,
                                 reg_site_path = reg_site_path,
                                 only_activatory = only_activatory)

      phosphosite_mech <- PsP %>% dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(ENTITY, '-', SEQUENCE),
                                                IDB = ID,
                                                INTERACTION = sign,
                                                MECHANISM = 'phosphorylation') %>%
        dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, INTERACTION, MECHANISM)

      phosphoscore_list[[length(phosphoscore_list) + 1]] <- phosphosite_mech
    }
  }

  # Dynamically bind the rows of all the entities
  phos_mech <- do.call("bind_rows", phosphoscore_list) %>% dplyr::distinct()

  # add activation row
  phos_mech$ACTIVATION <- ''

  phos_mech$ACTIVATION[(phos_mech$MECHANISM == 'dephosphorylation' & phos_mech$INTERACTION == '-1') |
                         (phos_mech$MECHANISM == 'phosphorylation' & phos_mech$INTERACTION == '1')] <- '1'

  phos_mech$ACTIVATION[(phos_mech$MECHANISM == 'phosphorylation' & phos_mech$INTERACTION == '-1') |
                         (phos_mech$MECHANISM == 'dephosphorylation' & phos_mech$INTERACTION == '1')] <- '-1'

  # select only phosphosites not controversial to use
  good_phos_df <- phos_mech %>%
    dplyr::distinct(PHOSPHO_KEY_GN_SEQ, ACTIVATION, .keep_all = T) %>%
    # grouping by same effect of phosphosite
    # DEPH -1 or PHOS 1,
    # DEPH 1 or PHOS -1
    dplyr::count(PHOSPHO_KEY_GN_SEQ) %>%
    dplyr::filter(n == 1) %>% # select only the ones with one effect
    dplyr::select(PHOSPHO_KEY_GN_SEQ)

  good_phos <- as.vector(unlist(good_phos_df))

  phosphoscore_table_final <- phos_mech %>% dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% good_phos) %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, IDB, ACTIVATION) %>% dplyr::rename(UNIPROT = IDB) %>%
    dplyr::distinct()

  return(phosphoscore_table_final)
}
