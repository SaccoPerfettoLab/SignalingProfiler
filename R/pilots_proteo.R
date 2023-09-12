

paste0('Patient', c(9,12,19,24,26,32,33,38,40,43)) -> tr_patients


# ==========================================================================
# Choose and process PKNs
# ==========================================================================
background_tr <- read_tsv('../../Lab/UseCase/processed_input/cleaning/transcriptome_logTPM_clean.tsv')
background_prot <- read_tsv('../../Lab/UseCase/processed_input/cleaning/proteomics.tsv')
background_phos <- read_tsv('../../Lab/UseCase/processed_input/cleaning/phosphoproteome_clean.tsv')

# * Choose indirect PKN ** #
PKN_table_ind <- choose_PKN(organism = 'human',
                            with_atlas = FALSE,
                            direct = FALSE)

PKN_expressed_ind <- preprocess_PKN(omics_data = list(background_tr, background_prot, background_phos),
                                    PKN_table = PKN_table_ind)

# * Choose direct PKN ** #
PKN_table_dir <- choose_PKN(organism = 'human',
                            with_atlas = FALSE,
                            direct = TRUE)

PKN_expressed_dir <- preprocess_PKN(omics_data = list(background_tr, background_prot, background_phos),
                                    PKN_table = PKN_table_dir)




# Define cplex paramteres
solver = 'cplex'
carnival_options = default_CARNIVAL_options(solver)

# ==========================================================================
# Creation of TRANSCRIPTIONAL NETWORK
# ==========================================================================

# Precheck of mutations
# * Starting proteins (mutations)
mutations <- read_tsv('../test_data/input/mutations_annotations/impact_matrix.tsv')

source_df <- mutations %>%
  #filter(patient_id == patient) %>%
  pivot_longer(cols = FLT3:`AML1-ETO`) %>%
  select(-patient_id, -value) %>%
  distinct() %>%
  rename('gene_name' = 'name') %>%
  #filter(!is.na(final_score)) %>%
  mutate(
    mf = 'rec', method = 'user',
    gene_name = str_to_upper(str_replace_all(gene_name, "[^[:alnum:]]", '_')))

coverage_of_inferred_proteins_in_db(source_df,
                                    organism = 'human',
                                    direct =  F,
                                    with_atlas = T,
                                    report = FALSE)

# * Choose patient

#patients_id <- c(7)

path_report <- '../test_data/result/report_transcriptional_naive_network_prot.txt'
file_conn <- file(path_report, open = "wt")
for(i in c(9,12,19,24,26,32,33,38,40,43)){

  patient <- patients[i]
  cat(paste0(patient, '\n'), file = file_conn, append = T)
  # * Target proteins
  target_df <- read_tsv(paste0(proteins_dir, patient, '_proteins.tsv'))

  # * Starting proteins (mutations)
  mutations <- read_tsv('../test_data/input/mutations_annotations/impact_matrix.tsv')

  source_df <- mutations %>%
    filter(patient_id == patient) %>%
    pivot_longer(cols = FLT3:`AML1-ETO`) %>%
    select(-patient_id) %>%
    rename('gene_name' = 'name', 'final_score' = 'value') %>%
    filter(!is.na(final_score)) %>%
    mutate(final_score = ifelse(final_score == 0, -1, 1),
           mf = 'rec', method = 'user',
           gene_name = str_to_upper(str_replace_all(gene_name, "[^[:alnum:]]", '_')))

  # Remove patients without mutations
  if(nrow(source_df) == 0){
    cat(paste0(patient, '\n', 'No mutations \n'), file = file_conn, append = T)
    next
  }

  # ============================================================================
  # NAIVE NETWORK BUILDING for TRANSCRIPTIONAL RUN
  # ============================================================================
  path_length = 3
  tr_naive <- one_layer_naive_network(starts_gn = source_df$gene_name,
                                      targets_gn = target_df$gene_name,
                                      PKN_table = PKN_expressed_ind, #or PKN_mouse
                                      max_length = path_length,
                                      rds_path = paste0(naive_dir, patient, '_transcriptional_naive_network.RDS'),
                                      sif_path = paste0(naive_dir, patient, '_transcriptional_naive_network.sif'))


  # ============================================================================
  # OPTIMIZATION RUN
  # ============================================================================

  naive_network <- read_tsv(paste0(naive_dir, patient, '_transcriptional_naive_network.sif'),
                            col_names = c('source', 'interaction', 'target'))

  source_df <- convert_gene_name_in_uniprotid(source_df, 'human')

  # * Target proteins
  target_df <- target_df %>% filter(!gene_name %in% source_df$gene_name)

  output_tr <- run_carnival_and_create_graph(source_df = source_df,
                                             target_df = target_df,
                                             naive_network = unique(naive_network),
                                             proteins_df = bind_rows(source_df, target_df),
                                             organism = 'human',
                                             topbottom = TRUE,
                                             carnival_options = carnival_options,
                                             files = TRUE,
                                             with_atlas = FALSE,
                                             direct = FALSE,
                                             path_sif = paste0(carnival_dir, patient, '_transcriptional_opt_network.sif'),
                                             path_rds = paste0(carnival_dir, patient, '_transcriptional_opt_object.RDS'))

  RCy3::createNetworkFromIgraph(output_tr$igraph_network, title = patient, collection = 'Transcriptional, proteoscore')

}

close(file_conn)


# ==========================================================================
# Creation of SIGNALLING NETWORK
# ==========================================================================

path_report <- '../test_data/result/report_signaling_naive_network_prot.txt'
file_conn <- file(path_report, open = "wt")

for(i in  c(9,12,19,24,26,32,33,38,40,43)){

  patient <- patients[i]
  cat(paste0(patient, '\n'), file = file_conn, append = T)

  # * Target proteins *
  target_df <- read_tsv(paste0(proteins_dir, patient, '_proteins.tsv'))# %>%
    #filter(method != 'proteoscore' | (method == 'proteoscore' & mf %in% c('kin', 'phos')))

  #   Divide proteins according to the molecular function
  kin_phos <- target_df %>%
    dplyr::filter(mf %in% c('kin', 'phos'))
  other <- target_df %>%
    dplyr::filter(mf == 'other')
  tfs <- target_df %>%
    dplyr::filter(mf == 'tf')

  if( length(c(kin_phos$gene_name, other$gene_name)) == 0){
    cat(paste0(patient, '\n', 'No kinases/phosphoproteins '), file = file_conn, append = T)
    next
  }

  # * Starting proteins (mutations) *
  mutations <- read_tsv('../test_data/input/mutations_annotations/impact_matrix.tsv')

  source_df <- mutations %>%
    filter(patient_id == patient) %>%
    pivot_longer(cols = FLT3:`AML1-ETO`) %>%
    select(-patient_id) %>%
    rename('gene_name' = 'name', 'final_score' = 'value') %>%
    filter(!is.na(final_score)) %>%
    mutate(final_score = ifelse(final_score == 0, -1, 1),
           mf = 'rec', method = 'user',
           gene_name = str_to_upper(str_replace_all(gene_name, "[^[:alnum:]]", '_')))

  if(nrow(source_df) == 0){
    cat(paste0(patient, '\n', 'No mutations \n'), file = file_conn, append = T)

    next
  }

  # ============================================================================
  # NAIVE NETWORK BUILDING for SIGNALLING RUN
  # ============================================================================

  signaling_naive <- three_layer_naive_network(starts_gn = source_df$gene_name,
                                               intermediate1_gn = kin_phos$gene_name,
                                               intermediate2_gn = unique(c(kin_phos$gene_name, other$gene_name)),
                                               targets_gn = tfs$gene_name,
                                               PKN_table = PKN_expressed_dir,
                                               both_intermediates = TRUE,
                                               keep_only_connected = TRUE,
                                               max_length_1 = 3,
                                               max_length_2 = 1,
                                               max_length_3 = 3,
                                               rds_path = paste0(naive_dir, patient, '_signaling_naive_network.RDS'),
                                               sif_path = paste0(naive_dir, patient, '_signaling_naive_network.sif'))

  # ============================================================================
  # OPTIMIZATION RUN
  # ============================================================================

  sig_network <- read_tsv(paste0(naive_dir, patient, '_signaling_naive_network.sif'),
                          col_names = c('source', 'interaction', 'target'))

  source_df <- convert_gene_name_in_uniprotid(source_df, 'human')

  if(!source_df$gene_name %in% unique(c(sig_network$source, sig_network$target))){
    cat(paste0('No', paste0(source_df$gene_name, collapse = '|'), 'in the network for ', patient), file = file_conn, append = T)
    next
  }

  # * Target proteins
  target_df <- target_df %>% filter(!gene_name %in% source_df$gene_name)

  output_sig <- run_carnival_and_create_graph(source_df = source_df,
                                              target_df = target_df,
                                              naive_network = unique(sig_network),
                                              proteins_df = bind_rows(source_df, target_df),
                                              organism = 'human',
                                              carnival_options = carnival_options,
                                              files = TRUE,
                                              with_atlas = FALSE,
                                              direct = TRUE,
                                              path_sif = paste0(carnival_dir, patient, '_signaling_opt_carnival.sif'),
                                              path_rds = paste0(carnival_dir, patient, '_signaling_opt_object.RDS'))

  if(is.null(output_sig)){
    cat(paste0('No signaling network for ', patient, 'with ', paste0(source_df$gene_name, collapse = '|')), file = file_conn, append = T)
    next
  }

  RCy3::createNetworkFromIgraph(output_sig$igraph_network, title = patient, collection = 'Signaling, proteoscore')

}
close(file_conn)

for( i in c(9,12,19,24,26,32,33,38,40,43)){

  patient <- patients[i]

  # * Starting proteins (mutations) *
  mutations <- read_tsv('../test_data/input/mutations_annotations/impact_matrix.tsv')

  source_df <- mutations %>%
    filter(patient_id == patient) %>%
    pivot_longer(cols = FLT3:`AML1-ETO`) %>%
    select(-patient_id) %>%
    rename('gene_name' = 'name', 'final_score' = 'value') %>%
    filter(!is.na(final_score)) %>%
    mutate(final_score = ifelse(final_score == 0, -1, 1),
           mf = 'rec', method = 'user',
           gene_name = str_to_upper(str_replace_all(gene_name, "[^[:alnum:]]", '_')))

  source_df <- convert_gene_name_in_uniprotid(source_df, 'human')

  # * Target proteins *
  target_df <- read_tsv(paste0(proteins_dir, patient, '_proteins.tsv')) %>%
    filter(method != 'proteoscore' & !gene_name %in% source_df$gene_name)

  proteins_df <- bind_rows(source_df, target_df)

  if(nrow(source_df) == 0){
    next
  }

  if(file.exists(paste0(carnival_dir, patient, '_signaling_opt_object.RDS')) &
     file.exists(paste0(carnival_dir, patient, '_transcriptional_opt_object.RDS'))){
    flag = ''

    # Union of graphs
    transc_graph = readRDS(paste0(carnival_dir, patient, '_transcriptional_opt_object.RDS'))$igraph_network
    sign_graph = readRDS(paste0(carnival_dir, patient, '_signaling_opt_object.RDS'))$igraph_network

    SP_object <- union_of_graphs(graph_1 = transc_graph,
                                 graph_2 = sign_graph,
                                 proteins_df = proteins_df,
                                 files = TRUE,
                                 path_sif = paste0(carnival_dir, patient, '_union_model.sif'),
                                 path_rds = paste0(carnival_dir, patient, '_union_object.RDS'))

  }else if(!file.exists(paste0(carnival_dir, patient, '_signaling_opt_object.RDS')) &
           file.exists(paste0(carnival_dir, patient, '_transcriptional_opt_object.RDS'))){

    SP_object <- readRDS(paste0(carnival_dir, patient, '_transcriptional_opt_object.RDS'))
    flag = '*'

  }

  phospho_df <- read_tsv(paste0(input_dir, 'phos/', patient, '_for_sp.tsv'))

  validated_SP_object <- expand_and_map_edges(optimized_object = SP_object,
                                              organism = 'human',
                                              phospho_df = phospho_df,
                                              files = TRUE,
                                              with_atlas = TRUE,
                                              direct = FALSE,
                                              path_sif = paste0(carnival_dir, patient, '_validated_graph.sif'),
                                              path_rds = paste0(carnival_dir, patient, '_validated_SP_object.RDS'))

  RCy3::createNetworkFromIgraph(validated_SP_object$igraph_network,
                                title = paste0(patient, ifelse(flag == '*', '*', '')),
                                collection = 'Union network, proteoscore')

}



