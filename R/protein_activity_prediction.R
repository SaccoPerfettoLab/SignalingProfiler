# scripts for the protein activity prediction module

###############################################################################
###                                VIPER                                    ###
###############################################################################

#' Create VIPER format
#'
#' @param omic_data phosphoproteomics or transcriptomics dataset
#' @param analysis type of analysis ksea or tfea
#' @param significance boolean value (TRUE or FALSE) if you want all analytes
#' or just significant ones
#'
#' @return dataset in VIPER format
#'
#' @examples
create_viper_format <- function(omic_data, analysis, significance){
  if(significance == TRUE){
    omic_filtered <- omic_data %>% dplyr::filter(significant == '+')

    if(analysis == 'tfea'){

      VIPER_format <- omic_filtered %>%
        dplyr::select(ID = gene_name, logFC = difference, t = difference, adj.P = logpval) %>%
        dplyr::mutate_at(c('logFC', 't'), as.numeric)

    }else if(analysis == 'ksea'){

      omic_filtered$ID <- paste0(omic_filtered$UNIPROT,
                                 '-',
                                 omic_filtered$aminoacid,
                                 '-',
                                 omic_filtered$position)

      VIPER_format <- omic_filtered %>%
        dplyr::select(ID, logFC = difference, t = difference, adj.P = logpval) %>%
        dplyr::mutate_at(c('logFC', 't'), as.numeric)

    }else{
      stop('please provide a valid analysis name')
    }

  }else{
    if(analysis == 'tfea'){

      VIPER_format <- omic_data %>%
        dplyr::select(ID = gene_name, logFC = difference, t = difference, adj.P = logpval) %>%
        dplyr::mutate_at(c('logFC', 't'), as.numeric)
    }else if(analysis == 'ksea'){

      omic_data$ID <- paste0(omic_data$UNIPROT,
                             '-',
                             omic_data$aminoacid,
                             '-',
                             omic_data$position)

      VIPER_format <- omic_data %>%
        dplyr::select(ID, logFC = difference, t = difference, adj.P = logpval) %>%
        dplyr::mutate_at(c('logFC', 't'), as.numeric)
    }else{
      stop('please provide a valid analysis name')
    }
  }
  return(VIPER_format)
}

#' Run VIPER analysis
#'
#' @param viper_format omic dataset in VIPER format
#' @param analysis string representing 'tfea' or 'ksea' analysis
#' @param organism string reporting the organism
#' @param reg_minsize viper function param: minimum regulon size to consider
#'
#' @return a list containing significantly enriched proteins and
#' all inferred proteins
#' @export
#'
#' @examples
run_viper <- function(viper_format, analysis, organism, reg_minsize){

  # viper_format <- v
  # analysis <- 'ksea'
  # organism <- 'mouse'
  # reg_minsize = 1

  # create the format needed for viper
  diff_matrix <- create_matrix_from_VIPER_format(viper_format)

  if(analysis == 'tfea'){
    if(organism == 'human'){
      regulons <- tfea_db_human
    }else if(organism == 'mouse'){
      regulons <- tfea_db_mouse
    }else{
      stop('please provide a valid organism name')
    }

  }else if(analysis == 'ksea'){
    if(organism == 'human'){
      regulons <- ksea_db_human
    }else if(organism == 'mouse'){
      regulons <- ksea_db_mouse
    }else{
      stop('please provide a valid organism name')
    }

  }else{
    stop('please provide a valid analysis name')
  }

  regulons <- regulons %>% dplyr::distinct()

  if(identical(intersect(rownames(diff_matrix), regulons$target), character(0))){
    stop('No measured analytes found in regulons, please check if regulons\' organism and experimental organism match')
    return(NULL)

  }else{
    viper_object <- viper::msviper(diff_matrix,
                                   dorothea::df2regulon(regulons),
                                   minsize = reg_minsize,
                                   #method = 'rank',
                                   ges.filter = FALSE,
                                   cores = 1,
                                   verbose = FALSE,
                                   pleiotropy = TRUE)

    activities_stat <- as.data.frame(viper_object$es$nes.bt) %>%
      dplyr::rename(NES = t) %>%
      dplyr::mutate(pvalues = viper_object$es$p.value) %>%
      tibble::rownames_to_column('gene_name')

    prot_sign <- activities_stat %>% dplyr::filter(pvalues < 0.05)
  }
  return(list(sign = prot_sign,
              all = activities_stat))
}

## -- EA VIPER CORRECTION -- ##

#' Correct viper output with hypergeometric test
#'
#' @param omic_data dataset containing experimental omic data
#' @param viper_output dataset containing viper enriched proteins
#' @param analysis type of analysis
#' @param organism organism
#'
#' @return table with weight for viper score correction
#' @export
#'
#' @examples
run_hypergeometric_test <- function(omic_data, viper_output,
                          analysis, organism){

  # omic_data <- tr_df
  # viper_output <- a$sign
  # analysis <- 'tfea'
  # organism <- 'mouse'

  # READING UNIVERSE INPUT AKA EXPERIMENTAL DATA #
  if(analysis == 'tfea'){
    uni_meas_v <- omic_data$gene_name
    uni_sign_f <- omic_data %>% dplyr::filter(significant == '+')
    uni_sign_v <- uni_sign_f$gene_name

  }else if(analysis == 'ksea'){

    uni_meas_f <- omic_data %>%
      dplyr::mutate(phosphositeID = paste0(toupper(UNIPROT),'-',
                                           aminoacid,'-',
                                           position))
    uni_meas_v <- uni_meas_f$phosphositeID

    uni_sign_f <- omic_data %>%
      dplyr::filter(significant == '+') %>%
      dplyr::mutate(phosphositeID = paste0(toupper(UNIPROT),'-',
                                           aminoacid,'-',
                                           position))
    uni_sign_v <- uni_sign_f$phosphositeID
  }else{
    stop('please provide a valid analysis type')
  }

  # REGULON INPUT
  #   # -- CREATE DATAFRAME WITH COUNT OF MEASURED AND
  #        SIGNIFICANT PHOSPHOSITES FOR EACH KIN/TF-- #

  if(analysis == 'tfea'){
    # set db
    if(organism == 'human'){
      df_regulons <- tfea_db_human
    }else if(organism == 'mouse'){
      df_regulons <- tfea_db_mouse
    }else{
      stop('please prove a valid organism')
    }

  }else if(analysis == 'ksea'){

    if(organism == 'human'){
      df_regulons <- ksea_db_human
    }else if(organism == 'mouse'){
      df_regulons <- ksea_db_mouse
    }else{
      stop('please provide a valid organism')
    }

  }else{stop('please provide a valid analysis type')}

  # from viper analysis get all measured analytes
  all_enriched_matrix <- create_matrix_from_VIPER_format(create_viper_format(omic_data, analysis, significance = FALSE))

  pho_in_reg <- df_regulons[df_regulons$target %in% (rownames(all_enriched_matrix)),] %>%
    dplyr::arrange(tf) %>%
    plyr::count('tf') %>%
    dplyr::rename(n = freq)

  # from viper analysis get all significant analytes
  all_significant_matrix <- create_matrix_from_VIPER_format(create_viper_format(omic_data, analysis, significance = TRUE))

  pho_in_reg_sign <- df_regulons[df_regulons$target %in% (rownames(all_significant_matrix)),] %>%
    dplyr::arrange(tf) %>%
    plyr::count('tf') %>%
    dplyr::rename(n = freq)


  pr_joined <- dplyr::left_join(viper_output, pho_in_reg, by = c('gene_name' = 'tf')) %>%
    dplyr::rename('Measured' = 'n')
  pr_joined <- dplyr::left_join(pr_joined, pho_in_reg_sign, by = c('gene_name' = 'tf')) %>%
    dplyr::rename('Significant' = 'n')


  # ENRICHMENT ANALYSIS WITH FISHER EXACT TEST
  proteins <- unlist(pr_joined$gene_name)
  pr_joined$pWeight <- NA

  for(protein in proteins){
    significant_members <- df_regulons[df_regulons$target %in% (rownames(all_significant_matrix)),] %>%
      dplyr::arrange(tf) %>%
      dplyr::filter(tf == protein) %>%
      dplyr::select(target) %>% unlist() %>% as.character

    regulon_members <- df_regulons[df_regulons$target %in% (rownames(all_enriched_matrix)),] %>%
      dplyr::arrange(tf) %>%
      dplyr::filter(tf == protein) %>%
      dplyr::select(target) %>% unlist() %>% as.character

    # creating Venn object
    VennObject <- RVenn::Venn(list(RegulonMeasured = regulon_members,
                                   RegulonSignificant = significant_members,
                                   UniverseMeasured = uni_meas_v,
                                   UniverseSignificant = uni_sign_v))

    # run enrichment test
    v <- RVenn::enrichment_test(VennObject, 'RegulonMeasured', 'RegulonSignificant', univ = uni_meas_v)

    pr_joined$pWeight[pr_joined$gene_name == protein] <- v$Significance
  }

  pr_joined <- pr_joined %>%
    dplyr::rename(reg_exp_meas = Measured,
                  reg_exp_sign = Significant)

  return(pr_joined)
}

#' Create a matrix from VIPER format dataset
#'
#' @param viper_format a dataset in VIPER format
#'
#' @return a matrix containing ID and t statistic
#'
#' @examples
create_matrix_from_VIPER_format <- function(viper_format){
  diff_matrix <- viper_format %>%
    dplyr::select(ID, t) %>%
    dplyr::filter(!is.na(t)) %>%
    tibble::column_to_rownames(var = "ID") %>%
    as.matrix()
  return(diff_matrix)
}

#
#' Title
#'
#' @param ea_output dataset representing ea output with viper and enrichment score
#'
#' @return dataset in input with weighted nes
#' @export
#'
#' @examples
weight_viper_score <- function(ea_output){

  #ea_output <- pr_joined
  #ea_output <- h

  ea_output_log <- ea_output %>%
    dplyr::filter(pWeight != 1) %>%
    dplyr::mutate(weight = -log(pWeight))

  ea_output_log$weight[ea_output_log$weight == Inf] <- 10
  ea_output_log$weight <- ea_output_log$weight/10

  ea_output_log <- ea_output_log %>%
    dplyr::mutate(weightedNES = weight * NES) %>%
    dplyr::arrange(gene_name)

  return(ea_output_log)
}

#' Filter VIPER output for molecular function
#'
#' @param inferred_proteins_mf dataset of inferred proteins with mf column
#' @param analysis string type of analysis
#'
#' @return dataset of inferred proteins with only the specific mf
#' @export
#'
#' @examples
filter_VIPER_output <- function(inferred_proteins_mf, analysis){
  if(analysis == 'tfea'){
    inferred_proteins_mf <- inferred_proteins_mf %>% dplyr::filter(mf == 'tf')
  }else if(analysis == 'ksea'){
    inferred_proteins_mf <- inferred_proteins_mf %>% dplyr::filter(mf == 'kin' | mf == 'phos')
  }

  return(inferred_proteins_mf)
}

#' run footprint based analysis
#'
#' @param omic_data dataset of experimental measured phosphosites or transcripts
#' @param analysis string, tfea or ksea
#' @param organism string, mouse or human
#' @param reg_minsize integer value, VIPER parameter for minsize of regulon
#' @param exp_sign boolean value, TRUE use only significant analytes,
#' FALSE: all measured analytes
#' @param hypergeom_corr boolean value, TRUE apply hypergeometric correction,
#' FALSE no correction
#' @param GO_annotation boolean value, TRUE perform GO molecular function annotaiton, FALSE default value
#' @return dataset of inferred proteins: transcription factor (tfea)
#' or kinases and phosphatases (ksea)
#' @export
#'
#' @examples
run_footprint_based_analysis <- function(omic_data, analysis, organism,
                                         reg_minsize, exp_sign,
                                         hypergeom_corr,
                                         GO_annotation = FALSE){

  message(' ** RUNNING FOOTPRINT BASED ANALYSIS ** ')
  message('Credits to Prof. Julio Saez-Rodriguez. For more information read this article: Dugourd A, Saez-Rodriguez J. Footprint-based functional analysis of multiomic data. Curr Opin Syst Biol. 2019 Jun;15:82-90. doi: 10.1016/j.coisb.2019.04.002. PMID: 32685770; PMCID: PMC7357600.')

  # omic_data <- readRDS('./data/JMD_phospho.RDS')
  # # omic_data <- phospho_toy_df
  # analysis <- 'ksea'
  # organism <- 'mouse'
  # reg_minsize <- 1
  # exp_sign <- FALSE
  # hypergeom_corr <- TRUE
  # GO_annotation = TRUE

  library(tidyverse)
  # run viper analysis
  viper_format <- create_viper_format(omic_data, analysis, significance = exp_sign)


  message('Starting VIPER analysis')
  output <- run_viper(viper_format, analysis, organism, reg_minsize)
  if(hypergeom_corr == TRUE){
    message('Starting hypergeometric test correction')
    output <- weight_viper_score(run_hypergeometric_test(omic_data,output$sign, analysis, organism))

    output_uniprot <- convert_gene_name_in_uniprotid(output, organism)
  }else{
    output_uniprot <- convert_gene_name_in_uniprotid(output$sign, organism)
  }

  #output_uniprot$UNIPROT

  if(GO_annotation == TRUE){
    message('GO molecular function annotation')

    output_uni_mf <- molecular_function_annotation(output_uniprot)
    output_final <- filter_VIPER_output(output_uni_mf, analysis) %>% dplyr::distinct()
    #footprint_output <- output_final
    return(output_final)
  }
  return(output_uniprot)


}

###############################################################################
###                             PHOSPHOSCORE                                ###
###############################################################################

#' Compute phosphoSCORE
#'
#' @param phosphoproteomic_data dataset of phosphoproteomics measurments
#' @param organism string human or mouse
#' @param activatory boolean value
#' @param path_fasta optional
#' @param local DEVELOPMENTAL PURPOSES
#' @param GO_annotation boolean value, TRUE perform GO molecular function annotaiton, FALSE default value
#' @param blastp_path optional, Windows path of blastp
#'
#' @return phosphoscore dataset with gene_name, inferred activity and
#' used phosphosites from experimental data
#' @export
#'
#' @examples
phosphoscore_computation <- function(phosphoproteomic_data,
                                    organism,
                                    activatory,
                                    path_fasta = './phospho.fasta',
                                    blastp_path = NULL,
                                    local = FALSE,
                                    GO_annotation = FALSE){
  message('** RUNNING PHOSPHOSCORE ANALYSIS **')

  # phosphoproteomic_data <- readRDS('./data/TKD_phospho.RDS') %>% mutate_at('gene_name', str_to_title)
  # path_fasta = './phospho.fasta'
  # organism = 'human'
  # activatory = TRUE
  # local = TRUE
  # GO_annotation = FALSE

  if(organism == 'mouse' | organism =='human'){
    phosphoscore_df_output <- map_experimental_on_regulatory_phosphosites(phosphoproteomic_data,
                                                                   organism, activatory, path_fasta)
    phosphoscore_df <- phosphoscore_df_output$phosphoscore_df
  }else if(organism == 'hybrid'){
    phosphoscore_df_output <- phospho_score_hybrid_computation(phosphoproteomic_data,
                                                        organism, activatory, blastp_path, path_fasta, local)

    phosphoscore_df <- phosphoscore_df_output$phosphoscore_df
  }else{
    stop('please provide a valid organism')
  }

  raw_output <- phosphoscore_df %>%
    #tidyr::separate(PHOSPHO_KEY_GN_SEQ, into = c('h_gene_name', 'phosphoseq'), sep = '-') %>%
    #dplyr::select(-c('h_gene_name')) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(phosphoscore = mean(inferred_activity))


  exp_fc_sub <- phosphoscore_df_output$used_exp_data %>%
    dplyr::mutate_at('PHOSPHO_KEY_GN_SEQ', toupper) %>%
    dplyr::select(gene_name, PHOSPHO_KEY_GN_SEQ, aminoacid, position) %>% dplyr::distinct()

  exp_fc_sub <- exp_fc_sub %>%
    dplyr::mutate(aa = paste0(aminoacid, position),
                  gene_name = (gene_name)) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n_sign_phos = dplyr::n(),
                     phos = paste0(aa, collapse = ';'))

  output <- dplyr::left_join(raw_output, exp_fc_sub, by = 'gene_name')

  if(organism == 'hybrid' & 'source_org' %in% colnames(output)){
    genes <- phosphoscore_df %>%
      dplyr::select(gene_name, source_org) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(source_org = paste0(source_org, collapse = ';'))
    output <- dplyr::left_join(output, genes, by = c('gene_name'))
  }

  # add uniprot and molecular function
  output <- convert_gene_name_in_uniprotid(bio_dataset = output, organism = organism)

  if(GO_annotation == TRUE){
    output <- molecular_function_annotation(output) %>%
      dplyr::distinct()
    return(output)
  }

  return(output)
}

#
#' Create fasta from phosphoproteomics files
#'
#' @param phospho_df dataset of phosphoproteomics
#' @param path path of fasta file
#'
#' @return nothing
#' @export
#'
#' @examples
create_fasta <- function(phospho_df, path){
  message('Creating fasta file from phosphoproteomics')
  phospho_df <- phospho_df %>%
    dplyr::mutate(gene_name = toupper(unlist(lapply(stringr::str_split(gene_name, ';'),
                                             function(x){x[1]})))) %>%
    dplyr::arrange(gene_name)

  # if it is not a 15mer then convert it
  if(nchar(phospho_df$sequence_window[1]) != 15){
    center <- (nchar(phospho_df$sequence_window[1])+1)/2
    phospho_df <- phospho_df %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
  }else{ #if it is rename the column just to have a homogeneous variable
    phospho_df <- phospho_df %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }

  fasta <- ''
  for(i in c(1:length(phospho_df$UNIPROT))){
    write(paste0('>', phospho_df$gene_name[i], '-',
                 phospho_df$sequence_window_sub[i], '-',
                 phospho_df$aminoacid[i], '-',
                 phospho_df$position[i], '-mouse\n',
                 gsub('_', '', phospho_df$sequence_window_sub[i])),
          path, append = TRUE)
  }
  return(NULL)
}

#
#' Title
#'
#' @param path_experimental_fasta_file string representing path of
#' phosphoproteomics data in fasta file
#' @param all boolean value representing if you want all alingments
#' or only the ones with same gene_name
#' #' @param local FOR DEVELOPMENTAL PURPOSES TO DELETE
#' @param blastp_path string specifying the blastp path on Windows
#' @param local boolean TRUE or FALSE, for developing purposes
#'
#' @return a list containing:
#' a dataframe of the sure mouse phosphopeptides aligned on human
#' a dataframe of the unsure human phosphopeptides aligned on human having
#' different gene_name, to check manually!
#' @export
#'
#' @examples
run_blast <- function(path_experimental_fasta_file, all = FALSE,
                      blastp_path = NULL, local = FALSE){

  message('Running blastp')
  #print(local)
  # local = TRUE
  if(local == TRUE){path_package <- './'
  }else{ path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')}

  if(is.null(blastp_path)){
    blastp <- paste0('blastp -query ', path_experimental_fasta_file,
                     ' -subject ', paste0(path_package, 'data/human_phosphosites_db.fasta '),
                     '-out map2.out -outfmt 7 -evalue 0.05')
  }else{
    blastp <- paste0("\"", blastp_path, "\"", ' -query ', path_experimental_fasta_file,
                     ' -subject ', paste0(path_package, 'data/human_phosphosites_db.fasta '),
                     '-out map2.out -outfmt 7 -evalue 0.05')
  }

  system(blastp)

  message('blastp finished')
  mapped <- readr::read_tsv('./map2.out',
                           comment = '#',
                           col_names = c('q_ID', 's_ID', 'pid', 'length',
                                         'mismatch', 'gapopen', 'qstart',
                                         'qend', 'sstart', 'ssend', 'evalue',
                                         'bitscore'))

  # selection of exact matching phosphosites
  mapped_exact <- mapped %>%
    dplyr::filter(qstart == 1 & qend == 15 & sstart == 1 & ssend == 15)

  mapped$position <- as.numeric(unlist(lapply(stringr::str_split(mapped$q_ID, '-'),
                                              function(x){x[4]})))

  mapped_exact <- dplyr::bind_rows(mapped %>%
                                     dplyr::filter(position < 8),mapped_exact)

  # create two columns for query protein (mouse gene name) and subject protein (human gene name)
  mapped_exact$qProt <- unlist(lapply(stringr::str_split(mapped_exact$q_ID, '-'), function(x){x[1]}))
  mapped_exact$sProt <- unlist(lapply(stringr::str_split(mapped_exact$s_ID, '-'), function(x){x[1]}))

  tocheck <- mapped_exact %>% dplyr::filter((qProt != sProt) & mismatch == 0)
  sure <- mapped_exact %>% dplyr::filter(qProt == sProt)
  system('rm ./map2.out')

  if(all == TRUE){
    return(list(mapped = sure, tocheck = tocheck))
  }
  else{
    return(list(mapped = sure))
  }
}

#' Generates hybrid regulatory db
#'
#' @param mh_alignment blastp output alignment among mouse and human
#' @param activatory boolean value
#'
#' @return hybrid database
#' @export
#'
#' @examples
generate_hybrid_db <- function(mh_alignment, activatory){

  if(activatory == TRUE){
    hreg_phos <- good_phos_df_human_act
  }else{
    hreg_phos <- good_phos_df_human_all
  }
  # regulatory role of phosphosites in human

  #mh_alignment = run_blast(path_fasta, local = local)$mapped

  good_phos_df_hybrid <- dplyr::inner_join(mh_alignment,
                                          hreg_phos,
                                          by = c('s_ID' = 'PHOSPHO_KEY_GN_SEQ')) %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ = q_ID, PHOSPHO_KEY_GN_SEQ_human = s_ID, UNIPROT, ACTIVATION) %>%
    # keeping only gene_name
    dplyr::mutate(PHOSPHO_KEY_GN_SEQ = sub('*-[T,Y,S]-\\d{1,5}-mouse', '', PHOSPHO_KEY_GN_SEQ)) %>%
    dplyr::arrange(PHOSPHO_KEY_GN_SEQ) %>%
    dplyr::distinct()

  return(good_phos_df_hybrid)
}


#' Title
#'
#' @param phosphoproteomic_data dataset of phosphoproteomics data
#' @param organism string, human mouse or hybrid
#' @param activatory boolean value, if relation activatory or all
#' @param organism string, human mouse or hybrid
#' @param local DEVELOPMENTAL PURPOSES TRUE OR FALSE
#'
#' @param path_fasta
#' @param blastp_path path of blastp in windows
#'
#' @return phosphoscore_df representing experimentally quantified phosphosites
#' associated to their fold-change
#' @export
#'
#' @examples
map_experimental_on_regulatory_phosphosites <- function(phosphoproteomic_data,
                                                        organism,
                                                        activatory,
                                                        path_fasta,
                                                        blastp_path = NULL,
                                                        local = FALSE){

  #phosphoproteomic_data <- readRDS('./data/TKD_phospho.RDS')
  # path_fasta = './phospho.fasta'
  # organism = 'mouse'
  # local = TRUE

  if(organism == 'human'){
    message('Mapping experimental phosphopeptides on human database of regulatory roles')
    if(activatory == TRUE){
      reg_phos_db <- good_phos_df_human_act
    }else{
      reg_phos_db <- good_phos_df_human_all
    }
  }else if(organism == 'mouse'){
    message('Mapping experimental phosphopeptides on mouse database of regulatory roles')
    if(activatory == TRUE){
      reg_phos_db <- good_phos_df_mouse_act
    }else{
      reg_phos_db <- good_phos_df_mouse_all
    }
  }else if(organism == 'hybrid'){
    message('Mapping mouse experimental phosphopeptides on human database of regulatory roles to enhance coverage')

    if(file.exists(path_fasta)){
     command <- paste0('rm ', path_fasta)
     system(command)
     create_fasta(phosphoproteomic_data, path_fasta)
    }else{
      create_fasta(phosphoproteomic_data, path_fasta)
    }

    reg_phos_db <- generate_hybrid_db(mh_alignment = run_blast(path_fasta, blastp_path = blastp_path, local = local)$mapped,
                                      activatory)
  }else{
    stop('please provide a valid organism')
  }

  if(nchar(phosphoproteomic_data$sequence_window[1]) != 15){
    center <- (nchar(phosphoproteomic_data$sequence_window[1])+1)/2
    phosphoproteomic_data <- phosphoproteomic_data %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
  }else{
    phosphoproteomic_data <- phosphoproteomic_data %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }


  #
  if(organism == 'hybrid'){
    exp_fc <- phosphoproteomic_data %>%
      dplyr::filter(significant == '+') %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(toupper(gene_name), '-', sequence_window_sub)) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% reg_phos_db$PHOSPHO_KEY_GN_SEQ)
  }else{
    exp_fc <- phosphoproteomic_data %>%
      dplyr::filter(significant == '+') %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0((gene_name), '-', sequence_window_sub)) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% reg_phos_db$PHOSPHO_KEY_GN_SEQ)
  }


  if(nrow(exp_fc) == 0){
    warning('No annotated regulatory phosphosites significantly modulated in your dataset')
  }else{
    phosphoscore_df <- dplyr::left_join(exp_fc, reg_phos_db, by = 'PHOSPHO_KEY_GN_SEQ') %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, ACTIVATION, difference) %>%
      dplyr::arrange(PHOSPHO_KEY_GN_SEQ) %>%
      dplyr::mutate(gene_name  = unlist(lapply(PHOSPHO_KEY_GN_SEQ,
                                               function(x){stringr::str_split(x, '-')[[1]][1]})),
                    inferred_activity = as.numeric(ACTIVATION) * as.numeric(difference)) %>%
      dplyr::distinct()

    return(list(used_exp_data = exp_fc,
                phosphoscore_df = phosphoscore_df))

  }
}

#' Title
#'
#' @param phosphoproteomic_data dataset containing experimental data
#' @param organism string specifying human, mouse or hybrid
#' @param organism boolean value, if activatory or all interactions
#' @param path_fasta optional, path of phosphoproteomic fasta file
#' @param activatory
#' @param blastp_path optional, path of Windows blastp exe file
#' @param local
#'
#' @return list containing used experimental data and phosphoscore dataframe
#'
#' @examples
phospho_score_hybrid_computation <- function(phosphoproteomic_data,
                                             organism,
                                             activatory,
                                             blastp_path = NULL,
                                             path_fasta = './phospho.fasta', local){

  # phosphoproteomic_data <- JMD_phospho_r
  # #   mutate_at('difference', as.numeric)
  # path_fasta = './phospho.fasta'
  # organism = 'hybrid'
  # local = TRUE
  # activatory = FALSE

  phosphoscore_df_mouse_output <- map_experimental_on_regulatory_phosphosites(phosphoproteomic_data, 'mouse',
                                                                              activatory = activatory, local)

  phosphoscore_df_hybrid_output <- map_experimental_on_regulatory_phosphosites(phosphoproteomic_data, 'hybrid',
                                                                               activatory = activatory,
                                                                               blastp_path = blastp_path, path_fasta = path_fasta, local = local)


  if(is.list(phosphoscore_df_mouse_output) & !is.list(phosphoscore_df_hybrid_output)){
    warning('No mapped mouse phosphosites on human, proceeding with mouse analysis only')
    return(phosphoscore_df_mouse_output)

  }else if(!is.list(phosphoscore_df_mouse_output) & is.list(phosphoscore_df_hybrid_output)){
    warning('No mouse phosphosites found, proceeding with mouse phosphosites analysis mapping on human')
    return(phosphoscore_df_hybrid_output)

  }else if(is.list(phosphoscore_df_mouse_output) & is.list(phosphoscore_df_hybrid_output)){
    phosphoscore_df_mouse <- phosphoscore_df_mouse_output$phosphoscore_df %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, inferred_activity, gene_name)

    phosphoscore_df_hybrid <- phosphoscore_df_hybrid_output$phosphoscore_df %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, inferred_activity, gene_name) %>%
      dplyr::mutate_at('gene_name', stringr::str_to_title)

    message('Computing Phosphoscore')
    joined_tables <- dplyr::full_join(phosphoscore_df_mouse, phosphoscore_df_hybrid,
                                      by = c('PHOSPHO_KEY_GN_SEQ', 'gene_name'),
                                      suffix = c('.m', '.h')) %>%
      dplyr::distinct() %>%
      dplyr::arrange(PHOSPHO_KEY_GN_SEQ)

    joined_tables$source_org <- 'both'
    joined_tables$source_org[is.na(joined_tables$inferred_activity.m)] <- 'human'
    joined_tables$source_org[is.na(joined_tables$inferred_activity.h)] <- 'mouse'

    joined_tables$inferred_activity <- NA
    joined_tables$inferred_activity[joined_tables$source_org == 'both'] <- joined_tables$inferred_activity.m[joined_tables$source_org == 'both']
    joined_tables$inferred_activity[joined_tables$source_org == 'human'] <- joined_tables$inferred_activity.h[joined_tables$source_org == 'human']
    joined_tables$inferred_activity[joined_tables$source_org == 'mouse'] <- joined_tables$inferred_activity.m[joined_tables$source_org == 'mouse']

    phosphoscore_df_flag <- joined_tables %>%
      dplyr::select(PHOSPHO_KEY_GN_SEQ, gene_name,source_org, inferred_activity)

    # experimental used data
    used_exp_data_both <- dplyr::bind_rows(phosphoscore_df_mouse_output$used_exp_data,
                                           phosphoscore_df_hybrid_output$used_exp_data) %>%
      dplyr::distinct()

    return(list(used_exp_data = used_exp_data_both,
                phosphoscore_df = phosphoscore_df_flag))
  }else{
    stop('No regulatory phosphosites found')
  }
}

#' Title
#'
#' @param footprint_output output of footprint based
#' @param phosphoscore_df output of phosphoscore based
#' @param analysis type of analysis
#'
#' @return dataset of inferred proteins with combined scores
#' @export
#'
#' @examples
combine_footprint_and_phosphoscore <- function(footprint_output, phosphoscore_df, analysis){

  # footprint_output
  # phosphoscore_df <- phosphoscore_output
  # analysis <- 'tfea'

  # if mf exists
  if('mf' %in% colnames(footprint_output) & 'mf' %in% colnames(phosphoscore_df)){
    if(analysis == 'tfea'){
      phosphoscore_df <- phosphoscore_df %>%
        dplyr::filter(mf == 'tf')
    }else if(analysis == 'ksea'){
      phosphoscore_df <- phosphoscore_df %>%
        dplyr::filter(mf == 'kin' | mf == 'phos')
    }else{stop('please provide valid analysis type')}

    comp <- dplyr::full_join(footprint_output, phosphoscore_df, by = c('gene_name', 'UNIPROT', 'mf'))
  }else{
    comp <- dplyr::full_join(footprint_output, phosphoscore_df, by = c('gene_name', 'UNIPROT'))
  }


  comp$final_score <- NA
  comp$method <- NA

  comp$final_score[!is.na(comp$weightedNES) & !is.na(comp$phosphoscore)] <- rowMeans(comp[!is.na(comp$weightedNES) & !is.na(comp$phosphoscore), c('weightedNES', 'phosphoscore')])
  comp$method[!is.na(comp$weightedNES) & !is.na(comp$phosphoscore)] <- 'VIPER+PhosphoScore'

  comp$final_score[!is.na(comp$weightedNES) & is.na(comp$phosphoscore)] <- comp$weightedNES[!is.na(comp$weightedNES) & is.na(comp$phosphoscore)]
  comp$method[!is.na(comp$weightedNES) & is.na(comp$phosphoscore)] <- 'VIPER'

  comp$final_score[is.na(comp$weightedNES) & !is.na(comp$phosphoscore)] <- comp$phosphoscore[is.na(comp$weightedNES) & !is.na(comp$phosphoscore)]
  comp$method[is.na(comp$weightedNES) & !is.na(comp$phosphoscore)] <- 'PhosphoScore'

  return(comp)
}


#' Title
#'
#' @param prot_df proteomic dataframe
#' @param organism string, 'human' or 'mouse'
#'
#' @return
#' @export
#'
#' @examples
activity_from_proteomics <- function(prot_df, organism){

  message('** COMPUTING PROTEOSCORE **')

  prot_df <- prot_df %>%
    dplyr::filter(significant == '+') %>%
    dplyr::select(gene_name, difference)

  prot_df <- convert_gene_name_in_uniprotid(prot_df, 'mouse')

  prot_df <- molecular_function_annotation(prot_df)

  return(prot_df)
}

#' Title
#'
#' @param activity_score df with activity measures derived from phosphoproteomics
#' @param proteo_score df with activity measures derived from proteomics
#'
#' @return df with combined activity measures
#' @export
#'
#' @examples
combine_activityscore_proteoscore <- function(activity_score, proteo_score){

  combined_score <- dplyr::full_join(activity_score,
                              proteo_score, by = c('gene_name', 'UNIPROT', 'mf')) %>%
    dplyr::rename(activity_score = final_score,
                  proteo_score = difference)

  combined_score$final_score <- combined_score$activity_score
  combined_score$final_score[is.na(combined_score$activity_score)] <- combined_score$proteo_score[is.na(combined_score$activity_score)]
  combined_score$final_score[is.na(combined_score$activity_score)] <- combined_score$proteo_score[is.na(combined_score$activity_score)]
  combined_score$method[is.na(combined_score$activity_score)] <- 'proteoscore'

  combined_score <- combined_score %>%
    dplyr::select(-c(activity_score, proteo_score)) %>%
    dplyr::arrange(gene_name)

  return(combined_score)
}


