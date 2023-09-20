
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

      # Separate analytes according to collapsed UNIPROT
      omic_filtered <- omic_filtered %>%
        dplyr::separate_rows(UNIPROT, sep = ';')

      # Create key for phosphosites
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
#' @param integrated_regulons boolean value, default FALSE; if TRUE, uses regulons derived from experimental data
#' @param collectri boolean, if TRUE uses collectri regulons for tfea
#'
#' @return a list containing significantly enriched proteins and
#' all inferred proteins
#' @export
#'
#' @examples
run_viper <- function(viper_format,
                      analysis,
                      organism,
                      reg_minsize,
                      integrated_regulons = FALSE,
                      collectri = FALSE){

  if(analysis == 'ksea' & collectri == TRUE){
    stop('collectri is only a \'tfea\' parameter')
  }

  # viper_format <- v
  # analysis <- 'ksea'
  # organism <- 'mouse'
  # reg_minsize = 1

  # create the format needed for viper
  diff_matrix <- create_matrix_from_VIPER_format(viper_format)

  if(analysis == 'tfea'){
    if(organism == 'human'){

      if(collectri == FALSE){
        regulons <- tfea_db_human
      }else{
        regulons <- tfea_db_human_collectri
      }

    }else if(organism == 'mouse'){
      regulons <- tfea_db_mouse
    }else{
      stop('please provide a valid organism name')
    }

  }else if(analysis == 'ksea'){
    if(organism == 'human'){
      if(integrated_regulons == TRUE){
        regulons <- ksea_db_human_atlas
      }else{
        regulons <- ksea_db_human
      }

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
#' @param collectri boolean, if TRUE uses collectri regulons
#' @param integrated_regulons boolean,  if TRUE uses Kinome Atlas regulons
#'
#' @return table with weight for viper score correction
#' @export
#'
#' @examples
run_hypergeometric_test <- function(omic_data, viper_output,
                                    analysis, organism, integrated_regulons, collectri){

  # omic_data <- tr_df
  # viper_output <- a$sign
  # analysis <- 'tfea'
  # organism <- 'mouse'
  #viper_output <- output$sign

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
      if(collectri == FALSE){
        df_regulons <- tfea_db_human
      }else{
        df_regulons <- tfea_db_human_collectri
      }
    }else if(organism == 'mouse'){
      df_regulons <- tfea_db_mouse
    }else{
      stop('please prove a valid organism')
    }

  }else if(analysis == 'ksea'){

    if(organism == 'human'){
      if(integrated_regulons == TRUE){
        df_regulons <- ksea_db_human_atlas
      }else{
        df_regulons <- ksea_db_human
      }
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
    dplyr::count(tf) %>%
    dplyr::rename('freq' = 'n')

  # from viper analysis get all significant analytes
  all_significant_matrix <- create_matrix_from_VIPER_format(create_viper_format(omic_data, analysis, significance = TRUE))

  pho_in_reg_sign <- df_regulons[df_regulons$target %in% (rownames(all_significant_matrix)),] %>%
    dplyr::arrange(tf) %>%
    dplyr::count(tf) %>%
    dplyr::rename('freq' = 'n')


  pr_joined <- dplyr::left_join(viper_output, pho_in_reg, by = c('gene_name' = 'tf')) %>%
    dplyr::rename('Measured' = 'freq')
  pr_joined <- dplyr::left_join(pr_joined, pho_in_reg_sign, by = c('gene_name' = 'tf')) %>%
    dplyr::rename('Significant' = 'freq')


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

    pvalue = phyper(q = length(significant_members),
                    m = length(uni_sign_v),
                    n = length(uni_meas_v) - length(regulon_members),
                    k = length(regulon_members), lower.tail = F)


    pr_joined$pWeight[pr_joined$gene_name == protein] <- pvalue

  }

  pr_joined <- pr_joined %>%
    dplyr::rename('reg_exp_meas' = 'Measured',
                  'reg_exp_sign' = 'Significant')


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
  # Transform hypergeomtric pvalue in -log to derive a raw_weight
  # if pWeight == 0 replace with min of pWeight otherwise the raw weight becomes Inf
  min <- min(ea_output$pWeight[ea_output$pWeight!= 0])

  ea_output_log <- ea_output %>%
    dplyr::mutate(pWeight = ifelse(pWeight == 0, min, pvalues)) %>%
    dplyr::mutate(raw_weight = -log(pWeight))

  # Create a weigth balanced on the max raw_weight
  ea_output_log <- ea_output_log %>%
    dplyr::mutate(weight = raw_weight/max(raw_weight))

  # Scale according to 4x weight the unsignificant genes
  # derived from the hypergeometric test
  ea_output_log <- ea_output_log %>%
    dplyr::mutate(weightedNES = ifelse(pWeight > 0.05,
                                       NES * weight*4,
                                       NES)) %>%
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
#' @param integrated_regulons boolean value, default FALSE; if TRUE uses regulons derived from experimental data
#' @param GO_annotation boolean value, TRUE perform GO molecular function annotaiton, FALSE default value
#' @param collectri boolean value, TRUE uses CollecTRI regulons
#' @param correct_proteomics boolean value, TRUE correct VIPER output accordin to protein abundance
#' @param prot_df dataframe, containing proteomics data in sp format
#'
#' @return dataset of inferred proteins: transcription factor (tfea)
#' or kinases and phosphatases (ksea)
#' @export
#'
#' @examples
run_footprint_based_analysis <- function(omic_data,
                                         analysis,
                                         organism,
                                         reg_minsize,
                                         exp_sign,
                                         integrated_regulons = FALSE,
                                         collectri = FALSE,
                                         hypergeom_corr,
                                         correct_proteomics = FALSE,
                                         prot_df = NULL,
                                         GO_annotation = FALSE){

  message(' ** RUNNING FOOTPRINT BASED ANALYSIS ** ')
  message('Credits to Prof. Julio Saez-Rodriguez. For more information read this article: Dugourd A, Saez-Rodriguez J. Footprint-based functional analysis of multiomic data. Curr Opin Syst Biol. 2019 Jun;15:82-90. doi: 10.1016/j.coisb.2019.04.002. PMID: 32685770; PMCID: PMC7357600.')

  if(collectri == TRUE & analysis == 'ksea'){
    stop('collectri is only a \'tfea\' parameter')

  }

  # Custom input
  # omic_data <- read_tsv('./input/transcriptomics.tsv')
  # exp_sign = FALSE
  # analysis = 'tfea'
  # organism = 'human'
  # reg_minsize = 5
  # integrated_regulons = FALSE
  # collectri = FALSE
  # prot_df <- read_tsv('./input/proteomics.tsv')
  # correct_proteomics = TRUE
  # hypergeom_corr = TRUE

  omic_data <- omic_data %>%
    dplyr::mutate(gene_name =  stringr::str_replace_all(gene_name, "[^[:alnum:]]", "_"))

  # run viper analysis
  viper_format <- create_viper_format(omic_data, analysis, significance = exp_sign)

  message('Starting VIPER analysis')

  output <- run_viper(viper_format = viper_format, analysis = analysis,
                      organism = organism, reg_minsize = reg_minsize,
                      integrated_regulons = integrated_regulons, collectri = collectri)

  # if no inferred protein from VIPER analysis
  if(nrow(output$sign) == 0){
    stop('No inferred protein found!')
  }

  # if you want to correct the analysis with proteomics
  if(correct_proteomics == TRUE){
    if(is.null(prot_df)){
      stop('Please provide proteomics')
    }else{
      viper_all <- output$all

      # Correct viper pvalues on proteomics
      # If a unsignificant enzyme in VIPER has a significant
      # modulation with same sign in proteomics, substitute the pvalue of VIPER
      # with the pvalue of proteomics, in order to make it significant

      viper_on_prot <- dplyr::left_join(viper_all, prot_df, by = c('gene_name')) %>%
        dplyr::mutate(pvalues = ifelse(!is.na(difference) & !is.na(logpval) & significant == '+' & #check on proteomics fields
                                         NES*difference > 0 & pvalues > 0.05, # check on VIPER fields
                                       10^(-logpval), # assign proteomic pvalue
                                       pvalues)) %>% # keep VIPER pvalue
        dplyr::select(gene_name, NES, pvalues)

      # Override VIPER output
      output$sign <- viper_on_prot %>% dplyr::filter(pvalues < 0.05)
      output$all <- viper_on_prot %>% dplyr::filter(pvalues > 0.05)
    }
  }

  if(hypergeom_corr == TRUE){
    message('Starting hypergeometric test correction')
    output <- weight_viper_score(run_hypergeometric_test(omic_data,output$sign,
                                                         analysis = analysis,
                                                         organism = organism,
                                                         collectri = collectri,
                                                         integrated_regulons = integrated_regulons))

    output_uniprot <- convert_gene_name_in_uniprotid(output, organism)
  }else{
    output_uniprot <- convert_gene_name_in_uniprotid(output$sign, organism)
  }

  if(GO_annotation == TRUE){
    message('GO molecular function annotation')

    output_uni_mf <- molecular_function_annotation(output_uniprot, organism)
    output_final <- filter_VIPER_output(output_uni_mf, analysis) %>% dplyr::distinct()
    #footprint_output <- output_final
    return(output_final)
  }
  return(output_uniprot)
}
