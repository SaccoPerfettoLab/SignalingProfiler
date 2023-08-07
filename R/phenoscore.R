
#' compute_phenoscore
#'
#' @param sp_output Signaling Profiler output list
#' @param desired_phenotypes NULL if you want all phenotypes, vector of phenotypes name in uppercase
#' @param compact TRUE if you want a compact result, FALSE otherwise
#' @param max_path_length integer, from 1 to 4, path length of
#' @param remove_cascade default TRUE
#' @param use_carnival_activity Boolean value, TRUE uses all proteins in network,
#' FALSE only experimentally derived proteins
#'
#' @return phenoSCORE table
#'
#' @examples
compute_phenoscore <- function(sp_output,
                               desired_phenotypes,
                               compact = TRUE,
                               max_path_length = 4,
                               remove_cascade = TRUE,
                               use_carnival_activity = FALSE){

  # sp_output <- CML_pheno
  # desired_phenotypes <- 'PROLIFERATION'
  # max_path_length = 4
  # remove_cascade = TRUE
  # use_carnival_activity = TRUE
#
#   sp_output <- CML_pheno_CTRL
#   desired_phenotypes = NULL
#   max_path_length = 4
#   remove_cascade = TRUE
#   compact = FALSE
#   use_carnival_activity = TRUE

  # filter out discordant nodes and the ones not coming from experimental data
  nodes_to_consider <- sp_output$nodes_df %>%
    #dplyr::filter(!is.na(final_score) & !discordant) %>%
    dplyr::mutate_at('gene_name', toupper)

  if(use_carnival_activity == TRUE){
    nodes_to_consider$final_score <- NULL
    nodes_to_consider <- nodes_to_consider %>%
      dplyr::rename(final_score = carnival_activity) %>%
      dplyr::filter(!(final_score == 0 & discordant == FALSE))
  }

  if(!is.null(desired_phenotypes)){
    paths_to_pheno <- paths_to_pheno %>%
      dplyr::filter(EndPathways %in% desired_phenotypes)
    if(nrow(paths_to_pheno) == 0){
      stop('Your phenotype cannot be connected to the proteins!')
    }
  }

  prot_to_pheno <- dplyr::inner_join(nodes_to_consider,
                                     paths_to_pheno %>% dplyr::filter(Path_Length <= max_path_length),
                                     by = c('gene_name' = 'QueryNode')) %>%
    dplyr::mutate(signed_exp_score = final_score * Final_Effect)


  if(remove_cascade == TRUE){
    sp_graph <- sp_output$igraph_network
    multiple_regulated_phen <- prot_to_pheno %>%
      dplyr::count(EndPathways) %>%
      dplyr::filter(n>1)

    # select phenotypes
    phenotypes <- multiple_regulated_phen$EndPathways

    # initialize empty list
    to_remove_list <- list()

    for(i_phen in c(1:length(phenotypes))){
      phenotype = phenotypes[i_phen]

      prot_to_phenotype <- prot_to_pheno %>%
        dplyr::filter(EndPathways == phenotype)

      combinatios <- combn(prot_to_phenotype$name, 2)
      dist_count <- matrix(nrow = 3, ncol = ncol(combinatios))

      # forward run
      for(i in ncol(combinatios)){
        #i = 1
        dist_count[1,i] <- unlist(igraph::distances(sp_graph, combinatios[1,i], to = combinatios[2,i], mode = "out"))
        if(dist_count[1,i] == Inf){
          dist_count[2,i] <- NA
        }else{
          dist_count[2,i] <- combinatios[2,i]
        }
      }
      # reverse run
      for(i in ncol(combinatios)){
        #i = 1
        dist_count[1,i] <- unlist(igraph::distances(sp_graph, combinatios[2,i], to = combinatios[1,i], mode = "out"))
        if(dist_count[1,i] == Inf){
          dist_count[3,i] <- NA
        }else{
          dist_count[3,i] <- combinatios[1,i]
        }
      }
      output <- dist_count[!(is.na(as.vector(dist_count)) | as.vector(dist_count) == Inf)]

      to_remove_list[[i_phen]] <- output
    }

    names(to_remove_list) <- phenotypes

    to_remove_list_clean <- to_remove_list[unlist(lapply(to_remove_list, function(x){!is.null(x)}))]

    to_remove_df <- tibble::tibble(phenotype = names(to_remove_list_clean),
                                   proteins = unlist(lapply(to_remove_list_clean, function(x){paste0(x, collapse = ',')})))
    to_remove_df <- to_remove_df %>% tidyr::separate_rows(proteins)

    prot_to_pheno <- dplyr::anti_join(prot_to_pheno,
                                            to_remove_df,
                                            by = c('EndPathways' = 'phenotype',
                                                   'name' = 'proteins'))
  }

  # compute phenotype score
  phenotype_score <- prot_to_pheno %>%
    dplyr::group_by(EndPathways) %>%
    dplyr::summarise(phenoSCORE = sum(signed_exp_score)/dplyr::n())

  extended_result <- dplyr::inner_join(prot_to_pheno, phenotype_score, by = 'EndPathways') %>%
    dplyr::select(name, gene_name, mf, final_score, EndPathways, phenoSCORE, Path_String, Path_Length, Final_Effect)

  compact_result <- extended_result %>%
    dplyr::group_by(EndPathways, phenoSCORE) %>% dplyr::summarise(Regulators = paste0(unique(gene_name), collapse = ';'),
                                                                  Regulators_Paths = paste0(Path_String, collapse = ';'),
                                                                  N_Paths = dplyr::n()) %>%
    dplyr::ungroup()

  if(compact == TRUE){
    return(compact_result)
  }else{
    return(extended_result)
  }
}

#' Title
#'
#' @param proteins_df dataframe containing up-reg or down-reg proteins
#' @param desired_phenotypes vector of interesting phenotypes
#' @param pvalue_threshold default 0.05
#' @param path_length default 4 (1 to 4 allowed)
#' @param zscore_threshold default -1.96
#' @param stat 'mean' or 'median'
#'
#' @return a list of a plot and result table
#'
#' @examples
phenoscore_computation_v1 <- function(proteins_df,
                                      desired_phenotypes,
                                      pvalue_threshold = 0.05,
                                      path_length = 4,
                                      stat = 'mean',
                                      zscore_threshold = -1.96){
  # proteins_df =  eval(parse(text = condition))
  # desired_phenotypes = desired_phenotypes
  # path_length <- 4
  # # zscore_threshold= -2.58 ###(3sd)
  # zscore_threshold= -1.96 ###(2sd)
  # pvalue_threshold = 0.05
  #

  # condition <- 'Ima_df_exp'
  # proteins_df <- eval(parse(text = condition))
  # desired_phenotypes = desired_phenotypes
  # path_length <- 4
  # # zscore_threshold= -2.58 ###(3sd)
  # zscore_threshold= -1.96 ###(2sd)
  # pvalue_threshold = 0.05



  ##############################################################################
  # PARAMETERS INPUT CHECK #
  ##############################################################################

  message('** PHENOSCORE ANALYSIS **')
  # check for path_length parameter
  if(path_length%%1 == 0){
    if(path_length <= 1 | path_length > 4){
      stop('path_length should be 2,3 or 4!')
    }
  }else{
    stop('Please provide an integer value between 1 and 4 for path_length!')
  }

  # check the zscore_threshold parameter
  if(!is.numeric(zscore_threshold)){
    stop('zscore_threshold should be numeric; choose between 1sd (-1.96) or 2sd (-2.58)')
  }

  # check pvalue_threshold
  if(!is.numeric(pvalue_threshold)){
    stop('pvalue_threshold should be numeric')
  }

  ##############################################################################
  # BUILD SIGNIFICANT PATHS TABLE #
  ##############################################################################
  message('Building significant paths to phenotypes table')
  ## filter by path_length (user defined)
  phenoscore_distances_table %>%
    dplyr::filter(Path_Length <= path_length) -> dist.glob.filtered

  output.table2 <- dist.glob.filtered[,1:8]

  output.table2 %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways,
                    .keep_all = TRUE) -> output.table2

  ##create effect name
  output.table2$Effect <- '-'
  output.table2$Effect[ output.table2$Final_Effect == 1] <- 'up-regulates'
  output.table2$Effect[ output.table2$Final_Effect == 0] <- '-'
  output.table2$Effect[ output.table2$Final_Effect == -1] <- 'down-regulates'

  ##create a Pathscore distrubution by pathway e calculate mean/median and standard deviation (sd)
  output.table2%>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(n = n(),
                     mean = mean(Path_Score),
                     median = median(Path_Score),
                     sd=sd(Path_Score)) -> summary.table.mean.sd.path.score ### here we can select between mean and median

  ##merge the path table with mean/median and standard deviation (sd)
  merge(x = output.table2,
        y = summary.table.mean.sd.path.score,
        by = c('EndPathways', 'Effect'), all.x = TRUE) -> output.table.mean.sd.path.score

  ##calculate z-score
  if(stat == 'mean'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$mean)/output.table.mean.sd.path.score$sd
  }else if(stat == 'median'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$median)/output.table.mean.sd.path.score$sd
  }else{
    stop('Please insert a valid stat: mean or median')
  }

  ## filter by z-score threshold (user defined)
  output.table.mean.sd.path.score %>%
    dplyr::filter(zscore <= zscore_threshold)  -> output.table.mean.sd.path.score.filtered

  ###remove duplicates
  output.table.mean.sd.path.score.filtered %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways, .keep_all = TRUE) -> output.table6

  ## analyse all the datasets:
  output.table6%>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(total = dplyr::n()) -> effect.total

  ##############################################################################
  # PHENOSCORE ANALYSIS #
  ##############################################################################


  ##create an empty dataframe to store results
  output.table.final.global <- data.frame(matrix (ncol = 12, nrow = 0))
  colnames(output.table.final.global) <- c( "EndPathways",
                                            "Effect",
                                            "hits_INPUT",
                                            "total_paths_impacting_phenotype",
                                            "Fraction",
                                            "mean(Frac_rand)",
                                            "sd",
                                            "min",
                                            "max",
                                            "t",
                                            "pvalue_raw")



  # parse proteins df
  proteins_vector <- proteins_df$gene_name
  INPUT_dataset <- data.frame(proteins_vector)
  colnames(INPUT_dataset)[1] <- "QueryNode"
  dim_datset = dim(unique(INPUT_dataset))[1]

  # ## start analysis to search for significant paths to up/down regulations of pathways/RLE
  output.table.final <- data.frame(matrix (ncol = 4,
                                           nrow = 0))
  colnames(output.table.final) <- c( "EndPathways",
                                     "Effect",
                                     "hits",
                                     "total")

  #### randomization

  message('Randomization')
  options(dplyr.summarise.inform = FALSE)
  for (i in c(1:n_random)){
    signor.random <- sample(background_phenoscore)[1:dim_datset]
    output.table6%>%
      dplyr::filter(QueryNode %in% signor.random) %>%
      dplyr::group_by(EndPathways, Effect)%>%
      dplyr::summarise(hits = dplyr::n()) -> output.table.randomized

    randomized <- merge(x= output.table.randomized,
                        y= effect.total,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)
    output.table.final<- rbind(output.table.final,randomized)
  }

  output.table.final$Frac_rand <- output.table.final$hits/output.table.final$total

  output.table.final %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(mean(Frac_rand),
                     sd=sd(Frac_rand),
                     min=min(Frac_rand),
                     max=max(Frac_rand)) -> aaa

  #### comparison with the input
  output.table6%>%
    dplyr::filter(QueryNode %in% INPUT_dataset$QueryNode) %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(hits = n()) -> output.table.INPUT

  output.table.INPUT <- merge(x= output.table.INPUT,
                              y= effect.total,
                              by = c('EndPathways', 'Effect'),
                              all.x = T)
  colnames(output.table.INPUT) <- c("EndPathways", "Effect", "hits_INPUT",
                                    "total_paths_impacting_phenotype" )
  output.table.INPUT$Fraction <- round(output.table.INPUT$hits_INPUT/output.table.INPUT$total_paths_impacting_phenotype, 2)
  output.table <- merge(x= output.table.INPUT,
                        y= aaa,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)

  #### t-score computation

  message('T-score computation between random and input list')
  output.table$t <- (output.table$Fraction-output.table$`mean(Frac_rand)`)/output.table$sd
  #
  output.table$pvalue_raw <- pnorm(output.table$t, lower.tail = FALSE)
  output.table.final.global <- rbind(output.table.final.global,output.table)
  output.table.final.global$pvalue <- p.adjust(output.table.final.global$pvalue_raw,
                                               method = 'BH',
                                               n = length(output.table.final.global$EndPathways))

  #filter output according to user thresholds
  output.table.final.global$sign <- ''
  output.table.final.global$sign[output.table.final.global$pvalue < 0.05]<- '*'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.001]<- '**'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.0001]<- '***'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.00001]<- '****'
  output.table.final.global$sign[output.table.final.global$pvalue >1 ]<- '!'

  output.table.final.global %>%
    dplyr::filter( pvalue < pvalue_threshold &
                     Effect != '-' &
                     EndPathways %in% desired_phenotypes)-> results.table

  if(nrow(results.table) == 0){
    message('No significant phenotype found! Try to change max_length parameters')
    return(NULL)
  }

  results.table$EndPathways <- sapply(strsplit(as.character(results.table$EndPathways ),
                                               "="),"[", 1)

  # generate plot
  message('Plot generation')
  ### display interactions from INPUT genes to phenotypes
  id_table <- data.frame(matrix (ncol = 1, nrow = 0))
  colnames(id_table) <- c("signor_ids")


  output.table6 <- output.table6 %>%
    dplyr::mutate_at(c('EndPathways', 'QueryNode'), as.character)

  # number of paths of each protein


  output.table6%>%
    dplyr::filter(QueryNode %in% proteins_vector)-> INPUT.Phen.Paths

  INPUT.Phen.Paths %>%
    dplyr::group_by(QueryNode, EndPathways, Effect) %>%
    dplyr::summarise(n = n()) -> count.query

  ProteinsPaths <- count.query %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(regulators = paste0((QueryNode),collapse = ';'),
                     node_idx = paste0((round(n/sum(n),2)), collapse = ';')) %>%
    dplyr::mutate_at('EndPathways', tolower) %>%
    dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))

  results.table
  # prepare data for plotting
  results.table$Log10_p_value_plot <- -log10(results.table$pvalue)

  results.table%>%
    dplyr::arrange(as.integer(Log10_p_value_plot))%>%
    dplyr::arrange(Effect)%>%
    dplyr::mutate(EndPathways = tolower(gsub('_', ' ', EndPathways)))%>%
    tibble::as_tibble()->results.table

  ggplot2::ggplot(results.table, ggplot2::aes(x = factor(EndPathways),
                                              y = Log10_p_value_plot,
                                              fill = Effect)) +
    ggplot2::geom_col(position = "dodge")+
    ggplot2::ggtitle(paste0('Significantly close phenotypes')) +
    ggplot2::ylab("-Log10(p-value)") + ggplot2::xlab('phenotype') +
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values = c("#9DD2C7","#0C838C"))+
    ggplot2::geom_text(ggplot2::aes(label = hits_INPUT),
                       colour = "black", size = 2, hjust= 1,
                       position = ggplot2::position_dodge(.9))+
    ggplot2::coord_flip()->barplot_phenotypes

  # join paths with regulators
  results.table_reg <- dplyr::left_join(results.table,
                                        ProteinsPaths,
                                        by = c('EndPathways', 'Effect'))


  return(list(barplot = barplot_phenotypes,
              table = results.table_reg))
}


#' Title
#'
#' @param proteins_df dataframe containing up-reg or down-reg proteins
#' @param desired_phenotypes vector of interesting phenotypes
#' @param pvalue_threshold default 0.05
#' @param path_length default 4 (1 to 4 allowed)
#' @param zscore_threshold default -1.96
#' @param sp_graph graph of signaling profiler output
#' @param remove_cascade Boolean value, TRUE or FALSE
#' @param n_random number of randomization run, 1000
#' @param stat 'mean' or 'median'
#' @param use_carnival_activity
#'
#' @return a list of a plot and result table
#'
#' @examples
phenoscore_computation_v2 <- function(proteins_df,
                                      desired_phenotypes,
                                      sp_graph,
                                      remove_cascade = TRUE,
                                      n_random = 1000,
                                      stat = 'mean',
                                      pvalue_threshold = 0.05,
                                      path_length = 4,
                                      zscore_threshold = -1.96){

  # i = 1
  # condition <- 'JMD_df'
  #
  # proteins_df <- eval(parse(text = condition))
  # sp_graph <- eval(parse(text = net_list[i]))
  # V(sp_graph)$gene_name <-   toupper(V(sp_graph)$gene_name)
  # remove_cascade = TRUE
  # pvalue_threshold = 0.05
  # path_length = 4
  # zscore_threshold = -1.96

  # CML
  # proteins_df
  # desired_phenotypes = desired_phenotypes
  # sp_graph = Ctrl$igraph_network
  # remove_cascade = TRUE
  # path_length = 4
  # pvalue_threshold = 0.05

  ##############################################################################
  # PARAMETERS INPUT CHECK #
  ##############################################################################

  message('** PHENOSCORE ANALYSIS **')
  # check for path_length parameter
  if(path_length%%1 == 0){
    if(path_length <= 1 | path_length > 4){
      stop('path_length should be 2,3 or 4!')
    }
  }else{
    stop('Please provide an integer value between 1 and 4 for path_length!')
  }

  # check the zscore_threshold parameter
  if(!is.numeric(zscore_threshold)){
    stop('zscore_threshold should be numeric; choose between 1sd (-1.96) or 2sd (-2.58)')
  }

  # check pvalue_threshold
  if(!is.numeric(pvalue_threshold)){
    stop('pvalue_threshold should be numeric')
  }

  ##############################################################################
  # BUILD SIGNIFICANT PATHS TABLE #
  ##############################################################################
  message('Building significant paths to phenotypes table')
  ## filter by path_length (user defined)
  phenoscore_distances_table %>%
    dplyr::filter(Path_Length <= path_length) -> dist.glob.filtered

  output.table2 <- dist.glob.filtered[,1:8]

  output.table2 %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways,
                    .keep_all = TRUE) -> output.table2

  ##create effect name
  output.table2$Effect <- '-'
  output.table2$Effect[ output.table2$Final_Effect == 1] <- 'up-regulates'
  output.table2$Effect[ output.table2$Final_Effect == 0] <- '-'
  output.table2$Effect[ output.table2$Final_Effect == -1] <- 'down-regulates'

  ##create a Pathscore distrubution by pathway e calculate mean/median and standard deviation (sd)
  output.table2%>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(n = n(),
                     mean = mean(Path_Score),
                     median = median(Path_Score),
                     sd=sd(Path_Score)) -> summary.table.mean.sd.path.score ### here we can select between mean and median

  ##merge the path table with mean/median and standard deviation (sd)
  merge(x = output.table2,
        y = summary.table.mean.sd.path.score,
        by = c('EndPathways', 'Effect'), all.x = TRUE) -> output.table.mean.sd.path.score

  ##calculate z-score
  if(stat == 'mean'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$mean)/output.table.mean.sd.path.score$sd
  }else if(stat == 'median'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$median)/output.table.mean.sd.path.score$sd
  }else{
    stop('Please insert a valid stat: mean or median')
  }

  ## filter by z-score threshold (user defined)
  output.table.mean.sd.path.score %>%
    dplyr::filter(zscore <= zscore_threshold)  -> output.table.mean.sd.path.score.filtered

  ###remove duplicates
  output.table.mean.sd.path.score.filtered %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways, .keep_all = TRUE) -> output.table6

  ## analyse all the datasets:
  output.table6%>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(total = dplyr::n()) -> effect.total

  ##############################################################################
  # PHENOSCORE ANALYSIS #
  ##############################################################################


  ##create an empty dataframe to store results
  output.table.final.global <- data.frame(matrix (ncol = 12, nrow = 0))
  colnames(output.table.final.global) <- c( "EndPathways",
                                            "Effect",
                                            "hits_INPUT",
                                            "total_paths_impacting_phenotype",
                                            "Fraction",
                                            "mean(Frac_rand)",
                                            "sd",
                                            "min",
                                            "max",
                                            "t",
                                            "pvalue_raw")

  # parse proteins df
  proteins_vector <- proteins_df$gene_name
  INPUT_dataset <- data.frame(proteins_vector)
  colnames(INPUT_dataset)[1] <- "QueryNode"
  dim_datset = dim(unique(INPUT_dataset))[1]

  # ## start analysis to search for significant paths to up/down regulations of pathways/RLE
  output.table.final <- data.frame(matrix (ncol = 4,
                                           nrow = 0))
  colnames(output.table.final) <- c( "EndPathways",
                                     "Effect",
                                     "hits",
                                     "total")

  #### randomization
  message('Randomization')
  options(dplyr.summarise.inform = FALSE)
  for (i in c(1:1000)){
    signor.random <- sample(background_phenoscore)[1:dim_datset]
    output.table6%>%
      dplyr::filter(QueryNode %in% signor.random) %>%
      dplyr::group_by(EndPathways, Effect)%>%
      dplyr::summarise(hits = dplyr::n()) -> output.table.randomized

    randomized <- merge(x= output.table.randomized,
                        y= effect.total,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)
    output.table.final<- rbind(output.table.final,randomized)
  }

  output.table.final$Frac_rand <- output.table.final$hits/output.table.final$total

  output.table.final %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(mean(Frac_rand),
                     sd=sd(Frac_rand),
                     min=min(Frac_rand),
                     max=max(Frac_rand)) -> aaa

  #### comparison with the input
  output.table6%>%
    dplyr::filter(QueryNode %in% INPUT_dataset$QueryNode) %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(hits = n()) -> output.table.INPUT

  output.table.INPUT <- merge(x= output.table.INPUT,
                              y= effect.total,
                              by = c('EndPathways', 'Effect'),
                              all.x = T)
  colnames(output.table.INPUT) <- c("EndPathways", "Effect", "hits_INPUT",
                                    "total_paths_impacting_phenotype" )
  output.table.INPUT$Fraction <- round(output.table.INPUT$hits_INPUT/output.table.INPUT$total_paths_impacting_phenotype, 2)
  output.table <- merge(x= output.table.INPUT,
                        y= aaa,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)

  #### t-score computation
  message('T-score computation between random and input list')
  output.table$t <- (output.table$Fraction-output.table$`mean(Frac_rand)`)/output.table$sd
  #
  output.table$pvalue_raw <- pnorm(output.table$t, lower.tail = FALSE)
  output.table.final.global <- rbind(output.table.final.global,output.table)
  output.table.final.global$pvalue <- p.adjust(output.table.final.global$pvalue_raw,
                                               method = 'BH',
                                               n = length(output.table.final.global$EndPathways))

  #filter output according to user thresholds
  output.table.final.global$sign <- ''
  output.table.final.global$sign[output.table.final.global$pvalue < 0.05]<- '*'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.001]<- '**'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.0001]<- '***'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.00001]<- '****'
  output.table.final.global$sign[output.table.final.global$pvalue >1 ]<- '!'

  output.table.final.global %>%
    dplyr::filter( pvalue < pvalue_threshold &
                     Effect != '-' &
                     EndPathways %in% desired_phenotypes)-> results.table

  if(nrow(results.table) == 0){
    message('No significant phenotype found! Try to change max_length parameters')
    return(NULL)
  }

  results.table$EndPathways <- sapply(strsplit(as.character(results.table$EndPathways ),
                                               "="),"[", 1)


  ### display interactions from INPUT genes to phenotypes
  id_table <- data.frame(matrix (ncol = 1, nrow = 0))
  colnames(id_table) <- c("signor_ids")


  output.table6 <- output.table6 %>%
    dplyr::mutate_at(c('EndPathways', 'QueryNode'), as.character)

  # number of paths of each protein

  # output.table6 <- output.table6 %>%
  #   dplyr::mutate_at('EndPathways', tolower) %>%
  #   dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))


  output.table7 <- inner_join(output.table6, results.table, by = c('EndPathways', 'Effect'))
  output.table7 %>%
    dplyr::filter(QueryNode %in% proteins_vector)-> INPUT.Phen.Paths

  # remove cascade
  message('Removing signaling cascade regulators')

  if(remove_cascade == TRUE){

    # filter phenotypes with effect multiply regulated

    INPUT.Phen.Paths_clean <- INPUT.Phen.Paths %>%
      dplyr::filter(Effect != '-')

    INPUT.Phen.Paths_clean <- INPUT.Phen.Paths_clean %>%
      mutate(key = ifelse(INPUT.Phen.Paths_clean$Effect == 'up-regulates',
                          paste0(INPUT.Phen.Paths_clean$EndPathways, '-up'),
                          paste0(INPUT.Phen.Paths_clean$EndPathways, '-down')))

    # keep only true multiple regulators
    prot_to_phenotypes <- INPUT.Phen.Paths_clean %>%
      dplyr::select(key, QueryNode) %>%
      dplyr::distinct()

    multiple_regulated_phen <- prot_to_phenotypes %>%
      dplyr::count(key) %>%
      dplyr::filter(n>1)

    # create a table with proteins on phenotypes multiply regulated
    prot_to_phenotypes_multiple <- prot_to_phenotypes %>%
      filter(key %in% multiple_regulated_phen$key)

    # select phenotypes
    phenotypes <- unique(prot_to_phenotypes_multiple$key)

    # initialize empty list
    to_remove_list <- list()

    for(i_phen in c(1:length(phenotypes))){
      # print('i phenotype')
      # print(i_phen)
      phenotype = phenotypes[i_phen]

      prot_to_phenotype <- prot_to_phenotypes_multiple %>%
        dplyr::filter(key == phenotype)

      combinatios <- combn(unique(prot_to_phenotype$QueryNode), 2)
      dist_count <- matrix(nrow = 3, ncol = ncol(combinatios))

      # forward run
      for(i in c(1:ncol(combinatios))){
        # print('i forward:')
        # print(i)

        # handling NA in gene_name column

        true_table_from <- igraph::V(sp_graph)$gene_name == combinatios[1,i]
        true_table_to <- igraph::V(sp_graph)$gene_name == combinatios[2,i]

        true_table_from[is.na(true_table_from)] <- FALSE
        true_table_to[is.na(true_table_to)] <- FALSE

        if(sum(true_table_from)==0 |  sum(true_table_to) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network')
          dist_count[1,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from]
          to = igraph::V(sp_graph)[true_table_to]
          dist_count[1,i] <- unlist(igraph::distances(sp_graph, v = from,
                                                      to = to,
                                                      mode = "out"))
          if(dist_count[1,i] == Inf){
            dist_count[2,i] <- NA
          }else{
            dist_count[2,i] <- combinatios[2,i]
          }
        }
      }
      # reverse run
      for(i in c(1:ncol(combinatios))){
        # print('i reverse')
        # print(i)

        true_table_from_rev <- igraph::V(sp_graph)$gene_name == combinatios[2,i]
        true_table_from_rev[is.na(true_table_from_rev)] <- FALSE

        true_table_to_rev <- igraph::V(sp_graph)$gene_name == combinatios[1,i]
        true_table_to_rev[is.na(true_table_to_rev)] <- FALSE

        if(sum(true_table_from_rev) == 0 | sum(true_table_to_rev) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network')
          dist_count[3,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from_rev]
          to = igraph::V(sp_graph)[true_table_to_rev]

          dist_count[1,i] <- unlist(igraph::distances(sp_graph, from, to = to, mode = "out"))

          if(dist_count[1,i] == Inf){
            dist_count[3,i] <- NA
          }else{
            dist_count[3,i] <- combinatios[1,i]
          }
        }
        output <- dist_count[!(is.na(as.vector(dist_count)) | as.vector(dist_count) == Inf)]
      }
      to_remove_list[[i_phen]] <- output
    }

    names(to_remove_list) <- phenotypes

    to_remove_list_clean <- to_remove_list[unlist(lapply(to_remove_list, function(x){!is.null(x)}))]

    to_remove_df <- tibble::tibble(phenotype = names(to_remove_list_clean),
                                   proteins = unlist(lapply(to_remove_list_clean, function(x){paste0(x, collapse = ',')})))
    to_remove_df <- to_remove_df %>% tidyr::separate_rows(proteins)

    #
    prot_to_phenotypes
    prot_to_phenotypes_clean <- dplyr::anti_join(prot_to_phenotypes,
                                                 to_remove_df,
                                                 by = c('key' = 'phenotype',
                                                        'QueryNode' = 'proteins'))

    INPUT.Phen.Paths <- inner_join(INPUT.Phen.Paths_clean,
                                   prot_to_phenotypes_clean,
                                   by = c('QueryNode', 'key'))

  }

  INPUT.Phen.Paths %>%
    dplyr::group_by(QueryNode, EndPathways, Effect) %>%
    dplyr::summarise(n = n()) -> count.query

  ProteinsPaths <- count.query %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(regulators = paste0((QueryNode),collapse = ';'),
                     node_idx = paste0((round(n/sum(n),2)), collapse = ';')) %>%
    dplyr::mutate_at('EndPathways', tolower) %>%
    dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))

  results.table
  # prepare data for plotting
  results.table$Log10_p_value_plot <- -log10(results.table$pvalue)

  results.table%>%
    dplyr::arrange(as.integer(Log10_p_value_plot))%>%
    dplyr::arrange(Effect)%>%
    dplyr::mutate(EndPathways = tolower(gsub('_', ' ', EndPathways)))%>%
    tibble::as_tibble()->results.table

  # generate plot
  message('Plot generation')
  ggplot2::ggplot(results.table, ggplot2::aes(x = factor(EndPathways),
                                              y = Log10_p_value_plot,
                                              fill = Effect)) +
    ggplot2::geom_col(position = "dodge")+
    ggplot2::ggtitle(paste0('Significantly close phenotypes')) +
    ggplot2::ylab("-Log10(p-value)") + ggplot2::xlab('phenotype') +
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values = c("#9DD2C7","#0C838C"))+
    ggplot2::geom_text(ggplot2::aes(label = hits_INPUT),
                       colour = "black", size = 2, hjust= 1,
                       position = ggplot2::position_dodge(.9))+
    ggplot2::coord_flip()->barplot_phenotypes

  # join paths with regulators
  results.table_reg <- dplyr::left_join(results.table,
                                        ProteinsPaths,
                                        by = c('EndPathways', 'Effect'))


  return(list(barplot = barplot_phenotypes,
              table = results.table_reg))
}

#' Title
#'
#' @param proteins_df dataframe containing up-reg or down-reg proteins
#' @param desired_phenotypes vector of interesting phenotypes
#' @param pvalue_threshold default 0.05
#' @param path_length default 4 (1 to 4 allowed)
#' @param zscore_threshold default -1.96
#' @param sp_graph graph of signaling profiler output
#' @param remove_cascade Boolean value, TRUE or FALSE
#' @param n_random number of randomization, 1000
#' @param stat 'mean' or 'median'
#' @param node_idx boolean
#' @param use_carnival_activity Boolean, TRUE or FALSE
#' @param phenoscore_distances_table df output of phenoscore network processing, in NULL take all
#'
#' @return a list of a plot and result table
#' @export
#'
#' @examples
phenoscore_computation <- function(proteins_df,
                                   desired_phenotypes = NULL,
                                   sp_graph,
                                   pheno_distances_table = NULL,
                                   path_length = 4,
                                   stat = 'mean',
                                   zscore_threshold = -1.96,
                                   n_random = 1000,
                                   pvalue_threshold = 0.05,
                                   remove_cascade = TRUE,
                                   node_idx = FALSE,
                                   use_carnival_activity = FALSE){

  # # TEST PARAMETERS #
  # proteins_df = Ctrl_df_exp
  # desired_phenotypes = desired_phenotypes
  # pheno_distances_table = pheno_table_distances
  # sp_graph = Ctrl$igraph_network
  # path_length = 4
  # stat = 'median'
  # zscore_threshold = -1.96
  # n_random = 1000
  # pvalue_threshold = 0.05
  # remove_cascade = TRUE
  # node_idx = TRUE
  # use_carnival_activity = FALSE


  # Distinguish between mouse and human
  if(sum(proteins_df$gene_name == stringr::str_to_title(proteins_df$gene_name)) == nrow(proteins_df)){
    flag = 'mouse'
    # if it is mouse organism, convert in uppercase
    proteins_df <- proteins_df %>%
      mutate(gene_name = stringr::str_to_upper(gene_name))
    igraph::V(sp_graph)$gene_name <- str_to_upper(igraph::V(sp_graph)$gene_name)
  }else{
    flag = 'human'
  }

  ##############################################################################
  # PARAMETERS INPUT CHECK #
  ##############################################################################

  message('** PHENOSCORE ANALYSIS **')
  # check for path_length parameter
  if(path_length%%1 == 0){
    if(path_length <= 1 | path_length > 4){
      stop('path_length should be 2,3 or 4!')
    }
  }else{
    stop('Please provide an integer value between 2 and 4 for path_length!')
  }

  # check the zscore_threshold parameter
  if(!is.numeric(zscore_threshold)){
    stop('zscore_threshold should be numeric; choose between 1sd (-1.96) or 2sd (-2.58)')
  }

  # check pvalue_threshold
  if(!is.numeric(pvalue_threshold)){
    stop('pvalue_threshold should be numeric')
  }

  ##############################################################################
  # BUILD SIGNIFICANT PATHS TABLE #
  ##############################################################################
  message('Building significant paths to phenotypes table')
  ## filter by path_length (user defined)

  # if you dont' want to prefilter, use build-in object
  if(is.null(pheno_distances_table)){
    pheno_distances_table <- phenoscore_distances_table
  }

  pheno_distances_table %>%
    dplyr::filter(Path_Length <= path_length) -> dist.glob.filtered

  output.table2 <- dist.glob.filtered[,1:8]

  output.table2 %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways,
                    .keep_all = TRUE) -> output.table2

  ##create effect name
  output.table2$Effect <- '-'
  output.table2$Effect[ output.table2$Final_Effect == 1] <- 'up-regulates'
  output.table2$Effect[ output.table2$Final_Effect == 0] <- '-'
  output.table2$Effect[ output.table2$Final_Effect == -1] <- 'down-regulates'

  ##create a Pathscore distrubution by pathway e calculate mean/median and standard deviation (sd)
  output.table2%>%
    dplyr::filter(Effect != '-') %>% # remove '-' effect
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(n = n(),
                     mean = mean(Path_Score),
                     median = median(Path_Score),
                     sd=sd(Path_Score)) -> summary.table.mean.sd.path.score ### here we can select between mean and median

  ##merge the path table with mean/median and standard deviation (sd)
  merge(x = output.table2 %>% dplyr::filter(Effect != '-'),
        y = summary.table.mean.sd.path.score,
        by = c('EndPathways', 'Effect'), all.x = TRUE) -> output.table.mean.sd.path.score

  ##calculate z-score
  if(stat == 'mean'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$mean)/output.table.mean.sd.path.score$sd
  }else if(stat == 'median'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$median)/output.table.mean.sd.path.score$sd
  }else{
    stop('Please insert a valid stat: mean or median')
  }

  ## filter by z-score threshold (user defined)
  output.table.mean.sd.path.score %>%
    dplyr::filter(zscore <= zscore_threshold)  -> output.table.mean.sd.path.score.filtered

  ###remove duplicates
  output.table.mean.sd.path.score.filtered %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways, .keep_all = TRUE) -> output.table6

  ## analyse all the datasets:
  output.table6%>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(total = dplyr::n()) -> effect.total

  ##############################################################################
  # PHENOSCORE ANALYSIS #
  ##############################################################################

  ##create an empty dataframe to store results
  output.table.final.global <- data.frame(matrix (ncol = 12, nrow = 0))
  colnames(output.table.final.global) <- c( "EndPathways",
                                            "Effect",
                                            "hits_INPUT",
                                            "total_paths_impacting_phenotype",
                                            "Fraction",
                                            "mean(Frac_rand)",
                                            "sd",
                                            "min",
                                            "max",
                                            "t",
                                            "pvalue_raw")



  # parse proteins df
  proteins_vector <- proteins_df$gene_name
  INPUT_dataset <- data.frame(proteins_vector)
  colnames(INPUT_dataset)[1] <- "QueryNode"
  dim_datset = dim(unique(INPUT_dataset))[1]

  # ## start analysis to search for significant paths to up/down regulations of pathways/RLE
  output.table.final <- data.frame(matrix (ncol = 4,
                                           nrow = 0))
  colnames(output.table.final) <- c( "EndPathways",
                                     "Effect",
                                     "hits",
                                     "total")

  #### randomization
  message('Randomization')
  options(dplyr.summarise.inform = FALSE)
  for (i in c(1:n_random)){
    signor.random <- sample(background_phenoscore)[1:dim_datset]
    output.table6%>%
      dplyr::filter(QueryNode %in% signor.random) %>%
      dplyr::group_by(EndPathways, Effect)%>%
      dplyr::summarise(hits = dplyr::n()) -> output.table.randomized

    randomized <- merge(x= output.table.randomized,
                        y= effect.total,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)
    output.table.final<- rbind(output.table.final,randomized)
  }

  output.table.final$Frac_rand <- output.table.final$hits/output.table.final$total

  output.table.final %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(mean(Frac_rand),
                     sd=sd(Frac_rand),
                     min=min(Frac_rand),
                     max=max(Frac_rand)) -> aaa

  #### comparison with the input
  output.table6%>%
    dplyr::filter(QueryNode %in% INPUT_dataset$QueryNode) %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(hits = n()) -> output.table.INPUT

  output.table.INPUT <- merge(x= output.table.INPUT,
                              y= effect.total,
                              by = c('EndPathways', 'Effect'),
                              all.x = T)
  colnames(output.table.INPUT) <- c("EndPathways", "Effect", "hits_INPUT",
                                    "total_paths_impacting_phenotype" )
  output.table.INPUT$Fraction <- round(output.table.INPUT$hits_INPUT/output.table.INPUT$total_paths_impacting_phenotype, 2)
  output.table <- merge(x= output.table.INPUT,
                        y= aaa,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)

  #### t-score computation
  message('T-score computation between random and input list')
  output.table$t <- (output.table$Fraction-output.table$`mean(Frac_rand)`)/output.table$sd
  #
  output.table$pvalue_raw <- pnorm(output.table$t, lower.tail = FALSE)
  output.table.final.global <- rbind(output.table.final.global,output.table)
  output.table.final.global$pvalue <- p.adjust(output.table.final.global$pvalue_raw,
                                               method = 'BH',
                                               n = length(output.table.final.global$EndPathways))

  #filter output according to user thresholds
  output.table.final.global$sign <- ''
  output.table.final.global$sign[output.table.final.global$pvalue < 0.05]<- '*'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.001]<- '**'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.0001]<- '***'
  output.table.final.global$sign[output.table.final.global$pvalue < 0.00001]<- '****'
  output.table.final.global$sign[output.table.final.global$pvalue >1 ]<- '!'

  if(is.null(desired_phenotypes)){
    output.table.final.global %>%
      dplyr::filter( pvalue < pvalue_threshold &
                       Effect != '-')-> results.table

  }else{
    output.table.final.global %>%
      dplyr::filter( pvalue < pvalue_threshold &
                       Effect != '-' &
                       EndPathways %in% desired_phenotypes)-> results.table
  }


  if(nrow(results.table) == 0){
    message('No significant phenotype found! Try to change max_length parameters')
    return(NULL)
  }

  results.table$EndPathways <- sapply(strsplit(as.character(results.table$EndPathways ),
                                               "="),"[", 1)


  output.table6 <- output.table6 %>%
    dplyr::mutate_at(c('EndPathways', 'QueryNode'), as.character)

  # number of paths of each protein

  # output.table6 <- output.table6 %>%
  #   dplyr::mutate_at('EndPathways', tolower) %>%
  #   dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))


  output.table7 <- inner_join(output.table6, results.table, by = c('EndPathways', 'Effect'))
  output.table7 %>%
    dplyr::filter(QueryNode %in% proteins_vector)-> INPUT.Phen.Paths

  # remove cascade
  message('Removing signaling cascade regulators')

  if(remove_cascade == TRUE){

    # filter phenotypes with effect multiply regulated

    INPUT.Phen.Paths_clean <- INPUT.Phen.Paths %>%
      dplyr::filter(Effect != '-')

    INPUT.Phen.Paths_clean <- INPUT.Phen.Paths_clean %>%
      mutate(key = ifelse(INPUT.Phen.Paths_clean$Effect == 'up-regulates',
                          paste0(INPUT.Phen.Paths_clean$EndPathways, '-up'),
                          paste0(INPUT.Phen.Paths_clean$EndPathways, '-down')))

    # keep only true multiple regulators
    prot_to_phenotypes <- INPUT.Phen.Paths_clean %>%
      dplyr::select(key, QueryNode) %>%
      dplyr::distinct()

    multiple_regulated_phen <- prot_to_phenotypes %>%
      dplyr::count(key) %>%
      dplyr::filter(n>1)

    # create a table with proteins on phenotypes multiply regulated
    prot_to_phenotypes_multiple <- prot_to_phenotypes %>%
      filter(key %in% multiple_regulated_phen$key)

    # select phenotypes
    phenotypes <- unique(prot_to_phenotypes_multiple$key)

    # initialize empty list
    to_remove_list <- list()

    for(i_phen in c(1:length(phenotypes))){
      #i_phen = 1
      # print('i phenotype')
      # print(i_phen)
      phenotype = phenotypes[i_phen]

      prot_to_phenotype <- prot_to_phenotypes_multiple %>%
        dplyr::filter(key == phenotype)

      combinatios <- combn(unique(prot_to_phenotype$QueryNode), 2)
      # create a matrix for the distance and the flag
      dist_count <- matrix(nrow = 3, ncol = ncol(combinatios))

      # forward run: check direction A --> B
      for(i in c(1:ncol(combinatios))){
        #i = 6
        # print('i forward:')
        # print(i)

        # handling NA in gene_name column

        true_table_from <- igraph::V(sp_graph)$gene_name == combinatios[1,i]
        true_table_to <- igraph::V(sp_graph)$gene_name == combinatios[2,i]

        true_table_from[is.na(true_table_from)] <- FALSE
        true_table_to[is.na(true_table_to)] <- FALSE

        if(sum(true_table_from)==0 |  sum(true_table_to) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network')
          dist_count[1,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from]
          to = igraph::V(sp_graph)[true_table_to]
          # compute the distance between two nodes
          # assign to the first row of the matrix the distance
          dist_count[1,i] <- unlist(igraph::distances(sp_graph, v = from,
                                                      to = to,
                                                      mode = "out"))

          # if the distance is infinite, NA node to remove
          if(dist_count[1,i] == Inf){
            dist_count[2,i] <- NA
          }else{ # if the nodes are connected, the second is more downstream and should be removed
            dist_count[2,i] <- combinatios[2,i] # oppposite
          }
        }
      }
      # reverse run: check the direction B --> A
      for(i in c(1:ncol(combinatios))){
        #i = 6
        # print('i reverse')
        # print(i)

        true_table_from_rev <- igraph::V(sp_graph)$gene_name == combinatios[2,i]
        true_table_from_rev[is.na(true_table_from_rev)] <- FALSE

        true_table_to_rev <- igraph::V(sp_graph)$gene_name == combinatios[1,i]
        true_table_to_rev[is.na(true_table_to_rev)] <- FALSE

        if(sum(true_table_from_rev) == 0 | sum(true_table_to_rev) == 0){
          warning('Some nodes aren\'t in the network, please check if you are using the right network')
          dist_count[3,i] <- NA
          next()
        }else{
          from = igraph::V(sp_graph)[true_table_from_rev]
          to = igraph::V(sp_graph)[true_table_to_rev]

          # override the previous distance

          distance_rev <- unlist(igraph::distances(sp_graph, from, to = to, mode = "out"))


          # if it is infinite, not connected no removal
          if(distance_rev == Inf){
            dist_count[3,i] <- NA
          }else{ # if it is finite, remove the first node because it is more downstream
            if(distance_rev < dist_count[1,i]){
              # but if the new distance is shorter than the forward one, remove this
              dist_count[3,i] <- combinatios[1,i] #opp
            }else{
              dist_count[3,i] <- combinatios[2,i] #opp
            }
          }
        }
      }
      output <- dist_count[!(is.na(as.vector(dist_count)) | as.vector(dist_count) == Inf)]

      # if i binary removed all the proteins


      to_remove_list[[i_phen]] <- output
    }

    names(to_remove_list) <- phenotypes

    to_remove_list_clean <- to_remove_list[unlist(lapply(to_remove_list, function(x){!is.null(x)}))]

    to_remove_df <- tibble::tibble(phenotype = names(to_remove_list_clean),
                                   proteins = unlist(lapply(to_remove_list_clean, function(x){paste0(x, collapse = ',')})))
    to_remove_df <- to_remove_df %>% tidyr::separate_rows(proteins)

    #
    prot_to_phenotypes
    prot_to_phenotypes_clean <- dplyr::anti_join(prot_to_phenotypes,
                                                 to_remove_df,
                                                 by = c('key' = 'phenotype',
                                                        'QueryNode' = 'proteins'))

    INPUT.Phen.Paths <- inner_join(INPUT.Phen.Paths_clean,
                                   prot_to_phenotypes_clean,
                                   by = c('QueryNode', 'key'))

  }

  INPUT.Phen.Paths %>%
    dplyr::group_by(QueryNode, EndPathways, Effect) %>%
    dplyr::summarise(n = n()) -> count.query

  # ProteinsPaths <- count.query %>%
  #   dplyr::group_by(EndPathways, Effect) %>%
  #   dplyr::summarise(regulators = paste0((QueryNode),collapse = ';'),
  #                    node_idx = paste0((round(n/sum(n),2)), collapse = ';')) %>%
  #   dplyr::mutate_at('EndPathways', tolower) %>%
  #   dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))

  ProteinsPaths <- count.query %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(regulators = paste0((QueryNode),collapse = ';'),
                     node_idx = paste0(round(n/sum(n),2), collapse = ';')) %>%
    dplyr::mutate_at('EndPathways', tolower) %>%
    dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))

  results.table
  # prepare data for plotting
  results.table$Log10_p_value_plot <- -log10(results.table$pvalue)

  results.table%>%
    dplyr::arrange(as.integer(Log10_p_value_plot))%>%
    dplyr::arrange(Effect)%>%
    dplyr::mutate(EndPathways = tolower(gsub('_', ' ', EndPathways)))%>%
    tibble::as_tibble()->results.table

  # prepare data for the plot
  message('Plot generation')

  # join paths with regulators
  results.table_reg <- dplyr::left_join(results.table,
                                        ProteinsPaths,
                                        by = c('EndPathways', 'Effect'))


  tidyr::separate_rows(results.table_reg,
                       regulators, node_idx, sep = ';') -> a_reg

  a_reg <- a_reg %>%
    dplyr::mutate_at('node_idx', as.numeric)

  only_one_reg_act <- dplyr::inner_join(a_reg,
                                        proteins_df,
                                        by = c('regulators' = 'gene_name'))  %>%
    dplyr::mutate(Sign = ifelse(Effect == 'up-regulates', 1, -1))

  if(node_idx == TRUE){
    if(use_carnival_activity == TRUE){
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Sign*Log10_p_value_plot*node_idx*10))

    }else{
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(final_score*Sign*Log10_p_value_plot*node_idx*10))
    }
  }else{
    if(use_carnival_activity == TRUE){
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Log10_p_value_plot*10))

    }else{
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(final_score*Sign*Log10_p_value_plot))
    }
  }

  phenoscore_df1 <- phenoscore_df %>%
    dplyr::filter(phenoscore != 0)

  if(nrow(phenoscore_df1) != nrow(phenoscore_df)){
    message('Some phenotypes are missing because had just one protein acting
            as both positive and negative regulator')
  }

  phenoscore_df1 <- phenoscore_df1 %>% dplyr::mutate(reg = ifelse(phenoscore < 0, 'down', 'up'))

  ggplot2::ggplot(phenoscore_df1,
                  ggplot2::aes(x = forcats::fct_reorder(EndPathways, phenoscore),
                                              y = phenoscore, fill = reg))+
    ggplot2::geom_bar(stat = 'identity', alpha = 0.8)+
    ggplot2::scale_fill_manual(values = c('#D46A6A', '#407F7F')) +
    ggplot2::ggtitle(paste0('Phenoscore')) +
    ggplot2::ylab("phenotype modulation") +
    ggplot2::xlab("") +
    ggplot2::theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          plot.title=element_text(hjust = 0.5),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text( size = 14, face = 'bold'),
          axis.text.x=element_text(size=10, angle=0, vjust=0.5, hjust=0.5))+
    ggplot2::coord_flip()->barplot_phenotypes

  if(flag == 'mouse'){
    results.table_reg$regulators <- stringr::str_to_title(results.table_reg$regulators)
  }

  return(list(barplot = barplot_phenotypes,
              table_regulators = results.table_reg,
              table_phenotypes = phenoscore_df))
}

# library(tidyverse)
# library(readxl)
# proteomics <- read_excel('../data-raw/livia/sara/data/proteomics_resCTRL_sensCTRL.xlsx')
# phospho <- read_excel('../data-raw/livia/sara/data/phosphoproteomics_resCTRL_sensCTRL.xlsx')

#' Title
#'
#' @param proteomics protemics dataframe processed for SP
#' @param phospho phosphoproteomics dataframe processed for SP
#'
#' @return SIGNOR cleaned according to expressed proteins
#' @export
#'
#' @examples
phenoscore_network_preprocessing <- function(proteomics, phospho, local = FALSE){

  if(local == FALSE){
    path_package <- './'
  }else{
    path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')
  }

  home_dir <- path.expand('~')

  write_tsv(proteomics, paste0(home_dir, '/proteomics.tsv'))
  write_tsv(phospho, paste0(home_dir, '/phosphoproteomics.tsv'))

  #reticulate::use_python("/usr/local/bin/python")
  reticulate::py_config()
  reticulate::py_run_file(paste0(path_package, "inst/python/script.py"))

  # Read Python processed file
  signor_filtered <- readr::read_tsv(paste0(home_dir, '/Global_result_final_table_minimized.txt'))

  file.remove(paste0(home_dir, '/Global_result_final_table_minimized.txt'))
  file.remove(paste0(home_dir, '/proteomics.tsv'))
  file.remove(paste0(home_dir, '/phosphoproteomics.tsv'))

  return(signor_filtered)
}
