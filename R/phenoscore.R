
#' phenoscore_network_preprocessing
#'
#' @param proteomics protemics dataframe processed for SP
#' @param phospho phosphoproteomics dataframe processed for SP
#' @param python_path specify python3 path if reticulate default doesn't work
#'
#' @return SIGNOR cleaned according to expressed proteins
#' @export
#'
#' @examples
phenoscore_network_preprocessing <- function(proteomics, phospho,
                                             local = FALSE,
                                             python_path = NULL){

  if(local == TRUE){
    path_package <- './inst/'
  }else{
    #path_package <- paste0(.libPaths()[1], '/SignalingProfiler/')
    path_package <- paste0(.libPaths(), '/SignalingProfiler/')
  }

  home_dir <- path.expand('~')

  write_tsv(proteomics, paste0(home_dir, '/proteomics.tsv'))
  write_tsv(phospho, paste0(home_dir, '/phosphoproteomics.tsv'))

  #reticulate::use_python("/usr/local/bin/python")
  if(!is.null(python_path)){
    reticulate::use_python(python_path)
  }

  reticulate::py_config()

  # Loop on all lib locations to find script.py,
  # if it doesn't work python3 location is the problem
  for(path in path_package){
    result <- tryCatch({
      reticulate::py_run_file(paste0(path, "/python/script.py"))
    }, error = function(e) {
      #message("An error occurred: ", e$message, 'with path ', path)
      NA  # Return NA if an error occurs
    })
  }

  if(file.exists(paste0(home_dir, '/Global_result_final_table_minimized.txt'))){
    signor_filtered <- readr::read_tsv(paste0(home_dir, '/Global_result_final_table_minimized.txt'))
    file.remove(paste0(home_dir, '/Global_result_final_table_minimized.txt'))
    file.remove(paste0(home_dir, '/proteomics.tsv'))
    file.remove(paste0(home_dir, '/phosphoproteomics.tsv'))
    return(signor_filtered)

  }else{
    stop('Try changing python path with \'python_path\' parameter for reticulate')
  }

}

#' phenoscore_computation
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
#' @param pheno_distances_table table of ProxPath phenotypes distances
#' @param create_pheno_network Boolean, TRUE or FALSE to add phenotypes to the model
#' @param use_carnival_activity Boolean, TRUE or FALSE
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
                                   use_carnival_activity = FALSE,
                                   create_pheno_network = TRUE){

  # # TEST PARAMETERS #

  # top20_phenoscore_anno[1,] -> parameter_row
  # proteins_df = nodes_df
  # desired_phenotypes =  str_replace_all(str_to_upper(gold_stadard_pheno$phenotypes), ' ', '_')
  # pheno_distances_table = pheno_table_distances
  # sp_graph = network_graph
  # path_length = parameter_row$path_length
  # stat = parameter_row$stat
  # zscore_threshold = -1.96
  # n_random = 1000
  # pvalue_threshold = 0.05
  # remove_cascade = parameter_row$remove_cascade
  # node_idx = parameter_row$node_idx
  # use_carnival_activity = parameter_row$use_carnival_activity
  # create_pheno_network = TRUE

  # proteins_df = nodes_df_aa
  # desired_phenotypes = desired_phenotypes
  # pheno_distances_table = pheno_table_distances
  # sp_graph = network_graph
  # path_length = 4
  # stat = 'mean'
  # zscore_threshold = -1.96
  # n_random = 1000
  # pvalue_threshold = 0.05
  # remove_cascade = TRUE
  # node_idx = FALSE
  # use_carnival_activity = FALSE
  # create_pheno_network = TRUE


  # Distinguish between mouse and human
  if(sum(proteins_df$gene_name == stringr::str_to_title(proteins_df$gene_name)) == nrow(proteins_df)){
    flag = 'mouse'
    # if it is mouse organism, convert in uppercase
    proteins_df <- proteins_df %>%
      mutate(gene_name = stringr::str_to_upper(gene_name))
    igraph::V(sp_graph)$name <- str_to_upper(igraph::V(sp_graph)$name)
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

  # if use carnival activity == FALSE, remove proteins without final_score
  if(!use_carnival_activity){
    message('Removing proteins without a final_score value')
    proteins_df <- proteins_df %>%
      dplyr::filter(!is.na(final_score))
  }

  if(is.null(desired_phenotypes)){
    message('No desired_phenotypes using all phenotypes!')
    desired_phenotypes <- unique(phenoscore_distances_table$EndPathways)
  }
  ##############################################################################
  # BUILD SIGNIFICANT PATHS TABLE #
  ##############################################################################
  message('Building significant paths to phenotypes table')
  ## filter by path_length (user defined)

  # if you dont' want to prefilter, use build-in object
  if(is.null(pheno_distances_table)){
    message('No custom pheno_distance_table using whole table without path filtering')
    pheno_distances_table <- phenoscore_distances_table
  }

  # Modify pheno_distance_table replacing all non numeric/characters as '_'
  pheno_distances_table <- pheno_distances_table %>%
    dplyr::mutate(QueryNode = stringr::str_replace_all(QueryNode, "[^[:alnum:]]", '_'))

  pheno_distances_table %>%
    dplyr::filter(Path_Length <= path_length &
                    EndPathways %in% desired_phenotypes) -> dist.glob.filtered

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

  ## Create a Pathscore distrubution by pathway and calculate mean/median and standard deviation (sd)
  output.table2%>%
    dplyr::filter(Effect != '-') %>% # remove '-' effect
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(n = n(),
                     mean = mean(Path_Score),
                     median = median(Path_Score),
                     sd=sd(Path_Score)) -> summary.table.mean.sd.path.score ### here we can select between mean and median

  ## Merge the path table with mean/median and standard deviation (sd)
  merge(x = output.table2 %>% dplyr::filter(Effect != '-'),
        y = summary.table.mean.sd.path.score,
        by = c('EndPathways', 'Effect'), all.x = TRUE) -> output.table.mean.sd.path.score

  ## Calculate z-score
  if(stat == 'mean'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$mean)/output.table.mean.sd.path.score$sd
  }else if(stat == 'median'){
    output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$median)/output.table.mean.sd.path.score$sd
  }else{
    stop('Please insert a valid stat: mean or median')
  }

  ## Filter by z-score threshold (user defined)
  output.table.mean.sd.path.score %>%
    dplyr::filter(zscore <= zscore_threshold)  -> output.table.mean.sd.path.score.filtered

  ### Remove duplicates
  output.table.mean.sd.path.score.filtered %>%
    dplyr::distinct(Path_String, QueryNode,
                    EndNode, EndPathways, .keep_all = TRUE) -> output.table6

  ## Analyse all the dataset:
  output.table6 %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::reframe(total = dplyr::n()) -> effect.total

  ##############################################################################
  # PHENOSCORE ANALYSIS #
  ##############################################################################

  ## Create an empty dataframe to store results
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



  # Parse proteins df
  proteins_vector <- proteins_df$gene_name
  INPUT_dataset <- data.frame(proteins_vector)
  colnames(INPUT_dataset)[1] <- "QueryNode"
  dim_datset = dim(unique(INPUT_dataset))[1]

  # Start analysis to search for significant paths to up/down regulations of pathways/RLE
  output.table.final <- data.frame(matrix (ncol = 4,
                                           nrow = 0))
  colnames(output.table.final) <- c( "EndPathways",
                                     "Effect",
                                     "hits",
                                     "total")

  #### Randomization
  message('Randomization')
  options(dplyr.summarise.inform = FALSE)

  # Create a base df with all 0s if it doesn't find a target
  base_df <- dplyr::bind_rows(tidyr::tibble(EndPathways = desired_phenotypes, Effect = 'down-regulates', hits = 0, total = 0),
                              tidyr::tibble(EndPathways = desired_phenotypes, Effect = 'up-regulates', hits = 0, total = 0))

  for (i in c(1:n_random)){
    signor.random <- sample(background_phenoscore)[1:dim_datset]
    output.table6 %>%
      dplyr::filter(QueryNode %in% signor.random) %>%
      dplyr::group_by(EndPathways, Effect)%>%
      dplyr::summarise(hits = dplyr::n()) -> output.table.randomized

    randomized <- merge(x= output.table.randomized,
                        y= effect.total,
                        by = c('EndPathways', 'Effect'),
                        all.x = T)

    base_df_anti <- anti_join(base_df, randomized, by = c('EndPathways', 'Effect'))

    randomized <- rbind(randomized, base_df_anti)

    output.table.final<- rbind(output.table.final, randomized)
  }

  output.table.final$Frac_rand <- output.table.final$hits/output.table.final$total
  output.table.final$Frac_rand[is.nan(output.table.final$Frac_rand)] <- 0

  output.table.final %>%
    dplyr::group_by(EndPathways, Effect)%>%
    dplyr::summarise(mean(Frac_rand),
                     sd=sd(Frac_rand),
                     min=min(Frac_rand),
                     max=max(Frac_rand)) -> aaa

  #### comparison with the input
  output.table6 %>%
    dplyr::filter(QueryNode %in% INPUT_dataset$QueryNode) %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::reframe(hits = n()) -> output.table.INPUT

  if(nrow(output.table.INPUT) == 0){
    message('No significant path found for your protein list; try to extend it!')
    return(NULL)
  }

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
  output.table$t <- (output.table$Fraction-output.table$`mean(Frac_rand)`)/(output.table$sd + 0.000000000000000001)
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
  output.table.final.global$sign[output.table.final.global$pvalue > 1 ]<- '!'

  if(is.null(desired_phenotypes)){
    output.table.final.global %>%
      dplyr::filter( pvalue < pvalue_threshold &
                       Effect != '-')-> results.table

  }else{
    output.table.final.global %>%
      dplyr::filter( pvalue < pvalue_threshold &
                       Effect != '-') -> results.table
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

  # Remove cascade
  message('Removing signaling cascade regulators')

  # only if there is more than one regulator!
  if(remove_cascade == TRUE & length(unique(INPUT.Phen.Paths$QueryNode)) != 1 ){

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

    # Remove cascade only if multiple_regulated_phen exist
    if(nrow(multiple_regulated_phen) != 0){

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

          true_table_from <- igraph::V(sp_graph)$name == combinatios[1,i]
          true_table_to <- igraph::V(sp_graph)$name == combinatios[2,i]

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

          true_table_from_rev <- igraph::V(sp_graph)$name == combinatios[2,i]
          true_table_from_rev[is.na(true_table_from_rev)] <- FALSE

          true_table_to_rev <- igraph::V(sp_graph)$name == combinatios[1,i]
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

  }

  INPUT.Phen.Paths %>%
    dplyr::group_by(QueryNode, EndPathways, Effect) %>%
    dplyr::summarise(n = n()) -> count.query

  ProteinsPaths <- count.query %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(regulators = paste0((QueryNode),collapse = ';'),
                     node_idx = paste0(round(n/sum(n),2), collapse = ';')) %>%
    dplyr::mutate_at('EndPathways', tolower) %>%
    dplyr::mutate(EndPathways = stringr::str_replace_all(EndPathways, '_', ' '))

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
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Sign*Log10_p_value_plot*node_idx))

    }else{
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(final_score*Sign*Log10_p_value_plot*node_idx))
    }
  }else{
    if(use_carnival_activity == TRUE){
      phenoscore_df <- only_one_reg_act %>%
        dplyr::group_by(EndPathways) %>%
        dplyr::summarise(phenoscore = mean(carnival_activity/100*Sign*Log10_p_value_plot))

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

  color_list <- list('down' = '#407F7F', 'up' = '#D46A6A')

  ggplot2::ggplot(phenoscore_df1,
                  ggplot2::aes(x = forcats::fct_reorder(EndPathways, phenoscore),
                               y = phenoscore, fill = reg))+
    ggplot2::geom_bar(stat = 'identity', alpha = 0.8)+
    ggplot2::scale_fill_manual(values = color_list, labels = names(color_list)) +
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
    igraph::V(sp_graph)$name <- stringr::str_to_title(igraph::V(sp_graph)$name)

  }

  # Create a network with the phenotypes linked
  if(create_pheno_network == TRUE){

    phenotype_regulators <- results.table_reg %>%
      dplyr::select(EndPathways, Effect, regulators) %>%
      tidyr::separate_rows(regulators, sep = ';') %>%
      dplyr::mutate(Effect = ifelse(Effect == 'down-regulates', -1, 1)) %>%
      dplyr::rename('source' = 'regulators' , 'sign' = 'Effect', 'target' = 'EndPathways')

    pheno_nodes <- tidyr::tibble(gene_name =  phenoscore_df$EndPathways,
                                 carnival_activity = NA,
                                 UNIPROT = NA,
                                 mf = 'phenotype',
                                 final_score = phenoscore_df$phenoscore,
                                 method = 'phenoscore',
                                 discordant = FALSE
    )

    pheno_nodes <- pheno_nodes %>%
      dplyr::mutate(carnival_activity = ifelse(final_score < 0, -100, 100),
                    gene_name = stringr::str_to_upper(gene_name))

    igraph::as_data_frame(sp_graph, what = 'vertices') -> nodes_df

    colnames(nodes_df)[1] <- 'gene_name'

    node_df_pheno <- dplyr::bind_rows(nodes_df, pheno_nodes)

    pheno_edges <- phenotype_regulators %>% dplyr::mutate(target = stringr::str_to_upper(target))
    pheno_edges <- tidyr::tibble(target = pheno_edges$target,
                                 sign = as.character(pheno_edges$sign),
                                 source = pheno_edges$source,
                                 carnival_weight = 100,
                                 direct = 'FALSE',
                                 aminoacid = '',
                                 is_quantified = 'FALSE',
                                 is_significant = 'FALSE',
                                 FC = '' )

    # generate edges table

    edges_df <- tidyr::as_tibble(igraph::as_data_frame(sp_graph, what = 'edges'))
    colnames(edges_df)[1:2] <- c('source', 'target')

    edges_df_pheno <- dplyr::bind_rows(edges_df, pheno_edges)

    pheno_graph <- igraph::graph_from_data_frame(edges_df_pheno,
                                                 vertices = node_df_pheno %>% dplyr::distinct() )

    pheno_graph_object <- list(igraph_network = pheno_graph,
                               nodes_df = node_df_pheno,
                               edges_df = edges_df_pheno)

    return(list(barplot = barplot_phenotypes,
                table_regulators = results.table_reg,
                table_phenotypes = phenoscore_df,
                sp_object_phenotypes = pheno_graph_object))

  }

  return(list(barplot = barplot_phenotypes,
              table_regulators = results.table_reg,
              table_phenotypes = phenoscore_df))
}


#' retrieve_functional_circuit
#'
#' @param SP_object sp_object list with network, nodes and edges
#' @param start_nodes character vector of nodes in the model
#' @param phenotype string reporting a phenotype
#' @param k integer, path length
#'
#' @return a subnetwork linking start nodes to phenotypes with a maximum of k length (shorter paths are included)
#' @export
#'
#' @examples
retrieve_functional_circuit <- function(SP_object, start_nodes, phenotype, k) {

  SP_graph <- SP_object$igraph_network

  for(i in c(1:length(start_nodes))){
    start_node <- start_nodes[i]

    # Find all paths of length k between start and end nodes
    all_paths <- igraph::all_simple_paths(SP_graph, from = start_node, to = phenotype, mode = "out", cutoff = k)
    all_paths_nodes <- unique(names(unlist(all_paths)))

    if(length(all_paths_nodes) == 0){
      warning(paste0('No path of length ', k, ' have been found for ', start_node))
    }
    if(i == 1){
      final_nodes <- all_paths_nodes
    }else{
      final_nodes <- c(final_nodes, all_paths_nodes)
    }
  }

  pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)

  igraph::as_data_frame(pheno_circuit, what = c('edges')) -> pheno_circuit_edges
  pheno_circuit_edges <- pheno_circuit_edges %>% dplyr::filter(!to %in% start_nodes)

  pheno_circuit_no_incoming <- igraph::graph_from_data_frame(d =pheno_circuit_edges,
                                                             vertices = igraph::as_data_frame(pheno_circuit, what = c('vertices')))

  return(pheno_circuit_no_incoming)
}

#' optimize_pheno_network
#'
#' @param sp_object sp_object list with network, nodes and edges
#' @param organism string, human or mouse
#' @param phospho_df tibble of phosphoproteomics data
#' @param carnival_options list of options returned by default_CARNIVAL_options
#' @param files boolean value, TRUE if you want output files
#' @param direct Boolean value, default FALSE, if TRUE uses only direct interactions
#' @param with_atlas Boolean value, default TRUE, if FALSE excludes Kinome Altas derived regulons
#' @param path_sif path of the sif output file of network
#' @param path_rds path of the rds output file of network
#'
#' @return SP object list with igraph object, optimized nodes df and optimized edges df
#' @export
#'
#' @examples
optimize_pheno_network <- function(sp_object,
                                   organism,
                                   phospho_df,
                                   carnival_options,
                                   files,
                                   direct = FALSE,
                                   with_atlas = FALSE,
                                   path_sif = './pheno_opt_graph.sif',
                                   path_rds = './pheno_opt_graph.rds'){

  # Create start and end nodes for CARNIVAL analysis
  start_df <- sp_object$sp_object_phenotypes$nodes_df %>%
    dplyr::filter(mf != 'phenotype') %>%
    dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf, method)

  pheno_df <- sp_object$sp_object_phenotypes$nodes_df %>%
    dplyr::filter(mf == 'phenotype') %>%
    dplyr::select(UNIPROT, gene_name, final_score = carnival_activity, mf, method)

  # Transform the optimized network edges in SIF format
  pheno_naive_df <- sp_object$sp_object_phenotypes$edges_df %>%
    dplyr::select(source, interaction = sign, target)

  pheno_out <- run_carnival_and_create_graph(source_df = start_df,
                                             target_df = pheno_df,
                                             naive_network = unique(pheno_naive_df),
                                             proteins_df = sp_object$sp_object_phenotypes$nodes_df %>%
                                               dplyr::select(UNIPROT, gene_name, final_score, mf, method),
                                             organism = organism,
                                             carnival_options = carnival_options,
                                             files = files,
                                             direct = direct,
                                             with_atlas = with_atlas,
                                             path_sif = path_sif,
                                             path_rds = path_rds)

  sp_pheno_out_validated <- expand_and_map_edges(optimized_object = pheno_out,
                                                 organism = organism,
                                                 phospho_df = phospho_df,
                                                 files = files,
                                                 direct = direct,
                                                 with_atlas = with_atlas,
                                                 path_sif = path_sif,
                                                 path_rds = path_rds)

  # Override old sp_object_phenotypes
  sp_object$sp_object_phenotypes <- sp_pheno_out_validated
  return(sp_object)
}

#' pheno_to_start_circuit
#'
#' @param SP_object sp_object list with network, nodes and edges
#' @param start_nodes character vector of nodes in the model
#' @param k integer, path length
#' @param phenotypes vector of phenotypes
#' @param start_to_top Boolean, to specify if you want to remove incoming edges on start nodes
#'
#' @return a subnetwork linking start nodes to phenotypes with a maximum of k length (shorter paths are included)
#' @export
#'
#' @examples
pheno_to_start_circuit <- function(SP_object, start_nodes, phenotypes, k, start_to_top = FALSE) {

  # SP_object = opt1$sp_object_phenotypes
  # start_nodes = start_nodes
  # phenotype = phenotypes[i_pheno]
  # k = k_vector[i_pheno]

  SP_graph <- SP_object$igraph_network

  for (i in c(1:length(phenotypes))) {
    phenotype <- phenotypes[i]
    to_nodes <- V(SP_graph)$name[V(SP_graph)$name %in% start_nodes]

    all_paths <- igraph::all_simple_paths(
      SP_graph,
      from = stringr::str_replace_all(string = phenotype, pattern = "[[:space:]\\\\/]", '_'),
      to = to_nodes,
      mode = "in",
      cutoff = k
    )

    all_paths_nodes <- unique(names(unlist(all_paths)))

    if (length(all_paths_nodes) == 0) {
      warning(paste0("No path of length ", k, " have been found for ",
                     phenotype))
    }

    if (i == 1) {
      final_nodes <- all_paths_nodes
    } else {
      final_nodes <- c(final_nodes, all_paths_nodes)
    }
  }

  pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)
  pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                               what = c("edges"))
  if(start_to_top == TRUE){
    pheno_circuit_edges <- pheno_circuit_edges %>% dplyr::filter(!to %in% start_nodes)
  }

  pheno_circuit_no_incoming <-
    igraph::graph_from_data_frame(d = pheno_circuit_edges,
                                  vertices = igraph::as_data_frame(pheno_circuit,
                                                                   what = c("vertices")))

  return(pheno_circuit_no_incoming)
}



