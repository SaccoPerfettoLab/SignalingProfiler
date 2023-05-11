
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
#' @export
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
    g <- sp_output$igraph_network
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
        dist_count[1,i] <- unlist(igraph::distances(g, combinatios[1,i], to = combinatios[2,i], mode = "out"))
        if(dist_count[1,i] == Inf){
          dist_count[2,i] <- NA
        }else{
          dist_count[2,i] <- combinatios[2,i]
        }
      }
      # reverse run
      for(i in ncol(combinatios)){
        #i = 1
        dist_count[1,i] <- unlist(igraph::distances(g, combinatios[2,i], to = combinatios[1,i], mode = "out"))
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
#'
#' @return a list of a plot and result table
#' @export
#'
#' @examples
phenoscore_computation_v1 <- function(proteins_df,
                                      desired_phenotypes,
                                      pvalue_threshold = 0.05,
                                      path_length = 4,
                                      zscore_threshold = -1.96){
  # proteins_df =  proteins_df
  # desired_phenotypes = desired_phenotypes
  # # ** TEST DATA ** #
  # proteins_df <- JMD_df[which(JMD_df$EAPred > 0),] %>%
  #   rename(gene_name = GeneName)
  #
  # phenotype_list_Sara_df <- read_xlsx('../data-raw/VERONICA/Phenotype_phenoscore.xlsx', sheet = 1, range = NULL, col_names = TRUE,
  #                                     col_types = NULL, na = "", trim_ws = TRUE, skip = 0,
  #                                     progress = readxl_progress(), .name_repair = "unique")
  # phenotype_list_Sara_df%>%
  #   filter(Choose == 'x')-> phenotype_list_Sara_df
  # phenotype_list_Sara <- toupper(phenotype_list_Sara_df$Phenotype_name)
  # phenotype_list_Sara <- gsub(' ', '_', phenotype_list_Sara)
  # desired_phenotypes <- phenotype_list_Sara
  #
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
    if(path_length <= 0 | path_length > 4){
      stop('path_length should be 1,2,3 or 4!')
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
                     sd=sd(Path_Score)) -> summary.table.mean.sd.path.score ### here we can select between mean and median

  ##merge the path table with mean/median and standard deviation (sd)
  merge(x = output.table2,
        y = summary.table.mean.sd.path.score,
        by = c('EndPathways', 'Effect'), all.x = TRUE) -> output.table.mean.sd.path.score

  ##calculate z-score
  output.table.mean.sd.path.score$zscore <- (output.table.mean.sd.path.score$Path_Score-output.table.mean.sd.path.score$mean)/output.table.mean.sd.path.score$sd

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

  results.table$EndPathways <- sapply(strsplit(as.character(results.table$EndPathways ),
                                               "="),"[", 1)

  # generate plot
  message('Plot generation')
  ### display interactions from INPUT genes to phenotypes
  id_table <- data.frame(matrix (ncol = 1, nrow = 0))
  colnames(id_table) <- c("signor_ids")


  output.table6 <- output.table6 %>%
    dplyr::mutate_at(c('EndPathways', 'QueryNode'), as.character)

  output.table6%>%
    dplyr::filter(QueryNode %in% proteins_vector)-> INPUT.Phen.Paths

  ProteinsPaths <- INPUT.Phen.Paths %>%
    dplyr::group_by(EndPathways, Effect) %>%
    dplyr::summarise(regulators = paste0(unique(QueryNode),collapse = ';')) %>%
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
