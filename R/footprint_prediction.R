#' Create VIPER-compatible dataset
#'
#' This function converts phosphoproteomics or transcriptomics data into a format
#' suitable for VIPER (Virtual Inference of Protein-activity by Enriched Regulon analysis)
#' (PMID:  27322546).
#'
#' @param omic_data A dataframe containing phosphoproteomics or transcriptomics data.
#' @param analysis A string indicating the type of analysis to perform. Options are:
#'   - "tfea": Transcription Factor Enrichment Analysis (TFEA).
#'   - "ksea": Kinase-Substrate Enrichment Analysis (KSEA).
#' @param significance A logical value indicating whether to filter for significant analytes.
#'   If `TRUE`, only significant analytes (where `significant == '+'`) are included.
#'
#' @return A data frame formatted for VIPER, containing columns:
#'   - `ID`: Gene name (TFEA) or phosphosite ID (KSEA).
#'   - `logFC`: Fold change value.
#'   - `t`: T-statistic for analysis.
#'   - `adj.P`: Adjusted p-value.
#'
#' @examples
#' # Example dataset with transcriptomics data
#' omic_data <- data.frame(
#'   gene_name = c("TP53", "MYC", "EGFR"),
#'   difference = c(1.5, -2.1, 0.8),
#'   logpval = c(0.01, 0.05, 0.2),
#'   significant = c("+", "+", NA)
#' )
#'
#' # Format for TFEA with only significant genes
#' formatted_data_sign <- create_viper_format(omic_data, analysis = "tfea", significance = TRUE)
#'
#' # Format for TFEA including all analytes
#' formatted_data_all <- create_viper_format(omic_data, analysis = "tfea", significance = FALSE)
#'
#' @export
create_viper_format <- function(omic_data, analysis, significance){
  if(significance){
    omic_data <- omic_data %>%
      dplyr::filter(significant == '+')
  }

  if(!analysis %in% c('tfea', 'ksea')){
    stop("Please provide a valid analysis name: 'tfea' or 'ksea'")
  }

  if(analysis == 'tfea'){
    VIPER_format <- omic_data %>%
      dplyr::select(ID = gene_name,
                    logFC = difference,
                    t = difference,
                    adj.P = logpval) %>%
      dplyr::mutate_at(c('logFC', 't'), as.numeric)
  }else{
    # Separate analytes according to collapsed UNIPROT
    omic_data <- omic_data %>%
      tidyr::separate_rows(UNIPROT, sep = ';') %>%
      dplyr::mutate(ID = paste0(UNIPROT, '-', aminoacid, '-', position))

    VIPER_format <- omic_data %>%
      dplyr::select(ID,
                    logFC = difference,
                    t = difference,
                    adj.P = logpval) %>%
      dplyr::mutate_at(c('logFC', 't'), as.numeric)
  }
  return(VIPER_format)
}

#' Convert VIPER Format to a Matrix
#'
#' This function converts a VIPER-formatted dataframe into a matrix
#' containing ID and t-statistic values.
#'
#' @param viper_format A dataframe formatted for VIPER.
#'
#' @return A matrix where rows correspond to gene IDs and columns contain t-statistic values.
#'
#' @examples
#' # Example dataset with transcriptomics data
#' omic_data <- data.frame(
#'   gene_name = c("TP53", "MYC", "EGFR"),
#'   difference = c(1.5, -2.1, 0.8),
#'   logpval = c(0.01, 0.05, 0.2),
#'   significant = c("+", "+", NA)
#' )
#' viper_format <- create_viper_format(omic_data, analysis = "tfea", significance = TRUE)
#' viper_matrix <- create_matrix_from_VIPER_format(viper_format)
#'
#' @export
create_matrix_from_VIPER_format <- function(viper_format) {
  diff_matrix <- viper_format %>%
    dplyr::select(ID, t) %>%
    dplyr::filter(!is.na(t)) %>%
    tibble::column_to_rownames(var = "ID") %>%
    as.matrix()
  return(diff_matrix)
}

#' Run VIPER Analysis
#'
#' This function runs the VIPER analysis for transcription factor or kinase activity inference.
#'
#' @param viper_format A dataframe in VIPER-compatible format.
#' @param analysis String, specifying the type of analysis (`"tfea"` or `"ksea"`).
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param reg_minsize Integer, minimum regulon size for VIPER.
#' @param integrated_regulons Logical, indicating whether to use regulons
#'   derived from experimental data (Ser/Thr and Tyr Kinome Atlas (PMIDs: 36631611, 38720073)). Default: `FALSE`.
#' @param collectri Logical, indicating whether to use CollecTRI regulons
#'      (only applicable to  `"tfea" `). Default: `FALSE`.
#' @param custom Logical indicating whether to use a custom regulon dataset. Default: `FALSE`.
#' @param custom_path String specifying the file path to a custom regulon dataset (if `custom = TRUE`). Default: `NULL`.
#'
#' @return A list containing:
#'   - `sign`: A dataframe of significantly enriched proteins.
#'   - `all`: A dataframe of all inferred proteins.
#'
#' @examples
#' # Run Transcription Factor Enrichment Analysis
#' viper_format <- create_viper_format(omic_data = tr_toy_df, analysis = 'tfea', significance = FALSE)
#' results <- run_viper(viper_format = viper_format,
#'                      analysis = "tfea",
#'                      organism = "human",
#'                      reg_minsize = 10,
#'                      integrated_regulons = FALSE)
#' results$sign
#' results$all
#'
#' @export
run_viper <- function(viper_format, analysis, organism, reg_minsize,
                      integrated_regulons = FALSE, collectri = FALSE,
                      custom = FALSE, custom_path = NULL) {

  # Validate inputs
  analysis <- match.arg(analysis, c("tfea", "ksea"))
  organism <- match.arg(organism, c("human", "mouse"))

  if (analysis == "ksea" & collectri) stop("Collectri is only a 'tfea' parameter.")
  if (custom & is.null(custom_path)) stop("Please provide a path to the custom regulon table.")

  # Convert VIPER format to matrix
  diff_matrix <- create_matrix_from_VIPER_format(viper_format)

  # Load regulons
  if (custom) {
    message("Reading custom regulons...")
    regulons <- readr::read_tsv(custom_path, show_col_types = FALSE)
  } else {
    regulons <- switch(
      analysis,
      "tfea" = {
        if (organism == "human") {
          if (collectri) get(data("tfea_db_human_collectri")) else get(data("tfea_db_human"))
        } else {
          get(data("tfea_db_mouse"))
        }
      },
      "ksea" = {
        if (organism == "human") {
          if (integrated_regulons) get(data("ksea_db_human_atlas")) else get(data("ksea_db_human"))
        } else {
          get(data("ksea_db_mouse"))
        }
      }
    )
  }

  # Ensure overlap between measured analytes and regulons
  measured_analytes <- intersect(rownames(diff_matrix), regulons$target)
  if (length(measured_analytes) == 0) {
    stop("No measured analytes found in regulons. Ensure that the regulons' organism and experimental organism match.")
  }

  # Run VIPER analysis
  viper_object <- viper::msviper(
    diff_matrix,
    dorothea::df2regulon(regulons),
    minsize = reg_minsize,
    ges.filter = FALSE,
    cores = 1,
    verbose = FALSE,
    pleiotropy = TRUE
  )

  # Extract results
  activities_stat <- as.data.frame(viper_object$es$nes.bt) %>%
    dplyr::rename(NES = t) %>%
    dplyr::mutate(pvalues = viper_object$es$p.value) %>%
    tibble::rownames_to_column("gene_name")

  prot_sign <- dplyr::filter(activities_stat, pvalues < 0.05)

  return(list(sign = prot_sign,
              all = activities_stat))
}

#' Run Hypergeometric Test for VIPER Correction
#'
#' This function corrects VIPER output using a hypergeometric test to weight
#' VIPER inferred protein activities  according to the number of
#' significantly modulated analytes (transcripts or phosphosites) in the regulon.
#'
#' @param omic_data A dataframe containing the experimental omics dataset.
#' @param viper_output A dataframe containing the VIPER-enriched proteins.
#' @param analysis String, specifying the type of analysis (`"tfea"` or `"ksea"`).
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param integrated_regulons Logical, indicating whether to use regulons
#'   derived from experimental data (Ser/Thr and Tyr Kinome Atlas (PMIDs: 36631611, 38720073)). Default: `FALSE`.
#' @param collectri Logical, indicating whether to use CollecTRI regulons
#'      (only applicable to  `"tfea" `). Default: `FALSE`.
#' @param custom Logical indicating whether to use a custom regulon dataset. Default: `FALSE`.
#' @param custom_path String specifying the file path to a custom regulon dataset (if `custom = TRUE`). Default: `NULL`.
#''
#' @return A dataframe containing the VIPER correction weights
#' based on the hypergeometric test.
#'
#' @examples
#' data('tr_toy_df')
#' viper_format <- create_viper_format(omic_data = tr_toy_df, analysis = 'tfea', significance = FALSE)
#' viper_output <- run_viper(viper_format = viper_format,
#'                      analysis = "tfea",
#'                      organism = "human",
#'                      reg_minsize = 10,
#'                      integrated_regulons = FALSE)
#' ea_output <- run_hypergeometric_test(omic_data = tr_toy_df,
#'                                      viper_output = viper_output$sign,
#'                                      analysis = "tfea",
#'                                      organism = "human")
#'
#' @export
#'
run_hypergeometric_test <- function(omic_data, viper_output,
                                    analysis, organism,
                                    integrated_regulons = FALSE,
                                    collectri = FALSE,
                                    custom = FALSE,
                                    custom_path = NULL) {
  # Validate inputs
  analysis <- match.arg(analysis, c("tfea", "ksea"))
  organism <- match.arg(organism, c("human", "mouse"))

  if (analysis == "ksea" & collectri) stop("Collectri is only a 'tfea' parameter.")
  if (custom & is.null(custom_path)) stop("Please provide a path to the custom regulon table.")

  # Prepare measured and significant analytes
  # Universe will be measured (uni_meas_v) and significant (uni_sign_v)
  # transcripts/phosphosites
  if (analysis == "tfea") {
    uni_meas_v <- omic_data$gene_name
    uni_sign_v <- omic_data %>%
      dplyr::filter(significant == "+") %>%
      dplyr::pull(gene_name)
  } else {
    omic_data <- omic_data %>%
      dplyr::mutate(phosphositeID = paste0(toupper(UNIPROT), "-", aminoacid, "-", position))
    uni_meas_v <- omic_data$phosphositeID
    uni_sign_v <- omic_data %>%
      dplyr::filter(significant == "+") %>%
      dplyr::pull(phosphositeID)
  }

  # Load regulons
  if (custom) {
    message("Reading custom regulons...")
    df_regulons <- readr::read_tsv(custom_path, show_col_types = FALSE)
  } else {
    df_regulons <- switch(
      analysis,
      "tfea" = if (organism == "human") {
        if (collectri) get(data("tfea_db_human_collectri")) else get(data("tfea_db_human"))
      } else {
        get(data("tfea_db_mouse"))
      },
      "ksea" = if (organism == "human") {
        if (integrated_regulons) access_remote_file(file = 'ksea_db_human_atlas.tsv', dir = 'PKN') else get(data("ksea_db_human"))
      } else {
        get(data("ksea_db_mouse"))
      }
    )
  }

  # Extract measured & significant analytes from regulons
  all_enriched_matrix <- create_matrix_from_VIPER_format(create_viper_format(omic_data, analysis, significance = FALSE))
  all_significant_matrix <- create_matrix_from_VIPER_format(create_viper_format(omic_data, analysis, significance = TRUE))

  an_in_reg <- df_regulons %>%
    dplyr::filter(target %in% rownames(all_enriched_matrix)) %>%
    dplyr::count(tf, name = "Measured")

  an_in_reg_sign <- df_regulons %>%
    dplyr::filter(target %in% rownames(all_significant_matrix)) %>%
    dplyr::count(tf, name = "Significant")

  # Merge enrichment data with VIPER output
  pr_joined <- viper_output %>%
    dplyr::left_join(an_in_reg, by = c("gene_name" = "tf")) %>%
    dplyr::left_join(an_in_reg_sign, by = c("gene_name" = "tf")) %>%
    dplyr::mutate(Measured = tidyr::replace_na(Measured, 0), Significant = tidyr::replace_na(Significant, 0))

  # Compute hypergeometric p-values for each gene_name
  pr_joined$pWeight <- purrr::map_dbl(pr_joined$gene_name, function(protein) {

    significant_members <- df_regulons %>%
      dplyr::filter(target %in% rownames(all_significant_matrix), tf == protein) %>%
      dplyr::pull(target)

    regulon_members <- df_regulons %>%
      dplyr::filter(target %in% rownames(all_enriched_matrix), tf == protein) %>%
      dplyr::pull(target)

    phyper(q = length(significant_members),
           m = length(uni_sign_v),
           n = length(uni_meas_v) - length(regulon_members),
           k = length(regulon_members),
           lower.tail = FALSE)
  })
  pr_joined <- pr_joined %>%
    dplyr::rename(reg_exp_meas = Measured, reg_exp_sign = Significant)

  return(pr_joined)
}

#' Adjust VIPER Output Scores with Hypergeometric Test
#'
#' This function adjusts the normalized enrichment scores (NES) from VIPER
#' using weights derived from a hypergeometric test.
#'
#' @param ea_output A dataframe containing VIPER output, including NES and enrichment scores.
#'
#' @return A dataframe with weighted NES scores for more refined protein activity inference.
#'
#' @examples
#'
#' data('tr_toy_df')
#' viper_format <- create_viper_format(omic_data = tr_toy_df, analysis = 'tfea', significance = FALSE)
#' viper_output <- run_viper(viper_format = viper_format,
#'                      analysis = "tfea",
#'                      organism = "human",
#'                      reg_minsize = 10,
#'                      integrated_regulons = FALSE)
#' ea_output <- run_hypergeometric_test(tr_toy_df,
#'                                      viper_output$sign,
#'                                      "tfea",
#'                                      "human")
#' weighted_results <- weight_viper_score(ea_output=ea_output)
#'
#' @export
#'
weight_viper_score <- function(ea_output) {

  # Ensure pWeight has no zeros to avoid log(0) errors
  # replacing 0 values with min of other pWeights
  min_pWeight <- min(ea_output$pWeight[ea_output$pWeight > 0], na.rm = TRUE)
  ea_output <- ea_output %>%
    dplyr::mutate(pWeight = dplyr::if_else(pWeight == 0, min_pWeight, pWeight)) %>%
    dplyr::mutate(
      raw_weight = -log(pWeight),  # Log-transform hypergeometric p-values
      weight = raw_weight / max(raw_weight, na.rm = TRUE),  # Normalize weights
      weightedNES = dplyr::if_else(pWeight < 0.05, NES * weight * 4, NES)  # Scale NES for significant proteins
    ) %>%
    dplyr::arrange(gene_name)
  return(ea_output)
}

#' Filter VIPER Output for Molecular Function
#'
#' This function filters inferred proteins based on their GO molecular function annotated
#' with [molecular_function_annotation] function of SignalingProfiler:
#' - Transcription Factors for TFEA;
#' - Kinases/Phosphatases for KSEA.
#'
#' @param inferred_proteins_mf A dataframe containing inferred proteins with molecular function annotations.
#' @param analysis A string specifying the analysis type ("tfea" or "ksea").
#'
#' @return A dataframe of inferred proteins filtered based on molecular function.
#'
#' @seealso [molecular_function_annotation]
#'
#' @examples
#' # Filter transcription factors in human
#' inferred_proteins <- data.frame(gene_name = c('TP53', 'CDK5'),
#'                                 mf = c('tf', 'kin'))
#' filtered_tfs <- filter_VIPER_output(inferred_proteins_mf = inferred_proteins,
#'                                      analysis = "tfea")
#'
#' @export
filter_VIPER_output <- function(inferred_proteins_mf, analysis) {
  if(analysis == 'tfea'){
    inferred_proteins_mf <- inferred_proteins_mf %>% dplyr::filter(mf == 'tf')
  }else if(analysis == 'ksea'){
    inferred_proteins_mf <- inferred_proteins_mf %>% dplyr::filter(mf == 'kin' | mf == 'phos')
  }
  return(inferred_proteins_mf)
}

#' Run Footprint-Based Protein Activity Analysis
#'
#' This function performs footprint-based analysis using VIPER to infer protein activities from omics data.
#' - `'analysis=='tfea'`: it performs transcription factor enrichment analysis (TFEA) from transcriptomics data;
#'   transcription factors-target genes pairs are from either Dorothea+SIGNOR or CollecTRI+SIGNOR.
#' - `'analysis=='ksea'`: it performs kinase substrates enrichment analysis (KSEA) from phosphoproteomics data;
#'   kinase/phosphatase-phosphosites pairs are from OmniPath+SIGNOR and can be integrated with pairs derived from
#'   experimental data collected in Serine/Threonine/Tyrosine Kinome Atlas (PMIDs: 36631611, 38720073).
#' If `'hypergeom_corr==TRUE`, VIPER activity score is weighted on the number of significant analytes within the regulon.
#' If `'correct_proteomics==TRUE`, VIPER output is turned significant considering proteomics fold-change of analytes.
#' If VIPER returns non-significant modulation but the same modulation is significant in proteomics,
#' includes the analyte in VIPER result.
#'
#' Custom set of regulons can be provided using `custom` and `custom_path` parameters.
#'
#' @param omic_data A dataframe of measured transcript (for `"tfea"`) or phosphosites (for `"ksea"`).
#' @param analysis A string specifying the type of analysis (`"tfea"` or `"ksea"`).
#' @param organism Character string, either `"human"` or `"mouse"`, specifying the organism.
#' @param reg_minsize Integer, minimum regulon size for VIPER.
#' @param exp_sign  Logical value indicating whether to filter for significant analytes.
#'   If `TRUE`, only significant analytes (where `significant == '+'`) are included.
#' @param integrated_regulons Logical, indicating whether to use regulons
#'   derived from experimental data (Ser/Thr and Tyr Kinome Atlas (PMIDs: 36631611, 38720073)). Default: `FALSE`.
#' @param collectri Logical, indicating whether to use CollecTRI regulons
#'      (only applicable to  `"tfea" `). Default: `FALSE`.
#' @param hypergeom_corr Logical indicating whether to apply hypergeometric correction.
#' @param correct_proteomics Logical indicating whether to adjust results using proteomics data.
#' @param prot_df A dataframe containing proteomics data (required if `correct_proteomics = TRUE`).
#' @param GO_annotation Logical indicating whether to perform GO molecular function annotation using [molecular_function_annotation].
#' @param custom Logical indicating whether to use a custom regulon dataset. Default: `FALSE`.
#' @param custom_path String specifying the file path to a custom regulon dataset (if `custom = TRUE`). Default: `NULL`.

#' @return A dataframe containing inferred transcription factors (TFEA) or kinases/phosphatases (KSEA).
#' The dataframe contains:
#' \describe{
#'  \item{`UNIPROT`}{The UNIPROT identifier.}
#'   \item{`gene_name`}{The gene symbol.}
#'   \item{`NES`}{The inferred activity score by VIPER.}
#'   \item{`pvalue`}{The inferred activity associated p-value by VIPER.}
#'   \item{`reg_exp_meas`}{The number of regulons targets experimentally measured in omics data.}
#'   \item{`reg_exp_sign`}{The number of regulons targets significantly modulated in omics data.}
#'   \item{`pWeight`}{The hypergeometric test derived  p-value.}
#'   \item{`raw_weight`}{The -logarithm of hypergeometric test p-value.}
#'   \item{`weight`}{The scaled logarithm of hypergeometric test p-value.}
#'   \item{`weightedNES`}{`NES` multiplied by `weight`, inferred activity weighted on the number of significantly modulated regulons members.}
#'   \item{`mf`}{GO Molecular Function tailored on signaling.}
#' }
#'
#' @examples
#'
#' # Transcription Factor Enrichment Analysis Example
#' data('tr_toy_df')
#' tfea_df <- run_footprint_based_analysis(omic_data = tr_toy_df,
#'                                         analysis = "tfea",
#'                                         organism = "human",
#'                                         reg_minsize = 10,
#'                                         hypergeom_corr = TRUE,
#'                                         collectri = FALSE,
#'                                         GO_annotation = TRUE)
#' # Kinase Enrichment Analysis Example
#' data('phospho_toy_df')
#' ksea_df <- run_footprint_based_analysis(omic_data = phospho_toy_df,
#'                                         analysis = "ksea",
#'                                         organism = "human",
#'                                         reg_minsize = 5,
#'                                         hypergeom_corr = TRUE,
#'                                         integrated_regulons = TRUE,
#'                                         GO_annotation = TRUE)
#'
#' @export
#'
run_footprint_based_analysis <- function(omic_data,
                                         analysis,
                                         organism,
                                         reg_minsize,
                                         exp_sign = FALSE,
                                         integrated_regulons = FALSE,
                                         collectri = FALSE,
                                         hypergeom_corr = TRUE,
                                         correct_proteomics = FALSE,
                                         prot_df = NULL,
                                         GO_annotation = FALSE,
                                         custom = FALSE,
                                         custom_path = NULL) {

  message(" ** RUNNING FOOTPRINT-BASED ANALYSIS ** ")
  message("Credits: Dugourd A, Saez-Rodriguez J. (2019) Curr Opin Syst Biol. doi: 10.1016/j.coisb.2019.04.002.")

  # Validate inputs
  analysis <- match.arg(analysis, c("tfea", "ksea"))
  organism <- match.arg(organism, c("human", "mouse"))

  if (analysis == "ksea" & collectri) stop("Collectri is only a 'tfea' parameter.")
  if (correct_proteomics & is.null(prot_df)) stop("Please provide proteomics data.")

  # Standardize gene names
  omic_data <- omic_data %>%
    dplyr::mutate(gene_name = stringr::str_replace_all(gene_name, "[^[:alnum:]]", "_"))

  # Run VIPER analysis
  message("Starting VIPER analysis...")
  viper_format <- create_viper_format(omic_data, analysis, significance = exp_sign)

  output <- run_viper(viper_format = viper_format,
                      analysis = analysis,
                      organism = organism,
                      reg_minsize = reg_minsize,
                      integrated_regulons = integrated_regulons,
                      collectri = collectri,
                      custom = custom,
                      custom_path = custom_path)

  if (nrow(output$sign) == 0) stop("No modulated proteins found!")

  # Correct VIPER output with proteomics if required
  if (correct_proteomics) {
    message("Applying proteomics correction...")
    viper_on_prot <- dplyr::left_join(output$all,
                                      prot_df,
                                      by = c('gene_name')) %>%
      dplyr::mutate(pvalues = ifelse(!is.na(difference) & !is.na(logpval) & significant == '+' & #check on proteomics fields
                                       NES*difference > 0 & pvalues > 0.05, # check on VIPER fields
                                     10^(-logpval), # assign proteomic pvalue
                                     pvalues)) %>% # keep VIPER pvalue
      dplyr::select(gene_name, NES, pvalues)

    # Override VIPER output
    output$sign <- viper_on_prot %>% dplyr::filter(pvalues < 0.05)
    output$all <- viper_on_prot %>% dplyr::filter(pvalues > 0.05)
  }

  # Apply hypergeometric test correction
  if (hypergeom_corr) {
    message("Starting hypergeometric test correction...")
    output <- weight_viper_score(run_hypergeometric_test(omic_data = omic_data,
                                                         viper_output = output$sign,
                                                         analysis = analysis,
                                                         organism = organism,
                                                         collectri = collectri,
                                                         integrated_regulons = integrated_regulons,
                                                         custom = custom,
                                                         custom_path = custom_path))
    output_uniprot <- convert_gene_name_in_uniprotid(output, organism)
  }else{
    output_uniprot <- convert_gene_name_in_uniprotid(output$sign, organism)
  }

  # Apply GO annotation
  if (GO_annotation) {
    message("Applying GO molecular function annotation...")
    output_uni_mf <- molecular_function_annotation(output_uniprot, organism)
    output_final <- filter_VIPER_output(output_uni_mf, analysis) %>% dplyr::distinct()
    return(output_final)
  }

  return(output_uniprot)
}

