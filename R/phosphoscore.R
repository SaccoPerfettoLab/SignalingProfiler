#' Infer protein actiity from phosphosite regulatory power (PhophoScore)
#'
#' This function infers protein activity by combining phosphosites' significant
#' modulation in phosphoproteomics data with their regulatory role (*activation* or *inhibition*)
#' available in SIGNOR and PhosphoSitePlus databases. This activity measure
#' is called PhosphoScore.
#'
#' @param phosphoproteomic_data A dataframe containing phosphoproteomics measurements.
#' @param organism Character string, either `"human"`, `"mouse"`, or `"hybrid"`, specifying the organism.
#' If  `"hybrid"` mouse regulatory information is remapped on human to enlarge coverage.
#' @param activatory Logical. If `TRUE`, considers only phosphosites regulating protein activity.
#' If `FALSE`, considers phosphosites regulating also protein abundance. Default: `TRUE`.
#' @param path_fasta A string specifying the path to the FASTA file for `blastp` Necessary when `organism='hybrid'` Default: `'./phospho.fasta'`.
#' @param blastp_path Optional. A string specifying the Windows path of `blastp`.
#' @param local Logical. For **developmental purposes only**. Default: `FALSE`.
#' @param GO_annotation Logical indicating whether to perform GO molecular function annotation using [molecular_function_annotation].
#' @param custom Logical. If `TRUE`, allows use of custom information on regulatory phosphosites.
#' @param custom_path A string specifying the path to a `.tsv` file containing custom  information on regulatory phosphosites.
#'
#' @return A dataframe containing:
#' \describe{
#'   \item{`gene_name`}{The gene symbol.}
#'   \item{`phosphoscore`}{The inferred activity score.}
#'   \item{`n_sign_phos`}{Number of significant phosphosites used in scoring.}
#'   \item{`phos`}{Concatenated phosphosites contributing to the score.}
#'   \item{`act_rol`}{Activation or inhibition roles of each phosphosite.}
#'   \item{`phosphosite_fc`}{Fold-change values for the phosphosites.}
#' }
#'
#' @examples
#' # Example phosphoproteomics dataset
#' phospho_data <- data.frame(
#'   gene_name = c("TP53", "MYC", "EGFR"),
#'   aminoacid = c("S", "T", "Y"),
#'   position = c(15, 22, 45),
#'   sequence_window = c("XXXXXXXXSXXXXXXX", "XXXXXXXTXXXXXXXX", "XXXXXXXYXXXXXXXX"),
#'   significant = c("+", "+", NA),
#'   difference = c(1.2, -0.8, 0.5)
#' )
#'
#' # Compute PhosphoSCORE
#' phosphoscore_results <- phosphoscore_computation(phospho_data, "human", TRUE)
#' head(phosphoscore_results)
#'
#' @export
#'
phosphoscore_computation <- function(phosphoproteomic_data,
                                     organism,
                                     activatory,
                                     path_fasta = './phospho.fasta',
                                     blastp_path = NULL,
                                     local = FALSE,
                                     GO_annotation = FALSE,
                                     custom = FALSE,
                                     custom_path = NULL) {

  message("** RUNNING PHOSPHOSCORE ANALYSIS **")

  # Validate inputs
  organism <- match.arg(organism, c("human", "mouse", "hybrid"))

  # Select regulatory phosphosites
  phosphoscore_df_output <- if (organism %in% c("human", "mouse")) {
    map_experimental_on_regulatory_phosphosites(phosphoproteomic_data, organism, activatory,
                                                path_fasta, custom, custom_path)
  } else {
    phospho_score_hybrid_computation(phosphoproteomic_data, organism, activatory, blastp_path, path_fasta, local)
  }
  phosphoscore_df <- phosphoscore_df_output$phosphoscore_df

  if (is.null(phosphoscore_df)) {
    warning("No proteins found from PhosphoScore computation!")
    return(NULL)
  }

  # Flag to deal with hybrid PhosphoScore computation
  if('source_org' %in% colnames(phosphoscore_df)){
    n_species <- length(unique(phosphoscore_df$source_org))
  }else{
    n_species <- 1
  }

  # Compute PhosphoScore (average the contribution of multiple phosphosites over a protein)
  raw_output <- phosphoscore_df %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(phosphoscore = mean(inferred_activity, na.rm = TRUE))

  # Extract and merge experimental phosphosites data
  exp_fc_sub <- phosphoscore_df_output$used_exp_data %>%
    dplyr::select(gene_name, PHOSPHO_KEY_GN_SEQ, aminoacid, position) %>%
    dplyr::distinct() %>%
    dplyr::left_join(phosphoscore_df %>%
                       dplyr::select(PHOSPHO_KEY_GN_SEQ, ACTIVATION, difference),
                     by = "PHOSPHO_KEY_GN_SEQ")

  # For mouse transform the mapping key to lowercase
  # in both experimental data and PhosphoScore output
  if(organism %in% c('hybrid', 'mouse')){
    exp_fc_sub <- exp_fc_sub %>%
      dplyr::mutate_at('PHOSPHO_KEY_GN_SEQ', toupper) %>%
      dplyr::distinct()

    exp_fc_sub <- exp_fc_sub %>% dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(gene_name, '-',
                                                                           sub('.*-', replacement = '', exp_fc_sub$PHOSPHO_KEY_GN_SEQ)))

    # change also raw output
    raw_output <- raw_output %>%
      dplyr::mutate_at('gene_name', str_to_title)
  }

  # Annotate for each protein the phosphosites used for the prediction
  # and merge it with PhosphoScore table
  exp_fc_sub <- exp_fc_sub %>%
    dplyr::mutate(
      aa = paste0(aminoacid, position),
      gene_name = (gene_name),
      difference = as.character(round(exp_fc_sub$difference,2))) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n_sign_phos = dplyr::n(),
                     phos = paste0(aa, collapse = ';'),
                     act_rol = paste0(ACTIVATION, collapse =';'),
                     phosphosite_fc = paste0(difference, collapse =';'))

  output <- dplyr::left_join(raw_output, exp_fc_sub, by = "gene_name")

  # If 'hybrid' organism, distinguish between human and mouse phosphosites
  if(organism == 'hybrid'){
    if(n_species == 2){
      genes <- phosphoscore_df %>%
        dplyr::select(gene_name, source_org) %>%
        dplyr::distinct() %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarise(source_org = paste0(source_org, collapse = ';'))
      output <- dplyr::left_join(output, genes, by = c('gene_name'))
    }
  }

  # Add UNIPROT ID
  output <- convert_gene_name_in_uniprotid(output, organism)

  # Apply GO annotation if required
  if (GO_annotation) {
    output <- molecular_function_annotation(output, organism) %>%
      dplyr::distinct()
  }

  return(output)
}

#' Infer protein actiity from phosphosite regulatory power (PhophoScore)
#' considering residue and position as mapping key
#'
#' The principle is the same as [compute_phosphoscore] but this function
#' should be used when in phosphoproteomics data no phosphopeptide is provided.
#'
#' @param phosphoproteomic_data A dataframe containing phosphoproteomics data.
#' @param organism A string specifying the organism (`"human"` or `"mouse"`).
#' @param activatory Logical. If `TRUE`, considers only activatory sites.
#' @param GO_annotation Logical. If `TRUE`, performs Gene Ontology (GO) annotation.
#' @param custom Logical. If `TRUE`, uses a custom regulatory database.
#' @param custom_path A string specifying the path to a `.tsv` file containing custom phosphosites.
#'
#' @return A dataframe containing inferred protein activities using amino acid and position-based mapping.
#'
#' @examples
#' phospho_data <- data.frame(
#'   gene_name = c("TP53", "MYC", "EGFR"),
#'   aminoacid = c("S", "T", "Y"),
#'   position = c(15, 22, 45),
#'   significant = c("+", "+", NA),
#'   difference = c(1.2, -0.8, 0.5)
#' )
#'
#' phosphoscore_aapos_results <- phosphoscore_computation_aapos(phospho_data, "human", TRUE)
#'
#' @export
#'
phosphoscore_computation_aapos <- function(phosphoproteomic_data,
                                           organism,
                                           activatory,
                                           GO_annotation = FALSE,
                                           custom = FALSE,
                                           custom_path = NULL) {

  message("** RUNNING PHOSPHOSCORE ANALYSIS WITH AMINO ACID AND POSITION **")

  # Validate inputs
  organism <- match.arg(organism, c("human", "mouse"))

  # Select regulatory phosphosites
  reg_phos_db <- if (custom) {
    if (is.null(custom_path)) stop("Please provide a path to the regulatory phosphosites table")
    readr::read_tsv(custom_path, show_col_types = FALSE)
  } else {
    if (activatory)
      if(organism == 'human') get(data("good_phos_df_human_act_aapos")) else get(data("good_phos_df_mouse_act_aapos"))
    else
      if(organism == 'human') get(data("good_phos_df_human_all_aapos")) else get(data("good_phos_df_mouse_all_aapos"))
  }

  # Filter experimental data
  exp_fc <- phosphoproteomic_data %>%
    dplyr::filter(significant == '+') %>%
    dplyr::mutate(residue_mapp = ifelse(aminoacid == 'S', 'Ser',
                                        ifelse(aminoacid == 'Y', 'Tyr',
                                               ifelse(aminoacid == 'T', 'Thr', aminoacid))))

  if(organism == 'human'){
    exp_fc <- exp_fc %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(toupper(gene_name), '-', residue_mapp, position)) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% reg_phos_db$PHOSPHO_KEY_GN_SEQ)
  }else{
    exp_fc <- exp_fc %>%
      dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(stringr::str_to_title(gene_name), '-', residue_mapp, position)) %>%
      dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% reg_phos_db$PHOSPHO_KEY_GN_SEQ)
  }

  if(nrow(exp_fc) == 0){
    warning('No annotated regulatory phosphosites significantly modulated in your dataset')
    return(NULL)
  }

  phosphoscore_df <- dplyr::left_join(exp_fc, reg_phos_db, by = 'PHOSPHO_KEY_GN_SEQ') %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, ACTIVATION, difference) %>%
    dplyr::arrange(PHOSPHO_KEY_GN_SEQ) %>%
    dplyr::mutate(gene_name  = unlist(lapply(PHOSPHO_KEY_GN_SEQ,
                                             function(x){stringr::str_split(x, '-')[[1]][1]})),
                  inferred_activity = as.numeric(ACTIVATION) * as.numeric(difference)) %>%
    dplyr::distinct()

  # Compute PhosphoScore (average the contribution of multiple phosphosites over a protein)
  raw_output <- phosphoscore_df %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(phosphoscore = mean(inferred_activity))

  exp_fc_sub <- exp_fc %>%
    dplyr::select(gene_name, PHOSPHO_KEY_GN_SEQ, aminoacid, position) %>%
    dplyr::distinct()

  exp_fc_sub <- dplyr::left_join(exp_fc_sub,
                                 phosphoscore_df %>%
                                   dplyr::select(PHOSPHO_KEY_GN_SEQ,
                                                 ACTIVATION,
                                                 difference),
                                 by = 'PHOSPHO_KEY_GN_SEQ')

  # Annotate for each protein the phosphosites used for the prediction
  # and merge it with PhosphoScore table
  exp_fc_sub <- exp_fc_sub %>%
    dplyr::mutate(aa = paste0(aminoacid, position),
                  gene_name = (gene_name),
                  difference = as.character(round(exp_fc_sub$difference,2))) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n_sign_phos = dplyr::n(),
                     phos = paste0(aa, collapse = ';'),
                     act_rol = paste0(ACTIVATION, collapse =';'),
                     phosphosite_fc = paste0(difference, collapse =';'))

  output <- dplyr::left_join(raw_output, exp_fc_sub, by = 'gene_name')

  # Add UNIPROT ID
  output <- convert_gene_name_in_uniprotid(output, organism)

  # Apply GO annotation if required
  if (GO_annotation) {
    output <- molecular_function_annotation(output, organism) %>%
      dplyr::distinct()
  }

  return(output)
}


#' Create a FASTA File from Phosphoproteomics Data
#'
#' This function generates a FASTA file from phosphoproteomics data, which can be used for sequence alignment.
#' The file has the following format:
#' >MAPK1-XXXXXXXXSXXXXXXX-S-15-mouse
#' XXXXXXXXSXXXXXXX
#'
#' @param phospho_df A dataframe containing phosphoproteomics data.
#' @param path A string specifying the output path of the FASTA file.
#'
#' @return This function writes a FASTA file to the specified path and returns `NULL`.
#'
#' @examples
#' # Example phosphoproteomics dataset
#' phospho_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   UNIPROT = c("P04637", "P01106"),
#'   sequence_window = c("XXXXXXXXSXXXXXXX", "XXXXXXXTXXXXXXXX"),
#'   aminoacid = c("S", "T"),
#'   position = c(15, 22)
#' )
#'
#' # Create a FASTA file
#' create_fasta(phospho_data, "output.fasta")
#'
#' @export
#'
create_fasta <- function(phospho_df, path) {
  message("Creating FASTA file from phosphoproteomics data...")

  # Ensure proper gene name formatting
  phospho_df <- phospho_df %>%
    dplyr::mutate(gene_name = toupper(stringr::word(gene_name, 1, sep = ";"))) %>%
    dplyr::arrange(gene_name)

  # Standardize sequence window length to 15-mer
  if(nchar(phospho_df$sequence_window[1]) != 15){
    center <- (nchar(phospho_df$sequence_window[1])+1)/2
    phospho_df <- phospho_df %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window,
                                                           center - 7,
                                                           center + 7))
  }else{
    phospho_df <- phospho_df %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }

  # Generate FASTA format
  fasta <- ''
  for(i in c(1:length(phospho_df$UNIPROT))){
    write(paste0('>', phospho_df$gene_name[i], '-',
                 phospho_df$sequence_window_sub[i], '-',
                 phospho_df$aminoacid[i], '-',
                 phospho_df$position[i], '-mouse\n',
                 gsub('_', '', phospho_df$sequence_window_sub[i])),
          path, append = TRUE)
  }

  message("FASTA file successfully created at: ", path)
  return(NULL)
}

#' Perform BLAST Alignment on Phosphoproteomics Data
#'
#' Runs BLAST to align experimental phosphoproteomics data against
#' a reference database of phosphopeptides.
#' In SignalingProfiler this function is used to map mouse phosphopeptides
#' on human phosphopeptides with annotated regulatory role
#'
#' @param path_experimental_fasta_file A string specifying the path to the experimental FASTA file for `blastp`.
#' @param all Logical. If `TRUE`, returns all alignments; otherwise, only those with the same gene name.
#' @param blastp_path Optional. A string specifying the Windows path of `blastp`.
#' @param local Logical. For **developmental purposes only**. Default: `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{`mapped`}{A dataframe of mouse phosphopeptides aligned to human phosphosites.}
#'   \item{`tocheck`}{A dataframe of ambiguous alignments (if `all = TRUE`).}
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Running BLAST with an experimental FASTA file
#' blast_results <- run_blast("experiment.fasta", all = FALSE)
#' head(blast_results$mapped)
#'}
#'
run_blast <- function(path_experimental_fasta_file, all = FALSE,
                      blastp_path = NULL, local = FALSE) {

  message("Running BLASTp alignment...")

  path_package <- if (local) "./" else paste0(.libPaths()[1], "/SignalingProfiler/")
  reference_db <- paste0(path_package, "extdata/human_phosphosites_db.fasta")

  command <- if (!is.null(blastp_path)) blastp_path else "blastp"
  
  args <- c(
    "-query", path_experimental_fasta_file,
    "-subject", reference_db,
    "-out", "map2.out",
    "-outfmt", "7",
    "-evalue", "0.05"
  )
  
  system2(command = command, args = args)
  
  message("BLASTp finished.")

  mapped <- readr::read_tsv("map2.out",
                            comment = "#",
                            col_names = c("q_ID", "s_ID", "pid", "length",
                                          "mismatch", "gapopen", "qstart",
                                          "qend", "sstart", "ssend", "evalue",
                                          "bitscore"))

  # Extract exact matches
  mapped_exact <- mapped %>%
    dplyr::filter(qstart == 1 & qend == 15 & sstart == 1 & ssend == 15)

  mapped$position <- as.numeric(stringr::word(mapped$q_ID, 4, sep = "-"))
  mapped_exact <- dplyr::bind_rows(mapped %>% dplyr::filter(position < 8), mapped_exact)

  # Identify gene names in query and subject sequences
  mapped_exact <- mapped_exact %>%
    dplyr::mutate(qProt = stringr::word(q_ID, 1, sep = "-"),
                  sProt = stringr::word(s_ID, 1, sep = "-"))

  tocheck <- mapped_exact %>% dplyr::filter(qProt != sProt & mismatch == 0)
  sure <- mapped_exact %>% dplyr::filter(qProt == sProt)

  system2("rm ./map2.out")

  return(if (all) list(mapped = sure, tocheck = tocheck) else list(mapped = sure))
}

#' Generate a Hybrid Regulatory Phosphosite Database
#'
#' This function generates a hybrid regulatory phosphosite database
#' by aligning mouse phosphosites to human counterparts.
#'
#' @param mh_alignment A dataframe of phosphosite alignments between mouse and human generated with [run_blast] function.
#' @param activatory Logical. If `TRUE`, considers only activatory sites.
#'
#' @seealso run_blast()
#' @return A dataframe containing hybrid regulatory phosphosites, mapped between species.
#'
#' @examples
#' # Example alignment dataset
#' alignment_data <- data.frame(
#'   q_ID = c("TP53-S15-mouse", "MYC-T22-mouse"),
#'   s_ID = c("TP53-S15-human", "MYC-T22-human"),
#'   UNIPROT = c("P04637", "P01106"),
#'   ACTIVATION = c('1', '-1')
#' )
#'
#' # Generate hybrid phosphosite database
#' hybrid_db <- generate_hybrid_db(alignment_data, activatory = TRUE)
#' head(hybrid_db)
#'
#' @export
generate_hybrid_db <- function(mh_alignment, activatory) {

  mh_alignment$UNIPROT <- NULL
  mh_alignment$ACTIVATION <- as.character(mh_alignment$ACTIVATION)
  # Select the regulatory phosphosite database
  reg_phos_db <- if (activatory) get(data("good_phos_df_human_act")) else get(data("good_phos_df_human_all"))

  hybrid_db <- mh_alignment %>%
    dplyr::inner_join(reg_phos_db, by = c("s_ID" = "PHOSPHO_KEY_GN_SEQ", 'ACTIVATION')) %>%
    dplyr::transmute(
      PHOSPHO_KEY_GN_SEQ = q_ID,
      PHOSPHO_KEY_GN_SEQ_human = s_ID,
      UNIPROT,
      ACTIVATION
    ) %>%
    dplyr::mutate(PHOSPHO_KEY_GN_SEQ = stringr::str_replace(PHOSPHO_KEY_GN_SEQ, "-[T,Y,S]-\\d+-mouse", "")) %>%
    dplyr::distinct()

  return(hybrid_db)
}

#' Map Experimental Phosphosites to Regulatory Phosphosites
#'
#' Maps experimentally identified phosphosites onto
#' known regulatory phosphosites for activity inference.
#'
#' @param phosphoproteomic_data A dataframe containing phosphoproteomics data.
#' @param organism Character string, either `"human"`, `"mouse"`, or `"hybrid"`, specifying the organism.
#' If  `"hybrid"` mouse regulatory information is remapped on human to enlarge coverage.
#' @param activatory Logical. If `TRUE`, considers only phosphosites regulating protein activity.
#' If `FALSE`, considers phosphosites regulating also protein abundance. Default: `TRUE`.
#' @param path_fasta A string specifying the path to the FASTA file for sequence mapping.
#' @param blastp_path Optional. A string specifying the Windows path of `blastp`.
#' @param local Logical. For **developmental purposes only**. Default: `FALSE`.
#' @param custom Logical. If `TRUE`, allows use of custom information on regulatory phosphosites.
#' @param custom_path A string specifying the path to a `.tsv` file containing custom  information on regulatory phosphosites.
#'
#' @return A list containing:
#' \describe{
#'   \item{`used_exp_data`}{A dataframe of experimental phosphosites that were successfully mapped to a regulatory phosphosite.}
#'   \item{`phosphoscore_df`}{A dataframe of inferred activity for each mapped phosphosite.}
#' }
#'
#' @examples
#' # Example phosphoproteomics dataset
#' phospho_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   sequence_window = c("XXXXXXXXSXXXXXXX", "XXXXXXXTXXXXXXXX"),
#'   aminoacid = c("S", "T"),
#'   position = c(15, 22),
#'   significant = c("+", "+"),
#'   difference = c(1.2, -0.8)
#' )
#'
#' # Map experimental phosphosites
#' mapped_results <- map_experimental_on_regulatory_phosphosites(phospho_data, "human", TRUE, "phospho.fasta")
#'
#' @export
map_experimental_on_regulatory_phosphosites <- function(phosphoproteomic_data,
                                                        organism,
                                                        activatory,
                                                        path_fasta,
                                                        blastp_path = NULL,
                                                        local = FALSE,
                                                        custom = FALSE,
                                                        custom_path = NULL) {

  message("** Mapping Experimental Phosphosites to Regulatory Phosphosites **")

  # Validate organism input
  organism <- match.arg(organism, c("human", "mouse", "hybrid"))

  # Load custom or predefined regulatory phosphosite database
  if (custom) {
    if (is.null(custom_path)) stop("Please provide a path to the custom regulatory phosphosites file.")
    message("Loading custom regulatory phosphosites...")
    reg_phos_db <- readr::read_tsv(custom_path, show_col_types = FALSE)
  } else {
    reg_phos_db <- switch(organism,
                          "human" = if (activatory) get(data('good_phos_df_human_act')) else  get(data('good_phos_df_human_all')),
                          "mouse" = if (activatory) get(data('good_phos_df_mouse_act')) else get(data('good_phos_df_mouse_all')),
                          "hybrid" = {
                            message("Generating hybrid regulatory phosphosite database...")
                            if (file.exists(path_fasta)) file.remove(path_fasta)
                            create_fasta(phosphoproteomic_data, path_fasta)
                            generate_hybrid_db(run_blast(path_fasta, blastp_path = blastp_path, local = local)$mapped, activatory)
                          }
    )
  }

  # Normalize sequence window length to 15-mer
  if(nchar(phosphoproteomic_data$sequence_window[1]) > 15){
    # if it is longer, subset the input sequence to 15-mer
    center <- (nchar(phosphoproteomic_data$sequence_window[1])+1)/2
    phosphoproteomic_data <- phosphoproteomic_data %>%
      dplyr::mutate(sequence_window_sub = stringr::str_sub(sequence_window, center - 7, center + 7))
  }else if(nchar(phosphoproteomic_data$sequence_window[1]) < 15){
    #if it is shorter, change the database
    offset <- (nchar(phosphoproteomic_data$sequence_window[1]) - 1) / 2
    reg_phos_db <- reg_phos_db %>%
      tidyr::separate(PHOSPHO_KEY_GN_SEQ, sep = '-', into = c('gene_name', 'seq')) %>%
      dplyr::mutate(seq = stringr::str_sub(seq, 7 - offset + 1, 7 + offset + 1)) %>%
      tidyr::unite('PHOSPHO_KEY_GN_SEQ', gene_name:seq, sep = '-')
    phosphoproteomic_data <- phosphoproteomic_data %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }else{
    phosphoproteomic_data <- phosphoproteomic_data %>%
      dplyr::rename(sequence_window_sub = sequence_window)
  }

  # Transform mouse gene_name in uppercase to match on human database
  if(organism == 'hybrid'){
    phosphoproteomic_data <- phosphoproteomic_data %>%
      mutate(gene_name = toupper(gene_name))
  }

  # Create mapping key and map experimental significantly modulated
  # phosphosites on regulatory database
  mapped_exp_data <- phosphoproteomic_data %>%
    dplyr::filter(significant == '+') %>%
    dplyr::mutate(PHOSPHO_KEY_GN_SEQ = paste0(gene_name, '-', sequence_window_sub)) %>%
    dplyr::filter(PHOSPHO_KEY_GN_SEQ %in% reg_phos_db$PHOSPHO_KEY_GN_SEQ)

  if (nrow(mapped_exp_data) == 0) {
    warning("No regulatory phosphosites found in the experimental dataset.")
    return(NULL)
  }

  # Compute for each phosphosite the regulatory power activity
  phosphoscore_df <- mapped_exp_data %>%
    dplyr::left_join(reg_phos_db, by = 'PHOSPHO_KEY_GN_SEQ') %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, ACTIVATION, difference) %>%
    dplyr::mutate(gene_name  = unlist(lapply(PHOSPHO_KEY_GN_SEQ,
                                             function(x){stringr::str_split(x, '-')[[1]][1]})),
                  inferred_activity = as.numeric(ACTIVATION) * as.numeric(difference)) %>%
    dplyr::distinct()

  return(list(used_exp_data = mapped_exp_data,
              phosphoscore_df = phosphoscore_df))
}

#' Compute PhosphoScore using Hybrid Reference Database
#'
#' This function computes the **PhosphoScore** for experimental phosphoproteomics data
#' by mapping mouse phoshosites on both mouse and human regulatory phosphosites
#' generating an hybrid functional annotation on phosphosites.
#' Mapping is performed using mouse and human 15-mers aligned with blastp.
#'
#' @param phosphoproteomic_data A dataframe containing phosphoproteomics data.
#' @param organism A string specifying the organism (`"human"`, `"mouse"`, or `"hybrid"`).
#' @param activatory Logical. If `TRUE`, considers only activatory phosphosites.
#' @param blastp_path Optional. A string specifying the path to the `blastp` executable (Windows users).
#' @param path_fasta A string specifying the path to the FASTA file for sequence mapping. Default: `"./phospho.fasta"`.
#' @param local Logical. If `TRUE`, runs locally (for development purposes). Default: `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{`used_exp_data`}{A dataframe of experimental phosphosites that were successfully mapped to a regulatory phosphosite.}
#'   \item{`phosphoscore_df`}{A dataframe of inferred activity for each mapped phosphosite, integrating mouse and hybrid phosphosites.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example phosphoproteomics dataset
#' phospho_data <- data.frame(
#'   gene_name = c("TP53", "MYC"),
#'   sequence_window = c("XXXXXXXXSXXXXXXX", "XXXXXXXTXXXXXXXX"),
#'   aminoacid = c("S", "T"),
#'   position = c(15, 22),
#'   significant = c("+", "+"),
#'   difference = c(1.2, -0.8)
#' )
#'
#' # Compute PhosphoScore for hybrid organism
#' phosphoscore_results <- phospho_score_hybrid_computation(phospho_data, "hybrid", TRUE, blastp_path = "path/to/blastp.exe")
#'}
phospho_score_hybrid_computation <- function(phosphoproteomic_data,
                                             organism,
                                             activatory,
                                             blastp_path = NULL,
                                             path_fasta = "./phospho.fasta",
                                             local = FALSE) {

  message("** Running PhosphoScore Hybrid Computation **")

  # Validate organism input
  organism <- match.arg(organism, c("human", "mouse", "hybrid"))

  # Compute phosphoscore for mouse-only mapping
  phosphoscore_mouse <- map_experimental_on_regulatory_phosphosites(
    phosphoproteomic_data = phosphoproteomic_data,
    organism = "mouse",
    activatory = activatory,
    path_fasta = path_fasta,
    local = local
  )

  # Compute phosphoscore for hybrid (mouse-to-human) mapping
  phosphoscore_hybrid <- map_experimental_on_regulatory_phosphosites(
    phosphoproteomic_data = phosphoproteomic_data,
    organism = "hybrid",
    activatory = activatory,
    path_fasta = path_fasta,
    blastp_path = blastp_path,
    local = local
  )

  # Handle cases where only one mapping is available
  if (is.null(phosphoscore_mouse) & is.null(phosphoscore_hybrid)) {
    stop("No regulatory phosphosites found in the experimental dataset.")
  } else if (is.null(phosphoscore_hybrid)) {
    warning("No mapped mouse phosphosites on human, proceeding with mouse-only analysis.")
    return(phosphoscore_mouse)
  } else if (is.null(phosphoscore_mouse)) {
    warning("No mouse phosphosites found, proceeding with human-mapped phosphosites.")
    return(phosphoscore_hybrid)
  }

  # Select relevant columns from both mappings
  phosphoscore_mouse_df <- phosphoscore_mouse$phosphoscore_df %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, inferred_activity, gene_name, ACTIVATION, difference)

  # Capitalize gene names to meet mouse Gene Symbol notation
  phosphoscore_hybrid_df <- phosphoscore_hybrid$phosphoscore_df %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, inferred_activity, gene_name, ACTIVATION, difference) %>%
    dplyr::mutate(gene_name = stringr::str_to_title(gene_name))

  message("Merging Mouse and Hybrid PhosphoScore Data...")

  # Merge mouse and hybrid phosphoscore results
  merged_phosphoscore <- dplyr::full_join(
    phosphoscore_mouse_df, phosphoscore_hybrid_df,
    by = c("PHOSPHO_KEY_GN_SEQ", "gene_name", "ACTIVATION", "difference"),
    suffix = c(".mouse", ".hybrid")
  ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(PHOSPHO_KEY_GN_SEQ)

  # Determine phosphosite origin
  merged_phosphoscore <- merged_phosphoscore %>%
    dplyr::mutate(
      source_org = dplyr::case_when(
        !is.na(inferred_activity.mouse) & !is.na(inferred_activity.hybrid) ~ "both",
        is.na(inferred_activity.mouse) ~ "human",
        is.na(inferred_activity.hybrid) ~ "mouse"
      ),
      inferred_activity = dplyr::coalesce(inferred_activity.mouse, inferred_activity.hybrid)
    ) %>%
    dplyr::select(PHOSPHO_KEY_GN_SEQ, gene_name, source_org, inferred_activity, ACTIVATION, difference)

  # Combine experimental data from both mappings
  used_exp_data <- dplyr::bind_rows(phosphoscore_mouse$used_exp_data,
                                    phosphoscore_hybrid$used_exp_data) %>%
    dplyr::distinct()

  return(list(
    used_exp_data = used_exp_data,
    phosphoscore_df = merged_phosphoscore
  ))
}



