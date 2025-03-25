#' Annotate Proteins with GO Molecular Function
#'
#' Annotates a table of proteins with Gene Ontology (GO) Molecular Function
#' categories based on a static reference table (`gomf_human` or `gomf_mouse`).
#' The annotation is performed using the `gene_name` column as the key.
#'
#' @param inferred_proteins_dataset Data frame, protein table containing a `gene_name` column.
#' @param organism Character, organism type. Must be one of `"human"` or `"mouse"`.
#'
#' @return A data frame identical to `inferred_proteins_dataset`, but with an
#' additional `mf` column containing GO Molecular Function annotations:
#'   - `"rec"` if the protein is user-defined.
#'   - `"kin"` if the protein is a kinase.
#'   - `"phos"` if the protein is a phosphatase.
#'   - `"tf"` if the protein is a transcription factor.
#'   - `"other"` for proteins without a matching annotation.
#'   - `"complex"` for proteins labeled with `SIGNOR-C`.
#'   - `"fusion protein"` for proteins labeled with `SIGNOR-F`.
#'
#' @details
#' The function joins the input dataset with preloaded molecular function
#' tables (`gomf_human` or `gomf_mouse`) using the `gene_name` column.
#' If an `mf` column is already present, it is removed before annotation.
#'
#' @export
#'
#' @examples
#' data('toy_prot_activity_df')
#' molecular_function_annotation(
#'   inferred_proteins_dataset = toy_prot_activity_df,
#'   organism = "human")
#'
molecular_function_annotation <-
  function(inferred_proteins_dataset, organism) {

    ## remove existing 'mf' column if present
    if ('mf' %in% colnames(inferred_proteins_dataset)) {
      inferred_proteins_dataset$mf <- NULL
    }

    ## select the correct molecular function table based on organism
    if (organism == 'human') {
      inferred_proteins_dataset <-
        dplyr::left_join(
          inferred_proteins_dataset,
          get(data('gomf_human')) %>% dplyr::select(gene_name, mf),
          by = c('gene_name')
        )

    } else if (organism == 'mouse' | organism == 'hybrid') {
      inferred_proteins_dataset <-
        dplyr::left_join(
          inferred_proteins_dataset,
          get(data('gomf_mouse')) %>% dplyr::select(gene_name, mf),
          by = c('gene_name')
        )

    } else{
      stop('please provide a valid organism')
    }

    ## manually annotate complex and fusion proteins
    inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf) &
                                   grepl('SIGNOR-C', inferred_proteins_dataset$UNIPROT)] <- 'complex'
    inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf) &
                                   grepl('SIGNOR-F', inferred_proteins_dataset$UNIPROT)] <-
      'fusion protein'

    ## assign 'other' to missing annotations
    inferred_proteins_dataset$mf[is.na(inferred_proteins_dataset$mf)] <-
      'other'

    ## assign 'rec' if the method is user-defined
    if ("method" %in% colnames(inferred_proteins_dataset)) {
      inferred_proteins_dataset$mf[inferred_proteins_dataset$method == 'user'] <-
        'rec'
    }

    return(inferred_proteins_dataset)
  }

