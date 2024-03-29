% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/protein_activity_prediction.R
\name{run_blast}
\alias{run_blast}
\title{Title}
\usage{
run_blast(
  path_experimental_fasta_file,
  all = FALSE,
  blastp_path = NULL,
  local = FALSE
)
}
\arguments{
\item{path_experimental_fasta_file}{string representing path of
phosphoproteomics data in fasta file}

\item{all}{boolean value representing if you want all alingments
or only the ones with same gene_name
#' @param local FOR DEVELOPMENTAL PURPOSES TO DELETE}

\item{blastp_path}{string specifying the blastp path on Windows}

\item{local}{boolean TRUE or FALSE, for developing purposes}
}
\value{
a list containing:
a dataframe of the sure mouse phosphopeptides aligned on human
a dataframe of the unsure human phosphopeptides aligned on human having
different gene_name, to check manually!
}
\description{
Title
}
