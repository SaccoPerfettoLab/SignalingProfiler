% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phosphoscore.R
\name{phosphoscore_computation_aapos}
\alias{phosphoscore_computation_aapos}
\title{Infer protein actiity from phosphosite regulatory power (PhophoScore)
considering residue and position as mapping key}
\usage{
phosphoscore_computation_aapos(
  phosphoproteomic_data,
  organism,
  activatory,
  GO_annotation = FALSE,
  custom = FALSE,
  custom_path = NULL
)
}
\arguments{
\item{phosphoproteomic_data}{A dataframe containing phosphoproteomics data.}

\item{organism}{A string specifying the organism (`"human"` or `"mouse"`).}

\item{activatory}{Logical. If `TRUE`, considers only activatory sites.}

\item{GO_annotation}{Logical. If `TRUE`, performs Gene Ontology (GO) annotation.}

\item{custom}{Logical. If `TRUE`, uses a custom regulatory database.}

\item{custom_path}{A string specifying the path to a `.tsv` file containing custom phosphosites.}
}
\value{
A dataframe containing inferred protein activities using amino acid and position-based mapping.
}
\description{
The principle is the same as [compute_phosphoscore] but this function
should be used when in phosphoproteomics data no phosphopeptide is provided.
}
\examples{
phospho_data <- data.frame(
  gene_name = c("TP53", "MYC", "EGFR"),
  aminoacid = c("S", "T", "Y"),
  position = c(15, 22, 45),
  significant = c("+", "+", NA),
  difference = c(1.2, -0.8, 0.5)
)

phosphoscore_aapos_results <- phosphoscore_computation_aapos(phospho_data, "human", TRUE)

}
