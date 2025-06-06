% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footprint_prediction.R
\name{create_matrix_from_VIPER_format}
\alias{create_matrix_from_VIPER_format}
\title{Convert VIPER Format to a Matrix}
\usage{
create_matrix_from_VIPER_format(viper_format)
}
\arguments{
\item{viper_format}{A dataframe formatted for VIPER.}
}
\value{
A matrix where rows correspond to gene IDs and columns contain t-statistic values.
}
\description{
This function converts a VIPER-formatted dataframe into a matrix
containing ID and t-statistic values.
}
\examples{
# Example dataset with transcriptomics data
omic_data <- data.frame(
  gene_name = c("TP53", "MYC", "EGFR"),
  difference = c(1.5, -2.1, 0.8),
  logpval = c(0.01, 0.05, 0.2),
  significant = c("+", "+", NA)
)
viper_format <- create_viper_format(omic_data, analysis = "tfea", significance = TRUE)
viper_matrix <- create_matrix_from_VIPER_format(viper_format)

}
