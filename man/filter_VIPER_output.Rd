% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footprint_prediction.R
\name{filter_VIPER_output}
\alias{filter_VIPER_output}
\title{Filter VIPER Output for Molecular Function}
\usage{
filter_VIPER_output(inferred_proteins_mf, analysis)
}
\arguments{
\item{inferred_proteins_mf}{A dataframe containing inferred proteins with molecular function annotations.}

\item{analysis}{A string specifying the analysis type ("tfea" or "ksea").}
}
\value{
A dataframe of inferred proteins filtered based on molecular function.
}
\description{
This function filters inferred proteins based on their GO molecular function annotated
with [molecular_function_annotation] function of SignalingProfiler:
- Transcription Factors for TFEA;
- Kinases/Phosphatases for KSEA.
}
\examples{
# Filter transcription factors in human
inferred_proteins <- data.frame(gene_name = c('TP53', 'CDK5'),
                                mf = c('tf', 'kin'))
filtered_tfs <- filter_VIPER_output(inferred_proteins_mf = inferred_proteins,
                                     analysis = "tfea")

}
\seealso{
[molecular_function_annotation]
}
