% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{proxpath_phenotypes}
\alias{proxpath_phenotypes}
\title{Summary of Phenotypes in ProxPath Results}
\format{
A data frame with the following columns:
\describe{
  \item{phenotype}{Character. Name of the phenotype or cellular process (originally from the EndPathways column).}
  \item{n_paths}{Integer. Number of signaling paths from SIGNOR proteins to the phenotype.}
}
}
\usage{
proxpath_phenotypes
}
\description{
A summary table listing all phenotypes (`EndPathways`) found in the ProxPath output,
along with the number of distinct signaling paths leading to each phenotype.
}
\details{
This object is derived from the `phenoscore_distances_table`. A data frame hosted on PerfettoLab server
annotating proximity scores computed between SIGNOR signaling proteins and phenotypic endpoints.
This table is the output of the ProxPath algorithm (38102483), which retrieves all directed paths from each SIGNOR protein
(QueryNode) to known phenotypes (EndNode) in the SIGNOR database. These paths are then scored based on edge weights
and annotated with their *final regulatory effect*.
The dataset is generated by running a Python script (`proxpath.py`) available upon request.
The script uses the SIGNOR network to extract relevant paths and exports the result to a minimized output file
which is then loaded into R and stored in the package.
by counting how many
times each phenotype appears as the target node of a signaling path. It can be used
to prioritize or filter phenotypes based on their connectivity within the SIGNOR
network.


This summary allows users to explore the coverage of phenotype annotations in the dataset,
and can guide selection of phenotypes for enrichment or proximity scoring analyses.
}
\examples{
data(proxpath_phenotypes)
head(proxpath_phenotypes)
}
\keyword{datasets}
