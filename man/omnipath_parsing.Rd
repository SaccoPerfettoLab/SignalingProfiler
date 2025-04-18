% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_databases.R
\name{omnipath_parsing}
\alias{omnipath_parsing}
\title{Parse OmniPath Interaction Data and Map to SIGNOR Format}
\usage{
omnipath_parsing(resources, file_path = NULL)
}
\arguments{
\item{resources}{Character vector. List of OmniPath resources to retrieve interactions from.}

\item{file_path}{Character. Optional. Path to save the output file. If `NULL`, results are not saved. Default: `NULL`.}
}
\value{
A `data.frame` containing OmniPath-derived protein interactions, mapped to the SIGNOR format.
}
\description{
This function retrieves, processes, and maps protein-protein interactions from OmniPath
to the SIGNOR format, handling 'protein complexes' misaligned notation and integrating
them into a Prior Knowledge Network (PKN).
}
\details{
The function:
- **Retrieves interactions** from OmniPath using `OmnipathR::omnipath_interactions()`.
- **Filters and formats interactions** into a SIGNOR-like structure.
- **Handles molecular complexes**, mapping them using SIGNOR internal data.
- **Resolves complex interactions** into individual protein-protein relationships.

If `file_path` is provided, the resulting table is saved as a `.tsv` file.
}
\examples{
# Retrieve and process OmniPath interactions
omni_data <- omnipath_parsing(resources = c("SIGNOR", "SignaLink3"))

}
\seealso{
[OmnipathR::omnipath_interactions]
}
