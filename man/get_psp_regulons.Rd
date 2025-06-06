% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regulons_creation.R
\name{get_psp_regulons}
\alias{get_psp_regulons}
\title{Retrieve PhosphoSitePlus Regulons}
\usage{
get_psp_regulons(reg_site_path, kin_sub_path, organism, with_atlas = FALSE)
}
\arguments{
\item{reg_site_path}{Character. Path to the PhosphoSitePlus "Regulatory_sites" file.}

\item{kin_sub_path}{Character. Path to the PhosphoSitePlus  "Kinases_Substrate_Dataset" file.}

\item{organism}{Character string, either `"human"` or `"mouse"`, specifying the organism.}

\item{with_atlas}{Logical. If `TRUE`, integrates Ser/Thr/Tyr Kinome Atlas interactions (human only). Default: `FALSE`.}
}
\value{
A tibble containing the kinase-substrate relationships data for `"ksea"`
with three columns:`tf` (kinase), `target` (substrate in `Gene-Residue-Position` format), and `mor` (mode of regulation).
`mor` is -1, 1, and is continous from 0 to 1 for Ser/Thr Kinome Atlas relationships.
}
\description{
This function extracts kinase-substrate interactions from the PsP user-provided files:
Kinase_Substrate_Dataset and Regulatory_sites.
}
\examples{
\dontrun{
psp_db <- psp_parsing(reg_site_path = './Regulatory_sites',
                      kin_sub_path = './Kinase_Substrate_Dataset',
                      organism = 'human',
                      with_atlas = FALSE,
                      local = FALSE)
}

}
