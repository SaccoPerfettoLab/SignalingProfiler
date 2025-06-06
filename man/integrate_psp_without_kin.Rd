% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reg_phosphosites_creation.R
\name{integrate_psp_without_kin}
\alias{integrate_psp_without_kin}
\title{Retrieve regulatory phosphosites without upstream kinases}
\usage{
integrate_psp_without_kin(
  organism,
  reg_site_path,
  psp_db,
  only_activatory = TRUE,
  local = FALSE
)
}
\arguments{
\item{organism}{Character string specifying the organism. Supported values are `"human"` and `"mouse"`.}

\item{reg_site_path}{Character string specifying the file path to the PhosphoSitePlus 'Regulatory Sites' data.
Required if `"PsP"` is included in `resources`.}

\item{psp_db}{PhosphoSitePlus database returned from the `psp_parsing` function.}

\item{only_activatory}{Logical, indicating whether to filter for only phosphosites that regulate protein activity. Default is `TRUE`.}

\item{local}{Logical, indicating whether to use a local or package-installed version of the SIGNOR dictionary. Default is `FALSE`.}
}
\value{
A tibble containing regulatory phosphosites with the following columns:
- `ENTITY`: Protein gene symbol.
- `IDB`: UniProt accession ID.
- `RESIDUE`: Modified residue (aminoacid and position).
- `INTERACTION`: Regulatory effect (`1` for activation, `-1` for inhibition).
}
\description{
This function keeps only phosphosites without an upstream kinase, hence
absent in PhosphoSitePlus prior knowledge network.
}
