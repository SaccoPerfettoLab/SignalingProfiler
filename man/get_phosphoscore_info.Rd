% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PKN_information.R
\name{get_phosphoscore_info}
\alias{get_phosphoscore_info}
\title{get_phosphoscore_info}
\usage{
get_phosphoscore_info(
  resources = c("SIGNOR", "PsP"),
  organism,
  psp_reg_site_path = NULL,
  only_activatory = TRUE,
  local = FALSE
)
}
\arguments{
\item{resources}{vector, source for regulatory phosphosites information}

\item{organism}{string, 'human' or 'mouse'}

\item{psp_reg_site_path}{path to the PhosphoSitePlus 'Regulatory Sites' file}

\item{only_activatory}{Boolean, if TRUE, only interactions regulating protein activity are included}
}
\description{
get_phosphoscore_info
}
