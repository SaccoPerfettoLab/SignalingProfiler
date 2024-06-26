% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PKN_information.R
\name{create_PKN}
\alias{create_PKN}
\title{create_PKN}
\usage{
create_PKN(
  database = c("SIGNOR", "Omnipath", "PsP", "SerThr_Atlas"),
  organism = "human",
  direct = FALSE,
  file_path = NULL,
  omnipath_resources = c("SignaLink3", "PhosphoSite", "SIGNOR"),
  psp_reg_site_path = NULL,
  psp_kin_sub_path = NULL
)
}
\arguments{
\item{database}{vector, resources can be SIGNOR, Omnipath, PsP, SerThrAtlas}

\item{organism}{string, 'human' or 'mouse'}

\item{direct}{Boolean, if TRUE only direct interactions are included}

\item{file_path}{string, path for file saving}

\item{omnipath_resources}{vector, resources included in OmniPath database (see OmniPath tutorial for the complete list)}

\item{reg_site_path}{path to the PhosphoSitePlus 'Regulatory Sites' file}

\item{kin_sub_path}{path to the PhosphoSitePlus 'Kinases-Substrates' file}
}
\value{
table of interactions derived from the user-specified resources
}
\description{
create_PKN
}
