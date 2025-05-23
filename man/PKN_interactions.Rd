% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PKN_interactions}
\alias{PKN_interactions}
\alias{PKN_human_dir}
\alias{PKN_human_ind}
\alias{PKN_human_atlas_dir}
\alias{PKN_human_atlas_ind}
\alias{PKN_mouse_dir}
\alias{PKN_mouse_ind}
\title{Prior Knowledge Network (PKN) Data Objects}
\format{
Each `PKN_*` object is a data frame with variable rows (depending on SignalingProfiler version) and 15 columns:
\describe{
  \item{ENTITYA}{Character, Source Biological Entity Name in uppercase with special characters as underscores}
  \item{ENTITYB}{Character, Target Biological Entity Name in uppercase with special characters as underscores}
  \item{INTERACTION}{Character, Sign of the Interaction, 1 Activating, -1 Inhibiting}
  \item{IDA}{Character, Source Biological Entity Identifier}
  \item{IDB}{Character, Target Biological Entity Identifier}
  \item{MECHANISM}{Character, Molecular Mechanism of the Interaction}
  \item{RESIDUE}{Character, if the `MECHANISM` is a phosphorylation, the modified residue in the interaction}
  \item{SEQUENCE}{Character, if the `MECHANISM` is a phosphorylation, a 15mer centered on the modified residue}
  \item{TYPEA}{Character, Source Biological Entity Molecular Type}
  \item{DATABASEA}{Character, Source Biological Entity Database}
  \item{TYPEB}{Character, Target Biological Entity Molecular Type}
  \item{DATABASEB}{Character, Target Biological Entity Database}
  \item{PHOSPHO_KEY_GN_SEQ}{Character composed of `ENTITYB` and `SEQUENCE` column}
  \item{SOURCE}{Character, Source Database for the interaction; either SIGNOR or PhosphoSitePlus or both}
  \item{DIRECT}{Logical, whether the interaction is direct}
}

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 34793 rows and 15 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 35624 rows and 15 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 62513 rows and 15 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 63344 rows and 15 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 26021 rows and 15 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 26793 rows and 15 columns.
}
\source{
- <https://signor.uniroma2.it>
- <https://www.phosphosite.org>
- <https://pubmed.ncbi.nlm.nih.gov/36631611/> (for Ser/Thr kinase-substrate interactions)
- <https://pubmed.ncbi.nlm.nih.gov/38720073/> (for Tyr kinase-substrate interactions)
}
\usage{
PKN_human_dir

PKN_human_ind

PKN_human_atlas_dir

PKN_human_atlas_ind

PKN_mouse_dir

PKN_mouse_ind
}
\description{
The Prior Knowledge Network (PKN) consists of signed and directed (causal)
interactomes that integrate curated signaling interactions from SIGNOR and PhosphoSitePlus.
Some versions also include experimentally derived kinase-substrate relationships
from Serine, Threonine, and Tyrosine Kinome Atlases (PMIDs: 36631611, 38720073).
}
\details{
Available datasets:
- `PKN_human_dir`: Human Prior Knowledge Network with *direct* interactions.
- `PKN_human_ind`: Human Prior Knowledge Network with *direct* and *indirect* interactions.
- `PKN_human_atlas_dir`: Human PKN with *direct* interactions and experimentally-derived kinase-substrate relationships.
- `PKN_human_atlas_ind`: Human PKN with *direct* and *indirect* interactions and experimentally-derived kinase-substrate relationships.
- `PKN_mouse_dir`: Mouse Prior Knowledge Network with *direct* interactions.
- `PKN_mouse_ind`: Mouse Prior Knowledge Network with *direct* and *indirect* interactions.
}
\keyword{datasets}
