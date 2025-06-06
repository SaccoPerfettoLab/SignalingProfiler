% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tfea_ksea_datasets}
\alias{tfea_ksea_datasets}
\alias{tfea_db_human}
\alias{tfea_db_human_collectri}
\alias{tfea_db_mouse}
\alias{tfea_db_mouse_collectri}
\alias{ksea_db_human}
\alias{ksea_db_human_atlas}
\alias{ksea_db_mouse}
\title{Transcription Factor Enrichment Analysis (TFEA) and Kinase-Substrate Enrichment Analysis (KSEA) Databases}
\format{
- **TFEA datasets** (`tfea_db_human`, `tfea_db_human_collectri`, `tfea_db_mouse`, `tfea_db_mouse_collectri`):
  Data frames with the following columns:
  \describe{
    \item{tf}{Character, Transcription Factor}
    \item{target}{Character, Target Protein Regulated by the TF}
    \item{mor}{Character, Mode of Regulation (Activation/Inhibition)}
    \item{confidence}{Character, Confidence Score for Interaction (Present only in `tfea_db_*` datasets)}
  }

- **KSEA datasets** (`ksea_db_human`, `ksea_db_human_atlas`, `ksea_db_mouse`):
  Data frames with the following columns:
  \describe{
    \item{tf}{Character, Kinase or Phosphatase}
    \item{target}{Character, Modified Target Protein}
    \item{mor}{Character, Mode of Regulation (Activation/Inhibition); for atlas, it is a continous value representing the probability of phosphorylation of the target}
  }

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 8972 rows and 4 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 44317 rows and 4 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 8340 rows and 4 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 8340 rows and 4 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 27850 rows and 3 columns.

An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 22317 rows and 3 columns.
}
\source{
- <https://signor.uniroma2.it>
- <https://www.phosphosite.org>
- <https://collectri.bioinf.uni-leipzig.de> (for CoLLECTRI-based datasets)
- <https://pubmed.ncbi.nlm.nih.gov/36631611/>
- <https://pubmed.ncbi.nlm.nih.gov/38720073/>
- <https://omnipathdb.org/>
}
\usage{
tfea_db_human

tfea_db_human_collectri

tfea_db_mouse

tfea_db_mouse_collectri

ksea_db_human

ksea_db_mouse
}
\description{
Processed databases for transcription factor enrichment analysis (TFEA) and kinase-substrate enrichment analysis (KSEA)
containing curated regulatory interactions between transcription factors (TFs) or kinases/phosphatases and their targets.
}
\details{
TFEA datasets optionally include confidence scores and are derived from SIGNOR, Dorothea, and CoLLECTRI.
KSEA datasets integrate interactions from SIGNOR, OmniPath, PhosphoSitePlus, and kinase-substrate atlases.


- `tfea_db_human`: A table of human transcription factor-target interactions from SIGNOR and Dorothea (confidence A).
- `tfea_db_human_collectri`: A curated table of human transcription factor-target interactions from SIGNOR and Collectri.
- `tfea_db_mouse`: A curated table of mouse transcription factor-target interactions from SIGNOR and Dorothea (confidence A)..
- `tfea_db_mouse_collectri`:  A curated table of mouse transcription factor-target interactions from SIGNOR and Collectri.
- `ksea_db_human`: A table of human kinase/phosphatase-substrate interactions from SIGNOR, Omnipath and PhosphoSitePlus.
- `ksea_db_human_atlas`: An extended human kinase-substrate dataset with additional interactions from Serine/Threonine and Tyrosine Kinome Atlases (PMIDs: 36631611, 38720073).
- `ksea_db_mouse`: A table of mouse kinase/phosphatase-substrate interactions from SIGNOR, Omnipath and PhosphoSitePlus.
}
\keyword{datasets}
