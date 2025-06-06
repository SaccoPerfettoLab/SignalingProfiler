% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_and_map_edges.R
\name{expand_and_map_edges}
\alias{expand_and_map_edges}
\title{Expand and Map Edges with Experimental Phosphoproteomics Data}
\usage{
expand_and_map_edges(
  optimized_object,
  organism,
  phospho_df,
  files,
  with_atlas = FALSE,
  direct = FALSE,
  path_sif,
  path_rds
)
}
\arguments{
\item{optimized_object}{A list containing the optimized CARNIVAL network, including:
\itemize{
  \item `nodes_df`: A dataframe with node attributes.
  \item `edges_df`: A dataframe with regulatory interactions.
  \item `igraph_network`: An igraph object representing the network.
}}

\item{organism}{Character string, either `"human"` or `"mouse"`, specifying the organism used.}

\item{phospho_df}{A tibble containing phosphoproteomics data with the following required columns:
\itemize{
  \item `gene_name` - Gene symbol of the phosphoprotein.
  \item `sequence_window` - The phosphopeptide sequence window surrounding the phosphorylated site.
  \item `aminoacid` - The phosphorylated amino acid (e.g., `"S"`, `"T"`, `"Y"`).
  \item `position` - The numeric position of the phosphosite in the protein sequence.
  \item `difference` - The fold-change in phosphorylation between experimental conditions.
  \item `significant` - Logical, indicating whether the phosphosite is significantly modulated.
}}

\item{files}{Logical (`TRUE` by default). If `TRUE`, saves the resulting expanded network to disk;
if `FALSE`, the function runs without writing files.}

\item{with_atlas}{Logical (`TRUE` by default). If `TRUE`, includes the Kinome Atlas-derived regulons
in the Prior Knowledge Network (PKN); if `FALSE`, excludes them.}

\item{direct}{Logical (`FALSE` by default). If `TRUE`, uses only direct interactions in the Prior
Knowledge Network; if `FALSE`, includes both direct and indirect interactions.}

\item{path_sif}{Character string specifying the output file path for the optimized network in
**SIF format**.}

\item{path_rds}{Character string specifying the output file path for the optimized **RDS object**.}
}
\value{
A list (*SP_object*) containing:
\item{igraph_network}{An igraph object representing the expanded network, incorporating phosphoproteomics data.}
\item{nodes_df}{A dataframe with node attributes, including regulatory activity scores and mapped metadata.}
\item{edges_df}{A dataframe with edge attributes, including:}
  \itemize{
    \item `source` - The source node of the regulatory interaction.
    \item `target` - The target node of the regulatory interaction.
    \item `sign` - The interaction sign (`+1` or `-1`).
    \item `carnival_weight` - The weight assigned to the interaction by CARNIVAL.
    \item `mechanism` - The molecular mechanism of the interaction annotated in SIGNOR.
    \item `direct` - Logical, indicating whether the interaction is direct.
    \item `residue` - The putative phosphosite modulated in the interaction.
    \item `is_quantified` - Logical, indicating whether the phosphosite was quantified in the dataset.
    \item `is_significant` - Logical, indicating whether the phosphosite is significantly modulated.
    \item `FC` - Character, fold-change in phosphorylation (starred values are significant).
  }
}
\description{
This function integrates CARNIVAL-optimized edges with experimental phosphoproteomics data.
It maps regulatory interactions to phosphorylation events incorporating quantitative phosphoproteomics information.
}
\examples{
# Load an optimized CARNIVAL network
data('toy_opt_network')
toy_opt_network$edges_df <- toy_opt_network$edges_df[, c('source', 'target', 'sign', 'carnival_weight')]

# Load phosphoproteomics dataset
data('phospho_toy_df')

expand_and_map_edges(
  optimized_object = toy_opt_network,
  organism = "human",
  phospho_df = phospho_toy_df,
  files = TRUE,
  with_atlas = TRUE,
  direct = FALSE,
  path_sif = "expanded_network.sif",
  path_rds = "expanded_network.rds"
)

}
\seealso{
[add_output_carnival_nodes_attributes()], [add_output_carnival_edges_attributes()], [molecular_function_annotation()]
}
