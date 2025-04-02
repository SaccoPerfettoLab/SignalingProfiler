# ---------------------------------------------------------------------- #
#                 Prior Knowledge Network (PKN) Data Objects            #
# ---------------------------------------------------------------------- #

#' Prior Knowledge Network (PKN) Data Objects
#'
#' The Prior Knowledge Network (PKN) consists of signed and directed (causal)
#' interactomes that integrate curated signaling interactions from SIGNOR and PhosphoSitePlus.
#' Some versions also include experimentally derived kinase-substrate relationships
#' from Serine, Threonine, and Tyrosine Kinome Atlases (PMIDs: 36631611, 38720073).
#'
#' @format Each `PKN_*` object is a data frame with variable rows (depending on SignalingProfiler version) and 15 columns:
#' \describe{
#'   \item{ENTITYA}{Character, Source Biological Entity Name in uppercase with special characters as underscores}
#'   \item{ENTITYB}{Character, Target Biological Entity Name in uppercase with special characters as underscores}
#'   \item{INTERACTION}{Character, Sign of the Interaction, 1 Activating, -1 Inhibiting}
#'   \item{IDA}{Character, Source Biological Entity Identifier}
#'   \item{IDB}{Character, Target Biological Entity Identifier}
#'   \item{MECHANISM}{Character, Molecular Mechanism of the Interaction}
#'   \item{RESIDUE}{Character, if the `MECHANISM` is a phosphorylation, the modified residue in the interaction}
#'   \item{SEQUENCE}{Character, if the `MECHANISM` is a phosphorylation, a 15mer centered on the modified residue}
#'   \item{TYPEA}{Character, Source Biological Entity Molecular Type}
#'   \item{DATABASEA}{Character, Source Biological Entity Database}
#'   \item{TYPEB}{Character, Target Biological Entity Molecular Type}
#'   \item{DATABASEB}{Character, Target Biological Entity Database}
#'   \item{PHOSPHO_KEY_GN_SEQ}{Character composed of `ENTITYB` and `SEQUENCE` column}
#'   \item{SOURCE}{Character, Source Database for the interaction; either SIGNOR or PhosphoSitePlus or both}
#'   \item{DIRECT}{Logical, whether the interaction is direct}
#' }
#'
#' @details Available datasets:
#' - `PKN_human_dir`: Human Prior Knowledge Network with *direct* interactions.
#' - `PKN_human_ind`: Human Prior Knowledge Network with *direct* and *indirect* interactions.
#' - `PKN_human_atlas_dir`: Human PKN with *direct* interactions and experimentally-derived kinase-substrate relationships.
#' - `PKN_human_atlas_ind`: Human PKN with *direct* and *indirect* interactions and experimentally-derived kinase-substrate relationships.
#' - `PKN_mouse_dir`: Mouse Prior Knowledge Network with *direct* interactions.
#' - `PKN_mouse_ind`: Mouse Prior Knowledge Network with *direct* and *indirect* interactions.
#'
#' @source
#' - <https://signor.uniroma2.it>
#' - <https://www.phosphosite.org>
#' - <https://pubmed.ncbi.nlm.nih.gov/36631611/> (for Ser/Thr kinase-substrate interactions)
#' - <https://pubmed.ncbi.nlm.nih.gov/38720073/> (for Tyr kinase-substrate interactions)
#'
#' @name PKN_interactions
#' @aliases PKN_human_dir PKN_human_ind PKN_human_atlas_dir PKN_human_atlas_ind
#'          PKN_mouse_dir PKN_mouse_ind
#'
#' @keywords datasets
NULL

#' @rdname PKN_interactions
"PKN_human_dir"

#' @rdname PKN_interactions
"PKN_human_ind"

#' @rdname PKN_interactions
"PKN_human_atlas_dir"

#' @rdname PKN_interactions
"PKN_human_atlas_ind"

#' @rdname PKN_interactions
"PKN_mouse_dir"

#' @rdname PKN_interactions
"PKN_mouse_ind"

# ---------------------------------------------------------------------- #
#                 Prior Knowledge Network Molecular Entities             #
# ---------------------------------------------------------------------- #

#' Molecular Entities in the Prior Knowledge Network (PKN)
#'
#' Molecular entities connected in the Prior Knowledge Network (PKN) encompassing
#' all entities involved in *direct* and *indirect* interactions.
#' All `ENTITY` values are formatted consistently: uppercase for human,
#' and start-case for mouse protein names. Complexes and other `TYPE` are
#' always uppercase with special characters converted as underscores.
#'
#' @format Each `PKN_proteins_*` object is a data frame with variable rows (depending on SignalingProfiler version) and 4 columns:
#' \describe{
#'   \item{ENTITY}{Character, Biological Entity Name}
#'   \item{ID}{Character, Biological Entity Identifier, from SIGNOR or UNIPROT}
#'   \item{TYPE}{Character, Biological Entity Molecular Type}
#'   \item{DATABASE}{Character, Database Source for Biological Entity Identifier}
#' }
#'
#' @details Available datasets:
#' - `PKN_proteins_human`: Molecular entities in the *human* Prior Knowledge Network.
#' - `PKN_proteins_human_atlas`: *Human* PKN entities, including kinase-substrate relationships.
#' - `PKN_proteins_mouse`: Molecular entities in the *mouse* Prior Knowledge Network.
#'
#' @source
#' - <https://signor.uniroma2.it>
#' - <https://www.phosphosite.org>
#' - <https://pubmed.ncbi.nlm.nih.gov/36631611/> (for Ser/Thr kinase-substrate interactions)
#' - <https://pubmed.ncbi.nlm.nih.gov/38720073/> (for Tyr kinase-substrate interactions)
#'
#' @name PKN_proteins
#' @aliases PKN_proteins_human PKN_proteins_human_atlas PKN_proteins_mouse
#'
#' @keywords datasets
#'
#' @rdname PKN_proteins
"PKN_proteins_human"
#' @rdname PKN_proteins
"PKN_proteins_human_atlas"
#' @rdname PKN_proteins
"PKN_proteins_mouse"

# ---------------------------------------------------------------------- #
#                    Prior Knowledge Networks as igraph                  #
# ---------------------------------------------------------------------- #

#' Prior Knowledge Networks (PKN) as igraph Objects
#'
#' Directed interactomes represented as `igraph` objects, reporting processed interactions
#' from SIGNOR and PhosphoSitePlus. Some versions include *experimentally-derived*
#' kinase-substrate relationships from Serine, Threonine, and Tyrosine Atlases (PMIDs: 36631611, 38720073).
#'
#' @format Each `db_*` object is an `igraph` object with:
#' - **Vertex attributes**:
#'   \describe{
#'     \item{name}{Character, Biological Entity Name, corresponds to `ENTITY` column of tabular form}
#'     \item{ID}{Character, Biological Entity Identifier}
#'     \item{TYPE}{Character, Biological Entity Molecular Type}
#'     \item{DATABASE}{Character, Database Source for Biological Entity Identifier}
#'   }
#' - **Edge attributes**:
#'   \describe{
#'     \item{INTERACTION}{Character, Sign of the Interaction, 1 Activating, -1 Inhibiting}
#'     \item{IDA}{Character, Source Biological Entity Identifier}
#'     \item{IDB}{Character, Target Biological Entity Identifier}
#'     \item{MECHANISM}{Character, Molecular Mechanism of the Interaction}
#'     \item{RESIDUE}{Character, if the MECHANISM is a phosphorylation, the modified residue in the interaction}
#'     \item{SEQUENCE}{Character, if the MECHANISM is a phosphorylation, a 15mer centered on the modified residue}
#'     \item{TYPEA}{Character, Source Biological Entity Molecular Type}
#'     \item{DATABASEA}{Character, Source Biological Entity Database}
#'     \item{TYPEB}{Character, Target Biological Entity Molecular Type}
#'     \item{DATABASEB}{Character, Target Biological Entity Database}
#'     \item{PHOSPHO_KEY_GN_SEQ}{Character composed of ENTITYB and SEQUENCE column}
#'     \item{SOURCE}{Character, Source Database for the interaction; either SIGNOR or PhosphoSitePlus or both}
#'     \item{DIRECT}{Logical, whether the interaction is direct}
#'   }
#'
#' @details Available datasets:
#' - `db_human_dir`: Human Prior Knowledge Network with *direct* interactions.
#' - `db_human_ind`: Human Prior Knowledge Network with *direct* and *indirect* interactions.
#' - `db_human_atlas_dir`: Human PKN with *direct* interactions and kinase-substrate relationships.
#' - `db_human_atlas_ind`: Human PKN with *direct* and *indirect* interactions and kinase-substrate relationships.
#' - `db_mouse_dir`: Mouse Prior Knowledge Network with *direct* interactions.
#' - `db_mouse_ind`: Mouse Prior Knowledge Network with *direct* and *indirect* interactions.
#'
#' @source
#' - <https://signor.uniroma2.it>
#' - <https://www.phosphosite.org>
#' - <https://pubmed.ncbi.nlm.nih.gov/36631611/> (for kinase-substrate interactions)
#' - <https://pubmed.ncbi.nlm.nih.gov/38720073/> (for kinase-substrate interactions)
#'
#' @name db_PKN
#' @aliases db_human_dir db_human_ind db_human_atlas_dir db_human_atlas_ind db_mouse_dir db_mouse_ind
#'
#' @keywords datasets

#' @rdname db_PKN
"db_human_dir"

#' @rdname db_PKN
"db_human_ind"

#' @rdname db_PKN
"db_human_atlas_dir"

#' @rdname db_PKN
"db_human_atlas_ind"

#' @rdname db_PKN
"db_mouse_dir"

#' @rdname db_PKN
"db_mouse_ind"

# ---------------------------------------------------------------------- #
#               TFEA/KSEA Regulons from multiple resources               #
# ---------------------------------------------------------------------- #

#' Transcription Factor Enrichment Analysis (TFEA) and Kinase-Substrate Enrichment Analysis (KSEA) Databases
#'
#' Processed databases for transcription factor enrichment analysis (TFEA) and kinase-substrate enrichment analysis (KSEA)
#' containing curated regulatory interactions between transcription factors (TFs) or kinases/phosphatases and their targets.
#'
#' TFEA datasets optionally include confidence scores and are derived from SIGNOR, Dorothea, and CoLLECTRI.
#' KSEA datasets integrate interactions from SIGNOR, OmniPath, PhosphoSitePlus, and kinase-substrate atlases.
#'
#' @format
#' - **TFEA datasets** (`tfea_db_human`, `tfea_db_human_collectri`, `tfea_db_mouse`, `tfea_db_mouse_collectri`):
#'   Data frames with the following columns:
#'   \describe{
#'     \item{tf}{Character, Transcription Factor}
#'     \item{target}{Character, Target Protein Regulated by the TF}
#'     \item{mor}{Character, Mode of Regulation (Activation/Inhibition)}
#'     \item{confidence}{Character, Confidence Score for Interaction (Present only in `tfea_db_*` datasets)}
#'   }
#'
#' - **KSEA datasets** (`ksea_db_human`, `ksea_db_human_atlas`, `ksea_db_mouse`):
#'   Data frames with the following columns:
#'   \describe{
#'     \item{tf}{Character, Kinase or Phosphatase}
#'     \item{target}{Character, Modified Target Protein}
#'     \item{mor}{Character, Mode of Regulation (Activation/Inhibition); for atlas, it is a continous value representing the probability of phosphorylation of the target}
#'   }
#'
#' @details
#' - `tfea_db_human`: A table of human transcription factor-target interactions from SIGNOR and Dorothea (confidence A).
#' - `tfea_db_human_collectri`: A curated table of human transcription factor-target interactions from SIGNOR and Collectri.
#' - `tfea_db_mouse`: A curated table of mouse transcription factor-target interactions from SIGNOR and Dorothea (confidence A)..
#' - `tfea_db_mouse_collectri`:  A curated table of mouse transcription factor-target interactions from SIGNOR and Collectri.
#' - `ksea_db_human`: A table of human kinase/phosphatase-substrate interactions from SIGNOR, Omnipath and PhosphoSitePlus.
#' - `ksea_db_human_atlas`: An extended human kinase-substrate dataset with additional interactions from Serine/Threonine and Tyrosine Kinome Atlases (PMIDs: 36631611, 38720073).
#' - `ksea_db_mouse`: A table of mouse kinase/phosphatase-substrate interactions from SIGNOR, Omnipath and PhosphoSitePlus.
#'
#' @source
#' - <https://signor.uniroma2.it>
#' - <https://www.phosphosite.org>
#' - <https://collectri.bioinf.uni-leipzig.de> (for CoLLECTRI-based datasets)
#' - <https://pubmed.ncbi.nlm.nih.gov/36631611/>
#' - <https://pubmed.ncbi.nlm.nih.gov/38720073/>
#' - <https://omnipathdb.org/>
#'
#' @name tfea_ksea_datasets
#' @aliases tfea_db_human tfea_db_human_collectri tfea_db_mouse tfea_db_mouse_collectri
#'          ksea_db_human ksea_db_human_atlas ksea_db_mouse
#' @keywords datasets
NULL

#' @rdname tfea_ksea_datasets
"tfea_db_human"

#' @rdname tfea_ksea_datasets
"tfea_db_human_collectri"

#' @rdname tfea_ksea_datasets
"tfea_db_mouse"

#' @rdname tfea_ksea_datasets
"tfea_db_mouse_collectri"

#' @rdname tfea_ksea_datasets
"ksea_db_human"

#' @rdname tfea_ksea_datasets
"ksea_db_mouse"


# ---------------------------------------------------------------------- #
#           PhosphoScore regulatory sites information                    #
# ---------------------------------------------------------------------- #

#' PhosphoScore - Regulatory Phosphosites Datasets
#'
#' Curated regulatory phosphosites for human and mouse proteins from SIGNOR and PhosphoSitePlus.
#' Each dataset reports whether the phosphosite has an activating or inhibitory effect on the
#' activity or abundance of the substrate protein. Versions are provided with or without
#' amino acid position data.
#'
#' @format A data frame with gene symbol, site ID, regulatory effect (activation or inhibition),
#' and optionally amino acid position.
#' Data frames with the following columns:
#'   \describe{
#'     \item{PHOSPHO_KEY_GN_SEQ}{Character, identifier of the phosphosite; `Gene-Sequence` or `Gene-Residue`}
#'     \item{UNIPROT}{Character, Target Protein Regulated by the TF}
#'     \item{ACTIVATION}{Character, Mode of Regulation (Activation/Inhibition)}
#'     \item{confidence}{Character, Confidence Score for Interaction (Present only in `tfea_db_*` datasets)}
#'   }
#'
#' @details These objects support the PhosphoScore pipeline for functional phosphoproteomics,
#' allowing integration of site-level regulatory knowledge. They differ by organism (human/mouse),
#' regulatory effect filtering (activity vs abundance and activity regulators), and presence of amino acid position.
#' @source
#' - <https://signor.uniroma2.it>
#' - <https://www.phosphosite.org>
#' @seealso \code{\link{get_phosphoscore_info}}
#'
#' @name phosphoscore_datasets
#' @aliases good_phos_df_human_act good_phos_df_human_act_aapos
#'          good_phos_df_human_all good_phos_df_human_all_aapos
#'          good_phos_df_mouse_act good_phos_df_mouse_act_aapos
#'          good_phos_df_mouse_all good_phos_df_mouse_all_aapos
#' @keywords datasets
NULL

# Individual object documentation
#' @rdname phosphoscore_datasets
"good_phos_df_human_act"

#' @rdname phosphoscore_datasets
"good_phos_df_human_act_aapos"

#' @rdname phosphoscore_datasets
"good_phos_df_human_all"

#' @rdname phosphoscore_datasets
"good_phos_df_human_all_aapos"

#' @rdname phosphoscore_datasets
"good_phos_df_mouse_act"

#' @rdname phosphoscore_datasets
"good_phos_df_mouse_act_aapos"

#' @rdname phosphoscore_datasets
"good_phos_df_mouse_all"

#' @rdname phosphoscore_datasets
"good_phos_df_mouse_all_aapos"

# ---------------------------------------------------------------------- #
#              GO Molecular Function Tables Annotations                  #
# ---------------------------------------------------------------------- #

#' Gene Ontology Molecular Function Annotations (Mouse and Human)
#'
#' Curated datasets linking gene symbols to their associated molecular functions (MF)
#' and UniProt identifiers, for both mouse and human. These annotations are
#' typically used to classify proteins into categories such as, trasncription
#' factors ("tf"), kinases ("kin"), phosphatases ("phos"), other proteins
#' ("other"), and phenotype ("pheno"). These annotations are tailored on
#' signaling key players.
#'
#' @format
#' - **gomf_human**: A data frame with:
#'   \describe{
#'     \item{gene_name}{Gene symbol (character).}
#'     \item{mf}{Molecular function category (character).}
#'     \item{UNIPROT}{UniProt protein accession identifier (character).}
#'   }
#'
#' - **gomf_mouse**: Same format as `gomf_mouse`, for human gene symbols.
#'
#' @details
#' - `\code{gomf_mouse}`: Mouse gene annotations with MF classes.
#' - `\code{gomf_human}`: Human gene annotations with MF classes.
#'
#' @source
#' Gene annotations were compiled from proteins GO Molecular Function annotations reported in UniProt database.
#'
#' @name gomf_datasets
#' @aliases gomf_mouse gomf_human
#' @keywords datasets
NULL

#' @rdname gomf_datasets
"gomf_mouse"

#' @rdname gomf_datasets
"gomf_human"

# ---------------------------------------------------------------------- #
#              Phenoscore Datasets Annotations                           #
# ---------------------------------------------------------------------- #

#' SIGNOR proteins as background for PhenoScore Randomization
#'
#' Vector of Gene Symbols from SIGNOR database, generated with [sp_parsing]
#' function; from this dataset a random set of proteins is sampled in
#' PhenoScore algorithm analysis.
"background_phenoscore"

#' Summary of Phenotypes in ProxPath Results
#'
#' A summary table listing all phenotypes (`EndPathways`) found in the ProxPath output,
#' along with the number of distinct signaling paths leading to each phenotype.
#'
#' This object is derived from the `phenoscore_distances_table`. A data frame hosted on PerfettoLab server
#' annotating proximity scores computed between SIGNOR signaling proteins and phenotypic endpoints.
#' This table is the output of the ProxPath algorithm (38102483), which retrieves all directed paths from each SIGNOR protein
#' (QueryNode) to known phenotypes (EndNode) in the SIGNOR database. These paths are then scored based on edge weights
#' and annotated with their *final regulatory effect*.
#' The dataset is generated by running a Python script (`proxpath.py`) available upon request.
#' The script uses the SIGNOR network to extract relevant paths and exports the result to a minimized output file
#' which is then loaded into R and stored in the package.
#' by counting how many
#' times each phenotype appears as the target node of a signaling path. It can be used
#' to prioritize or filter phenotypes based on their connectivity within the SIGNOR
#' network.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{phenotype}{Character. Name of the phenotype or cellular process (originally from the EndPathways column).}
#'   \item{n_paths}{Integer. Number of signaling paths from SIGNOR proteins to the phenotype.}
#' }
#'
#' @details
#' This summary allows users to explore the coverage of phenotype annotations in the dataset,
#' and can guide selection of phenotypes for enrichment or proximity scoring analyses.
#'
#' @seealso \code{\link{phenoscore_distances_table}}
#'
#' @examples
#' data(proxpath_phenotypes)
#' head(proxpath_phenotypes)
"proxpath_phenotypes"

# ---------------------------------------------------------------------- #
#              Toy datasets to run SignalingProfiler                     #
# ---------------------------------------------------------------------- #

#' Transcriptomics toy dataset
#'
#' A toy dataset representing transcriptomic differential expression results
#' from MCF7 cell lines upon metformin treatment (metformin vs control) (PMID: 27135362).
#' This dataset is structured to be compatible with SignalingProfiler.
#'
#' @format A data frame with one row per gene, and the following columns:
#' \describe{
#'   \item{gene_name}{Character. Gene symbol used as key in SignalingProfiler.}
#'   \item{difference}{Numeric. Fold-change in gene expression between the two conditions.}
#'   \item{logpval}{Numeric. -log10(p-value) of the statistical test assessing differential expression.}
#'   \item{significant}{Character. A '+' if the gene is significantly modulated, 'NA' otherwise.}
#' }
#'
#' @examples
#' data('tr_toy_df')
"tr_toy_df"

#' Proteomics toy dataset
#'
#' A toy dataset representing proteomic differential expression results
#' from MCF7 cell lines upon metformin treatment (metformin vs control) (PMID: 27135362).
#' This dataset is structured to be compatible with SignalingProfiler.
#'
#' @format A data frame with one row per gene, and the following columns:
#' \describe{
#'   \item{gene_name}{Character. Gene symbol used as key in SignalingProfiler.}
#'   \item{UNIPROT}{Character. UniProt protein accession.}
#'   \item{difference}{Numeric. Fold-change in gene expression between the two conditions.}
#'   \item{logpval}{Numeric. -log10(p-value) of the statistical test assessing differential expression.}
#'   \item{significant}{Character. A '+' if the gene is significantly modulated, 'NA' otherwise.}
#' }
#'
#' @examples
#' data('prot_toy_df')
"prot_toy_df"

#' Phosphoproteomics toy dataset
#'
#' A toy dataset representing phosphoproteomics differential expression results
#' from MCF7 cell lines upon metformin treatment (metformin vs control) (PMID: 27135362).
#' This dataset is structured to be compatible with SignalingProfiler.
#'
#' @format A data frame with one row per gene, and the following columns:
#' \describe{
#'  \item{UNIPROT}{Character. UniProt protein accession number.}
#'  \item{gene_name}{Character. Gene symbol corresponding to the protein.}
#'  \item{aminoacid}{Character. Single-letter code of the phosphorylated residue (e.g., 'S', 'T', 'Y').}
#'  \item{position}{Integer. Position of the phosphorylated amino acid within the protein sequence.}
#'  \item{sequence_window}{Character. 15-mer amino acid sequence centered on the modified site.}
#'  \item{difference}{Numeric. Log fold-change in phosphorylation level between two conditions.}
#'  \item{logpval}{Numeric. -Log10(p-value) of the significance test for the phosphosite change.}
#'  \item{significant}{Character. A '+' if the gene is significantly modulated, 'NA' otherwise.}
#' }
#'
#' @examples
#' data('phospho_toy_df')
"phospho_toy_df"

# ---------------------------------------------------------------------- #
#                    Toy outputs of SignalingProfiler                    #
# ---------------------------------------------------------------------- #

#' Toy Output of Signaling Proteins Activity Inference (Step 1)
#'
#' A toy example of the result of protein activity prediction with \code{SignalingProfiler}.
#' This dataset contains 92 modulated proteins upon metformin treatment
#' inferred using `run_footprint_based_analysis` and `compute_phosphoscore`
#' from transcriptomics, proteomics, and phosphoproteomics datasets
#' of MCF7 cell lines upon metformin treatment (PMID:  27135362).
#'
#' @format A data frame of 92 proteins and 5 columns:
#' \describe{
#'   \item{`gene_name`}{The Gene Symbol of the inferred protein.}
#'   \item{`final_score`}{The combined activity score from both methods.}
#'   \item{`method`}{Indicates whether the score is derived from `"VIPER"`, `"PhosphoScore"`, or `"VIPER+PhosphoScore"`.}
#'   \item{`mf`}{Molecular function (`"tf"` for transcription factors, `"kin"` or `"phos"` for kinases and phosphatases).}
#'   \item{`UNIPROT`}{The UniProt ID of the inferred protein.}
#' }
#' @examples
#' head(toy_prot_activity_df)
"toy_prot_activity_df"

#' Toy Output of Signaling Network construction (Step 2)
#'
#' A toy example of the signaling rewiring network generated as part of the phenotype scoring step
#' in the \code{SignalingProfiler} pipeline. This object contains a signaling network
#' composed of molecular entities (tf, kin, phos, etc) connected through causal interactions (signed and directed).
#'
#' The toy network was computed using \code{two_layer_naive_network()}, \code{run_carnival_and_create_graph()} and
#' \code{expand_and_map_edges()} functions on a multi-omic dataset
#' from MCF7 breast cancer cells treated with metformin.
#'  It includes the graph structure, node metadata, and edge annotations using phosphoproteomics data.
#'
#' @format A named list with 3 components:
#' \describe{
#'   \item{igraph_network}{An \code{igraph} object representing the signaling network (99 nodes, 219 edges).}
#'
#'   \item{nodes_df}{A data frame with metadata about each node in the network. Each row corresponds to a protein, complex, or fusion protein. Columns include:
#'     \describe{
#'       \item{gene_name}{Gene symbol or entity name (e.g., protein, complex, fusion protein).}
#'       \item{carnival_activity}{Numeric score indicating activity state from CARNIVAL or related tools; typically +100 (activated) or -100 (inhibited).}
#'       \item{UNIPROT}{UniProt accession for the entity; complexes, fusion or phenotypes have SIGNOR IDs (e.g., "SIGNOR-C150").}
#'       \item{final_score}{Quantitative score from integrated omic layers (e.g., VIPER, user input).}
#'       \item{method}{Origin of the activity inference: "VIPER", "CARNIVAL", "user", etc.}
#'       \item{discordant}{Logical value indicating whether predictions across sources are inconsistent.}
#'       \item{mf}{Molecular function class: e.g., "kin" (kinase), "tf" (transcription factor), "rec" (receptor), "complex", etc.}
#'     }
#'   }
#'
#'   \item{edges_df}{A data frame describing causal relationships between nodes. Each row corresponds to a signed and directed interaction. Columns include:
#'     \describe{
#'       \item{source}{Source node of the interaction.}
#'       \item{target}{Target node of the interaction.}
#'       \item{sign}{Character or numeric sign of the effect: "1" for activation, "-1" for inhibition.}
#'       \item{carnival_weight}{Weight of the edge as inferred by CARNIVAL, reflecting its relevance in logic modeling.}
#'       \item{direct}{Logical indicating whether the interaction is direct (curated from sources like SIGNOR).}
#'       \item{mechanism}{Interaction type, e.g., "phosphorylation", "relocalization"; multiple mechanisms are semicolon-separated.}
#'       \item{aminoacid}{Modified residue (e.g., "S473") when applicable.}
#'       \item{is_quantified}{Logical or character indicating whether the interaction is backed by omics quantification.}
#'       \item{is_significant}{Logical or character indicating whether the quantified evidence is statistically significant.}
#'       \item{FC}{Fold change of the associated signal (e.g., phosphosite level); empty if not applicable.}
#'     }
#'   }
#' }
#'
#' @details
#' This network enables users to visualize and analyze how perturbations—such as drug treatments—affect
#' signaling flow withing the cells modelled through causal signaling paths. To functionally
#' interpret this proteic network use it as input for #' \code{phenoscore_computation()}.
#'
#' @seealso \code{\link{run_carnival_and_create_graph}}, \code{\link{two_layer_naive_network}}
#'
#' @examples
#' data(toy_opt_network)
#' head(toy_opt_network$nodes_df)
#' head(toy_opt_network$edges_df)
"toy_opt_network"

#' Toy Output of PhenoScore analysis (Step 3)
#'
#' A toy example of the result of phenotypes activity prediction with \code{SignalingProfiler}.
#' This object was obtained by applying `phenoscore_computation` function
#' on a multi-omic dataset derived from  MCF7 breast cancer cells treated with metformin.
#' The output includes graphical and tabular summaries
#' of predicted phenotype modulation based on signaling network analysis, and
#' (optionally) a network as *igraph object*.
#'
#' @format A named list of 4 elements:
#' \describe{
#'   \item{barplot}{A \code{ggplot2} object representing metformin-induced phenotypes modulation as a barplot.}
#'   \item{table_regulators}{A data frame of 19 rows (10 phenotypes) with their predicted proteic activators and/or inhibitors in the network. Each phenotype can occur twice, one reporting activators and one regulators}
#'   \item{table_phenotypes}{A data frame of 10 rows (10 phenotypes) with their activity modulation `phenoscore` upon metformin treatment.}
#'   \item{sp_object_phenotypes}{A list reporting the metformin-induced signaling rewiring network (109 nodes and 309 edges), as *igraph_network*, *nodes_df*, and *edges_df*.}
#' }
#'
#' @details
#' This object illustrates the structure and content of the typical output generated by the \code{phenoscore_computation()} function in the SignalingProfiler workflow. It enables downstream visualization, interpretation, and pathway analysis of phenotype modulation in response to a perturbation—in this case, metformin treatment in MCF7 cells.
#'
#' @seealso \code{\link{phenoscore_computation}}, \code{\link{toy_sp_output}}
#'
#' @examples
#' data(toy_phenoscore_output)
#' names(toy_phenoscore_output)
"toy_phenoscore_output"

#' Toy Signaling Network Output from PhenoScore Analysis
#'
#' A toy example of the signaling rewiring network generated as part of the phenotype scoring step
#' in the \code{SignalingProfiler} pipeline. This object corresponds to the \code{sp_object_phenotypes}
#' element of \code{toy_phenoscore_output}, and represents the signaling network composed of
#' molecular entities (tf, kin, phos, etc) to phenotypic endpoints (cellular functional traits).
#'
#' The toy network was computed using the \code{phenoscore_computation()} function on a multi-omic dataset
#' from MCF7 breast cancer cells treated with metformin. It includes the graph structure, node metadata,
#' and edge annotations.
#'
#' @format A named list with 3 components:
#' \describe{
#'   \item{igraph_network}{An \code{igraph} object representing the signaling rewiring network (109 nodes, 309 edges).}
#'
#'   \item{nodes_df}{A data frame with metadata about each node in the network. Each row corresponds to a protein, complex, or phenotype. Columns include:
#'     \describe{
#'       \item{gene_name}{Gene symbol or entity name (e.g., protein, complex, or phenotype).}
#'       \item{carnival_activity}{Numeric score indicating activity state from CARNIVAL or related tools; typically +100 (activated) or -100 (inhibited).}
#'       \item{UNIPROT}{UniProt accession for the entity; complexes, fusion or phenotypes have SIGNOR IDs (e.g., "SIGNOR-C150").}
#'       \item{final_score}{Quantitative score from integrated omic layers (e.g., VIPER, user input).}
#'       \item{method}{Origin of the activity inference: "VIPER", "CARNIVAL", "user", etc.}
#'       \item{discordant}{Logical value indicating whether predictions across sources are inconsistent.}
#'       \item{mf}{Molecular function class: e.g., "kin" (kinase), "tf" (transcription factor), "rec" (receptor), "complex", etc.}
#'     }
#'   }
#'
#'   \item{edges_df}{A data frame describing causal relationships between nodes. Each row corresponds to a directed interaction. Columns include:
#'     \describe{
#'       \item{source}{Source node of the interaction.}
#'       \item{target}{Target node of the interaction.}
#'       \item{sign}{Character or numeric sign of the effect: "1" for activation, "-1" for inhibition.}
#'       \item{carnival_weight}{Weight of the edge as inferred by CARNIVAL, reflecting its relevance in logic modeling.}
#'       \item{direct}{Logical indicating whether the interaction is direct (curated from sources like SIGNOR).}
#'       \item{mechanism}{Interaction type, e.g., "phosphorylation", "relocalization"; multiple mechanisms are semicolon-separated.}
#'       \item{aminoacid}{Modified residue (e.g., "S473") when applicable.}
#'       \item{is_quantified}{Logical or character indicating whether the interaction is backed by omics quantification.}
#'       \item{is_significant}{Logical or character indicating whether the quantified evidence is statistically significant.}
#'       \item{FC}{Fold change of the associated signal (e.g., phosphosite level); empty if not applicable.}
#'     }
#'   }
#' }
#'
#' @details
#' This network enables users to visualize and analyze how perturbations—such as drug treatments—affect
#' cellular phenotypes through causal signaling paths. The toy network includes transcription factors, kinases, phosphatases,
#' intermediate nodes, and phenotypic endpoints, and serves as a reference for interpreting output from the
#' \code{phenoscore_computation()}
#'
#' @seealso \code{\link{toy_phenoscore_output}}, \code{\link{phenoscore_computation}}
#'
#' @examples
#' data(toy_sp_output)
#' head(toy_sp_output$nodes_df)
#' head(toy_sp_output$edges_df)
"toy_sp_output"




