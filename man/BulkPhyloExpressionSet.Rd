% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloset_bulk.R
\name{BulkPhyloExpressionSet}
\alias{BulkPhyloExpressionSet}
\title{Bulk PhyloExpressionSet Class}
\usage{
BulkPhyloExpressionSet(
  strata = stop("@strata is required"),
  strata_values = stop("@strata_values is required"),
  gene_ids = stop("@gene_ids is required"),
  name = "Phylo Expression Set",
  species = character(0),
  index_type = "TXI",
  identities_label = "Identities",
  null_conservation_sample_size = 5000L,
  precomputed_null_conservation_txis = NULL,
  .expression = stop("@.expression is required"),
  .groups = stop("@.groups is required"),
  precomputed_bootstrapped_txis = NULL
)
}
\arguments{
\item{strata}{Factor vector of phylostratum assignments for each gene}

\item{strata_values}{Numeric vector of phylostratum values used in TXI calculations}

\item{gene_ids}{Character vector of gene identifiers}

\item{name}{Character string naming the dataset (default: "Phylo Expression Set")}

\item{species}{Character string specifying the species (default: NULL)}

\item{index_type}{Character string specifying the transcriptomic index type (default: "TXI")}

\item{identities_label}{Character string labeling the identities (default: "Identities")}

\item{null_conservation_sample_size}{Numeric value for null conservation sample size (default: 5000)}

\item{precomputed_null_conservation_txis}{Precomputed null conservation TXI values (default: NULL)}

\item{.expression}{Matrix of expression counts with genes as rows and samples as columns}

\item{.groups}{Factor vector indicating which identity each sample belongs to}

\item{precomputed_bootstrapped_txis}{Precomputed bootstrapped TXI values (default: NULL)}
}
\description{
S7 class for bulk RNA-seq phylotranscriptomic expression data.
This class handles expression data with biological replicates.
}
