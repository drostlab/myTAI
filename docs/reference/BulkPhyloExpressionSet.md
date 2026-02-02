# Bulk PhyloExpressionSet Class

S7 class for bulk RNA-seq phylotranscriptomic expression data. This
class handles expression data with biological replicates and provides
bootstrapping functionality for statistical analysis.

## Usage

``` r
BulkPhyloExpressionSet(
  strata = stop("@strata is required"),
  strata_values = stop("@strata_values is required"),
  expression = stop("@expression is required"),
  groups = stop("@groups is required"),
  name = "Phylo Expression Set",
  species = character(0),
  index_type = "TXI",
  identities_label = "Identities",
  gene_ids = character(0),
  null_conservation_sample_size = 5000L,
  .null_conservation_txis = NULL,
  .bootstrapped_txis = NULL
)
```

## Arguments

- strata:

  Factor vector of phylostratum assignments for each gene

- strata_values:

  Numeric vector of phylostratum values used in TXI calculations

- expression:

  Matrix of expression counts with genes as rows and samples as columns

- groups:

  Factor vector indicating which identity each sample belongs to

- name:

  Character string naming the dataset (default: "Phylo Expression Set")

- species:

  Character string specifying the species (default: NULL)

- index_type:

  Character string specifying the transcriptomic index type (default:
  "TXI")

- identities_label:

  Character string labeling the identities (default: "Stages")

- gene_ids:

  Character vector of gene identifiers (default: character(0),
  auto-generated from expression rownames if not provided)

- null_conservation_sample_size:

  Numeric value for null conservation sample size (default: 5000)

- .null_conservation_txis:

  Precomputed null conservation TXI values (default: NULL)

- .bootstrapped_txis:

  Precomputed bootstrapped TXI values (default: NULL)

## Value

A BulkPhyloExpressionSet object

## Details

The BulkPhyloExpressionSet class is designed for bulk RNA-seq data with
biological replicates. It extends the base PhyloExpressionSetBase class
with bulk-specific functionality.

**Replicate Handling:** Expression data across biological replicates is
collapsed by taking row means within each experimental condition or
developmental stage.

**Computed Properties:** In addition to inherited computed properties
from the base class, this class provides:

- `expression_collapsed` - Matrix of expression data collapsed across
  replicates (genes x identities)

- `bootstrapped_txis` - Matrix of bootstrapped TXI values for
  statistical inference (500 bootstrap samples x identities)

Inherited computed properties from PhyloExpressionSetBase include:

- `gene_ids` - Character vector of gene identifiers

- `identities` - Character vector of identity labels

- `sample_names` - Character vector of sample names

- `num_identities` - Integer count of unique identities

- `num_samples` - Integer count of total samples

- `num_genes` - Integer count of genes

- `num_strata` - Integer count of phylostrata

- `index_full_name` - Full name of the transcriptomic index type

- `group_map` - List mapping identity names to sample names

- `TXI` - Numeric vector of TXI values for each identity

- `TXI_sample` - Numeric vector of TXI values for each sample

- `null_conservation_txis` - Matrix of null conservation TXI values for
  statistical testing

**Statistical Analysis:** The class supports confidence interval
estimation and standard deviation calculation through bootstrapped TXI
values, enabling robust statistical analysis of developmental or
experimental patterns.
