# PhyloExpressionSet Base Class

Abstract S7 base class for storing and manipulating phylotranscriptomic
expression data. This class provides the common interface for both bulk
and single-cell phylotranscriptomic data.

## Usage

``` r
PhyloExpressionSetBase(
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
  .null_conservation_txis = NULL
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

  Character string labeling the identities (default: "Identities")

- gene_ids:

  Character vector of gene identifiers (default: character(0),
  auto-generated from expression rownames if not provided)

- null_conservation_sample_size:

  Numeric value for null conservation sample size (default: 5000)

- .null_conservation_txis:

  Precomputed null conservation TXI values (default: NULL)

## Value

A PhyloExpressionSetBase object

## Details

The PhyloExpressionSetBase class serves as the foundation for
phylotranscriptomic analysis, providing shared functionality for both
bulk and single-cell data types.

**Abstract Properties:** Subclasses must implement the
`expression_collapsed` property to define how expression data should be
collapsed across replicates or cells.

**Computed Properties:** Several properties are computed automatically
when accessed:

- `gene_ids` - Character vector of gene identifiers (rownames of
  expression matrix)

- `identities` - Character vector of identity labels (colnames of
  collapsed expression)

- `sample_names` - Character vector of sample names (colnames of
  expression matrix)

- `num_identities` - Integer count of unique identities

- `num_samples` - Integer count of total samples

- `num_genes` - Integer count of genes

- `num_strata` - Integer count of phylostrata

- `index_full_name` - Full name of the transcriptomic index type

- `group_map` - List mapping identity names to sample names

- `TXI` - Numeric vector of TXI values for each identity (computed from
  collapsed expression)

- `TXI_sample` - Numeric vector of TXI values for each sample (computed
  from raw expression)

- `null_conservation_txis` - Matrix of null conservation TXI values for
  statistical testing

**Validation:** The class ensures consistency between expression data,
phylostratum assignments, and groupings. All gene-level vectors must
have matching lengths, and sample groupings must be consistent.
