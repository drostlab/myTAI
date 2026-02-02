# Single-Cell PhyloExpressionSet Class

S7 class for single-cell phylotranscriptomic expression data. This class
stores expression matrices and metadata, with support for dimensional
reductions and pseudobulking functionality.

## Usage

``` r
ScPhyloExpressionSet(
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
  .pseudobulk_cache = list(),
  .TXI_sample = numeric(0),
  metadata = NULL,
  selected_idents = character(0),
  idents_colours = list(),
  reductions = list()
)
```

## Arguments

- strata:

  Factor vector of phylostratum assignments for each gene

- strata_values:

  Numeric vector of phylostratum values used in TXI calculations

- expression:

  Sparse or dense matrix of expression counts with genes as rows and
  cells as columns

- groups:

  Factor vector indicating which identity each cell belongs to (derived
  from selected_idents column in metadata)

- name:

  Character string naming the dataset (default: "Phylo Expression Set")

- species:

  Character string specifying the species (default: NULL)

- index_type:

  Character string specifying the transcriptomic index type (default:
  "TXI")

- identities_label:

  Character string labeling the identities (default: "Cell Type")

- gene_ids:

  Character vector of gene identifiers (default: character(0),
  auto-generated from expression rownames if not provided)

- null_conservation_sample_size:

  Numeric value for null conservation sample size (default: 5000)

- .null_conservation_txis:

  Precomputed null conservation TXI values (default: NULL)

- .pseudobulk_cache:

  Internal cache for pseudobulked expression matrices by different
  groupings

- .TXI_sample:

  Internal storage for computed TXI values

- metadata:

  Data frame with cell metadata, where rownames correspond to cell IDs
  and columns contain cell attributes

- selected_idents:

  Character string specifying which metadata column is currently used
  for grouping cells

- idents_colours:

  List of named character vectors specifying colors for each identity
  level, organized by metadata column name

- reductions:

  List of dimensional reduction matrices (PCA, UMAP, etc.) with cells as
  rows and dimensions as columns

## Value

A ScPhyloExpressionSet object

## Details

The ScPhyloExpressionSet class provides a comprehensive framework for
single-cell phylotranscriptomic analysis. Key features include:

**Identity Management:** The `selected_idents` property determines which
metadata column is used for grouping cells. When changed, it
automatically updates the `groups` property and invalidates cached
pseudobulk data to ensure consistency.

**Dimensional Reductions:** The `reductions` property stores
pre-computed dimensional reductions (PCA, UMAP, etc.). If not provided
during construction from Seurat objects, basic PCA and UMAP are computed
automatically.

**Color Management:** `idents_colours` allows custom color schemes for
different metadata columns, ensuring consistent visualization across
plots.

**Computed Properties:** Several properties are computed automatically
when accessed:

- `available_idents` - Character vector of factor columns in metadata
  that can be used for grouping (automatically detected from metadata)

- `expression_collapsed` - Matrix of pseudobulked expression data (genes
  x identities), created by summing expression within each identity
  group

- `TXI_sample` - Named numeric vector of TXI (Transcriptomic Age Index)
  values for each cell, computed using efficient C++ implementation

Inherited computed properties from PhyloExpressionSetBase include:

- `gene_ids` - Character vector of gene identifiers

- `identities` - Character vector of identity labels

- `sample_names` - Character vector of sample names (cell IDs)

- `num_identities` - Integer count of unique cell types/identities

- `num_samples` - Integer count of total cells

- `num_genes` - Integer count of genes

- `num_strata` - Integer count of phylostrata

- `index_full_name` - Full name of the transcriptomic index type

- `group_map` - List mapping identity names to cell IDs

- `TXI` - Numeric vector of TXI values for each identity (computed from
  pseudobulked expression)

- `null_conservation_txis` - Matrix of null conservation TXI values for
  statistical testing

These properties use lazy evaluation and caching for optimal
performance.

## Examples

``` r
# Create from Seurat object
data(example_phyex_set_sc)
sc_set <- example_phyex_set_sc

# Switch to different cell grouping
sc_set@selected_idents <- "day"
#> Loading required package: Matrix

# Access pseudobulked data (computed au  tomatically)
pseudobulk <- sc_set@expression_collapsed

# Access TXI values for each cell
txi_values <- sc_set@TXI_sample
```
