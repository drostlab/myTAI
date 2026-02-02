# Select Top Variable Genes

Select genes with the highest variance across samples.

## Usage

``` r
genes_top_variance(phyex_set, top_p = 0.99, top_k = NULL)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- top_p:

  Quantile threshold for gene selection (default: 0.99). Ignored if
  top_k is specified.

- top_k:

  Absolute number of top genes to select (default: NULL). Takes
  precedence over top_p.

## Value

Character vector of gene IDs with variance \>= top_p quantile or top
top_k genes

## Details

This function identifies genes with the highest variance across samples,
which are often the most informative for downstream analyses.

## Examples

``` r
# Select top 1% most variable genes
high_var_genes <- genes_top_variance(example_phyex_set, top_p = 0.99)

# Select top 500 most variable genes
top_500_var_genes <- genes_top_variance(example_phyex_set, top_k = 500)
```
