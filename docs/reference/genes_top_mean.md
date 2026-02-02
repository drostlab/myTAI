# Select Top Mean Expressed Genes

Select genes with the highest mean expression across samples.

## Usage

``` r
genes_top_mean(phyex_set, top_p = 0.99, top_k = NULL)
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

Character vector of gene IDs with mean expression \>= top_p quantile or
top top_k genes

## Details

This function identifies genes with the highest mean expression levels,
which are often the most reliably detected and functionally important.

## Examples

``` r
# Select top 1% most expressed genes by mean
high_expr_genes <- genes_top_mean(example_phyex_set, top_p = 0.99)

# Select top 1000 most expressed genes
top_1000_genes <- genes_top_mean(example_phyex_set, top_k = 1000)
```
