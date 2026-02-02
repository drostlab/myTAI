# Gene Expression Filtering Functions

Collection of functions for filtering genes based on expression patterns
in PhyloExpressionSet objects.

Generic function to select genes with the highest values for a given
expression metric.

## Usage

``` r
genes_top_expr(phyex_set, FUN = rowMeans, top_p = 0.99, top_k = NULL, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- FUN:

  Function to calculate gene-wise expression metric (default: rowMeans)

- top_p:

  Quantile threshold for gene selection (default: 0.99). Ignored if
  top_k is specified.

- top_k:

  Absolute number of top genes to select (default: NULL). Takes
  precedence over top_p.

- ...:

  Additional arguments passed to FUN

## Value

Character vector of gene IDs with metric values \>= top_p quantile or
top top_k genes

## Details

This function applies the specified function to calculate a metric for
each gene across samples, then selects genes above the specified
quantile threshold or the top k genes by absolute count. If both top_p
and top_k are specified, top_k takes precedence.

## Examples

``` r
# Select top 1% most expressed genes by mean
high_expr_genes <- genes_top_expr(example_phyex_set, function(x) apply(x, 1, mean), top_p = 0.99)

# Select top 100 most expressed genes
top_100_genes <- genes_top_expr(example_phyex_set, function(x) apply(x, 1, mean), top_k = 100)
```
