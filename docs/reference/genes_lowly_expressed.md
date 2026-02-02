# Select Lowly Expressed Genes

Select genes with mean expression below a specified threshold.

## Usage

``` r
genes_lowly_expressed(phyex_set, threshold = 1)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- threshold:

  Mean expression threshold (default: 1)

## Value

Character vector of gene IDs with mean expression \<= threshold

## Details

This function identifies genes with low mean expression levels, which
might be candidates for filtering or separate analysis.

## Examples

``` r
# Select genes with mean expression <= 1
low_expr_genes <- genes_lowly_expressed(example_phyex_set, threshold = 1)
```
