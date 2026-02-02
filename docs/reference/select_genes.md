# Select Genes from PhyloExpressionSet

Extract a subset of genes from a PhyloExpressionSet object.

## Usage

``` r
select_genes(phyex_set, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- ...:

  Additional arguments passed to methods (typically includes 'genes'
  parameter)

## Value

A PhyloExpressionSet object containing only the selected genes

## Examples

``` r
# Select specific genes
selected_set <- select_genes(example_phyex_set, example_phyex_set@gene_ids[1:10])
```
