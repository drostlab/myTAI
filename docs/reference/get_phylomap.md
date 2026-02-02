# Extract Phylomap from PhyloExpressionSet

Extract a phylomap (tibble of strata-to-gene mappings) from a
PhyloExpressionSet object.

## Usage

``` r
get_phylomap(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

## Value

A tibble with two columns: Stratum (gene age value) and GeneID (gene
identifiers)

## Examples

``` r
# Extract phylomap
phylomap <- get_phylomap(example_phyex_set)
```
