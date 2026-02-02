# Remove Genes from PhyloExpressionSet

Remove specified genes from a PhyloExpressionSet object.

## Usage

``` r
remove_genes(
  phyex_set,
  genes,
  new_name = paste(phyex_set@name, "perturbed"),
  reuse_null_txis = TRUE
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- genes:

  Character vector of gene IDs to remove

- new_name:

  Character string for the new dataset name (default: auto-generated)

- reuse_null_txis:

  Logical indicating whether to reuse precomputed null conservation TXIs
  (default: TRUE)

## Value

A PhyloExpressionSet object with the specified genes removed

## Examples

``` r
# Remove specific genes
filtered_set <- remove_genes(example_phyex_set, example_phyex_set@gene_ids[1:5], 
                            new_name = "Filtered Dataset")
```
