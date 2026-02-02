# Calculate Phylostratum Enrichment

Calculate log2(observed/expected) enrichment ratios for phylostrata in a
selected gene set compared to the background distribution.

## Usage

``` r
strata_enrichment(strata, selected_gene_ids)
```

## Arguments

- strata:

  Named factor vector of phylostratum assignments (names are gene IDs)

- selected_gene_ids:

  Character vector of gene IDs to test for enrichment

## Value

A data frame with columns:

- Stratum:

  Phylostratum factor levels

- log_obs_exp:

  Log2 ratio of observed vs expected proportions

## Details

This function calculates enrichment or depletion of phylostrata in a
gene set by comparing the observed proportion of each stratum in the
selected genes to the expected proportion based on the background
distribution in all genes.

Positive values indicate enrichment (more genes than expected), while
negative values indicate depletion (fewer genes than expected).

## Examples

``` r
# Calculate enrichment for a gene set
enrichment <- strata_enrichment(example_phyex_set@strata, example_phyex_set@gene_ids[1:30])
print(enrichment)
#>                 Stratum log_obs_exp
#> 1    Cellular Organisms  0.08837558
#> 2             Eukaryota  0.20650847
#> 3         Viridiplantae  0.87596724
#> 4          Streptophyta  0.87596724
#> 5        Streptophytina        -Inf
#> 6           Embryophyta        -Inf
#> 7           Traceophyta        -Inf
#> 8         Euphyllophyta        -Inf
#> 9         Spermatophyta        -Inf
#> 10        Magnoliopsida        -Inf
#> 11      Mesangiospermae        -Inf
#> 12         Pentapetalae        -Inf
#> 13               Rosids        -Inf
#> 14              Malvids        -Inf
#> 15          Brassicales        -Inf
#> 16         Brassicaceae        -Inf
#> 17           Camelineae        -Inf
#> 18          Arabidopsis  2.19789534
#> 19 Arabidopsis thaliana        -Inf
```
