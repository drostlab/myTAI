# Plot Signature Under Different Transformations

Compare transcriptomic signatures under various data transformations to
assess the robustness of phylotranscriptomic patterns.

## Usage

``` r
plot_signature_transformed(phyex_set, transformations = COUNT_TRANSFORMS, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- transformations:

  Named list of transformation functions (default: COUNT_TRANSFORMS)

- ...:

  Additional arguments passed to plot_signature_multiple

## Value

A ggplot2 object showing signatures under different transformations

## Details

This function applies different transformations to the same dataset and
compares the resulting transcriptomic signatures. This is useful for
assessing whether phylotranscriptomic patterns are robust to different
data processing approaches or are artifacts of specific transformations.

The analysis works with both bulk and single-cell data, helping to
determine whether phylotranscriptomic patterns are consistent across
different normalization and transformation methods.

## Examples

``` r
# Single-cell data with custom transformations

phyex_set <- example_phyex_set

custom_transforms <- list(raw = identity, log = log1p)
p <- plot_signature_transformed(phyex_set, transformations = custom_transforms)
#> 
Computing: [========================================] 100% (done)                         
```
