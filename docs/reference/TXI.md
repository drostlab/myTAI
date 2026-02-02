# Calculate Transcriptomic Index (TXI)

Calculate the transcriptomic index values for a PhyloExpressionSet. This
function provides backward compatibility with the old TXI() function.

## Usage

``` r
TXI(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

## Value

Numeric vector of TXI values for each identity

## Examples

``` r
# Calculate TXI values
txi_values <- TXI(example_phyex_set)
```
