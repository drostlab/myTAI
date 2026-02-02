# Calculate Transcriptomic Age Index (TAI)

Calculate the transcriptomic age index values for a PhyloExpressionSet.
This function provides backward compatibility with the old TAI()
function.

## Usage

``` r
TAI(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

## Value

Numeric vector of TAI values for each identity

## Examples

``` r
# Calculate TAI values
tai_values <- TAI(example_phyex_set)
```
