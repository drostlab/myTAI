# Set Identity Colours for ScPhyloExpressionSet

Set custom colours for one or more identity columns in the
ScPhyloExpressionSet object.

## Usage

``` r
set_identity_colours(phyex_set, identity_name, colours)
```

## Arguments

- phyex_set:

  A ScPhyloExpressionSet object

- identity_name:

  Character, name of the metadata column to set colours for

- colours:

  Named character vector of colours where names correspond to identity
  values

## Value

ScPhyloExpressionSet object with updated identity colours

## Examples

``` r
# Set colours for specific identity
# colours <- c("TypeA" = "red", "TypeB" = "blue", "TypeC" = "green")
# sc_set <- set_identity_colours(sc_set, "groups", colours)
```
