# Calculate Gene Expression Angles

Calculate developmental stage angles for genes based on their expression
patterns.

## Usage

``` r
get_angles(e)
```

## Arguments

- e:

  Matrix of expression values with genes as rows and stages as columns

## Value

Numeric vector of angles representing each gene's expression pattern

## Details

This function uses PCA to project genes and ideal expression patterns
into 2D space, then calculates the angle of each gene relative to ideal
patterns. The angles represent the peak expression timing during
development.
