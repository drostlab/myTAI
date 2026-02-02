# Standardise Expression Data

Standardise gene expression data by centring and scaling each gene.

## Usage

``` r
.to_std_expr(e)
```

## Arguments

- e:

  Matrix of expression values with genes as rows and samples as columns

## Value

Matrix with standardised expression values (mean=0, sd=1 for each gene)

## Details

This function standardises each gene's expression profile by subtracting
the mean and dividing by the standard deviation. Genes with zero or
undefined variance are set to zero. This is useful for comparing
expression patterns across genes with different absolute expression
levels. \# Standardise expression data \# std_expr \<-
.to_std_expr(expression_matrix)
