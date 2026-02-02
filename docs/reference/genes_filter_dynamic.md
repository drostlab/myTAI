# Filter Dynamic Expression Genes

Filter genes based on expression variance to select the most dynamically
expressed genes.

## Usage

``` r
genes_filter_dynamic(e, thr = 0.9)
```

## Arguments

- e:

  Matrix of expression values with genes as rows and samples as columns

- thr:

  Threshold quantile for variance filtering (default: 0.9)

## Value

Matrix containing only genes above the variance threshold

## Details

This function calculates the variance for each gene across samples and
retains only genes with variance above the specified quantile threshold.
This helps focus analysis on genes that show significant expression
changes. \# Filter top 10 \# filtered_expr \<-
genes_filter_dynamic(expression_matrix, thr = 0.9)
