# Get Expression Matrix from Seurat Object

Extract expression matrix from Seurat object, preserving sparse format
when possible.

## Usage

``` r
.get_expression_from_seurat(seurat, layer)
```

## Arguments

- seurat:

  A Seurat object

- layer:

  Character string specifying which layer to use

## Value

Expression matrix (sparse or dense)
