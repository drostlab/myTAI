# Plot Sample Space Visualization

Create a dimensional reduction plot to visualize sample relationships in
gene expression space using PCA or UMAP.

## Usage

``` r
plot_sample_space(
  phyex_set,
  method = c("PCA", "UMAP"),
  colour_by = c("identity", "TXI"),
  seed = 42,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- method:

  Character string specifying the dimensionality reduction method: "PCA"
  or "UMAP" (default: "PCA")

- colour_by:

  Character string specifying what to colour by: "identity" (default),
  "TXI"

- seed:

  Integer seed for reproducible UMAP results (default: 42)

- ...:

  Additional arguments passed to specific methods

## Value

A ggplot2 object showing the sample space visualisation

## Details

This function performs log1p transformation on expression data, removes
genes with zero variance, and applies the specified dimensionality
reduction method. Samples are coloured by their group assignments or TAI
values.

## Examples

``` r
# Create PCA plot coloured by identity
pca_plot <- plot_sample_space(example_phyex_set, method = "PCA", colour_by = "identity")

# Create UMAP plot coloured by TXI
if (requireNamespace("uwot", quietly = TRUE)) {
    umap_plot <- plot_sample_space(example_phyex_set, method = "UMAP", colour_by = "TXI")
}
```
