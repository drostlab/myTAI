# Plot Multiple Transcriptomic Signatures

Create a plot comparing multiple transcriptomic signatures on the same
axes, with options for statistical testing and transformations.

## Usage

``` r
plot_signature_multiple(
  phyex_sets,
  legend_title = "Phylo Expression Set",
  show_p_val = TRUE,
  conservation_test = stat_flatline_test,
  transformation = NULL,
  colours = NULL,
  ...
)
```

## Arguments

- phyex_sets:

  A vector of PhyloExpressionSet objects to compare
  (BulkPhyloExpressionSet or ScPhyloExpressionSet)

- legend_title:

  Title for the legend (default: "Phylo Expression Set")

- show_p_val:

  Logical indicating whether to show p-values (default: TRUE)

- conservation_test:

  Function to use for conservation testing (default: stat_flatline_test)

- transformation:

  Optional transformation function to apply to all datasets (default:
  NULL)

- colours:

  Optional vector of colors for each dataset (default: NULL)

- ...:

  Additional arguments passed to plot_signature

## Value

A ggplot2 object showing multiple transcriptomic signatures

## Details

This function allows comparison of multiple transcriptomic signatures by
overlaying them on the same plot. Each signature is colored differently
and can be tested for conservation patterns (bulk data only). If a
transformation is provided, it's applied to all datasets before
plotting.

The function automatically adapts to the data type: - \*\*Bulk data\*\*:
Line plots with optional statistical testing - \*\*Single-cell data\*\*:
Violin plots showing distributions

All datasets must use the same axis labels (developmental stages or cell
types).

## Examples

``` r
# Compare multiple bulk datasets
phyex_set <- example_phyex_set

bulk_list <- c(phyex_set, 
  phyex_set |> remove_genes(phyex_set@gene_ids[1:5]))
p <- plot_signature_multiple(bulk_list, legend_title = "Dataset")


```
