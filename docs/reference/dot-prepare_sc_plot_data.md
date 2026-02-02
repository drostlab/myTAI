# Prepare Single-Cell Plot Data

Prepare data frame for single-cell plotting with flexible identities

## Usage

``` r
.prepare_sc_plot_data(
  phyex_set,
  primary_identity = NULL,
  secondary_identity = NULL
)
```

## Arguments

- phyex_set:

  A ScPhyloExpressionSet object

- primary_identity:

  Primary identity column name

- secondary_identity:

  Secondary identity column name (optional)

## Value

List with sample data and aggregated data
