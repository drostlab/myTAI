# Format P-Value for Scientific Notation

Format p-values in scientific notation for plot annotations.

## Usage

``` r
exp_p(p, sci_thresh = 4)
```

## Arguments

- p:

  Numeric p-value

- sci_thresh:

  Numeric threshold for using scientific notation (number of decimal
  places)

## Value

Expression object for use in plot annotations

## Details

This function formats p-values in scientific notation using the format
"p = a × 10^b" which is suitable for ggplot2 annotations and maintains
proper mathematical formatting.

## Examples

``` r
# Format p-value for plotting
expr <- exp_p(0.001)
```
