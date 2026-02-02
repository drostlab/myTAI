# Plot Null TXI Sample Distribution

Create a plot showing the null TXI distribution sample compared to the
observed test TXI values across developmental stages.

## Usage

``` r
plot_null_txi_sample(test_result)
```

## Arguments

- test_result:

  A ConservationTestResult object containing null TXI distributions

## Value

A ggplot2 object showing null samples as gray lines and test TXI as
colored line

## Details

This function creates a visualization of the null hypothesis testing by
plotting: - Gray lines representing individual null TXI samples from
permutations - A horizontal line showing the mean of null
distributions - A colored line showing the observed test TXI values

The plot helps visualize how the observed TXI pattern compares to what
would be expected under the null hypothesis of no conservation signal.
