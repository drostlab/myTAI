# Plot Cullen-Frey Diagram for Distribution Assessment

Create a Cullen-Frey diagram to assess which distribution family best
fits the null sample data.

## Usage

``` r
plot_cullen_frey(test_result)
```

## Arguments

- test_result:

  A TestResult object

## Value

A Cullen-Frey plot from the fitdistrplus package

## Details

The Cullen-Frey diagram plots skewness vs. kurtosis to help identify
appropriate distribution families for the null sample data.
