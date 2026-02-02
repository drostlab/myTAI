# Perform Permutation Tests Under Different Transformations

*tf_stability* statistically evaluates the stability of
phylotranscriptomics permutation tests (e.g., `stat_flatline_test`,
`stat_reductive_hourglass_test`, etc.) under different data
transformations using a `PhyloExpressionSet`.

## Usage

``` r
tf_stability(
  phyex_set,
  conservation_test = stat_flatline_test,
  transforms = COUNT_TRANSFORMS
)
```

## Arguments

- phyex_set:

  a `PhyloExpressionSet`.

- conservation_test:

  a conservation test function (e.g. `stat_flatline_test`,
  `stat_reductive_hourglass_test`, etc.)

- transforms:

  named list of transformation functions (default: `COUNT_TRANSFORMS`)

## Value

Named numeric vector of p-values for each transformation.

## Details

Assesses the stability of data transforms on the permutation test of
choice. See [`tf`](https://drostlab.github.io/myTAI/reference/tf.md),
[`stat_flatline_test`](https://drostlab.github.io/myTAI/reference/stat_flatline_test.md),
[`stat_reductive_hourglass_test`](https://drostlab.github.io/myTAI/reference/stat_reductive_hourglass_test.md),
etc.

## References

Lotharukpong JS et al. (2023) (unpublished)

## Author

Jaruwatana Sodai Lotharukpong
