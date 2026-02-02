# Transform Phylostratum Values

This function performs transformation of phylostratum values.

## Usage

``` r
tf_PS(phyex_set, transform = "qr")
```

## Arguments

- phyex_set:

  a `PhyloExpressionSet` S7 object.

- transform:

  a character vector of any valid function that transforms PS values.
  Possible values can be:

  - `transform` = `"qr"` (or `"quantilerank"`) : quantile rank
    transformation analogous to Julia function `StatsBase.quantilerank`
    using `method = :tied`.

## Value

a `PhyloExpressionSet` object storing transformed Phylostratum levels.

## Details

This function transforms the phylostratum assignment. The return value
of this function is a PhyloExpressionSet object with transformed
phylostratum `tfPhylostratum` as the first column, satisfying
[`PhyloExpressionSetBase`](https://drostlab.github.io/myTAI/reference/PhyloExpressionSetBase.md).
Note that the input `transform` must be an available function, currently
limited to only `"qr"` (or `"quantilerank"`).

## See also

[`tf`](https://drostlab.github.io/myTAI/reference/tf.md)

## Author

Jaruwatana Sodai Lotharukpong and Lukas Maischak

## Examples

``` r
# get the relative expression profiles for each phylostratum
tfPES <- tf_PS(example_phyex_set, transform = "qr")
```
