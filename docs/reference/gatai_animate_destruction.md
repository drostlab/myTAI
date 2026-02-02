# Animate GATAI Destruction Process

Create an animation showing how the transcriptomic signature changes
during the GATAI gene removal process across generations.

## Usage

``` r
gatai_animate_destruction(
  phyex_set,
  save_file = NULL,
  fps = 10,
  width = 1000,
  height = 800,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object containing the original gene expression
  data.

- save_file:

  Optional file path to save the animation (default: NULL, returns
  animation object).

- fps:

  Frames per second for the animation (default: 20).

- width:

  Width of the animation in pixels (default: 1000).

- height:

  Height of the animation in pixels (default: 800).

- ...:

  Additional arguments passed to
  [`gataiR::gatai`](https://rdrr.io/pkg/gataiR/man/gatai.html).

## Value

If save_file is NULL, returns a gganimate animation object. If save_file
is specified, saves the animation and returns the file path invisibly.

## Details

This function runs GATAI for a single run while saving intermediate TAI
values at each generation. It then creates an animated plot showing how
the transcriptomic signature evolves as genes are progressively removed.
The animation shows: - The original signature (generation 0) -
Progressive changes through each generation - Final signature after
convergence

The intermediate file format contains generation numbers in the first
column and TAI values for each developmental stage in subsequent
columns.
