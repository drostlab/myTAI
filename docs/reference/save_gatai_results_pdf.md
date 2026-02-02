# Save GATAI Analysis Results to PDF

Save removed gene IDs and all GATAI analysis plots to a PDF file.

## Usage

``` r
save_gatai_results_pdf(
  phyex_set,
  gatai_result,
  analysis_dir = "gatai_analysis",
  prefix = "report",
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object containing the original gene expression
  data.

- gatai_result:

  Result list from
  [`destroy_pattern()`](https://drostlab.github.io/myTAI/reference/destroy_pattern.md),
  containing GATAI analysis output.

- analysis_dir:

  Directory to save the PDF file.

- prefix:

  Optional prefix for the PDF filename (default: "report").

- ...:

  Additional arguments passed to
  [`plot_gatai_results()`](https://drostlab.github.io/myTAI/reference/plot_gatai_results.md).

## Value

Invisibly returns the path to the saved PDF.
