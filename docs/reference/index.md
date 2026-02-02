# Package index

## Visualisation

Visualise and visually inspect evolutionary transcriptomics patterns as
well as summary statistics.

- [`plot_signature()`](https://drostlab.github.io/myTAI/reference/plot_signature.md)
  : Plot Transcriptomic Signature
- [`plot_signature_gene_quantiles()`](https://drostlab.github.io/myTAI/reference/plot_signature_gene_quantiles.md)
  : Plot Signature Across Gene Expression Quantiles
- [`plot_signature_multiple()`](https://drostlab.github.io/myTAI/reference/plot_signature_multiple.md)
  : Plot Multiple Transcriptomic Signatures
- [`plot_signature_transformed()`](https://drostlab.github.io/myTAI/reference/plot_signature_transformed.md)
  : Plot Signature Under Different Transformations
- [`plot_contribution()`](https://drostlab.github.io/myTAI/reference/plot_contribution.md)
  : Plot Phylostratum Contribution to Transcriptomic Index
- [`plot_cullen_frey()`](https://drostlab.github.io/myTAI/reference/plot_cullen_frey.md)
  : Plot Cullen-Frey Diagram for Distribution Assessment
- [`plot_distribution_expression()`](https://drostlab.github.io/myTAI/reference/plot_distribution_expression.md)
  : Comparing expression levels distributions across developmental
  stages
- [`plot_distribution_pTAI()`](https://drostlab.github.io/myTAI/reference/plot_distribution_pTAI.md)
  : Partial TAI Distribution Plotting Functions
- [`plot_distribution_pTAI_qqplot()`](https://drostlab.github.io/myTAI/reference/plot_distribution_pTAI_qqplot.md)
  : QQ plot comparing partial TAI distributions across developmental
  stages against a reference stage
- [`plot_distribution_strata()`](https://drostlab.github.io/myTAI/reference/plot_distribution_strata.md)
  : Plot Distribution of Genes Across Phylostrata
- [`plot_gatai_results()`](https://drostlab.github.io/myTAI/reference/plot_gatai_results.md)
  : Plot Comprehensive GATAI Results
- [`plot_gene_heatmap()`](https://drostlab.github.io/myTAI/reference/plot_gene_heatmap.md)
  : Plot Gene Expression Heatmap
- [`plot_gene_profiles()`](https://drostlab.github.io/myTAI/reference/plot_gene_profiles.md)
  : Plot Individual Gene Expression Profiles
- [`plot_gene_space()`](https://drostlab.github.io/myTAI/reference/plot_gene_space.md)
  : Plot Gene Space Using PCA
- [`plot_mean_var()`](https://drostlab.github.io/myTAI/reference/plot_mean_var.md)
  : Plot Mean-Variance Relationship
- [`plot_null_txi_sample()`](https://drostlab.github.io/myTAI/reference/plot_null_txi_sample.md)
  : Plot Null TXI Sample Distribution
- [`plot_relative_expression_bar()`](https://drostlab.github.io/myTAI/reference/plot_relative_expression_bar.md)
  : Plot Mean Relative Expression Levels as Barplot
- [`plot_relative_expression_line()`](https://drostlab.github.io/myTAI/reference/plot_relative_expression_line.md)
  : Plot Relative Expression Profiles (Line Plot)
- [`plot_sample_space()`](https://drostlab.github.io/myTAI/reference/plot_sample_space.md)
  : Plot Sample Space Visualization
- [`plot_strata_expression()`](https://drostlab.github.io/myTAI/reference/plot_strata_expression.md)
  : Plot Expression Levels by Phylostratum

## Statistical tests

Functions for statistical testing of evolutionary transcriptomics
patterns (with visualisation).

- [`stat_early_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_early_conservation_test.md)
  : Early Conservation Test
- [`stat_flatline_test()`](https://drostlab.github.io/myTAI/reference/stat_flatline_test.md)
  : Flat Line Test for Conservation Pattern
- [`stat_generic_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md)
  : Generic Conservation Test Framework
- [`stat_late_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_late_conservation_test.md)
  : Late Conservation Test
- [`stat_pairwise_test()`](https://drostlab.github.io/myTAI/reference/stat_pairwise_test.md)
  : Pairwise Conservation Test
- [`stat_reductive_hourglass_test()`](https://drostlab.github.io/myTAI/reference/stat_reductive_hourglass_test.md)
  : Reductive Hourglass Test
- [`stat_reverse_hourglass_test()`](https://drostlab.github.io/myTAI/reference/stat_reverse_hourglass_test.md)
  : Reverse Hourglass Test

## Data structures

Core S7 data classes used for evolutionary transcriptomics analysis.

- [`BulkPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/BulkPhyloExpressionSet.md)
  : Bulk PhyloExpressionSet Class
- [`ConservationTestResult()`](https://drostlab.github.io/myTAI/reference/ConservationTestResult.md)
  : Conservation Test Result S7 Class
- [`PhyloExpressionSetBase()`](https://drostlab.github.io/myTAI/reference/PhyloExpressionSetBase.md)
  : PhyloExpressionSet Base Class
- [`ScPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/ScPhyloExpressionSet.md)
  : Single-Cell PhyloExpressionSet Class
- [`TestResult()`](https://drostlab.github.io/myTAI/reference/TestResult.md)
  : Test Result S7 Class
- [`Distribution()`](https://drostlab.github.io/myTAI/reference/Distribution.md)
  : Distribution S7 Class

## Data structures-related

Functions to convert between dataframe and `BulkPhyloExpressionSet` or
`ScPhyloExpressionSet`.

- [`as_BulkPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/as_BulkPhyloExpressionSet.md)
  : as_BulkPhyloExpressionSet
- [`BulkPhyloExpressionSet_from_df()`](https://drostlab.github.io/myTAI/reference/BulkPhyloExpressionSet_from_df.md)
  : Convert Data to BulkPhyloExpressionSet
- [`ScPhyloExpressionSet_from_matrix()`](https://drostlab.github.io/myTAI/reference/ScPhyloExpressionSet_from_matrix.md)
  : Create Single-Cell PhyloExpressionSet from Expression Matrix
- [`ScPhyloExpressionSet_from_seurat()`](https://drostlab.github.io/myTAI/reference/ScPhyloExpressionSet_from_seurat.md)
  : Convert Seurat Object to Single-Cell PhyloExpressionSet
- [`as_data_frame()`](https://drostlab.github.io/myTAI/reference/as_data_frame.md)
  : Convert BulkPhyloExpressionSet to Data Frame
- [`check_PhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/check_PhyloExpressionSet.md)
  : Check if object is a PhyloExpressionSet
- [`check_BulkPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/check_BulkPhyloExpressionSet.md)
  : Check if object is a BulkPhyloExpressionSet
- [`check_ScPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/check_ScPhyloExpressionSet.md)
  : Check if object is a ScPhyloExpressionSet

## Single-cell specific functions

Functions specifically designed for single-cell RNA-seq data analysis.

- [`downsample_expression()`](https://drostlab.github.io/myTAI/reference/downsample_expression.md)
  : Downsample Expression Matrix by Groups
- [`downsample()`](https://drostlab.github.io/myTAI/reference/downsample.md)
  : Downsample ScPhyloExpressionSet
- [`match_map_sc_matrix()`](https://drostlab.github.io/myTAI/reference/match_map_sc_matrix.md)
  : Match Expression Matrix with Phylostratum Map
- [`match_map_sc_seurat()`](https://drostlab.github.io/myTAI/reference/match_map_sc_seurat.md)
  : Match Single-Cell Expression Data with Phylostratum Map (Seurat)

## Taxonomic information retrieval

Functions to compute transcriptome indices.

- [`taxid()`](https://drostlab.github.io/myTAI/reference/taxid.md) :
  Retrieve taxonomy categories from NCBI Taxonomy

## Transcriptomic index calculation

Compute full and partial evolutionary transcriptomics measures.

- [`TAI()`](https://drostlab.github.io/myTAI/reference/TAI.md) :
  Calculate Transcriptomic Age Index (TAI)
- [`TDI()`](https://drostlab.github.io/myTAI/reference/TDI.md) :
  Calculate Transcriptomic Divergence Index (TDI)
- [`TEI()`](https://drostlab.github.io/myTAI/reference/TEI.md) :
  Calculate Transcriptomic Evolutionary Index (TEI)
- [`TPI()`](https://drostlab.github.io/myTAI/reference/TPI.md) :
  Calculate Transcriptomic Polymorphism Index (TPI)
- [`TXI()`](https://drostlab.github.io/myTAI/reference/TXI.md) :
  Calculate Transcriptomic Index (TXI)
- [`sTXI()`](https://drostlab.github.io/myTAI/reference/sTXI.md) :
  Calculate Stratum-Specific Transcriptomic Index
- [`pTXI()`](https://drostlab.github.io/myTAI/reference/pTXI.md) :
  Calculate Phylostratum-Specific Transcriptomic Index

## Expression or index data-related

Functions to transform, normalise, and manipulate expression datasets.

- [`age.apply()`](https://drostlab.github.io/myTAI/reference/age.apply.md)
  : Age Category Specific apply Function
- [`collapse()`](https://drostlab.github.io/myTAI/reference/collapse.md)
  : Collapse PhyloExpressionSet Replicates
- [`destroy_pattern()`](https://drostlab.github.io/myTAI/reference/destroy_pattern.md)
  : Destroy Phylotranscriptomic Pattern Using GATAI
- [`diagnose_test_robustness()`](https://drostlab.github.io/myTAI/reference/diagnose_test_robustness.md)
  : Diagnose Test Robustness
- [`distributions`](https://drostlab.github.io/myTAI/reference/distributions.md)
  : Predefined Distribution Objects
- [`exp_p()`](https://drostlab.github.io/myTAI/reference/exp_p.md) :
  Format P-Value for Scientific Notation
- [`genes_lowly_expressed()`](https://drostlab.github.io/myTAI/reference/genes_lowly_expressed.md)
  : Select Lowly Expressed Genes
- [`genes_top_expr()`](https://drostlab.github.io/myTAI/reference/genes_top_expr.md)
  : Gene Expression Filtering Functions
- [`genes_top_mean()`](https://drostlab.github.io/myTAI/reference/genes_top_mean.md)
  : Select Top Mean Expressed Genes
- [`genes_top_variance()`](https://drostlab.github.io/myTAI/reference/genes_top_variance.md)
  : Select Top Variable Genes
- [`get_phylomap()`](https://drostlab.github.io/myTAI/reference/get_phylomap.md)
  : Extract Phylomap from PhyloExpressionSet
- [`get_strata_legend()`](https://drostlab.github.io/myTAI/reference/get_strata_legend.md)
  : Get Strata Legend from PhyloExpressionSet
- [`match_map()`](https://drostlab.github.io/myTAI/reference/match_map.md)
  : Match Gene Expression Data with Phylostratum Map
- [`normalise_stage_expression()`](https://drostlab.github.io/myTAI/reference/normalise_stage_expression.md)
  : Normalise Stage Expression Data
- [`omit_matrix()`](https://drostlab.github.io/myTAI/reference/omit_matrix.md)
  : Compute TXI Profiles Omitting Each Gene
- [`permute_PS()`](https://drostlab.github.io/myTAI/reference/permute_PS.md)
  : Permute Strata in PhyloExpressionSet
- [`quantile_rank()`](https://drostlab.github.io/myTAI/reference/quantile_rank.md)
  : Calculate Quantile Ranks
- [`rel_exp_matrix()`](https://drostlab.github.io/myTAI/reference/rel_exp_matrix.md)
  : Compute Relative Expression Matrix for PhyloExpressionSet
- [`relative_expression()`](https://drostlab.github.io/myTAI/reference/relative_expression.md)
  : Relative Expression Functions
- [`remove_genes()`](https://drostlab.github.io/myTAI/reference/remove_genes.md)
  : Remove Genes from PhyloExpressionSet
- [`select_genes()`](https://drostlab.github.io/myTAI/reference/select_genes.md)
  : Select Genes from PhyloExpressionSet
- [`set_expression()`](https://drostlab.github.io/myTAI/reference/set_expression.md)
  : Gene Expression Transformation Functions
- [`strata_enrichment()`](https://drostlab.github.io/myTAI/reference/strata_enrichment.md)
  : Calculate Phylostratum Enrichment
- [`transform_counts()`](https://drostlab.github.io/myTAI/reference/transform_counts.md)
  : Transform Expression Counts in PhyloExpressionSet
- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) : Short
  Alias for Transform Counts
- [`tf_PS()`](https://drostlab.github.io/myTAI/reference/tf_PS.md) :
  Transform Phylostratum Values
- [`tf_stability()`](https://drostlab.github.io/myTAI/reference/tf_stability.md)
  : Perform Permutation Tests Under Different Transformations

## Utilities/Metadata

Helper functions and metadata constants.

- [`COUNT_TRANSFORMS`](https://drostlab.github.io/myTAI/reference/COUNT_TRANSFORMS.md)
  : Count Transformation Functions
- [`PS_colours()`](https://drostlab.github.io/myTAI/reference/PS_colours.md)
  : Generate Phylostratum Colors
- [`save_gatai_results_pdf()`](https://drostlab.github.io/myTAI/reference/save_gatai_results_pdf.md)
  : Save GATAI Analysis Results to PDF
- [`gatai_animate_destruction()`](https://drostlab.github.io/myTAI/reference/gatai_animate_destruction.md)
  : Animate GATAI Destruction Process
- [`print.TestResult`](https://drostlab.github.io/myTAI/reference/print.TestResult.md)
  : Print Method for TestResult
- [`rename_phyex_set()`](https://drostlab.github.io/myTAI/reference/rename_phyex_set.md)
  : Rename a PhyloExpressionSet
- [`TI_map`](https://drostlab.github.io/myTAI/reference/TI_map.md) :
  Transcriptomic Index Name Mapping
- [`TXI_conf_int()`](https://drostlab.github.io/myTAI/reference/TXI_conf_int.md)
  : Confidence Intervals for Transcriptomic Index (TXI)
- [`TXI_std_dev()`](https://drostlab.github.io/myTAI/reference/TXI_std_dev.md)
  : Standard Deviation for TXI

## Example datasets

Example datasets for demonstrating and testing package functions.

- [`example_phyex_set`](https://drostlab.github.io/myTAI/reference/example_phyex_set.md)
  : Example phyex set
- [`example_phyex_set_old`](https://drostlab.github.io/myTAI/reference/example_phyex_set_old.md)
  : Example phyex set old
- [`example_phyex_set_sc`](https://drostlab.github.io/myTAI/reference/example_phyex_set_sc.md)
  : Load Example Single-Cell PhyloExpressionSet
