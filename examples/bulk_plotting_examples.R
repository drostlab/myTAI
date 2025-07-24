#' Bulk PhyloExpressionSet Plotting Examples
#' 
#' This script demonstrates the main plotting functions available for 
#' BulkPhyloExpressionSet objects using the example_phyex_set dataset.
#' 
#' The example_phyex_set contains Arabidopsis thaliana embryogenesis data
#' with 8 developmental stages (3 replicates each) and 27,520 genes.

library(myTAI)

# Load the example bulk dataset
data("example_phyex_set")
example_phyex_set


# Signature plot - Show overall transcriptomic index patterns
plot_signature(example_phyex_set, show_reps = TRUE, show_p_val = TRUE)


plot_gene_heatmap(example_phyex_set, top_p = 0.01, cluster_rows=TRUE)

plot_gene_profiles(example_phyex_set, max_genes = 10)


plot_sample_space(example_phyex_set, method = "PCA")


plot_gene_space(example_phyex_set)



plot_distribution_expression(example_phyex_set |> tf(log1p), show_strata=TRUE)



plot_contribution(example_phyex_set)


plot_strata_expression(example_phyex_set |> tf(log1p))

plot_mean_var(example_phyex_set)
