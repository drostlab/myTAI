#' Single-cell PhyloExpressionSet Plotting Examples
#' 
#' This script demonstrates the main plotting functions available for 
#' ScPhyloExpressionSet objects using the example_phyex_set_sc dataset.
#' 
#' The example_phyex_set_sc is a synthetic single-cell dataset designed
#' for testing single-cell phylotranscriptomic functionality.

library(myTAI)

# Load the example single-cell dataset
example_phyex_set_sc <- load_example_phyex_set_sc()
example_phyex_set_sc


# Check dataset structure
cat("Number of cells:", example_phyex_set_sc@num_samples, "\n")
cat("Number of genes:", length(example_phyex_set_sc@gene_ids), "\n")
cat("Sample groups:", unique(example_phyex_set_sc@groups), "\n")


# SIGNATURE PLOTS - Transcriptomic index patterns in single cells
plot_signature(example_phyex_set_sc, show_reps = TRUE, show_p_val=TRUE)


plot_gene_heatmap(example_phyex_set_sc, top_p = 0.01)

plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, top_p = 0.005)



plot_gene_profiles(example_phyex_set_sc, max_genes = 8, show_reps = FALSE)




plot_sample_space(example_phyex_set_sc, method = "UMAP", colour_by="TXI")


plot_gene_space(example_phyex_set_sc)


plot_distribution_expression(example_phyex_set_sc)

plot_distribution_expression(example_phyex_set_sc, show_strata = TRUE)

plot_contribution(example_phyex_set_sc)

plot_strata_expression(example_phyex_set_sc)

plot_mean_var(example_phyex_set_sc, colour_by = "strata")

