## ----message = FALSE, warning = FALSE-----------------------------------------
library(myTAI); library(S7); library(ggplot2); library(patchwork)
data("example_phyex_set")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output with stat_flatline_test", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature(example_phyex_set, 
                      show_p_val = TRUE, 
                      conservation_test = stat_flatline_test,
                      colour = "lavender") +
  # as the plots are ggplot2 objects, we can simply modify them using ggplot2
  ggplot2::labs(title = "Developmental stages of A. thaliana")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output with stat_reductive_hourglass_test", dev.args = list(bg = 'transparent'), fig.align='center'----
module_info <- list(early = 1:3, mid = 4:6, late = 7:8)
myTAI::plot_signature(example_phyex_set,
                      show_p_val = TRUE,
                      conservation_test = stat_reductive_hourglass_test,
                      modules = module_info,
                      colour = "lavender")

## ----message = FALSE, fig.height=5, fig.width=8, fig.alt="plot_signature_transformed function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature_transformed(
  example_phyex_set)

## ----message = FALSE, fig.height=5, fig.width=8, fig.alt="plot_signature_gene_quantiles function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature_gene_quantiles(
  example_phyex_set)

## ----message = FALSE, warning = FALSE, fig.height=3, fig.width=5, fig.alt="stat_flatline_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::stat_flatline_test(
  example_phyex_set, plot_result = TRUE)

## ----message = FALSE, warning = FALSE-----------------------------------------
res_flt <- myTAI::stat_flatline_test(example_phyex_set, plot_result = FALSE)

## ----message = FALSE, warning = FALSE, fig.height=5, fig.width=5, fig.alt="plot_cullen_frey function output for stat_flatline_test", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_cullen_frey(res_flt)

## ----message = FALSE, warning = FALSE, fig.height=5, fig.width=5, fig.alt="plot_null_txi_sample function output for stat_flatline_test", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_null_txi_sample(res_flt) +
  ggplot2::guides(x =  guide_axis(angle = 90))

## ----message = FALSE, warning = FALSE, fig.height=3, fig.width=5, fig.alt="stat_reductive_hourglass_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
module_info <- list(early = 1:3, mid = 4:6, late = 7:8)
myTAI::stat_reductive_hourglass_test(
  example_phyex_set, plot_result = TRUE,
  modules = module_info)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_strata_expression function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_strata_expression(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_strata_expression function output ggplot2", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_strata_expression(example_phyex_set) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "Expression aggregated by mean (log-scaled)")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=8, fig.alt="plot_strata_expression function output", dev.args = list(bg = 'transparent'), fig.align='center'----
library(patchwork)
p1 <- myTAI::plot_strata_expression(example_phyex_set |> myTAI::tf(log1p))

# equivalent to 
p2 <- example_phyex_set |> myTAI::tf(log1p) |> myTAI::plot_strata_expression() 

p1+p2

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_contribution function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_contribution(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=8, fig.alt="plot_distribution_expression function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_distribution_expression(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=6, fig.width=9, fig.alt="plot_distribution_pTAI_qqplot function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_distribution_pTAI_qqplot(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_distribution_strata function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_distribution_strata(example_phyex_set@strata) /
myTAI::plot_distribution_strata(
  example_phyex_set@strata,
  selected_gene_ids = myTAI::genes_top_variance(example_phyex_set, top_p = 0.95),
  as_log_obs_exp = TRUE
) + plot_annotation(title = "Distribution of gene ages (top), Observed vs Expected plot of top 5% variance genes (bottom)")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=12, fig.width=12, fig.alt="plot_gene_heatmap function output default", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_gene_heatmap(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=12, fig.width=12, fig.alt="plot_gene_heatmap function output clustered", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_gene_heatmap(example_phyex_set, cluster_rows = TRUE, show_reps=TRUE, show_gene_ids=TRUE, top_p=0.005)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=12, fig.width=10, fig.alt="plot_gene_heatmap function output nonstd", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_gene_heatmap(example_phyex_set, cluster_rows = TRUE, show_reps=TRUE, top_p=0.005, std=FALSE, show_gene_ids=TRUE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_gene_space function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_gene_space(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_gene_space function output by strata", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_gene_space(example_phyex_set,colour_by = "strata")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=12, fig.alt="plot_sample_space function output by TXI", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_sample_space(example_phyex_set) | myTAI::plot_sample_space(example_phyex_set, colour_by = "TXI")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_sample_space function output by TXI", dev.args = list(bg = 'transparent'), fig.align='center', eval = requireNamespace("uwot", quietly = TRUE)----
# we can even do a UMAP
myTAI::plot_sample_space(example_phyex_set, method = "UMAP")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=3, fig.width=8, fig.alt="plot_mean_var function output simple vs highlighted", dev.args = list(bg = 'transparent'), fig.align='center'----

# highlighting top variance genes
top_var_genes <- myTAI::genes_top_variance(example_phyex_set, top_p = 0.9995)
p1 <- myTAI::plot_mean_var(example_phyex_set)
p2 <- myTAI::plot_mean_var(example_phyex_set, 
                     highlight_genes = top_var_genes)

p1 + p2 + plot_annotation(title = "Mean-variance: simple vs. highlighted top variance genes")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=3, fig.width=6, fig.alt="plot_gene_space function output log transform coloured by strata", dev.args = list(bg = 'transparent'), fig.align='center'----
# with log transform and colouring by phylostratum
myTAI::plot_mean_var(example_phyex_set |> myTAI::tf(log1p), 
                     colour_by = "strata") +
  ggplot2::guides(colour = guide_legend(ncol=2))

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=12, fig.alt="plot_gene_profiles function output manual vs strata coloring", dev.args = list(bg = 'transparent'), fig.align='center'----
# side by side: manual coloring vs strata coloring
p1 <- myTAI::plot_gene_profiles(example_phyex_set, max_genes = 10, colour_by = "manual")
p2 <- myTAI::plot_gene_profiles(example_phyex_set, max_genes = 10, colour_by = "strata")

p1 + p2 + plot_annotation(title = "Gene profiles: manual vs. strata coloring")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_gene_profiles function output stage std_log transformation", dev.args = list(bg = 'transparent'), fig.align='center'----
# stage colouring with standardized log transformation
myTAI::plot_gene_profiles(example_phyex_set, max_genes = 10, 
                          transformation = "std_log", colour_by = "stage")

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=8, fig.width=12, fig.alt="plot_gene_profiles function output faceted by strata", dev.args = list(bg = 'transparent'), fig.align='center'----
# faceted by phylostratum
myTAI::plot_gene_profiles(example_phyex_set, max_genes = 1000, 
                          colour_by = "strata", facet_by_strata = TRUE, show_set_mean = TRUE,
                          show_labels = FALSE)

## ----message = FALSE, warning = FALSE-----------------------------------------
# Load example single-cell data
data(example_phyex_set_sc)

## ----message = FALSE, warning = FALSE-----------------------------------------
example_phyex_set_sc

## ----message = FALSE, warning = FALSE-----------------------------------------
# Check available identities
cat("Available identities for plotting:\n")
print(example_phyex_set_sc@available_idents)

## ----message = FALSE, warning = FALSE-----------------------------------------
# Set up custom color schemes for better visualization
day_colors <- c("Day1" = "#3498db", "Day3" = "#2980b9", "Day5" = "#1f4e79", "Day7" = "#0d2a42")
condition_colors <- c("Control" = "#27ae60", "Treatment" = "#e74c3c")
group_colors <- c("TypeA" = "#e74c3c", "TypeB" = "#f39c12", "TypeC" = "#9b59b6")

example_phyex_set_sc@idents_colours[["day"]] <- day_colors
example_phyex_set_sc@idents_colours[["condition"]] <- condition_colors
example_phyex_set_sc@idents_colours[["groups"]] <- group_colors

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature single-cell basic", dev.args = list(bg = 'transparent'), fig.align='center'----
# Basic signature plot showing TXI distribution across cell types
myTAI::plot_signature(example_phyex_set_sc)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature single-cell without individual cells", dev.args = list(bg = 'transparent'), fig.align='center'----
# Plot without showing individual cells (just means)
myTAI::plot_signature(example_phyex_set_sc, show_reps = FALSE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature single-cell by day", dev.args = list(bg = 'transparent'), fig.align='center'----
# Plot TXI distribution by developmental day instead of cell type
myTAI::plot_signature(example_phyex_set_sc, primary_identity = "day", show_p_val = FALSE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature single-cell by condition", dev.args = list(bg = 'transparent'), fig.align='center'----
# Plot TXI distribution by experimental condition
myTAI::plot_signature(example_phyex_set_sc, primary_identity = "condition", show_p_val=FALSE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=4, fig.width=8, fig.alt="plot_signature single-cell with secondary coloring", dev.args = list(bg = 'transparent'), fig.align='center'----
# Plot by day, colored by condition
myTAI::plot_signature(example_phyex_set_sc, 
                     primary_identity = "day", 
                     secondary_identity = "condition",
                     show_p_val=FALSE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=6, fig.width=10, fig.alt="plot_signature single-cell with faceting", dev.args = list(bg = 'transparent'), fig.align='center'----
# Plot by day, faceted by condition
myTAI::plot_signature(example_phyex_set_sc, 
                     primary_identity = "day", 
                     secondary_identity = "batch",
                     facet_by_secondary = TRUE,
                     show_p_val = FALSE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=8, fig.width=10, fig.alt="plot_gene_heatmap single-cell", dev.args = list(bg = 'transparent'), fig.align='center'----
# Gene heatmap for single-cell data (aggregated by cell type)
myTAI::plot_gene_heatmap(example_phyex_set_sc, top_p = 0.1, cluster_rows=TRUE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=8, fig.width=12, fig.alt="plot_gene_heatmap single-cell with individual cells", dev.args = list(bg = 'transparent'), fig.align='center'----
# Gene heatmap showing individual cells (subsampled)
myTAI::plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, max_cells_per_type = 10, top_p = 0.05, cluster_rows=TRUE)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=8, fig.width=12, fig.alt="plot_gene_heatmap single-cell grouped by day", dev.args = list(bg = 'transparent'), fig.align='center'----
# Change identity to "day" and plot heatmap grouped by developmental time
example_sc_by_day <- example_phyex_set_sc
example_sc_by_day@selected_idents <- "day"
myTAI::plot_gene_heatmap(example_sc_by_day, show_reps = TRUE, max_cells_per_type = 8, top_p = 0.05, cluster_rows=TRUE, show_gene_ids=TRUE, std=FALSE)

