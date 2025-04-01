
#' @import ggplot2
plot_gene_profiles <- function(phyex_set,
                               genes=phyex_set@gene_ids) {
    phyex_set <- select_genes(phyex_set, genes)
    df <- as_tibble(phyex_set@count_matrix) |>
        add_column(Phylostratum = phyex_set@strata_vector,
                   GeneID = phyex_set@gene_ids) |>
        tidyr::pivot_longer(cols=phyex_set@conditions, names_to="Condition", values_to="Value")
        
    ggplot(data=df, 
           aes(x=factor(Condition, levels=unique(Condition)), 
               y=Value, 
               colour=Phylostratum, 
               group=GeneID)) + 
        geom_line(alpha=0.6, lwd=0.8) +
        labs(x=phyex_set@conditions_label) +
        theme_minimal()

}