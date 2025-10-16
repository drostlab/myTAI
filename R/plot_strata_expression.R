
#' @title Plot Expression Levels by Phylostratum
#' @description Create a boxplot showing the distribution of expression levels
#' for each phylostratum.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param aggregate_FUN Function to aggregate expression across identities (default: mean)
#' 
#' @return A ggplot2 object showing expression distributions by phylostratum
#' 
#' @details
#' This function creates a boxplot visualization showing how expression levels
#' vary across different phylostrata. Each point represents a gene, and the
#' boxes show the distribution of expression levels within each phylostratum.
#' 
#' For bulk data, expression is aggregated across developmental stages.
#' For single-cell data, expression is aggregated across cell types.
#' 
#' @examples
#' # Plot expression by strata using mean aggregation for bulk data
#' p1 <- plot_strata_expression(example_phyex_set, aggregate_FUN = mean)
#' 
#' # Plot using median aggregation for single-cell data
#' p2 <- plot_strata_expression(example_phyex_set_sc, aggregate_FUN = median)
#' 
#' @import ggplot2 
#' @export
plot_strata_expression <- function(phyex_set,
                                   aggregate_FUN = mean) {
    check_PhyloExpressionSet(phyex_set)

    agg_name <- deparse(substitute(aggregate_FUN))
    aggregate_FUN <- match.fun(aggregate_FUN)
    
    # Create data frame with gene info and aggregated expression
    df <- data.frame(
        GeneID = phyex_set@gene_ids,
        Stratum = phyex_set@strata,
        Expression = apply(phyex_set@expression_collapsed, 1, aggregate_FUN)
    )
    p <- ggplot(df, aes(x=Stratum, y=Expression, colour=Stratum, fill=Stratum)) +
        geom_jitter(size=0.2, width=0.2, alpha=0.3) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3, size=0.1, colour="black") +

        scale_colour_manual(values = PS_colours(phyex_set@num_strata)) +
        scale_fill_manual(values = PS_colours(phyex_set@num_strata)) +
        labs(y=paste("Expresssion aggregated by", agg_name)) + 
        guides(colour="none", fill="none") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
}
