
#' @title Plot Expression Levels by Phylostratum
#' @description Create a boxplot showing the distribution of expression levels
#' for each phylostratum.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param aggregate_FUN Function to aggregate expression across conditions (default: mean)
#' 
#' @return A ggplot2 object showing expression distributions by phylostratum
#' 
#' @details
#' This function creates a boxplot visualization showing how expression levels
#' vary across different phylostrata. Each point represents a gene, and the
#' boxes show the distribution of expression levels within each phylostratum.
#' 
#' @examples
#' # Plot expression by strata using mean aggregation
#' # p1 <- plot_strata_expression(phyex_set, aggregate_FUN = mean)
#' 
#' # Plot using median aggregation
#' # p2 <- plot_strata_expression(phyex_set, aggregate_FUN = median)
#' 
#' @import ggplot2 
#' @export
plot_strata_expression <- function(phyex_set,
                                   aggregate_FUN = mean) {

    agg_name <- deparse(substitute(aggregate_FUN))
    aggregate_FUN <- match.fun(aggregate_FUN)
    df <- phyex_set@data |>
        mutate(Expression = apply(phyex_set@counts_collapsed, 1, aggregate_FUN))
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

#' @title Plot Expression Rank by Phylostratum
#' @description Create a boxplot showing the distribution of expression ranks
#' for each phylostratum.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param aggregate_FUN Function to aggregate expression across conditions (default: mean)
#' 
#' @return A ggplot2 object showing expression rank distributions by phylostratum
#' 
#' @details
#' This function creates a boxplot visualization showing how expression ranks
#' vary across different phylostrata. Expression values are first aggregated
#' and then ranked, providing a normalized view of expression patterns.
#' 
#' @examples
#' # Plot expression ranks by strata
#' # p <- plot_strata_expression_rank(phyex_set)
#' 
#' @import ggplot2
#' @export 
plot_strata_expression_rank <- function(phyex_set,
                                   aggregate_FUN = mean) {
    
    agg_name <- deparse(substitute(aggregate_FUN))
    aggregate_FUN <- match.fun(aggregate_FUN)
    df <- phyex_set@data |>
        mutate(Expression = apply(phyex_set@counts_collapsed, 1, aggregate_FUN)) |>
        arrange(Expression) |>
        mutate(Rank = row_number()) |>
        mutate(PS_num = as.numeric(Stratum)) |>
        mutate(AvgPS = zoo::rollapply(PS_num, width=20, FUN=mean, align= "center", fill="extend"))
    
    p <- ggplot(df, aes(x=Rank, y=Expression, colour=AvgPS), alpha=0.3) +
        geom_point() + 
        scale_colour_gradient(low = "white", high = "blue") +
        labs(y=paste("Expresssion aggregated by", agg_name)) + 
        theme_minimal()
    
    return(p)
}