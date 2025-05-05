
#' @import ggplot2 
plot_strata_expression <- function(phyex_set,
                                   aggregate_FUN = mean) {

    agg_name <- deparse(substitute(aggregate_FUN))
    aggregate_FUN <- match.fun(aggregate_FUN)
    df <- phyex_set@data |>
        mutate(Expression = apply(phyex_set@count_matrix, 1, aggregate_FUN))
    
    p <- ggplot(df, aes(x=Stratum, y=Expression, colour=Stratum)) +
        geom_boxplot() + 
        scale_colour_manual(values = PS_colours(phyex_set@num_strata)) +
        labs(y=paste("Expresssion aggregated by", agg_name)) + 
        guides(colour="none") +
        theme_minimal()
    
    return(p)
}