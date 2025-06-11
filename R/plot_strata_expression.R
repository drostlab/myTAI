
#' @import ggplot2 
plot_strata_expression <- function(phyex_set,
                                   aggregate_FUN = mean) {

    agg_name <- deparse(substitute(aggregate_FUN))
    aggregate_FUN <- match.fun(aggregate_FUN)
    df <- phyex_set@data |>
        mutate(Expression = apply(phyex_set@counts_collapsed, 1, aggregate_FUN))
    
    p <- ggplot(df, aes(x=Stratum, y=Expression, colour=Stratum)) +
        geom_jitter() + 
        scale_colour_manual(values = PS_colours(phyex_set@num_strata)) +
        labs(y=paste("Expresssion aggregated by", agg_name)) + 
        guides(colour="none") +
        theme_minimal()
    
    return(p)
}

#' @import ggplot2 
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