#' @title Plot the Mean Gene Expression Profiles with Optional Standard Deviation and Labels
#' @description
#' This function visualizes the mean gene expression profiles for the stages defined in the input \code{ExpressionSet}.
#' It can also display the standard deviation as a shaded area around the mean, add labels with the expression values at each stage, and customize various plot parameters.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object containing gene expression data.
#' @param y_ticks a numeric value specifying the number of ticks to be drawn on the y-axis. Default is 6.
#' @param digits a numeric value specifying the number of digits to display in the expression value labels. Default is 2.
#' @param color a string specifying the color of the line and points representing the mean expression. Default is \code{"#009999"}.
#' @param line_width a numeric value specifying the line width for the mean expression line. Default is 2.
#' @param point_size a numeric value specifying the size of the points representing the mean expression values. Default is 3.
#' @param add_values a logical value indicating whether expression values should be displayed on the plot. Default is \code{TRUE}.
#' @param xlab a string specifying the label for the x-axis. Default is \code{"Ontogeny"}.
#' @param ylab a string specifying the label for the y-axis. Default is \code{"Mean Expression Level"}.
#' @param add_sd a logical value indicating whether the standard deviation should be shown as a shaded ribbon around the mean expression. Default is \code{FALSE}.
#' @param yaxis_range a numeric value specifying the extra range added to the y-axis when standard deviation is shown. Default is 4.
#' @param shadow_color a string specifying the color of the shaded area representing the standard deviation. Default is \code{"grey"}.
#' 
#' @examples
#' 
#' #In case you want the expression values show together with the standard deviation
#' plot_gene_set_mean(ExpressionSet = PhyloExpressionSetExample,
#'                   add_sd = TRUE)
#'                   
#' #In case you would like for the expression values and the standard deviation don't show up
#' plot_gene_set_mean(ExpressionSet = PhyloExpressionSetExample,
#'                   add_values = FALSE)
#' 
#' @author Filipa Martins Costa
#' @export

plot_gene_set_mean <- function(ExpressionSet,
                               gene_set,
                               y_ticks = 6,
                               digits = 2,
                               color = "#009999",
                               line_width  = 2,
                               point_size = 3,
                               add_expression_values = T,
                               xlab = "Ontogeny",
                               ylab = "Mean Expression Level",
                               add_sd = F,
                               yaxis_range = 4,
                               shadow_color = "grey"){
    
    is.ExpressionSet(ExpressionSet)
    
    if (!is.vector(gene_set))
        stop("Please provide a valid character vector as input gene_ids.", call. = FALSE)
    
    GeneSubSet_indixes <- stats::na.omit(
        match(
            tolower(gene_set),
            dplyr::pull(
                dplyr::filter(
                    dplyr::mutate(ExpressionSet, GeneID = tolower(GeneID)),
                    GeneID %in% tolower(gene_set)
                ),
                GeneID
            )
        )
    )
    
    if (length(GeneSubSet_indixes) == 0)
        stop ("None of your input gene ids could be found in the ExpressionSet.", call. = FALSE)
    
    if (length(GeneSubSet_indixes) != length(gene_set))
        warning ("Only ",length(GeneSubSet_indixes), " out of your ", length(gene_set), " gene ids could be found in the ExpressionSet.", call. = FALSE)
    
    GeneStats <- ExpressionSet[ExpressionSet$GeneID %in% gene_set,2:ncol(ExpressionSet)] |>
        dplyr::select(-GeneID) |>
        tidyr::gather(key = "Stage", value = "Expression") |>
        dplyr::group_by(Stage) |>
        dplyr::summarize(
            MeanExpression = mean(Expression, na.rm = TRUE), 
            SDExpression = sd(Expression, na.rm = TRUE)
        ) |>
        dplyr::ungroup()
    
    GeneStats$Stage <- factor(GeneStats$Stage, levels = colnames(ExpressionSet)[3:ncol(ExpressionSet)])
    
    yaxis_min = min(GeneStats$MeanExpression - GeneStats$SDExpression)
    yaxis_max = max(GeneStats$MeanExpression + GeneStats$SDExpression)
    
    p <- ggplot2::ggplot(GeneStats, ggplot2::aes(x = Stage, y = MeanExpression), group = 1)
    
    if (add_sd){
        p <- p + ggplot2::geom_ribbon(
            ggplot2::aes(ymin = MeanExpression - SDExpression, ymax = MeanExpression + SDExpression),
            fill = shadow_color, group = 1
        )
    }
    
    p <- p + ggplot2::geom_line(ggplot2::aes(group = 1), color = color, linewidth = line_width) +   
        ggplot2::geom_point(size = point_size, color = color) +                              
        ggplot2::labs(x = xlab, y = ylab) + 
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks))
    
    if (add_sd){
        p <- p + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks), 
                                             limits = c(yaxis_min-yaxis_range, yaxis_max+yaxis_range))
    }
    
    if (add_expression_values) {
        p <- p + ggplot2::geom_text( 
            ggplot2::aes(label = round(MeanExpression, digits = digits)),  
            size = 3, 
            vjust = -0.5, 
            hjust = 0.5,   
            nudge_y = 0.05
        )
    }
    
    print(p)
}

