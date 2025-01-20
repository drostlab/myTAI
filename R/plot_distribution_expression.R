#' @title Comparing expression levels distributions across the different developmental stages
#' @description \emph{plot_distribution_partialTAI} generates 2 plots that help to compare the distribution
#' of expression levels through various developmental stages, highlighting each stage with 
#' distinct colors.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param stages a numeric vector specifying the indices of the stages to compare. Each index 
#' corresponds to a stage in the ExpressionSet. Starts in one.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param seed defines the colors for the different developmetal stages
#' @section Recomendation - Apply a square root transformation to enhance the visualization of differences 
#' in the distributions: plot_distribution_partialTAI(tf(ExpressionSet, sqrt))
#' @author Filipa Martins Costa
#' @export

plot_distribution_expression <- function(ExpressionSet,
                                         stages = 1:(ncol(ExpressionSet)-2),
                                         xlab = "Expression",
                                         ylab = "Density",
                                         main = "Density Distribution of Expression by Developmental Stage",
                                         seed = 123){
  
  is.ExpressionSet(ExpressionSet)
  
  if (any(stages > ncol(ExpressionSet))) {
    stop("Some indices in 'stages' exceed the number of columns in 'ExpressionSet'.")
  }
  
  expression_long <- tidyr::pivot_longer(
    ExpressionSet[,c(2,stages+2)],
    cols = -GeneID,
    names_to = "Stage",
    values_to = "Expression"
  )
  
  qual_col_pals <-  RColorBrewer::brewer.pal.info[ RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply( RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(seed)
  colors <- sample(col_vector, ncol(ExpressionSet)-1)
  
  P1 <- ggplot2::ggplot(expression_long, ggplot2::aes(x = Expression, fill = Stage)) +
    ggplot2::geom_density(alpha = 0.7, color = "black") +
    ggplot2::labs(
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_fill_manual(values = colors)
  
  P2 <- ggplot2::ggplot(expression_long, ggplot2::aes(x = Expression, y = Stage, fill = Stage)) +
    ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
    ggplot2::labs(
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_fill_manual(values = colors)

  cowplot::plot_grid(P1, P2, labels = main)
}
