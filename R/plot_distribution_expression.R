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
#' @import ggplot2 tibble
#' @export

plot_distribution_expression <- function(phyex_set,
                                         conditions = phyex_set@conditions,
                                         xlab = "Expression",
                                         ylab = "Density",
                                         main = "Density Distribution of Expression by Developmental Stage",
                                         seed = 123){
  if(!all(conditions %in% phyex_set@conditions))
      stop("Some elements in `conditions` do not occur in the `phyex_set`")
  
  expression_long <- tidyr::pivot_longer(
    tibble(GeneID=phyex_set@strata_vector, as_tibble(phyex_set@count_matrix)),
    cols = -GeneID,
    names_to = "Stage",
    values_to = "Expression"
  )
  
  qual_col_pals <-  RColorBrewer::brewer.pal.info[ RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply( RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(seed)
  colors <- sample(col_vector, ncol(phyex_set@data)-1)
  
  P1 <- ggplot(expression_long, 
               aes(x = Expression, fill = Stage)) +
      geom_density(alpha = 0.7, color = "black") +
      labs(
          x = xlab,
          y = ylab) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
    scale_fill_manual(values = colors)
  
  P2 <- ggplot(expression_long, 
               aes(x = Expression, y = Stage, fill = Stage)) +
      ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
      labs(
          x = xlab,
          y = ylab) +
      theme_minimal() +
      theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 12)) +
    scale_fill_manual(values = colors)
  
  cowplot::plot_grid(P1, P2, labels = main)
}
