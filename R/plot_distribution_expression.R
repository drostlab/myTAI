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
#' @import ggplot2 tibble patchwork
#' @export

plot_distribution_expression <- function(phyex_set, 
                                         show_conditions=T,
                                         show_strata=F,
                                         seed=123){
  
  df <- tidyr::pivot_longer(
    phyex_set@data_collapsed,
    cols = -c(Stratum,GeneID),
    names_to = "Condition",
    values_to = "Expression"
  )
  
  qual_col_pals <-  RColorBrewer::brewer.pal.info[ RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply( RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(seed)
  colors <- sample(col_vector, ncol(phyex_set@data)-1)
  
  p <- ggplot(df, 
               aes(x = Expression)) +
      geom_density(alpha = 0.7, color = "black", fill="darkgray") +
      labs(
          x = "Density",
          y = "Expression") +
    theme_minimal()
  
  if (show_conditions) {
      p_cond <- ggplot(df, 
                   aes(x = Expression, y = factor(Condition, levels=unique(Condition)), fill = factor(Condition, levels=unique(Condition)))) +
          ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
          labs(
              x = "Density",
              y = "Expression",
              fill = phyex_set@conditions_label) +
          theme_minimal() +
          scale_fill_manual(values = colors)
      
      p <- p + p_cond
  }
  if (show_strata) {
      p_strata <- ggplot(df, 
                       aes(x = Expression, y = Stratum, fill = Stratum)) +
          ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
          labs(
              x = "Density",
              y = "Expression") +
          theme_minimal() +
          scale_fill_manual(values = PS_colours(phyex_set@num_strata))
      
      p <- p + p_strata
  }
  p + plot_layout(nrow=1)
}
