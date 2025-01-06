#' @title Comparing expression/partial TAI distributions across the different developmental stages
#' @description \emph{plot_distribution_expression_partialTAI} generates 2 plots that help to compare the distribution
#' of the quotient of expression by partial TAI through various developmental stages, highlighting each stage with 
#' distinct colors.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param stages a numeric vector specifying the indices of the stages to compare. Each index 
#' corresponds to a stage in the ExpressionSet.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @author Filipa Martins Costa
#' @export

plot_distribution_expression_partialTAI <- function(ExpressionSet,
                                                    stages = 1:ncol(ExpressionSet),
                                                    xlab = "Expression / Partial TAI",
                                                    ylab = "Density",
                                                    main = "Density Distribution of Expression / Partial TAI by Developmental Stage"){
  
  if (any(stages > ncol(ExpressionSet))) {
    stop("Some indices in 'stages' exceed the number of columns in 'ExpressionSet'.")
  }
  
  partial_TAI_matrix <- pMatrix(ExpressionSet[stages])

  partial_TAI_df <- tibble::rownames_to_column(as.data.frame(partial_TAI_matrix), var = "GeneID")
  partial_TAI_long <- tidyr::pivot_longer(
    partial_TAI_df,
    cols = -GeneID,
    names_to = "Stage",
    values_to = "PartialTAI"
  )
  
  expression_long <- tidyr::pivot_longer(
    ExpressionSet[, 2:length(stages)],
    cols = -GeneID,
    names_to = "Stage",
    values_to = "Expression"
  )
  
  combined_long <- partial_TAI_long
  combined_long <- dplyr::inner_join(combined_long, expression_long, by = c("GeneID", "Stage"))

  #combined_long <- combined_long%>% dplyr::mutate(Ratio = Expression / PartialTAI)
  combined_long <- dplyr::mutate(
    combined_long,
    Ratio = Expression / PartialTAI
  )

  combined_long$Stage <- factor(combined_long$Stage, levels = colnames(partial_TAI_matrix))
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(123)
  colors <- sample(col_vector, ncol(partial_TAI_df))
  
  P1 <- ggplot2::ggplot(combined_long, ggplot2::aes(x = Ratio, fill = Stage)) +
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
  
  P2 <- ggplot2::ggplot(combined_long, ggplot2::aes(x = Ratio, y = Stage, fill = Stage)) +
    ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
    ggplot2::labs(
      x = xlab,
      y =ylab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_fill_manual(values = colors)
  
  
  cowplot::plot_grid(P1, P2, labels = main)
  
}
