#' @title QQ plot comparing partial TAI distributions across developmental stages against a reference stage
#' @description \emph{plot_distribution_partialTAI_qqplot} generates a QQ plot to compare the partial TAI 
#' distributions of various developmental stages against a reference stage.
#' It visualizes quantile differences between the reference and other stages,
#' highlights each stage with distinct colors, and annotates the plot with the p-values
#' from the nonparametric \code{\link{ks.test}} to indicate the significance of distribution differences.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param reference_stage_index an integer specifying the index of the reference developmental stage 
#' of the ExpressionSet. The partial TAI distribution of this stage will be used as the reference 
#' for comparisons with other stages (as default stage index 1).
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param alpha transparency of the point.
#' @param size size of the points.
#' @author Filipa Martins Costa
#' @export

plot_distribution_partialTAI_qqplot <- function(ExpressionSet, 
                               reference_stage_index = 1,
                               xlab = "Quantiles of Reference Stage",
                               ylab = "Quantiles of Other Stages",
                               main = "QQ Plot: Developmental Stages vs Reference Stage (p-values from KS test)",
                               alpha = 0.7,
                               size = 1.2) {
  
  partial_TAI_matrix <- pMatrix(ExpressionSet)
  
  if (reference_stage_index > ncol(partial_TAI_matrix)){
    stop("The specified reference stage index exceeds the number of stages in the ExpressionSet. Please provide a valid index.")
  }
  
  # define the reference stage
  reference_stage <- partial_TAI_matrix[, reference_stage_index]
  
  # function to compute QQ data and KS p-value for a given stage
  compute_qq_ks <- function(stage_index) {
    stage_values <- partial_TAI_matrix[, stage_index]
    
    # compute qq plot data
    qq <- qqplot(reference_stage, stage_values, plot.it = FALSE)
    
    # perform KS test
    ks_p_value <- ks.test(reference_stage, stage_values)$p.value
    
    # return a list containing qq data and KS p-value
    list(
      qq_data = data.frame(
        reference = qq$x,  # quantiles of the reference stage
        stage_quantiles = qq$y,  # quantiles of the current stage
        stage = colnames(partial_TAI_matrix)[stage_index]  # stage name
      ),
      ks_result = data.frame(
        stage = colnames(partial_TAI_matrix)[stage_index],
        p_value = ks_p_value
      )
    )
  }
  
  # apply the function to all stages except the reference
  results <- purrr::map(setdiff(1:ncol(partial_TAI_matrix), reference_stage_index), compute_qq_ks)
  
  # combine QQ plot data
  qq_data <- dplyr::bind_rows(purrr::map(results, "qq_data"))
  
  # combine KS test results
  ks_results <- dplyr::bind_rows(purrr::map(results, "ks_result"))
  
  # add p-values to stage labels
  ks_results <- dplyr::mutate(
    ks_results,
    stage_label = paste(stage, "(p-value =", format(p_value, digits = 3), ")")
  )
  
  # merge stage labels into the qq_data frame
  qq_data <- dplyr::left_join(
    qq_data,
    ks_results[, c("stage", "stage_label")],
    by = "stage"
  )
  qq_data <- dplyr::mutate(
    qq_data,
    stage_label = factor(stage_label, levels = unique(ks_results$stage_label))
  )
  
  ggplot2::ggplot(qq_data, ggplot2::aes(x = reference, y = stage_quantiles, color = stage_label)) +
    ggplot2::geom_point(alpha = alpha, size = size) +  
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
    ggplot2::labs(
      title = main,
      x = xlab,
      y = ylab,
      color = "Developmental Stage (p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_color_brewer(palette = "Set1")
}
