#' @title Plot Transcriptome Index using bootstrapping and confidence intervals
#' @description 
#' Function to plot and compare the confidence intervals of Transcriptome Index between transformed and non-transformed expression data 
#' by using bootstrapping appraoches instead of permutation tests used in \code{\link{PlotSignature}}.
#' @details 
#' This function can be used to check potential outliers (e.g. a few exramly highly expressed genes) in transcriptome. 
#' Since Transcriptome Index is weigthed mean, it could be easily affectd by outliers. So, we have to check potential outliers in the transcriptome data.
#' Because log or sqrt trandformation can alleviate the effect of outliers, if there are some outliers, we could see the confidence intervals (genetated by bootstrapping) 
#' from non-trandformed expression data are much higher and more variable than from log or sqrt trandformed expression data.
#' In order to compare the range of confidence intervals in the same scale, we plotted the ratio of upper to lower confidence interval boundary across development.
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
#' @param measure type of transcriptome index that shall be computed. E.g. measure = "TAI" (Transcriptome Age Index), measure = "TDI" (Transcriptome Divergence Index), measure = "TPI" (Transcriptome Polymorphism Index).
#' @param nbootstraps number of independent bootstraps.
#' @param plot_type a string defining specifying the type of visualization for the results. If:
#'  \itemize{
#'    \item \code{plot_type = "ratio_and_CI"} plots of both CI ratio and the span of the CI are shown.
#'    \item \code{plot_type = "only_ratio"} only plot of the CI is shown. 
#'  }
#' @author Jialin Liu and Filipa Martins Costa
#' @examples 
#' data("PhyloExpressionSetExample")
#' PlotCIRatio(PhyloExpressionSetExample,"TAI",5)
#' @seealso \code{\link{PlotSignature}}
#' @export

PlotCIRatio <- function(ExpressionSet, measure, nbootstraps, plot_type = "ratio_and_CI") {
  if (!is.element(measure, c("TAI", "TDI", "TPI"))) {
    stop(
      "Measure '", measure, "' is not available for this function. Please specify a measure supported by this function.",
      call. = FALSE
    )
  }
  
  if (!plot_type %in% c("only_ratio", "ratio_and_CI")) {
    stop("Plot type '", plot_type, "' is not available for this function. Please specify a plot type supported by this function.", call. = FALSE)
  }
  
  ExpressionSet <- as.data.frame(ExpressionSet)
  names(ExpressionSet)[2] <- "GeneID"
  bootID <- replicate(nbootstraps, sample(ExpressionSet$GeneID, replace = TRUE))
  
  rawIndex <- c()
  logIndex <- c()
  sqrtIndex <- c()
  
  for (i in 1:nbootstraps) {
    tempID <- data.frame(bootID[ , i])
    names(tempID) <- "GeneID"
    tempData <- merge(ExpressionSet, tempID, by = "GeneID")
    tempData <- tempData[c(2, 1, 3:ncol(ExpressionSet))]
    
    if (measure == "TAI") {
      ## raw expression value
      rawIndex <- rbind(rawIndex, TAI(tempData))
      ## log2 expression value
      tempData_log2 <- tf(tempData, function(x) log2(x + 1))
      logIndex <- rbind(logIndex, TAI(tempData_log2))
      ## square root expression value
      tempData_sqrt <- tf(tempData, sqrt)
      sqrtIndex <- rbind(sqrtIndex, TAI(tempData_sqrt))
    }
    if (measure == "TDI") {
      ## raw expression value
      rawIndex <- rbind(rawIndex, TDI(tempData))
      ## log2 expression value
      tempData_log2 <- tf(tempData, function(x) log2(x + 1))
      logIndex <- rbind(logIndex, TDI(tempData_log2))
      ## square root expression value
      tempData_sqrt <- tf(tempData, sqrt)
      sqrtIndex <- rbind(sqrtIndex, TDI(tempData_sqrt))
    }
  }
  
  # calculate confidence interval ratios
  compute_ratio_and_ci <- function(index_matrix) {
    ci_upper <- apply(index_matrix, 2, function(x) stats::quantile(x, probs = 0.975))
    ci_lower <- apply(index_matrix, 2, function(x) stats::quantile(x, probs = 0.025))
    ci_ratio <- ci_upper / ci_lower
    data.frame(Ratio = ci_ratio, CI_Upper = ci_upper, CI_Lower = ci_lower)
  }
  
  rawData <- compute_ratio_and_ci(rawIndex)
  logData <- compute_ratio_and_ci(logIndex)
  sqrtData <- compute_ratio_and_ci(sqrtIndex)
  
  # create a data frame for ggplot
  devNames <- colnames(ExpressionSet)[c(-1, -2)]
  
  ratioData <- rbind(
    data.frame(CI_Ratio = rawData$Ratio, CI_Upper = rawData$CI_Upper, CI_Lower = rawData$CI_Lower,
               Transformation = "Without transformation", DevelopmentStage = colnames(ExpressionSet)[3:ncol(ExpressionSet)]),
    data.frame(CI_Ratio = logData$Ratio, CI_Upper = logData$CI_Upper, CI_Lower = logData$CI_Lower,
               Transformation = "Log transformation", DevelopmentStage = colnames(ExpressionSet)[3:ncol(ExpressionSet)]),
    data.frame(CI_Ratio = sqrtData$Ratio, CI_Upper = sqrtData$CI_Upper, CI_Lower = sqrtData$CI_Lower,
               Transformation = "Square root transformation", DevelopmentStage = colnames(ExpressionSet)[3:ncol(ExpressionSet)])
  )
  
  ratioData$Transformation <- factor(ratioData$Transformation, 
                                     levels = unique(ratioData$Transformation))
  
  ratioData$DevelopmentStage <- factor(ratioData$DevelopmentStage, 
                                       levels = devNames)
  
  # plots
  p_ratio <- ggplot2::ggplot(ratioData, ggplot2::aes(x = DevelopmentStage, y = CI_Ratio, color = Transformation, group = Transformation)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      x = "Developmental Stages",
      y = "CI Boundary Ratio",
      title = paste("Confidence Interval Ratio for", measure)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "top"
    ) +
    ggplot2::scale_x_discrete(labels = devNames)
  
  if(plot_type == "only_ratio"){
    return(p_ratio)
  }
  
  p_CI <- ggplot2::ggplot(ratioData, ggplot2::aes(x = DevelopmentStage, y = CI_Ratio, color = Transformation, group = Transformation)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_Lower, ymax = CI_Upper, fill = Transformation), alpha = 0.5, color = NA) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      x = "Developmental Stages",
      y = "CI Boundary Ratio",
      title = paste("Confidence Intervals for", measure)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "top"
    ) +
    ggplot2::scale_x_discrete(labels = devNames)
  
  plots <- cowplot::plot_grid(p_ratio, p_CI)
  print(plots)

}
