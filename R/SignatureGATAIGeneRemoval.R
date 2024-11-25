#' @title Removed genes after GATAI gene removal and plots of evolutionary signatures across transcriptomes before, after GATAI gene removal and with randomly removed genes
#' @description Main function to return the removed genes after GATAI gene removal and plots transcriptome indices pattern before, after GATAI gene removal and by removing random genes (by default same number as GATAI removed).
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
#' @param plot_type a string defining specifying the type of visualization for the results. If:
#'  \itemize{
#'    \item \code{plot_type = "separate"} individual plots are generated for the original dataset, the dataset with genes removed by GATAI, and the dataset with randomly removed genes.
#'    \item \code{plot_type = "combined"} all datasets are visualized together in a single combined plot. 
#'  }
#' @param measure type of transcriptome index that shall be computed. E.g.
#' \itemize{
#' \item \code{measure = "TAI"} (Transcriptome Age Index)
#' \item \code{measure = "TDI"} (Transcriptome Divergence Index)
#' \item \code{measure = "TPI"} (Transcriptome Polymorphism Index)
#' }
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be:
#'  \itemize{
#' \item \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line
#' \item \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern)
#' \item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a early conservation pattern (low-high-high pattern)
#' \item \code{TestStatistic} = \code{"LateConservationTest"} : Statistical test for the existence of a late conservation pattern (high-high-low pattern)
#' \item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
#' }
#' @param modules a list storing three elements for the \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}},
#' or \code{\link{ReverseHourglassTest}}: early, mid, and late.
#' Each element expects a numeric vector specifying the developmental stages
#' or experiments that correspond to each module. For example:
#' \itemize{
#' \item \code{module} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} divides a dataset storing seven developmental stages into 3 modules.
#' }
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}}, \code{\link{ReductiveHourglassTest}} or \code{\link{ReverseHourglassTest}}.
#' @param lillie.test a boolean value specifying whether the Lilliefors Kolmogorov-Smirnov Test shall be performed.
#' @param p.value a boolean value specifying whether the p-value of the test statistic shall be printed as a subtitle.
#' @param shaded.area a boolean value specifying whether a shaded area shall
#' be drawn for the developmental stages defined to be the presumptive phylotypic period.
#' @param custom.perm.matrix a custom \code{\link{bootMatrix}} (permutation matrix) to perform the underlying test statistic visualized by \code{PlotSignature}. Default is \code{custom.perm.matrix = NULL}.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param lwd line width.
#' @param alpha transparency of the shaded area (between [0,1]). Default is \code{alpha = 0.1}.
#' @param y.ticks number of ticks on the y-axis. Default is \code{ticks = 10}.
#' @param n_random_removal number of randomly removed genes for the third plot. Default is \code{"gatai"} (same number of genes removed by GATAI).
#' @author Filipa Martins Costa
#' @export
SignatureGATAIGeneRemoval <- function(ExpressionSet,
                                      plot_type = "separate",
                                      measure = "TAI",
                                      TestStatistic = "FlatLineTest",
                                      modules = NULL,
                                      permutations = 1000,
                                      lillie.test = FALSE,
                                      p.value = TRUE,
                                      shaded.area = FALSE,
                                      custom.perm.matrix = NULL,
                                      xlab = "Ontogeny",
                                      ylab = "Transcriptome Index",
                                      main = "",
                                      lwd = 4,
                                      alpha = 0.1,
                                      y.ticks = 10,
                                      n_random_removal = "gatai",
                                      ...) {
  
  if (!plot_type %in% c("separate", "combined")) {
    stop("Plot type '", plot_type, "' is not available for this function. Please specify a plot type supported by this function.", call. = FALSE)
  }
  
  removed_gene_list <- GATAI(ExpressionSet, ...)
  
  if (is.character(n_random_removal) && n_random_removal == "gatai") {
    n_random_genes <- length(removed_gene_list[[1]]) 
  } else if (is.numeric(n_random_removal)) {
    n_random_genes <- n_random_removal  
  } else {
    message("Invalid value for 'n_random_removal'. It should be 'gatai' or a numeric value.")
    message("Removing the same number as GATAI.")
    n_random_genes <- length(removed_gene_list[[1]])
  }
  
  random_removed_genes <- sample(ExpressionSet$GeneID, n_random_genes)

  if (plot_type == "separate"){
    #  original sample of genes
    P1 <-  myTAI::PlotSignature(ExpressionSet,
                                measure = measure,
                                TestStatistic = TestStatistic,
                                modules = modules,
                                permutations = permutations,
                                lillie.test = lillie.test,
                                p.value = p.value,
                                shaded.area = shaded.area,
                                custom.perm.matrix = custom.perm.matrix,
                                xlab = xlab,
                                ylab = ylab,
                                main = main,
                                lwd = lwd,
                                alpha = alpha,
                                y.ticks = y.ticks)
    
    # removed genes
    P2 <-  myTAI::PlotSignature(dplyr::filter(ExpressionSet, !GeneID %in% as.character(removed_gene_list[[1]])),
                                measure = measure,
                                TestStatistic = TestStatistic,
                                modules = modules,
                                permutations = permutations,
                                lillie.test = lillie.test,
                                p.value = p.value,
                                shaded.area = shaded.area,
                                custom.perm.matrix = custom.perm.matrix,
                                xlab = xlab,
                                ylab = ylab,
                                main = main,
                                lwd = lwd,
                                alpha = alpha,
                                y.ticks = y.ticks)
    
    P3 <-  myTAI::PlotSignature(dplyr::filter(ExpressionSet, !GeneID %in% random_removed_genes),
                                measure = measure,
                                TestStatistic = TestStatistic,
                                modules = modules,
                                permutations = permutations,
                                lillie.test = lillie.test,
                                p.value = p.value,
                                shaded.area = shaded.area,
                                custom.perm.matrix = custom.perm.matrix,
                                xlab = xlab,
                                ylab = ylab,
                                lwd = lwd,
                                alpha = alpha,
                                y.ticks = y.ticks)
    
    # save the data from both plots
    data_p1 <- ggplot2::ggplot_build(P1)$data[[1]]
    data_p2 <- ggplot2::ggplot_build(P2)$data[[1]]
    data_p3 <- ggplot2::ggplot_build(P3)$data[[1]]
    
    # deciding the y index range by looking at the minimum and maximum values of the line and shade from both plots
    min_y_p1 <- min(data_p1$y)
    max_y_p1 <- max(data_p1$y)
    min_y_p2 <- min(data_p2$y)
    max_y_p2 <- max(data_p2$y)
    min_y_p3 <- min(data_p3$y)
    max_y_p3 <- max(data_p3$y)
    min_y_p1_shading <- min(data_p1$ymin)
    max_y_p1_shading <- max(data_p1$ymax)
    min_y_p2_shading <- min(data_p2$ymin)
    max_y_p2_shading <- max(data_p2$ymax)
    min_y_p3_shading <- min(data_p3$ymin)
    max_y_p3_shading <- max(data_p3$ymax)
    min_y <- min(min_y_p1, min_y_p2, min_y_p3, min_y_p1_shading, min_y_p2_shading, min_y_p3_shading)
    max_y <- max(max_y_p1, max_y_p2, max_y_p3, max_y_p1_shading, max_y_p2_shading, max_y_p3_shading)
    
    P1 <- P1 + ggplot2::scale_y_continuous(limits = c(min_y-0.05, ymax = max_y+0.05), 
                                                      breaks = scales::breaks_pretty(n = y.ticks))

    P2 <- P2 + ggplot2::scale_y_continuous(limits = c(min_y-0.05, ymax = max_y+0.05), 
                                           breaks = scales::breaks_pretty(n = y.ticks))
    
    P3 <- P3 + ggplot2::scale_y_continuous(limits = c(min_y-0.05, ymax = max_y+0.05), 
                                           breaks = scales::breaks_pretty(n = y.ticks))
    
    
    plots <- cowplot::plot_grid(P1, P2, P3, labels = c(paste("Original sample of genes:", nrow(ExpressionSet), "genes"), 
                                                   paste("GATAI:", nrow(removed_gene_list), "Removed genes"), paste(n_random_genes, "Randomly Removed Genes")))
    print(plots)
  }
  
  if(plot_type == "combined"){
    expression_sets = list(ExpressionSet, dplyr::filter(ExpressionSet, !GeneID %in% as.character(removed_gene_list[[1]])),
                           dplyr::filter(ExpressionSet, !GeneID %in% random_removed_genes))

    set_labels = c(paste("Original sample of genes:", nrow(ExpressionSet), "genes"), 
                   paste("GATAI:", nrow(removed_gene_list), "Removed genes"), paste(n_random_genes, "Randomly Removed Genes"))
    
    print(myTAI::PlotSignatureMultiple(ExpressionSets = expression_sets,
                          set.labels = set_labels,
                          measure = measure,
                          TestStatistic=TestStatistic,
                          main = "Comparison of TAI Patterns",
                          y.tick = y.ticks))
  }
  
  return(removed_gene_list)
  
}
