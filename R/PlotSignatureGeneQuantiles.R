
#' @title Plot evolutionary signatures across transcriptomes, removing top level genes
#' @description \emph{PlotSignatureGeneQuantiles} is used to investigate
#' the robustness of a transcriptomic index signal to removing the most highly expressed genes
#' from the given expression set. The resulting perturbed signal is plotted for 
#' different quantile probability thresholds defining the set of top level genes.
#'
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
#' @param quantiles a numeric vector of quantile probabilities (between [0,1])
#' @param gene.selection.criterion a string defining the criterion by which genes should be selected
#' Possible values can be:
#' \itemize{
#' \item \code{gene.selection.criterion} = \code{"mean"} : Sort the genes by the mean of their gene expression (equivalently, by their total expression counts across stages)
#' \item \code{gene.selection.criterion} = \code{"variance"} : Sort the genes by the variance of their gene expression across stages
#' }
#' @param measure type of transcriptome index that shall be computed. E.g.
#' \itemize{
#' \item \code{measure = "TAI"} (Transcriptome Age Index)
#' \item \code{measure = "TDI"} (Transcriptome Divergence Index)
#' \item \code{measure = "TPI"} (Transcriptome Polymorphism Index)
#' }
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be:
#' \itemize{
#' \item \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line
#' \item \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern)
#' \item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
#' \item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a early conservation pattern (low-high-high pattern)
#' \item \code{TestStatistic} = \code{"LateConservationTest"} : Statistical test for the existence of a late conservation pattern (high-high-low pattern)
#' }
#' @param modules a list storing three elements for the \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}},
#' or \code{\link{ReverseHourglassTest}}: early, mid, and late.
#' Each element expects a numeric vector specifying the developmental stages
#' or experiments that correspond to each module. For example:
#' \itemize{
#' \item \code{modules} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} divides a dataset storing seven developmental stages into 3 modules.
#' }
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}}, \code{\link{ReductiveHourglassTest}} or \code{\link{ReverseHourglassTest}}.
#' @param p.value a boolean value specifying whether the p-value of the test statistic shall be printed within the legend, for each expression set.
#' @param shaded.area a boolean value specifying whether a shaded area shall
#' be drawn for the developmental stages defined to be the presumptive phylotypic period.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param lwd line width.
#' @param alpha transparency of the shaded area and error ribbon (between [0,1]). Default is \code{alpha = 0.1}.
#' @param y.ticks number of ticks on the y-axis. Default is \code{ticks = 3}.
#'
#' @return a ggplot object visualising the transcriptome index of the
#' expression set, together with its standard deviation per stage,
#' obtained by permuting the gene ages.
#' For each quantile probability threshold, the resulting perturbed signal is shown
#' as a separate profile.
#' The profiles are shown on the same axes, so that they can be readily compared.
#' Optionally, the p-value of each profile, with respect to the choice of statistic,
#' is shown.
#' 
#' @author Stefan Manolache
#' 
#' @export
#'
#' @examples
#' data(PhyloExpressionSetExample)
#'                       
#' # Flat line test, select top expressed genes
#' PlotSignatureGeneQuantiles(ExpressionSet = PhyloExpressionSetExample,
#'                            quantiles=c(1.0, 0.99, 0.95, 0.90, 0.80),
#'                            main="Excluding top level genes by total 
#'                            expression using different thresholds",
#'                            gene.selection.criterion="mean")
#' # Flat line test, select top genes by variance of expression
#' PlotSignatureGeneQuantiles(ExpressionSet = PhyloExpressionSetExample,
#'                            main="Excluding top level genes by variance of 
#'                            expression using different thresholds",
#'                            quantiles=c(1.0, 0.99, 0.95, 0.90, 0.80),
#'                            gene.selection.criterion="variance")


PlotSignatureGeneQuantiles <-
    function(ExpressionSet,
             quantiles=c(1.0, 0.99, 0.95, 0.90, 0.80),
             gene.selection.criterion="mean",
             measure = "TAI",
             TestStatistic = "FlatLineTest",
             modules = NULL,
             permutations = 1000,
             p.value = TRUE,
             shaded.area = FALSE,
             xlab = "Ontogeny",
             ylab = "Transcriptome Index",
             main = "",
             lwd = 4,
             alpha = 0.1,
             y.ticks = 3) {
        
        TopGenes <- switch(gene.selection.criterion,
                           "mean" = TopExpressionGenes,
                           "variance" = TopVarianceGenes,
                           stop("Gene selection criterion '",gene.selection.criterion, "' is not available for this function. Please specify a selection criterion supported by this function.",
                                call. = FALSE)
                           )
        
        expression_sets <- quantiles |>
                           lapply(\(p) TopGenes(ExpressionSet, p=p)) |>
                           lapply(\(genes) subset(ExpressionSet, !(ExpressionSet[[2]] %in% genes)))
        
        set_labels <- quantiles
        
        plot <- PlotSignatureMultiple(expression_sets,
                                      set_labels,
                                      measure = measure,
                                      TestStatistic = TestStatistic,
                                      modules = modules,
                                      permutations = permutations,
                                      p.value = p.value,
                                      shaded.area = shaded.area,
                                      xlab = xlab,
                                      ylab = ylab,
                                      main = main,
                                      legend.title = "Quantile",
                                      lwd = lwd,
                                      alpha = alpha,
                                      y.ticks = y.ticks)
        
        return(plot)
    }