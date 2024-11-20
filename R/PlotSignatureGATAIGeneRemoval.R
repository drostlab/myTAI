#' @title Plot evolutionary signatures across transcriptomes after GATAI gene removal
#' @description Main function to visualize transcriptome indices after GATAI gene removal  .
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
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
#' @param \dots parameters passed on to \code{\link{GATAI}}.
#' @author Filipa Martins Costa and Hajk-Georg Drost
#' @export
PlotSignatureGATAIGeneRemoval <- function(ExpressionSet,
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
                                          y.ticks = 10, ...) {
  
  removed_gene_list <- GATAI(ExpressionSet, ...)
 
  p_final <- myTAI::PlotSignature(dplyr::filter(ExpressionSet, !ExpressionSet[[2]] %in% as.character(removed_gene_list[[1]])),
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
                      main = paste(nrow(removed_gene_list), "Removed genes"),
                      lwd = lwd,
                      alpha = alpha,
                      y.ticks = y.ticks)
  
  # plot_final <- cowplot::plot_grid(p1, p2, labels = c('A', 'B'), label_size = 18)
  
  print(p_final)
  
  return(p_final)
}
