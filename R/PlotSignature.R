#' @title Plot evolutionary signatures across transcriptomes  
#' @description E
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be: \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line.
#' \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern).
#' \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a earlyconservation pattern (low-high-high pattern).
#' @param modules a list storing three elements for the \code{\link{ReductiveHourglassTest}} or \code{\link{EarlyConservationTest}}: early, mid, and late. 
#' Each element expects a numeric vector specifying the developmental stages 
#' or experiments that correspond to each module. For example, 
#' \code{module} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} devides a dataset storing seven developmental stages into 3 modules.
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}} or \code{\link{ReductiveHourglassTest}}.
#' @param xlab
#' @param ylab
#' @param main
#' @param lwd
#' @param alpha
#' @param y.ticks
#' @author Hajk-Georg Drost
#' @export

PlotSignature <-
        function(ExpressionSet,
                 measure = "TAI",
                 TestStatistic = "FlatLineTest",
                 modules = NULL,
                 permutations = 1000,
                 xlab = "Ontogeny",
                 ylab = "Transcriptome Index",
                 main = "",
                 lwd = 4,
                 alpha = 0.1,
                 y.ticks = 10) {
                
        
        if (!is.element(measure, c("TAI", "TDI", "TPI")))
                stop("Measure '",measure,"' is not available for this function. Please specify a measure supporting by this function.", call.= FALSE)
        
        # store transcriptome index in tibble
        if (measure == "TAI") {
                TI <-
                        tibble::tibble(Stage = names(TAI(ExpressionSet)),
                                       TI = TAI(ExpressionSet))
        }
        
        if (measure == "TDI") {
                TI <-
                        tibble::tibble(Stage = names(TDI(ExpressionSet)),
                                       TI = TDI(ExpressionSet))
        }
        
        # generate bootMatrix
        bm <- tibble::as_tibble(bootMatrix(ExpressionSet = ExpressionSet, 
                                           permutations  = permutations))
        
        TI.ggplot <- ggplot2::ggplot(TI, ggplot2::aes(
                x = factor(Stage, levels = unique(Stage)),
                y = TI,
                group = 1
        ))  + ggplot2::geom_ribbon(ggplot2::aes(
                ymin = TI - apply(bm, 2, stats::sd),
                ymax = TI + apply(bm, 2, stats::sd)
        ), alpha = alpha) +
                ggplot2::geom_line(lwd = lwd) +
                ggplot2::theme_minimal() +
                ggplot2::labs(x = xlab, y = ylab, title = main) +
                ggplot2::theme(
                        title            = ggplot2::element_text(size = 18, face = "bold"),
                        legend.title     = ggplot2::element_text(size = 18, face = "bold"),
                        legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                        axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                        axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                        axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
                        panel.background = ggplot2::element_blank(),
                        strip.text.x     = ggplot2::element_text(
                                size           = 18,
                                colour         = "black",
                                face           = "bold"
                        )
                ) +
                ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))
        
        
        return (TI.ggplot)
}


