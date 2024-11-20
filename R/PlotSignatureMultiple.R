
#' @title Plot evolutionary signatures across transcriptomes, for multiple expression sets
#' @description \emph{PlotSignatureMultiple} is used to compare multiple 
#' transcriptomic index profiles (e.g. TAI) over a common process of interest
#' (e.g. a developmental process). 
#' The main use case is visualising how removing different subsets of genes from
#' the expression set can perturb the transcriptome index signal.
#'
#' @param ExpressionSets a list of PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet objects.
#' @param set.labels a character vector of labels, one for each given expression set
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
#' @return a ggplot object visualising the transcriptome index of each given
#' expression set, together with its standard deviation per stage,
#' obtained by permuting the gene ages.
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
#' # remove top 1% expressed genes
#' genes.top_expression <- TopExpressionGenes(PhyloExpressionSetExample, p=.99)
#' PhyloExpressionSetExample.top_removed <- subset(PhyloExpressionSetExample, 
#'                                         !(GeneID %in% genes.top_expression))
#' 
#' expression_sets = list(PhyloExpressionSetExample, PhyloExpressionSetExample.top_removed)
#' set_labels = c("100%", "99%")
#' 
#' # Flat line test
#' PlotSignatureMultiple(ExpressionSets = expression_sets,
#'                       set.labels = set_labels,
#'                       TestStatistic="FlatLineTest",
#'                       main = "A. thaliana embryogenesis",
#'                       legend.title = "Quantile")
#'                       
#' # Reductive hourglass test
#' PlotSignatureMultiple(ExpressionSets = expression_sets,
#'                       set.labels = set_labels,
#'                       TestStatistic="ReductiveHourglassTest",
#'                       main = "A. thaliana embryogenesis",
#'                       legend.title = "Quantile",
#'                       modules=list(early=1:2, mid=3:5, late=6:7),
#'                       shaded.area=TRUE)                                                        
PlotSignatureMultiple <-
    function(ExpressionSets,
             set.labels,
             measure = "TAI",
             TestStatistic = "FlatLineTest",
             modules = NULL,
             permutations = 1000,
             p.value = TRUE,
             shaded.area = FALSE,
             xlab = "Ontogeny",
             ylab = "Transcriptome Index",
             main = "",
             legend.title = "Expression Sets",
             lwd = 4,
             alpha = 0.1,
             y.ticks = 3) {
        
        ## Parameter validation
        TXI <- switch(
            measure,
            "TAI" = TAI,
            "TDI" = TDI,
            "TPI" = TPI,
            stop("Measure '",measure, "' is not available for this function. Please specify a measure supporting by this function.",
                call. = FALSE
            )
        )
        Test <- switch(
            TestStatistic,
            "FlatLineTest" = FlatLineTest,
            "ReductiveHourglassTest" = ReductiveHourglassTest,
            "ReverseHourglassTest" = ReverseHourglassTest,
            "EarlyConservationTest" = EarlyConservationTest,
            "LateConservationTest" = LateConservationTest,
            stop("Please select the available test: 'FlatLineTest', 'ReductiveHourglassTest', 'ReverseHourglassTest', 'EarlyConservationTest' or 'LateConservationTest' using the argument test = 'FlatLineTest'",
                call. = FALSE
            )
        )
        
        if (TestStatistic %in% c(
            "ReductiveHourglassTest",
            "ReverseHourglassTest",
            "EarlyConservationTest",
            "LateConservationTest"
        ) & is.null(modules))
            stop(
                "Please specify the three modules: early, mid, and late using the argument 'module = list(early = ..., mid = ..., late = ...)'.",
                call. = FALSE
            )
        
        
        if (length(ExpressionSets) != length(set.labels))
            stop("Please specify an equal number of set labels and expression sets.", call. = FALSE)
        
        # The expression sets should have the same domains
        same_domain <- ExpressionSets |>
                       lapply(function(x) names(x[, -c(1:2)])) |>
                       unique() |>
                       (\(x) length(x) == 1)()
        
        if (!same_domain) 
            stop("Please ensure the expression sets have the same domain of expression stages.")
        
        
        
        # [TODO] Usually the expression sets will share a common superset of genes
        # Would it be worth it then to compute only one shared bootstrap matrix?
        if (TestStatistic == "FlatLineTest")
            test_results <- lapply(ExpressionSets, \(x) Test(x, permutations=permutations))
        else
            test_results <- lapply(ExpressionSets, \(x) Test(x, permutations=permutations, modules=modules))
        
        TIs <- lapply(ExpressionSets, TXI)
        
        plot <- ggplot2::ggplot()
        p_label <- switch(
            TestStatistic,
            "FlatLineTest" = "p_flt",
            "ReductiveHourglassTest" = "p_rht",
            "ReverseHourglassTest" = "p_reverse_hourglass",
            "EarlyConservationTest" = "p_ect",
            "LateConservationTest" = "p_lct",
            "p_test"
        )
        labels <- vector()
        # Plot TXI curves and errors in a shared domain
        for (i in 1:length(ExpressionSets)) {
            
            p_value <- test_results[[i]]$p.value
            std_devs <- test_results[[i]]$std.dev
            
            label <- paste0("<br><span style='font-size:16pt'>**",
                           set.labels[[i]],
                           "**</span>"
            )
            if (p.value)
                label <- paste0(label,
                               "<br><br>**",
                               p_label,
                               "**: ",
                               as.character(signif(p_value,digits=3)),
                               "<br>")
            
            labels <- c(labels, label)
            
            TI <- tibble::tibble(Stage = names(TIs[[i]]), 
                                 TI = TIs[[i]], 
                                 Group = label
                                 )
            
            
            plot <- plot +
                    ggplot2::geom_line(
                        data=TI, 
                        ggplot2::aes(
                            x=factor(Stage, levels=unique(Stage)), 
                            y=TI, 
                            group = Group,
                            color = Group
                        ), 
                        lwd = lwd,
                    ) +
                    ggplot2::geom_ribbon(
                        data=TI, 
                        ggplot2::aes(
                            x = factor(Stage, levels=unique(Stage)),
                            ymin = TI - std_devs, 
                            ymax = TI + std_devs,
                            group=Group,
                            fill=Group
                        ), 
                        alpha=alpha
                    )
        }
        
        # Plot shaded area
        if (shaded.area) 
            plot <- plot + ggplot2::geom_rect(ggplot2::aes(
                xmin = modules[[2]][1] - 0.3,
                xmax = modules[[2]][length(modules[[2]])] + 0.3,
                ymin = -Inf,
                ymax = Inf), fill = "#4d004b", alpha = alpha)  
        
            
        # Labels and style
        plot <- plot + 
            ggplot2::scale_y_continuous(
                breaks = scales::pretty_breaks(n = y.ticks)
            ) +
            ggplot2::scale_colour_discrete(breaks=labels) + 
            ggplot2::scale_fill_discrete(breaks=labels) + 
            ggplot2::labs(
                x = xlab, 
                y = ylab,
                title=main,
                colour = legend.title,
                fill = legend.title
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                title            = ggplot2::element_text(face = "bold", size=18),
                legend.title     = ggplot2::element_text(face = "bold", size=18),
                legend.text      = ggtext::element_markdown(size=14),
                axis.title       = ggplot2::element_text(face = "bold", size=20),
                axis.text.y      = ggplot2::element_text(face = "bold", size=18),
                axis.text.x      = ggplot2::element_text(face = "bold", size=18),
                panel.background = ggplot2::element_blank(),
                strip.text.x     = ggplot2::element_text(
                    colour         = "black",
                    face           = "bold",
                    size           = 18
                )
            )
                
        return(plot)
    }