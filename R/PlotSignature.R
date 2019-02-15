#' @title Plot evolutionary signatures across transcriptomes  
#' @description Main function to visualize transcriptome indices.
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
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
#' \item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a earlyconservation pattern (low-high-high pattern)
#' \item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
#' }
#' @param modules a list storing three elements for the \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}, or \code{\link{ReverseHourglassTest}}: early, mid, and late. 
#' Each element expects a numeric vector specifying the developmental stages 
#' or experiments that correspond to each module. For example:
#' \itemize{
#' \item \code{module} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} devides a dataset storing seven developmental stages into 3 modules.
#' } 
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}}, \code{\link{ReductiveHourglassTest}} or \code{\link{ReverseHourglassTest}}.
#' @param lillie.test a boolean value specifying whether the Lilliefors Kolmogorov-Smirnov Test shall be performed.
#' @param p.value a boolean value specifying whether the p-value of the test statistic shall be printed within the plot area.
#' @param shaded.area a boolean value specifying whether a shaded area shall 
#' be drawn for the developmental stages defined to be the presumptive phylotypic period.
#' @param custom.perm.matrix a custom \code{\link{bootMatrix}} (permutation matrix) to perform the underlying test statistic visualized by \code{PlotSignature}. Default is \code{custom.perm.matrix = NULL}.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param lwd line width.
#' @param alpha transparency of the shaded area (between [0,1]). Default is \code{alpha = 0.1}.
#' @param y.ticks number of ticks on the y-axis. Default is \code{ticks = 10}.
#' @author Hajk-Georg Drost
#' @details 
#' This function substitutes the functionality of the \code{\link{PlotPattern}} function
#' and is based on ggplot2 insead of base R graphics.
#' 
#' The following transcriptome indices can be computed and visualized with this function:
#' \itemize{
#' \item Transcriptome Age Index (\code{\link{TAI}})
#' \item Transcriptome Divergence Index (\code{\link{TDI}})
#' \item Transcriptome Polymorphism Index (\code{\link{TPI}})
#' }
#' 
#' @examples 
#' data(PhyloExpressionSetExample)
#' 
#' # plot TAI pattern and perform flat line test
#' PlotSignature(PhyloExpressionSetExample, 
#'               measure       = "TAI", 
#'               permutations  = 1000,
#'               TestStatistic = "FlatLineTest",
#'               ylab = "Transcriptome Age Index")
#'               
#' @export

PlotSignature <-
        function(ExpressionSet,
                 measure = "TAI",
                 TestStatistic = "FlatLineTest",
                 modules = NULL,
                 permutations = 1000,
                 lillie.test  = FALSE,
                 p.value = TRUE,
                 shaded.area  = FALSE,
                 custom.perm.matrix = NULL,
                 xlab = "Ontogeny",
                 ylab = "Transcriptome Index",
                 main = "",
                 lwd = 4,
                 alpha = 0.1,
                 y.ticks = 10) {
                
        
        if (!is.element(measure, c("TAI", "TDI", "TPI")))
                stop(
                        "Measure '",
                        measure,
                        "' is not available for this function. Please specify a measure supporting by this function.",
                        call. = FALSE
                )
        
        if (!is.element(TestStatistic, c("FlatLineTest", "ReductiveHourglassTest", "EarlyConservationTest", "ReverseHourglassTest")))
            stop("Please choose a 'TestStatistic' that is supported by this function. E.g. TestStatistic = 'FlatLineTest', TestStatistic = 'ReductiveHourglassTest', TestStatistic = 'EarlyConservationTest', TestStatistic = 'ReverseHourglassTest'.", call. = FALSE)
            
        cat("Plot signature: '",measure, "' and test statistic: '",TestStatistic,"' running ", permutations, " permutations." )
        cat("\n")
                
        Stage <- NULL        
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
        
        if (measure == "TPI") {
                TI <-
                        tibble::tibble(Stage = names(TPI(ExpressionSet)),
                                       TI = TPI(ExpressionSet))
        }
        
        # generate bootMatrix
        bm <- tibble::as_tibble(bootMatrix(ExpressionSet = ExpressionSet, 
                                           permutations  = permutations))
        
        
        if (TestStatistic == "FlatLineTest") {
                if (is.null(custom.perm.matrix)) {
                        resList <- FlatLineTest(ExpressionSet = ExpressionSet,
                                                permutations  = permutations)
                }
                
                else if (!is.null(custom.perm.matrix)) {
                        resList <- FlatLineTest(ExpressionSet      = ExpressionSet,
                                                custom.perm.matrix = custom.perm.matrix)
                        
                }
        }
        
        if (TestStatistic == "ReductiveHourglassTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReductiveHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReductiveHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReductiveHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReductiveHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
        }
        
        if (TestStatistic == "ReverseHourglassTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReverseHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReverseHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReverseHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReverseHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
        }
        
        if (TestStatistic == "EarlyConservationTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- EarlyConservationTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        EarlyConservationTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- EarlyConservationTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        EarlyConservationTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                        
                }
        }
        
        # get p-value and standard deviation values from the test statistic
        pval <- resList$p.value
        pval <- format(pval,digits = 3)
        
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
        
        if ((TestStatistic == "FlatLineTest") && p.value){
                TI.ggplot <-
                        TI.ggplot + ggplot2::annotate(
                                "text",
                                x = 2,
                                y = max(TI$TI) + (max(TI$TI) / 30),
                                label = paste0("p_flt = ", pval),
                                size = 6
                        )
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
                cat("\n")
        }
                
        if (TestStatistic == "ReductiveHourglassTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::annotate(
                                        "text",
                                        x = 2,
                                        y = max(TI$TI) + (max(TI$TI) / 30),
                                        label = paste0("p_rht = ", pval),
                                        size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }   
        
        if (TestStatistic == "ReverseHourglassTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::annotate(
                                        "text",
                                        x = 2,
                                        y = max(TI$TI) + (max(TI$TI) / 30),
                                        label = paste0("p_reverse_hourglass = ", pval),
                                        size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }  
        
        if (TestStatistic == "EarlyConservationTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::annotate(
                                        "text",
                                        x = 2,
                                        y = max(TI$TI) + (max(TI$TI) / 30),
                                        label = paste0("p_ect = ", pval),
                                        size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha * 0.5)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }
        
        TI.ggplot <- TI.ggplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
        return (TI.ggplot)
}


