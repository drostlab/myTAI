#' @title Plot the Transcriptome Age Index or Transcriptome Divergence Index
#' @description Function to plot the TAI or TDI of a given PhyloExpressionSet or DivergenceExpressionSet object.
#' This function plot the \code{\link{TAI}} or \code{\link{TDI}} of a given PhyloExpressionSet or DivergenceExpressionSet object. 
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
#' @param lillie.test a boolean value specifying whether the Lilliefors Kolmogorov-Smirnov Test shall be performed.
#' @param digits.ylab a numeric value specifying the number of digits shown for the TAI or TDI values on the y-axis.
#' @param p.value a boolean value specifying whether the p-value of the test statistic shall be printed within the plot area.
#' @param shaded.area a boolean value specifying whether a shaded area shall 
#' be drawn for the developmental stages defined to be the presumptive phylotypic period.
#' @param y.ticks a numeric value specifying the number of ticks to be drawn on the y-axis.
#' @param custom.perm.matrix a custom \code{\link{bootMatrix}} (permutation matrix) to perform the underlying test statistic visualized by \code{PlotPattern}. Default is \code{custom.perm.matrix = NULL}.
#' @param \dots default plot parameters.
#' @details 
#' 
#' This function computes a permutation test quantifying the statistical significance of the prensent phylotranscriptomics pattern. 
#' The user can choose between the \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}, or \code{\link{EarlyConservationTest}}. 
#' The \code{\link{FlatLineTest}} tests for any significant deviation from a flat line. 
#' Each period or stage that significantly deviates from a flat line, might be governed by stronger selective pressure (in terms of natural selection) compared to other stages or periods of development.
#' The \code{\link{ReductiveHourglassTest}} specificly tests for the statistical significance of a molecular hourglass pattern (high-low-high pattern) with prior biological knowlegde.
#' The corresponding p-value that is printed within the plot (by default) specifies the statistical significance of the chosen test statistic.
#' 
#' The \code{\link{EarlyConservationTest}} specificly tests for the statistical significance of a low-high-high pattern (monotonically increasing function)
#' with prior biological knowlegde.
#' The corresponding p-value that is printed within the plot (by default) specifies the statistical significance of the chosen test statistic.
#' 
#' The x-axis denotes the developmental series (time course / experiments / ontogeny) of the input ExpressionSet and the y-axis 
#' denotes the corresponding mean transcriptome age value (\code{\link{TAI}} or \code{\link{TDI}}) of the given ExpressionSet. 
#' 
#' Furthermore, the grey lines above and below the actual phylotranscriptomics pattern denotes the standard deviations of \code{\link{TAI}} or \code{\link{TDI}} 
#' values that have been estimated from the \code{\link{bootMatrix}}.
#' A low mean transcriptome age value denotes an evolutionary older transcriptome being active during the corresponding periods, 
#' whereas a high mean transcriptome age value denotes an evolutionary younger transcriptome being active during the corresponding periods.
#' For mean transcriptome divergence, a low mean transcriptome divergence value denotes a more conserved transcriptome 
#' being active (between two species), whereas a high mean transcriptome divergence value denotes a more divergent transcriptome 
#' being active (between two species) - in terms of protein-sequence substitution rates.
#' 
#' This function is useful to fastly plot the \code{\link{TAI}} or \code{\link{TDI}} 
#' profile of a given PhyloExpressionSet or DivergenceExpressionSet object and 
#' the statistical significance of the corresponding pattern.
#' Internally the function calls several graphics functions, 
#' such as \code{\link{plot}}, \code{\link{axis}}, and \code{\link{legend}}. 
#' For the ellipsis argument \code{\link{...}} all graphics specific arguments can be defined. 
#' Internally the function specific arguments for e.g. \code{\link{plot}}, \code{\link{axis}}, 
#' and \code{\link{legend}} will be detected and are passed to the corresponding function.
#' 
#' Hence, when calling the function \code{PlotPattern}, one can specify arguments 
#' for \code{\link{plot}} and \code{\link{axis}} and \code{\link{legend}} 
#' as \code{\link{...}}. 
#' 
#' In case prior biological knowledge is present for a specific period of development,
#' the \code{shaded.area} argument can be set to \code{TRUE} and the function will use
#' the values stored in the \code{mid} argument to draw a shaded area within the corresponding period of development.
#' @return a plot visualizing the phylotranscriptomic pattern of a given PhyloExpressionSet or DivergenceExpressionSet object.
#' 
#' The corresponding \emph{p-value} of the test statistic is named as follows:
#' 
#' \code{p_flt} = p-value of the corresponding \code{\link{FlatLineTest}}
#' 
#' \code{p_rht} = p-value of the corresponding \code{\link{ReductiveHourglassTest}}
#' 
#' \code{p_ect} = p-value of the corresponding \code{\link{EarlyConservationTest}}
#' 
#' @references 
#' Domazet-Loso T and Tautz D. (2010). \emph{A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012). \emph{A transcriptomic hourglass in plant embryogenesis}. Nature (490): 98-101.
#' 
#' Drost HG et al. (2015) \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{FlatLineTest}},
#'  \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}
#' @examples
#' 
#' # load PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # only visualize the TAI profile without any test statistics...
#' # this is equavalent to performing: plot(TAI(PhyloExpressionSetExample), type = "l", lwd = 6)
#' PlotPattern(ExpressionSet = PhyloExpressionSetExample,
#'             TestStatistic = NULL,
#'             type          = "l",
#'             xlab          = "Ontogeny",
#'             ylab          = "TAI",
#'             lwd           = 9)
#'
#' # the simplest example of plotting the TAI profile of a given PhyloExpressionSet:
#' # In this case (default) the FlatLineTest will be performed to quantify
#' # the statistical significance of the present TAI pattern and will be drawn as 'p = ... '
#' # in the plot
#'
#' PlotPattern(ExpressionSet = PhyloExpressionSetExample, 
#'             TestStatistic = "FlatLineTest",
#'             permutations  = 100, 
#'             type          = "l", 
#'             xlab          = "Ontogeny", 
#'             ylab          = "TAI", 
#'             lwd           = 9)
#' 
#' # an example performing the ReductiveHourglassTest and printing the p-value
#' # and shaded area of the presumptive phylotypic period into the plot
#' # Here the 'p = ...' denotes the p-value that is returned by the ReductiveHourglassTest
#'
#' PlotPattern(
#'             ExpressionSet = PhyloExpressionSetExample,
#'             TestStatistic = "ReductiveHourglassTest",
#'             modules       = list(early = 1:2,mid = 3:5,late = 6:7), 
#'             permutations  = 100, 
#'             p.value       = TRUE, 
#'             shaded.area   = TRUE, 
#'             xlab          = "Ontogeny", 
#'             ylab          = "TAI", 
#'             type          = "l", 
#'             lwd           = 9)
#'
#' # testing for early conservation model 
#' PlotPattern( ExpressionSet = PhyloExpressionSetExample,
#'              TestStatistic = "EarlyConservationTest",
#'             modules        = list(early = 1:2,mid = 3:5,late = 6:7), 
#'             permutations   = 100,
#'             p.value        = TRUE, 
#'             shaded.area    = TRUE, 
#'             xlab           = "Ontogeny", 
#'             ylab           = "TAI", 
#'             type           = "l", 
#'             lwd            = 9)
#'             
#' 
#' # use your own permutation matrix
#' custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)
#' 
#' PlotPattern(ExpressionSet      = PhyloExpressionSetExample, 
#'             TestStatistic      = "FlatLineTest",
#'             custom.perm.matrix = custom_perm_matrix, 
#'             type               = "l", 
#'             xlab               = "Ontogeny", 
#'             ylab               = "TAI", 
#'             lwd                = 9)
#' 
#' 
#' @export

PlotPattern <- function(ExpressionSet, 
                        TestStatistic      = "FlatLineTest",
                        modules            = NULL, 
                        permutations       = 1000,
                        lillie.test        = FALSE,
                        digits.ylab        = 4, 
                        p.value            = TRUE, 
                        shaded.area        = FALSE, 
                        y.ticks            = 4,
                        custom.perm.matrix = NULL, ...)
{
        
        is.ExpressionSet(ExpressionSet)
        
        if (is.null(TestStatistic)){
                
                graphics::plot(TAI(ExpressionSet), xaxt = "n", ...)
                graphics::axis(1,1:(ncol(ExpressionSet)-2),names(ExpressionSet)[3:ncol(ExpressionSet)])
        } else {
        
        if (!(is.element(TestStatistic, c("FlatLineTest","ReductiveHourglassTest","EarlyConservationTest")))){
                stop("Please enter a correct string for the test statistic: 'FlatLineTest', 'EarlyConservationTest' or 'ReductiveHourglassTest'.")
        }
        
        if ((is.element(TestStatistic,c("ReductiveHourglassTest","EarlyConservationTest"))) & is.null(modules))
                stop ("Please specify the modules for the ReductiveHourglassTest or EarlyConservationTest: modules = list(early = ..., mid = ..., late = ...).")
        
        if ((!is.null(modules)) & (TestStatistic == "FlatLineTest"))
                warning ("You don't need to specify the modules argument for the FlatLineTest.")
        
        nCols <- dim(ExpressionSet)[2]
        resList <- vector("list", length = 2)
        age <- vector(mode = "numeric",length = (nCols-2))
        age <- cpp_TAI(as.matrix(ExpressionSet[ , 3:nCols]),as.vector(ExpressionSet[ , 1]))
        

        if (TestStatistic == "FlatLineTest"){
                
                if (is.null(custom.perm.matrix)){
                        
                        resList <- FlatLineTest( ExpressionSet = ExpressionSet,
                                                 permutations  = permutations )
                }
                
                else if (!is.null(custom.perm.matrix)){
                        
                        resList <- FlatLineTest( ExpressionSet      = ExpressionSet,
                                                 custom.perm.matrix = custom.perm.matrix)
                        
                }
        }
        
        if (TestStatistic == "ReductiveHourglassTest"){
                
                if (lillie.test){
                        
                        if (is.null(custom.perm.matrix)){
                                resList <- ReductiveHourglassTest( ExpressionSet = ExpressionSet,
                                                                   modules       = modules, 
                                                                   permutations  = permutations, 
                                                                   lillie.test   = TRUE )
                        }
                        
                        else if (!is.null(custom.perm.matrix)){
                                resList <- ReductiveHourglassTest( ExpressionSet      = ExpressionSet,
                                                                   modules            = modules, 
                                                                   lillie.test        = TRUE,
                                                                   custom.perm.matrix = custom.perm.matrix)
                        }
                }
                
                if(!lillie.test){
                        
                        if (is.null(custom.perm.matrix)){
                                resList <- ReductiveHourglassTest( ExpressionSet = ExpressionSet,
                                                                   modules       = modules, 
                                                                   permutations  = permutations, 
                                                                   lillie.test   = FALSE )
                        }
                        
                        else if (!is.null(custom.perm.matrix)){
                                resList <- ReductiveHourglassTest( ExpressionSet      = ExpressionSet,
                                                                   modules            = modules, 
                                                                   lillie.test        = FALSE,
                                                                   custom.perm.matrix = custom.perm.matrix)
                        }
                }
        }
        
        if(TestStatistic == "EarlyConservationTest"){
                
                if(lillie.test){
                        
                        if (is.null(custom.perm.matrix)){
                                resList <- EarlyConservationTest( ExpressionSet = ExpressionSet,
                                                                  modules       = modules, 
                                                                  permutations  = permutations, 
                                                                  lillie.test   = TRUE )
                        }
                        
                        else if (!is.null(custom.perm.matrix)){
                                resList <- EarlyConservationTest( ExpressionSet      = ExpressionSet,
                                                                  modules            = modules, 
                                                                  lillie.test        = TRUE,
                                                                  custom.perm.matrix = custom.perm.matrix)
                        }
                }
                
                if(!lillie.test){
                        
                        if (is.null(custom.perm.matrix)){
                                resList <- EarlyConservationTest( ExpressionSet = ExpressionSet,
                                                                  modules       = modules, 
                                                                  permutations  = permutations, 
                                                                  lillie.test   = FALSE )
                        }
                        
                        else if (!is.null(custom.perm.matrix)){
                                resList <- EarlyConservationTest( ExpressionSet      = ExpressionSet,
                                                                  modules            = modules, 
                                                                  lillie.test        = FALSE,
                                                                  custom.perm.matrix = custom.perm.matrix)
                        }
                        
                }
        }
        
        # get p-value and standard deviation values from the test statistic 
        pval <- resList$p.value
        sd_vals <- resList$std.dev
        #random_mean <- resList$mean
        max_sd <- max(sd_vals)
        
        # plot the age pattern surrounded by the standard deviation
        
        # computing the standard error of the TAI/TDI pattern using bootstrap analyses
        x <- min(age)
        y <- max(age)
        
        # adjust the plotting margins 
        aspect_ratio <- ((max_sd + y) - (x - max_sd)) / 20
        ylim.range <- range((x - max_sd - aspect_ratio),(y + max_sd + aspect_ratio))  
        
        # define arguments for different graphics functions
        plot.args <- c("type","lwd","col","cex.lab","main","xlab","ylab")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        #
        #   plot phylotranscriptomic age
        # 
        
        if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) | (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                
                do.call(graphics::plot,c(list(x = age,ylim = c(ylim.range[1],ylim.range[2]),axes = FALSE), 
                                         dots[!is.element(names(dots),c(axis.args,legend.args))]))
        }
        
        # default: xlab = "Ontogeny" and ylab = "Age Index"
        else{
                do.call(graphics::plot,c(list(x = age,ylim = c(ylim.range[1],ylim.range[2]),axes = FALSE, 
                                              xlab = "Ontogeny",ylab = "Age Index" ), dots[!is.element(names(dots),c(axis.args,legend.args))]))  
        }
        
        do.call(graphics::axis,c(list(side = 1,at = 1:(nCols-2),
                                      labels = names(ExpressionSet)[3:nCols]),dots[!is.element(names(dots),c(plot.args,legend.args))]))
        
        do.call(graphics::axis,c(list(side = 2,at = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab),
                                      labels = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab)), 
                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
        
        
        # age + std.err
        graphics::lines(age + sd_vals,lwd = 2,col = "darkgrey")
        # age - std.err
        graphics::lines(age - sd_vals,lwd = 2,col = "darkgrey")
        
        if(p.value == TRUE){
                
                if(TestStatistic == "FlatLineTest"){
                        do.call(graphics::legend,c(list(x = "topleft",bty = "n",legend = paste("p_flt = ",format(pval,digits = 3),sep = ""),
                                           dots[!is.element(names(dots),c(plot.args,axis.args))])))
                }
                
                if(TestStatistic == "ReductiveHourglassTest"){
                        do.call(graphics::legend,c(list(x = "topleft",bty = "n",legend = paste("p_rht = ",format(pval,digits = 3),sep = ""),
                                                   dots[!is.element(names(dots),c(plot.args,axis.args))])))
                }
                
                if(TestStatistic == "EarlyConservationTest"){
                        do.call(graphics::legend,c(list(x = "topleft",bty = "n",legend = paste("p_ect = ",format(pval,digits = 3),sep = ""),
                                                   dots[!is.element(names(dots),c(plot.args,axis.args))])))
                }
                
        }
        
        # draw a dotted line passing the global minimum value of age
        #abline(h=min(age),col="black",lty="dotted");
        
        if(shaded.area == TRUE){
                #abline(v=c(mid[1],mid[length(mid)]),col="black",lty="dotted")
                col2alpha <- function (color, alpha = 0.2){
                        rgbCode <- grDevices::col2rgb(color)[,1]
                        grDevices::rgb(rgbCode[1], rgbCode[2], rgbCode[3], 255 * alpha, maxColorValue = 255)
                }
                usr <- graphics::par('usr')
                graphics::rect(modules[[2]][1], usr[3], modules[[2]][length(modules[[2]])], usr[4], col = col2alpha("midnightblue",alpha = 0.2)) 
        }
        }
}














