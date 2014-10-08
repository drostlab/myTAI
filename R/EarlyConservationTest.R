#' @title Function to perform the Early Conservation Test.
#' @description The \emph{Early Conservation Test} has been developed to statistically evaluate the
#' existence of a monotonically increasing phylotranscriptomic pattern based on \code{\link{TAI}} or \code{\link{TDI}} computations.
#' The corresponding p-value quantifies the probability that a given TAI or TDI pattern (or any phylotranscriptomics pattern) 
#' does not follow an early conservation like pattern. A p-value < 0.05 indicates that the corresponding phylotranscriptomics pattern does
#' indeed follow an early conservation (low-high-high) shape.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param modules a list storing three elements: early, mid, and late. Each element expects a numeric
#' vector specifying the developmental stages or experiments that correspond to each module. 
#' For example, \code{module} = list(early = 1:2, mid = 3:5, late = 6:7) devides a dataset 
#' storing seven developmental stages into 3 modules.
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{ReductiveHourglassTest}.
#' @param lillie.test a boolean value specifying whether the Lilliefors Kolmogorov-Smirnov Test shall be performed to quantify the goodness of fit.
#' @param plotHistogram a boolean value specifying whether a \emph{Lillifor's Kolmogorov-Smirnov-Test} 
#' shall be performed to test the goodness of fit of the approximated distribution, as well as additional plots quantifying the significance
#' of the observed phylotranscriptomic pattern.
#' @param parallel a boolean value specifying whether goodness of fit computations 
#' shall be performed in parallel using the \pkg{doMC} package.
#' @param runs specify the number of runs to be performed for goodness of fit computations, in case \code{plotHistogram} = \code{TRUE}.
#' In most cases \code{runs} = 100 is a reasonable choice. Default is \code{runs} = 10 (because it takes less computation time for demonstration purposes).
#' @details The \emph{Early Conservation Test}
#' @return a list object containing the list elements:
#' 
#' p.value : the p-value quantifying the statistical significance (low-high-high pattern) of the given phylotranscriptomics pattern.
#'
#' std.dev : the standard deviation of the N sampled phylotranscriptomics patterns for each developmental stage S.
#' 
#' lillie.test : a boolean value specifying whether the \emph{Lillifors KS-Test} returned a p-value > 0.05, 
#' which indicates that fitting the permuted scores with a normal distribution seems plausible.
#' @references 
#' 
#' Drost et al. (2014), Active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis.
#'
#' Quint M et al. (2012). "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' 
#' Piasecka B, Lichocki P, Moretti S, et al. (2013) The hourglass and the early conservation models--co-existing patterns of developmental constraints in vertebrates. PLoS Genet. 9(4): e1003476.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{gpScore}}, \code{\link{bootMatrix}}, \code{\link{FlatLineTest}}, \code{\link{PlotPattern}}
#' @examples \dontrun{
#' 
#' 
#' 
#' }
#' @export

EarlyConservationTest <- function(ExpressionSet,modules = NULL,
                                  permutations = 1000, lillie.test = FALSE, 
                                  plotHistogram = FALSE,parallel = FALSE,
                                  runs = 10){
        
        is.ExpressionSet(ExpressionSet)
        
        if(length(modules) != 3)
                stop("Please specify three modules: early, mid, and late to perform the ReductiveHourglassTest.")
        
        if(length(unlist(modules)) != (dim(ExpressionSet)[2] - 2))
                stop("The number of stages classified into the three modules does not match the total number of stages stored in the given ExpressionSet.")
        
        if(is.null(modules))
                stop("Please specify the three modules: early, mid, and late using the argument 'module = list(early = ..., mid = ..., late = ...)'.")
        
        if(any(table(unlist(modules)) > 1))
                stop("Intersecting modules are not defined for the ReductiveHourglassTest.")
        
        nCols <- dim(ExpressionSet)[2]
        score_vector <- vector(mode = "numeric",length = permutations)
        resMatrix <- matrix(NA_real_, permutations,(nCols-2))
        real_age <- vector(mode = "numeric",length = nCols-2)
        real_age <- TAI(ExpressionSet)
        
        
        
        
        
}