#' @title Plot the Quality of Biological Replicates
#' @description This function performs several quality checks to validate the
#' biological variation between replicates and stages (experiments).
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param FUN a function that should be applied to quantify the variablity or quality of replicates.
#' The default function is the log(var(x)) quantifying the variance between replicates.
#' @param legend.pos the position of the legend, e.g. 'topleft', or 'topright' (see \code{\link{legend}}).
#' @param stage.names a character vector specifying the stage names.
#' @param ... additional graphics parameters.
#' @author Hajk-Georg Drost
#' @details The following quality checks can be performed:
#' \itemize{
#' \item Quantification of variability between replicates as density function. 
#' }
#' 
#' @examples 
#' data(PhyloExpressionSetExample)
#' 
#' # visualize log(var(x)) between replicates for each gene and developmental stage 
#' PlotReplicateQuality(ExpressionSet = PhyloExpressionSetExample[1:5000 , 1:8],
#'                      nrep          = 2,
#'                      legend.pos   = "topright",
#'                      ylim          = c(0,0.2),
#'                      lwd           = 6)
#'                      
#' 
#' @export

PlotReplicateQuality <- function(ExpressionSet,
                                 nrep,
                                 FUN         = function(x) log(stats::var(x)),
                                 legend.pos  = "topleft",
                                 stage.names = NULL, ...){
        
        
        if (!all(sapply(nrep,function(x) x > 1, simplify = TRUE)))
                stop("Please insert at least 2 replicates per stage.")
        
        ncols <- ncol(ExpressionSet)
        
        if (length(nrep) == 1){
                if ((ncols - 2) %% nrep != 0)
                        stop("The number of stages and the number of replicates do not match.")
                # in case nrep = 2
                nStages <- (ncols - 2) / nrep
                # get all combinations of stages to perform
                # t-test computations
        }
        
        else if (length(nrep) > 1){
                if (!((ncols - 2) == sum(nrep)))
                        stop("The number of stages and the number of replicates do not match.")
                nStages <- length(nrep)
        }
        
        stage.cols <- re.colors(nStages)
        
        custom.FUN <- match.fun(FUN)
        # receive variance distribution of replicates
        CollapsedExpressionSet <- CollapseReplicates(ExpressionSet = ExpressionSet,
                                                     nrep          = nrep,
                                                     FUN           = custom.FUN)
        
        col.index <- 1
        graphics::plot(stats::density(CollapsedExpressionSet[ , 3]), col = stage.cols[1],main = "Distributions of replicate log variances", ...)
        apply(CollapsedExpressionSet[ , 4:(3 + nStages - 1)], 2 ,function(x) {
                
                col.index <<- col.index + 1
                graphics::lines(stats::density(x),col = stage.cols[col.index], ...)
                
                
        })
        
        if (is.null(stage.names))
                graphics::legend(legend.pos, bty = "n", legend = paste0("S",1:nStages), fill = stage.cols, ncol = ifelse(nStages <= 4, 1, floor(nStages/2)))
        
        if (!is.null(stage.names))
                graphics::legend(legend.pos, bty = "n", legend = stage.names, fill = stage.cols, ncol = ifelse(nStages <= 4, 1, floor(nStages/2)))
        
        
}


