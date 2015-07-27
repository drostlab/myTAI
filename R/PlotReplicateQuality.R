#' @title Plot the Quality of Biological Replicates
#' @description This function performs several quality checks to validate the
#' biological variation between replicates and stages (experiments).
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param stage.names a character vector specifying the new names of collapsed stages.
#' @param ... additional graphics parameters.
#' @author Hajk-Georg Drost
#' @details The following quality checks can be performed:
#' \itemize{
#' \item
#' \item
#' }
#' 
#' @export
PlotReplicateQuality <- function(ExpressionSet, nrep, stage.names, ...){
        
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
        
        # receive variance distribution of replicates
        CollapsedExpressionSet <- CollapseReplicates(ExpressionSet = ExpressionSet,
                                                     nrep          = nrep,
                                                     FUN           = function(x) log(var(x)),
                                                     stage.names   = stage.names)
        
        
        graphics::plot(density(CollapsedExpressionSet[ , 3]), col = stage.cols[1],main = "Distributions of replicate log variances", ...)
        apply(CollapsedExpressionSet[ , 4:(3 + nStages - 1)], 2 ,function(x) {
                
                graphics::lines(density(x),col = stage.cols[col.index], ...)
                col.index <- col.index + 1
                
        })
        
        graphics::legend("topleft", bty = "n", legend = stage.names, fill = stage.cols)
}


