#' @title Combine Replicates in an ExpressionSet
#' @description 
#' This function takes an \code{ExpressionSet} object storing either a constant or variable number
#' of biological or technical replicates per stage and collapses replicate expression levels using a defined \code{FUN} (window function).
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param FUN a window function (e.g., \code{\link{mean}}, \code{\link{median}}, \code{\link{max}}, \code{\link{min}}, etc.) specifying how replicate expression levels should be collapsed.
#' @param stage.names a character vector specifying the new names of collapsed stages.
#' @author Hajk-Georg Drost
#' @examples 
#' data(PhyloExpressionSetExample)
#' 
#' # combine the expression levels of the 2 replicates (const) per stage
#' # using mean as window function and rename new stages: "S1","S2","S3"
#' CollapseReplicates(ExpressionSet = PhyloExpressionSetExample[1:5,1:8], 
#'                    nrep          = 2, 
#'                    FUN           = mean, 
#'                    stage.names   = c("S1","S2","S3"))
#' 
#' # combine the expression levels of the 2 replicates (stage one), 2 replicates (stage two),
#' # and 3 replicates (stage three) using mean as window function 
#' # and rename new stages: "S1","S2","S3"
#' CollapseReplicates(ExpressionSet = PhyloExpressionSetExample[1:5,1:9], 
#'                    nrep          = c(2,2,3), 
#'                    FUN           = mean, 
#'                    stage.names   = c("S1","S2","S3"))
#' 
#' 
#' @export 

CollapseReplicates <- function(ExpressionSet, nrep, FUN, stage.names = NULL){
        
        is.ExpressionSet(ExpressionSet)
        ncols <- dim(ExpressionSet)[2]
        
        # in case a constant number of replicates per stage is given
        if(length(nrep) == 1){
                if((ncols - 2) %% nrep != 0)
                        stop("The number of stages and the number of replicates do not match.")
                
                # automatically rename replicate stages to 1.1, 1.2, ... , n.1, n.2
                # in case nrep = 2
                nStages <- (ncols - 2) / nrep
                replicate.names <- lapply(lapply(1:nStages,rep,times = nrep),paste0,paste0(".",1:nrep))
                stage.names.repl <- as.vector(unlist(replicate.names))
                colnames(ExpressionSet)[3:ncols] <- stage.names.repl
                
                # indices for further computations
                IndexOne <- seq(3, ncol(ExpressionSet), nrep)
                IndexTwo <- seq(3 + nrep - 1, ncol(ExpressionSet), nrep)
                
                DatasetWithCombinedReplicates <- matrix(NA,nrow(ExpressionSet), nStages)
                
                if(length(IndexOne) == length(IndexTwo)){
                        for(i in 1:length(IndexOne)){
                                DatasetWithCombinedReplicates[ , i] <- CollapseFromTo(ExpressionSet,
                                                                                      from = IndexOne[i], 
                                                                                      to   = IndexTwo[i], 
                                                                                      FUN  = match.fun(FUN))
                        }
                } else {
                        stop("Something went wrong with the constant number of replicates per stage. Are you sure that each stage has the same exact number of replicates?")
                }
                # concatenating the geneIDs to the new dataset at the first position 
                FinalDataset <- data.frame(ExpressionSet[,1:2],DatasetWithCombinedReplicates);
                
                if(!is.null(stage.names))
                        names(FinalDataset) <- c("Phylostratum","GeneID",stage.names)
                if(is.null(stage.names))
                        names(FinalDataset) <- c("Phylostratum","GeneID",colnames(FinalDataset)[3:ncol(FinalDataset)])
                
                return(FinalDataset)
        }
        
        # in case a variable number of replicates per stage is given
        if(length(nrep) > 1){
                nStages <- length(nrep)
                
                if(sum(nrep) != (ncols - 2))
                        stop("The number of stages and the number of replicates do not match.")
                
                # automatically rename replicate stages to 1.1, 1.2, ... , n.1, n.2 taking
                # variable replicate numbers per stage into account in case nrep = c(2,2,3)
                # -> 2 replicates at stage 1, 2 replicates at stage 2, 3 replicates at stage 3
                # and so on ...
                f <- function(x){ 
                        unlist(lapply(lapply(x,rep,times = nrep[x]),paste0,paste0(".",1:nrep[x])))
                }
                
                replicate.names <- lapply(1:nStages,f)
                stage.names.repl <- as.vector(unlist(replicate.names))
                colnames(ExpressionSet)[3:ncols] <- stage.names.repl
                DatasetWithCombinedReplicates <- matrix(NA,nrow(ExpressionSet), nStages)
                i <- 3
                for(j in 1:length(nrep)){
                        
                DatasetWithCombinedReplicates[ , j] <-  CollapseFromTo(ExpressionSet,
                                                                       from = i, 
                                                                       to   = i + (nrep[j]-1), 
                                                                       FUN  = match.fun(FUN))
                i <- i + nrep[j]
                        
                }
                
                FinalDataset <- data.frame(ExpressionSet[ ,1:2],DatasetWithCombinedReplicates);
                
                if(!is.null(stage.names))
                        names(FinalDataset) <- c("Phylostratum","GeneID",stage.names)
                if(is.null(stage.names))
                        names(FinalDataset) <- c("Phylostratum","GeneID",colnames(FinalDataset)[3:ncol(FinalDataset)])
                
                return(FinalDataset)
        }
}






