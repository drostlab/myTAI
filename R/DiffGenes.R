#' @title Differential Gene Expression Analysis
#' @description 
#' Detect differentially expressed genes (DEGs) in a standard \code{ExpressionSet} object.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param method method to detect differentially expressed genes.
#' @param p.adjust.method p value correction method.
#' @param stage.names a character vector specifying the new names of collapsed stages.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Available methods for the detection of differentially expressed genes:
#' 
#' \itemize{
#' \item \code{method = "foldchange"}: ratio of replicate means between developmental stages.
#' \item \code{method = "log-foldchange"}: difference of replicate log-means between developmental stages.
#' \item \code{method = "t.test"}: Welch t.test between replicate expression levels between two samples.
#' }
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' 
#' # Detection of DEGs using the fold-change measure
#' DEGs <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[1:5,1:8],
#'                   nrep          = 2,
#'                   method        = "foldchange",
#'                   stage.names   = c("S1","S2","S3"))
#' 
#' 
#' head(DEGs)
#' @export

DiffGenes <- function(ExpressionSet,
                      nrep,
                      method          = "foldchange",
                      p.adjust.method = NULL,
                      stage.names     = NULL){
        
        is.ExpressionSet(ExpressionSet)
        
        if (!is.element(method,c("foldchange","log-foldchange","t.test")))
                stop("Please enter a method to detect differentially expressed genes that is implemented in DiffGenes().")
        
        ncols <- ncol(ExpressionSet)
        
        if (is.element(method,c("foldchange","log-foldchange"))){
                CollapsedExpressionSet <- CollapseReplicates(ExpressionSet = ExpressionSet,
                                                             nrep          = nrep,
                                                             FUN           = mean,
                                                             stage.names   = stage.names)
                
                nStages <- ncol(CollapsedExpressionSet) - 2
                
                # get all combinations of stages to perform
                # foldchange computations
                combin.stages <- expand.grid(1:nStages,1:nStages)
                
                test_combin_func <- function(x){
                        if (x[1] == x[2])
                                return (FALSE)
                        
                        if (x[1] != x[2])
                                return (TRUE)
                }
                
                # delete all comparisons: 1->1, 2->2, 3->3, ...
                false_comb <- which(!apply(combin.stages,1,test_combin_func))
                combin.stages <- as.data.frame(combin.stages[-false_comb, ])
                
                idx <- vector("numeric",2)
                DEGMatrix <- matrix(NA_real_,nrow = nrow(CollapsedExpressionSet),ncol = nrow(combin.stages))
                
                for (i in 1:nrow(combin.stages)){
                        idx <- as.numeric(combin.stages[i, ])
                        
                        if (method == "foldchange"){
                                DEGMatrix[ , i] <- CollapsedExpressionSet[ , idx[1] + 2] / CollapsedExpressionSet[ , idx[2] + 2]
                        }
                        
                        if (method == "log-foldchange"){
                                DEGMatrix[ , i] <- log2(CollapsedExpressionSet[ , idx[1] + 2]) - log2(CollapsedExpressionSet[ , idx[2] + 2])
                        }
                }
        }
        
        else if (is.element(method,c("t.test"))){
                
                # in case a constant number of replicates per stage is given
                if(length(nrep) == 1){
                        if((ncols - 2) %% nrep != 0)
                                stop("The number of stages and the number of replicates do not match.")
                        
                        # automatically rename replicate stages to 1.1, 1.2, ... , n.1, n.2
                        # in case nrep = 2
                        nStages <- (ncols - 2) / nrep
                        # get all combinations of stages to perform
                        # foldchange computations
                        combin.stages <- expand.grid(1:nStages,1:nStages)
                        
                        test_combin_func <- function(x){
                                if (x[1] == x[2])
                                        return (FALSE)
                                
                                if (x[1] != x[2])
                                        return (TRUE)
                        }
                        # delete all comparisons: 1->1, 2->2, 3->3, ...
                        false_comb <- which(!apply(combin.stages,1,test_combin_func))
                        combin.stages <- as.data.frame(combin.stages[-false_comb, ])
                        
                        idx <- vector("numeric",2)
                        DEGMatrix <- matrix(NA_real_,nrow = nrow(ExpressionSet),ncol = nrow(combin.stages))
                        # indices for further computations
                        IndexOne <- seq(1, ncol(ExpressionSet)-2, nrep)
                        IndexTwo <- seq(1 + nrep - 1, ncol(ExpressionSet)-2, nrep)
                        
                        for (k in 1:nrow(combin.stages)){
                                idx <- as.numeric(combin.stages[k, ])
                                # perform Welch t-test
                                DEGMatrix[ , k] <-  apply(ExpressionSet[, 3:ncol(ExpressionSet)], 1, function(x){
                                        t.test(x[seq(IndexOne[idx[1]],IndexTwo[idx[1]])],
                                               x[seq(IndexOne[idx[2]],IndexTwo[idx[2]])],
                                               alternative = "two.sided",
                                               var.equal   = FALSE)$p.value
                                })
                        }
                        
                        if(!is.null(p.adjust.method)){
                                
                                DEGMatrix <- t(apply(DEGMatrix,1,p.adjust,method = p.adjust.method)) 
                        }
                        
                        
                        } else {
                                stop("Something went wrong with the constant number of replicates per stage.
                                     Are you sure that each stage has the same exact number of replicates?")
                                }
                }
                
                if(!is.null(stage.names))
                        names(ExpressionSet)[3:ncol(ExpressionSet)] <- stage.names
        
                DEG.ExpressionSet <- data.frame(ExpressionSet[ , 1:2], DEGMatrix) 
                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                              apply(combin.stages,1,function(x) paste0(names(ExpressionSet)[x[1] + 2],"->",names(ExpressionSet)[x[2] + 2])))
                return(DEG.ExpressionSet)
}








