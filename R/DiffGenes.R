#' @title Differential Gene Expression Analysis
#' @description 
#' Detect differentially expressed genes (DEGs) in a standard \code{ExpressionSet} object.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param method method to detect differentially expressed genes.
#' @param comparison
#' @param alpha a numeric value specifying the cut-off value above which Genes fulfilling the corresponding fold-change, log-fold-change, or p-value should be retained and returned by \code{DiffGenes}.
#' @param filter.method a method how to \code{alpha} values in multiple stages. Options are \code{"const"}, \code{"min-set"}, and \code{"n-set"}.
#' @param n a numeric value for \code{method = "n-set"}.
#' @param p.adjust.method p value correction method.
#' @param stage.names a character vector specifying the new names of collapsed stages.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Available methods for the detection of differentially expressed genes:
#' 
#' \itemize{
#' \item \code{method = "foldchange"}: ratio of replicate means between developmental stages. Here, the \code{DiffGenes} functions assumes that absolute expression levels are stored in your input \code{ExpresisonSet}.
#' \item \code{method = "log-foldchange"}: difference of replicate log-means between developmental stages.Here, the \code{DiffGenes} functions assumes that \code{log2} transformed expression levels are stored in your input \code{ExpresisonSet}.
#' \item \code{method = "t.test"}: Welch t.test between replicate expression levels between two samples.
#' }
#' 
#' 
#' Exclude non differentially expressed genes from the result dataset:
#' 
#' When specifying the \code{alpha} argument you furthermore, need to specify the \code{filter.method} to decide how non differentially expressed genes should be classified in multiple sample comparisons and which genes should be retained in the final dataset returned by \code{DiffGenes}. In other words, all genes < \code{alpha} based on the following \code{filter.method} are removed from the result dataset.
#' 
#' Following extraction criteria are implemented in this function: 
#' 
#' \itemize{
#' \item \code{const}: all genes that have at least one sample comparison that undercuts or exceeds the \code{alpha} value \code{cut.off} will be excluded from the \code{ExpressionSet}. Hence, for a 7 stage \code{ExpressionSet} genes passing the \code{alpha} threshold in 6 stages will be retained in the \code{ExpressionSet}.
#' \item \code{min-set}: genes passing the \code{alpha} value in \code{ceiling(n/2)} stages will be retained in the \code{ExpressionSet}, where \emph{n} is the number of stages in the \code{ExpressionSet}.
#' \item \code{n-set}: genes passing the \code{alpha} value in \code{n} stages will be retained in the \code{ExpressionSet}. Here, the argument \code{n} needs to be specified.
#' }
#' 
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' 
#' # Detection of DEGs using the fold-change measure
#' DEGs <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[ ,1:8],
#'                   nrep          = 2,
#'                   method        = "foldchange",
#'                   stage.names   = c("S1","S2","S3"))
#' 
#' 
#' head(DEGs)
#' 
#' 
#' # Detection of DEGs using the log-fold-change measure
#' # when choosing method = "log-foldchange" it is assumed that
#' # your input expression matrix stores log2 expression levels 
#' log.DEGs <- DiffGenes(ExpressionSet = tf(PhyloExpressionSetExample[1:5,1:8],log2),
#'                       nrep          = 2,
#'                       method        = "log-foldchange",
#'                       stage.names   = c("S1","S2","S3"))
#' 
#' 
#' head(log.DEGs)
#' 
#' 
#' # Remove fold-change values < 2 from the dataset:
#' 
#' ## first have a look at the range of fold-change values of all genes 
#' apply(DEGs[ , 3:8],2,range)
#' 
#' # now remove genes undercutting the alpha = 2 threshold
#' # hence, remove genes having p-values <= 0.05 in at
#' # least one sample comparison
#' DEGs.alpha <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[1:250 ,1:8],
#'                         nrep          = 2,
#'                         method        = "t.test",
#'                         alpha         = 0.05,
#'                         comparison    = "above",
#'                         filter.method = "n-set",
#'                         n             = 1,
#'                         stage.names   = c("S1","S2","S3"))
#' 
#' # now again have a look at the range and find
#' # that fold-change values of 2 are the min value
#' apply(DEGs.alpha[ , 3:8],2,range)
#' 
#' # now check whether each example has at least one stage with a p-value <= 0.05
#' head(DEGs.alpha)
#' 
#' 
#' @seealso \code{\link{Expressed}}
#' @export

DiffGenes <- function(ExpressionSet,
                      nrep,
                      method          = "foldchange",
                      comparison      = NULL,
                      alpha           = NULL,
                      filter.method   = NULL,
                      n               = NULL,
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
                #combin.stages <- expand.grid(1:nStages,1:nStages)
                combin.stages <- data.frame(Var1 = as.vector(sapply(1:nStages,function(x) rep(x,nStages))),
                           Var2 = rep(1:nStages,nStages))
                
                test_combin_func <- function(x){
                        ifelse(x[1] == x[2],FALSE,TRUE) 
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
                                DEGMatrix[ , i] <- CollapsedExpressionSet[ , idx[1] + 2] - CollapsedExpressionSet[ , idx[2] + 2]
                        }
                }
        }
        
        else if (is.element(method,c("t.test"))){
                
                # in case a constant number of replicates per stage is given
                if (length(nrep) == 1){
                        if ((ncols - 2) %% nrep != 0)
                                stop("The number of stages and the number of replicates do not match.")
                        
                        # automatically rename replicate stages to 1.1, 1.2, ... , n.1, n.2
                        # in case nrep = 2
                        nStages <- (ncols - 2) / nrep
                        # get all combinations of stages to perform
                        # foldchange computations
                        #combin.stages <- expand.grid(1:nStages,1:nStages)
                        combin.stages <- data.frame(Var1 = as.vector(sapply(1:nStages,function(x) rep(x,nStages))),
                                                    Var2 = rep(1:nStages,nStages))
                        
                        test_combin_func <- function(x){
                                ifelse(x[1] == x[2],FALSE,TRUE) 
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
                        
                        if (!is.null(p.adjust.method)){
                                
                                DEGMatrix <- t(apply(DEGMatrix,1,p.adjust,method = p.adjust.method)) 
                        }
                        
                        
                        } else {
                                stop("Something went wrong with the constant number of replicates per stage.
                                     Are you sure that each stage has the same exact number of replicates?")
                                }
                }
                
                if (!is.null(stage.names))
                        names(ExpressionSet)[3:ncol(ExpressionSet)] <- stage.names
        
                DEG.ExpressionSet <- data.frame(ExpressionSet[ , 1:2], DEGMatrix) 
                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                              apply(combin.stages,1,function(x) paste0(names(ExpressionSet)[x[1] + 2],"->",names(ExpressionSet)[x[2] + 2])))
                
                if (!is.null(alpha)){
                        
                        if (!dplyr::between(alpha,
                                          min(DEG.ExpressionSet[ , 3:ncol(DEG.ExpressionSet)]),
                                          max(DEG.ExpressionSet[ , 3:ncol(DEG.ExpressionSet)]))){
                                stop("Please specify a value for alpha that lies within the range of fold-change or p-values.")
                        }
                        
                        if (any(c(is.null(alpha),is.null(filter.method),is.null(comparison))))
                                stop("Arguments alpha, comparison, and filter.method neet to be specified to remove non expressed genes.")
                        
                        return( Expressed(ExpressionSet    = DEG.ExpressionSet,
                                                cut.off    = alpha, 
                                                method     = filter.method,
                                                comparison = comparison,
                                                n          = n) )
                        
                } else {
                        return(DEG.ExpressionSet)
                }
}








