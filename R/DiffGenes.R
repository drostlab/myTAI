#' @title Differential Gene Expression Analysis
#' @description 
#' Detect differentially expressed genes (DEGs) in a standard \code{ExpressionSet} object.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param nrep either a numeric value specifying the constant number of replicates per stage or a numeric vector specifying the variable number of replicates for each stage position.
#' @param method method to detect differentially expressed genes.
#' @param lib.size the library sizes to equalize library sizes by quantile-to-quantile normalization (see \code{\link[edgeR]{equalizeLibSizes}}).
#' @param p.adjust.method p value correction method.
#' @param comparison a character string specifying whether genes having fold-change or p-values
#'  below, above, or below AND above (both) the \code{alpha} value should be excluded from the dataset.
#'  In case \code{comparison = "both"} is chosen, the \code{cut.off} argument must be a two dimensional vector defining the lower \code{alpha} value at the first position and the upper \code{alpha} value
#' at the second position. 
#' @param alpha a numeric value specifying the cut-off value above which Genes fulfilling the corresponding fold-change, log-fold-change, or p-value should be retained and returned by \code{DiffGenes}.
#' @param filter.method a method how to \code{alpha} values in multiple stages. Options are \code{"const"}, \code{"min-set"}, and \code{"n-set"}.
#' @param n a numeric value for \code{method = "n-set"}.
#' @param stage.names a character vector specifying the new names of collapsed stages.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' All methods to perform dection of differentially expressed genes assume that your input
#' dataset has been normalized before passing it to \emph{DiffGenes}. For RNA-Seq data
#' \emph{DiffGenes} assumes that the libraries have been normalized to have the same size, i.e.,
#' to have the same expected column sum under the null hypothesis. If this isn't the case
#' please run \code{\link[edgeR]{equalizeLibSizes}} before calling \emph{DiffGenes}. 
#' 
#' Available methods for the detection of differentially expressed genes:
#' 
#' \itemize{
#' \item \code{method = "foldchange"}: ratio of replicate geometric means between developmental stages.
#'  Here, the \emph{DiffGenes} functions assumes that absolute expression levels are stored in your input \code{ExpresisonSet}.
#' \item \code{method = "log-foldchange"}: difference of replicate arithmetic means between developmental stages. Here, the \emph{DiffGenes} functions assumes that \emph{log2a}
#' transformed expression levels are stored in your input \code{ExpresisonSet}.
#' \item \code{method = "t.test"}: Welch t.test between replicate expression levels of two samples.
#' \item \code{method = "wilcox.test"}: Wilcoxon Rank Sum Test between replicate expression levels of two samples.
#' \item \code{method = "doubletail"}: Computes two-sided p-values by doubling the smaller tail probability (see \code{\link[edgeR]{exactTestDoubleTail}} for details).
#' \item \code{method = "smallp"}: Performs the method of small probabilities as proposed by Robinson and Smyth (2008) (see \code{\link[edgeR]{exactTestBySmallP}} for details).
#' \item \code{method = "deviance"}: Uses the deviance goodness of fit statistics to define the rejection region, and is therefore equivalent to a conditional likelihood ratio test (see \code{\link[edgeR]{exactTestByDeviance}} for details).
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
#'                   comparison    = "below",
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
#'                       comparison    = "below",
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
#' apply(DEGs.alpha[ , 3:5],2,range)
#' 
#' # now check whether each example has at least one stage with a p-value <= 0.05
#' head(DEGs.alpha)
#' 
#' @note In case input \code{ExpressionSet} objects store 0 values, internally all expression levels are 
#' shifted by \code{+1} to allow sufficient fold-change and p-value computations. Additionally, a warning
#' is printed to the console in case expression levels have been automatically shifted.
#' @seealso \code{\link{Expressed}}
#' @export

DiffGenes <- function(ExpressionSet,
                      nrep,
                      method          = "foldchange",
                      lib.size        = NULL,
                      p.adjust.method = NULL,
                      comparison      = NULL,
                      alpha           = NULL,
                      filter.method   = NULL,
                      n               = NULL,
                      stage.names     = NULL){
        
        is.ExpressionSet(ExpressionSet)
        
        if (!is.element(method,c("foldchange","log-foldchange","t.test",
                                 "wilcox.test","doubletail","smallp","deviance")))
                stop ("Please enter a method to detect differentially expressed genes that is implemented in DiffGenes().", call. = FALSE)
        
        ncols <- ncol(ExpressionSet)
        
        # test whether or not input data stores 0 values
        # if yes -> shift expression levels by +1
        if (any(apply(ExpressionSet[ , 3:ncols],2, function(x) x == 0))){
                
                ExpressionSet <- tf(ExpressionSet, function(x) x + 1)
                warning ("Your input ExpressionSet stores 0 values. Therefore, all expression levels have been shifted by +1 to allow sufficient fold-change or p-value computations.")
        }
        
        if (is.element(method,c("foldchange","log-foldchange"))){
                
                # assume absolute expression levels: so accumulation (window) function = geom.mean
                if ( method == "foldchange"){
                        CollapsedExpressionSet <- CollapseReplicates(ExpressionSet = ExpressionSet,
                                                                     nrep          = nrep,
                                                                     FUN           = geom.mean,
                                                                     stage.names   = stage.names)
                        
                } 
                
                # assume log expression levels: so accumulation (window) function = mean
                if (method == "log-foldchange"){
                        CollapsedExpressionSet <- CollapseReplicates(ExpressionSet = ExpressionSet,
                                                                     nrep          = nrep,
                                                                     FUN           = mean,
                                                                     stage.names   = stage.names)
                }
                
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
        
        else if (is.element(method,c("t.test","wilcox.test","doubletail","smallp","deviance"))){
                
                # in case a constant number of replicates per stage is given
                if (length(nrep) > 0){
                        
                        # automatically rename replicate stages to 1.1, 1.2, ... , n.1, n.2
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
                        
                        # determine all possible combinations of stagewise comparisons
                        combin.tbl <- t(utils::combn(1:nStages,m = 2))
                        combin.stages <- data.frame(Var1 = combin.tbl[ , 1],
                                                    Var2 = combin.tbl[ , 2])

                        idx <- vector("numeric",2)
                        DEGMatrix <- matrix(NA_real_,nrow = nrow(ExpressionSet),ncol = nrow(combin.stages))
                        # determine stage indices for nrep = const
                        if (length(nrep) == 1){
                                # indices for further computations
                                IndexOne <- seq(1, ncol(ExpressionSet)-2, nrep)
                                IndexTwo <- seq(1 + nrep - 1, ncol(ExpressionSet)-2, nrep)
                        }
                        
                        else if (length(nrep) > 1){
                                
                                IndexOne <- GetColumnIndexFromTo(nrep)[ , "From"]
                                IndexTwo <- GetColumnIndexFromTo(nrep)[ , "To"]
                        }
                        
                        for (k in 1:nrow(combin.stages)){
                                idx <- as.numeric(combin.stages[k, ])
                                # perform Welch t-test
                                
                                if (method == "t.test"){
                                        DEGMatrix[ , k] <-  apply(ExpressionSet[, 3:ncol(ExpressionSet)], 1, function(x){
                                                stats::t.test(x[seq(IndexOne[idx[1]],IndexTwo[idx[1]])],
                                                              x[seq(IndexOne[idx[2]],IndexTwo[idx[2]])],
                                                              alternative = "two.sided",
                                                              var.equal   = FALSE)$p.value
                                        })
                                }
                                        
                                if (method == "wilcox.test"){
                                        DEGMatrix[ , k] <-  apply(ExpressionSet[, 3:ncol(ExpressionSet)], 1, function(x){
                                                                           stats::wilcox.test(x[seq(IndexOne[idx[1]],IndexTwo[idx[1]])],
                                                                           x[seq(IndexOne[idx[2]],IndexTwo[idx[2]])],
                                                                           alternative = "two.sided")$p.value
                                                                })
                                } 
                                
                                if (is.element(method, c("doubletail","smallp","deviance"))){
                                                
                                        if (is.null(lib.size))
                                                stop ("Please specify a library size.", call. = FALSE)
                                                        
                                        # combine stage replicates to a common replicate matrix
                                        # fulfilling the edgeR exactTest() function specification
                                        compStageMatrix <- cbind(ExpressionSet[ , 2 + seq(IndexOne[idx[1]],IndexTwo[idx[1]])],
                                                                 ExpressionSet[ , 2 + seq(IndexOne[idx[2]],IndexTwo[idx[2]])])
                                        
                                        # number of replicates in each stage to compare
                                        nrep.1 <- length(seq(IndexOne[idx[1]],IndexTwo[idx[1]]))
                                        nrep.2 <- length(seq(IndexOne[idx[2]],IndexTwo[idx[2]]))
                                        
                                        exactTestObject <- edgeR::DGEList(counts   = compStageMatrix,
                                                                          group    = c(rep(1,nrep.1),rep(2,nrep.2)),
                                                                          lib.size = rep(lib.size, nrep.1 + nrep.2))
                                        
                                        DEGMatrix[ , k] <- edgeR::exactTest(object           = exactTestObject,
                                                                            pair             = 1:2,
                                                                            dispersion       = 0.2,
                                                                            rejection.region = method,
                                                                            big.count        = 900,
                                                                            prior.count      = 0.125)$table[ , "PValue"]
                                                        
#                                       DEGMatrix[ , k] <-  edgeR::exactTestDoubleTail(y1 = ExpressionSet[ , 2 + seq(IndexOne[idx[1]],IndexTwo[idx[1]])],
#                                                                                 y2 = ExpressionSet[ , 2 + seq(IndexOne[idx[2]],IndexTwo[idx[2]])],
#                                                                                 dispersion = 0.2,
#                                                                                 big.count = 900)
                                                        
                                        }
                                
                        }
                        
                        if (!is.null(p.adjust.method)){
                                
                                DEGMatrix <- apply(DEGMatrix,2,stats::p.adjust,method = p.adjust.method) 
                        }
                } else {
                        stop("Something went wrong with the number of replicates per stage.
                                     Are you sure that each stage has the correct number of replicates?", call. = FALSE)
                                }
                }

                DEG.ExpressionSet <- data.frame(ExpressionSet[ , 1:2], DEGMatrix) 
                
                if (!is.null(stage.names)){
                        
                        DefaultStageNames <- stage.names
                        
                        if (is.element(method,c("foldchange","log-foldchange"))){
                                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                                              apply(combin.stages,1,function(x) paste0(DefaultStageNames[x[1]],"->",DefaultStageNames[x[2]])))
                        } else {
                                
                                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                                              apply(combin.stages,1,function(x) paste0(DefaultStageNames[x[1]],"<->",DefaultStageNames[x[2]])))
                        }
                        
                        
                } else {
                        
                        if (is.element(method,c("foldchange","log-foldchange"))){
                                
                                DefaultStageNames <- paste0("S",1:nStages)
                                
                                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                                              apply(combin.stages,1,function(x) paste0(DefaultStageNames[x[1]],"->",DefaultStageNames[x[2]])))
                        } else {
                                
                                DefaultStageNames <- paste0("S",1:nStages)
                                
                                names(DEG.ExpressionSet) <- c(names(ExpressionSet)[1:2],
                                                              apply(combin.stages,1,function(x) paste0(DefaultStageNames[x[1]],"<->",DefaultStageNames[x[2]])))
                                
                        }
                        
                }
                
                
                if (!is.null(alpha)){
                        
                        if (is.element(method, c("foldchange","log-foldchange"))){
                                if (!dplyr::between(alpha,
                                                    min(DEG.ExpressionSet[ , 3:ncol(DEG.ExpressionSet)]),
                                                    max(DEG.ExpressionSet[ , 3:ncol(DEG.ExpressionSet)]))){
                                        stop("Please specify a value for alpha that lies within the range of fold-change or p-values.", call. = FALSE)
                                }       
                        }
                        
                        else if (is.element(method, c("t.test"))){
                                
                                if (!dplyr::between(alpha,0,1))
                                        stop("Please specify a value for alpha that lies within the range of fold-change or p-values.", call. = FALSE)
                        }
                        
                        
                        if (any(c(is.null(alpha),is.null(filter.method),is.null(comparison))))
                                stop ("Arguments alpha, comparison, and filter.method neet to be specified to remove non expressed genes.", call. = FALSE)
                        
                        return( Expressed(ExpressionSet    = DEG.ExpressionSet,
                                                cut.off    = alpha, 
                                                method     = filter.method,
                                                comparison = comparison,
                                                n          = n) )
                        
                } else {
                        return(DEG.ExpressionSet)
                }
}








