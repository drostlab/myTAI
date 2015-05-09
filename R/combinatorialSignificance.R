#' @title Compute the Statistical Significance of Each Replicate Combination
#' @description In case a PhyloExpressionSet or DivergenceExpressionSet stores replicates for each
#' developmental stage or experiment, this function allows to compute the p-values quantifying
#' the statistical significance of the underlying pattern for all combinations of replicates.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param replicates a numeric vector storing the number of replicates within each developmental stage or experiment.
#' In case replicate stores only one value, then the function assumes that each developmental stage or experiment
#' stores the same number of replicates.
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Default is \code{TestStatistic} = \code{"FlatLineTest"}.
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}.
#' @param parallel a boolean value specifying whether parallel processing (multicore processing) shall be performed.
#' @details 
#' 
#' The intention of this analysis is to validate that there exists no sequence of replicates 
#' (for all possible combination of replicates) that results in a non-significant pattern,
#' when the initial pattern with combined replicates was shown to be significant.
#' 
#' A small Example: 
#' 
#'      
#' Assume PhyloExpressionSet stores 3 developmental stages with 3 replicates measured for each stage.
#' The 9 replicates in total are denoted as: 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3. Now the function computes the
#' statistical significance of each pattern derived by the corresponding combination of replicates, e.g.
#'
#' \itemize{
#' \item 1.1, 2.1, 3.1 -> p-value for combination 1
#'
#' \item 1.1, 2.2, 3.1  -> p-value for combination 2
#'
#' \item 1.1, 2.3, 3.1 -> p-value for combination 3
#'
#' \item 1.2, 2.1, 3.1 -> p-value for combination 4
#'
#' \item 1.2, 2.1, 3.1 -> p-value for combination 5
#'
#' \item 1.2, 2.1, 3.1 -> p-value for combination 6
#'
#' \item 1.3, 2.1, 3.1 -> p-value for combination 7
#'
#' \item 1.3, 2.2, 3.1 -> p-value for combination 8
#'
#' \item 1.3, 2.3, 3.1 -> p-value for combination 9
#' 
#' \item \dots
#' }
#' This procedure yields 27 p-values for the \eqn{3^3} (\eqn{n_stages^n_replicates}) replicate combinations.
#' 
#' Note, that in case you have a large amount of stages/experiments and a large amount of replicates
#' the computation time will increase by \eqn{n_stages^n_replicates}. For 11 stages and 4 replicates, 4^11 = 4194304 p-values have to be computed. 
#' Each p-value computation itself is based on a permutation test running with 1000 or more permutations. Be aware that this might take some time.
#'
#' The p-value vector returned by this function can then be used to plot the p-values to see
#' whether an critical value \eqn{\alpha} is exeeded or not (e.g. \eqn{\alpha = 0.05}).
#' 
#' 
#' The function receives a standard PhyloExpressionSet or DivergenceExpressionSet object and a vector storing the number of replicates present in each stage or experiment. Based on these arguments the function computes all possible replicate combinations using the \code{\link{expand.grid}} function and performs a permutation test (either a \code{\link{FlatLineTest}} for each replicate combination. The \emph{permutation} parameter of this function specifies the number of permutations that shall be performed for each permutation test. When all p-values are computed, a numeric vector storing the corresponding p-values for each replicate combination is returned.  
#' 
#' In other words, for each replicate combination present in the PhyloExpressionSet or DivergenceExpressionSet object, the TAI or TDI pattern of the corresponding replicate combination is tested for its statistical significance based on the underlying test statistic.
#' 
#' This function is also able to perform all computations in parallel using multicore processing. The underlying statistical tests are written in C++ and optimized for fast computations.
#' 
#' @return a numeric vector storing the p-values returned by the underlying test statistic for all possible replicate combinations.
#' @references Drost HG et al. (2015). \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{expand.grid}}, \code{\link{FlatLineTest}}
#' @examples
#' 
#' # load a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # we assume that the PhyloExpressionSetExample 
#' # consists of 3 developmental stages 
#' # and 2 replicates for stage 1, 3 replicates for stage 2, 
#' # and 2 replicates for stage 3
#' # FOR REAL ANALYSES PLEASE USE: permutations = 1000 or 10000
#' # BUT NOTE THAT THIS TAKES MUCH MORE COMPUTATION TIME
#' p.vector <- CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample, 
#'                                       replicates    = c(2,3,2), 
#'                                       TestStatistic = "FlatLineTest", 
#'                                       permutations  = 10, 
#'                                       parallel      = FALSE)
#'
#'
#'
#'
#' @import foreach
#' @export
CombinatorialSignificance <- function(ExpressionSet,
                                      replicates,
                                      TestStatistic = "FlatLineTest", 
                                      permutations  = 1000, 
                                      parallel      = FALSE)
{
  
  is.ExpressionSet(ExpressionSet)
  
  if(!is.element(TestStatistic, c("FlatLineTest"))){
    stop("Please enter a correct string for the test statistic: 'FlatLineTest'.")
  }
  
  ncols <- dim(ExpressionSet)[2]
  
  # in case all stages have the exact same number of replicates
  if(length(replicates) == 1){
  
     if((ncols - 2) %% replicates != 0)
        stop("The number of stages and the number of replicates do not match.")
  
     nStages <- (ncols - 2) / replicates
     replicateName.List <- lapply(lapply(1:nStages,rep,times = replicates),paste0,paste0(".",1:replicates))
     stageNames <- as.vector(unlist(replicateName.List))
  
     colnames(ExpressionSet)[3:ncols] <- stageNames
  
     # compute all possible combinations
     combinatorialMatrix <- expand.grid(replicateName.List, stringsAsFactors = FALSE)
  }
  
  # in case stages have variable number of replicates per stage
  if(length(replicates) > 1){
    
    nStages <- length(replicates)
    
    if(sum(replicates) != (ncols - 2))
      stop("The number of stages and the number of replicates do not match.")
    
    f <- function(x){ 
      unlist(lapply(lapply(x,rep,times = replicates[x]),paste0,paste0(".",1:replicates[x])))
    }
    
    
    replicateName.List <- lapply(1:nStages,f)
    stageNames <- as.vector(unlist(replicateName.List))
    colnames(ExpressionSet)[3:ncols] <- stageNames
    
    combinatorialMatrix <- expand.grid(replicateName.List, stringsAsFactors = FALSE)
    
  }
  
  nCombinations <- dim(combinatorialMatrix)[1]
  p.vals <- vector(mode = "numeric",length = nCombinations)
  first_cols_names <- as.character(colnames(ExpressionSet)[1:2])
  
  if(parallel){
          
    
          
    ### Parallellizing the sampling process using the 'doParallel' and 'parallel' package
    ### register all given cores for parallelization
    par_cores <- parallel::makeForkCluster(parallel::detectCores())
    doParallel::registerDoParallel(par_cores)
    
    # perform the sampling process in parallel
    p.vals <- as.vector(foreach::foreach(i              = 1:nCombinations,
                                         .combine       = "c",
                                         .errorhandling = "stop") %dopar% {
      
      FlatLineTest(as.data.frame(ExpressionSet[c(first_cols_names,as.character(combinatorialMatrix[i , ]))]), permutations = permutations)$p.value
      
    })
    
    parallel::stopCluster(par_cores)
    
  }
  
  
  if(!parallel){
    # sequential computations of p-values 
#     if(nCombinations > 10){
#       # initializing the progress bar
#       progressBar <- txtProgressBar(min = 1,max = nCombinations,style = 3)
#       
#     }
    
    for(i in 1:nCombinations){
      
      p.vals[i] <- FlatLineTest(as.data.frame(ExpressionSet[c(first_cols_names,as.character(combinatorialMatrix[i , ]))]), permutations = permutations)$p.value
      
#       if(nCombinations > 10){
#         # printing out the progress
#         setTxtProgressBar(progressBar,i)
#       }
      
    }
  }
  
  cat("\n")
  return(p.vals)
}

