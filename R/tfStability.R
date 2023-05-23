#' @title Perform Permutation Tests Under Different Transformations
#' @description \emph{tfStability} aims to statistically evaluate the
#' stability of \code{\link{ReductiveHourglassTest}}, \code{\link{FlatLineTest}}, 
#' \code{\link{ReverseHourglassTest}}, \code{\link{EarlyConservationTest}}, or 
#' \code{\link{LateConservationTest}} 
#' (all based on \code{\link{TAI}} or \code{\link{TDI}} computations) against different
#' data transformations.
#' The corresponding p-value quantifies the probability that a given TAI or TDI pattern (or any phylotranscriptomics pattern) 
#' does not support the chosen test. A p-value < 0.05 indicates that the corresponding phylotranscriptomics pattern does
#' indeed support the chosen test.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param modules a list storing three elements: early, mid, and late. Each element expects a numeric
#' vector specifying the developmental stages or experiments that correspond to each module. 
#' For example, \code{module} = list(early = 1:2, mid = 3:5, late = 6:7) devides a dataset 
#' storing seven developmental stages into 3 modules.
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{ReductiveHourglassTest}.
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be: 
#' \itemize{
#' \item \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line
#' \item \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern)
#' \item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
#' \item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a early conservation pattern (low-high-high pattern)
#' \item \code{TestStatistic} = \code{"LateConservationTest"} : Statistical test for the existence of a late conservation pattern (high-high-low pattern)
#' }
#' @param transforms a character vector of any valid function that transforms gene expression levels.
#' @param pseudocount any valid number to the expression matrix prior to transformation.
#' @details 
#' An assessment for the stability of a data transforms on the results of the permutation test of choice.
#' The reductive hourglass test is a permutation test based on the following test statistic. 
#' For details, please consult \code{\link{ReductiveHourglassTest}}, \code{\link{FlatLineTest}}, 
#' \code{\link{ReverseHourglassTest}} or \code{\link{EarlyConservationTest}}
#'
#' @return a vector object containing the vector elements:
#' 
#' \code{p.value} : the p-value quantifying the statistical significance (depending on the chosen test) of the given phylotranscriptomics pattern under the given data transformation(s).
#' @references 
#' 
#' Lotharukpong JS et al. (2023) bioRxiv
#'
#'   
#' @author Jaruwatana Sodai Lotharukpong
#' @seealso \code{\link{rhScore}}, \code{\link{bootMatrix}}, \code{\link{FlatLineTest}},
#' \code{\link{ReverseHourglassTest}}, \code{\link{EarlyConservationTest}}, 
#' \code{\link{ReductiveHourglassTest}}, \code{\link{PlotSignature}}, \code{\link{LateConservationTest}}
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#'
#' # perform the reductive hourglass test for a PhyloExpressionSet
#' # here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
#' # stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
#' # module 3 = late
#' tfStability(ExpressionSet = PhyloExpressionSetExample,
#'                      TestStatistic = "ReductiveHourglassTest",
#'                      transforms = c("log1p", "sqrt", "none"),
#'                      modules = list(early = 1:2, mid = 3:5, late = 6:7))
#'
#'
#' @export

# returns p value only
tfStability <- function(ExpressionSet,
                        TestStatistic      = "FlatLineTest",
                        transforms         = c("log1p", "sqrt", "none"),
                        modules            = NULL,
                        permutations       = 1000,
                        pseudocount        = 0)
{
  
  myTAI::is.ExpressionSet(ExpressionSet)
  
  if(!TestStatistic %in% c("FlatLineTest", "ReductiveHourglassTest", "ReverseHourglassTest", "EarlyConservationTest", "LateConservationTest"))
    stop("Please select the availagetble test: 'FlatLineTest', 'ReductiveHourglassTest', 'ReverseHourglassTest', 'EarlyConservationTest' or 'LateConservationTest' using the argument test = 'FlatLineTest'", call. = FALSE)
  
  if(TestStatistic %in% c("ReductiveHourglassTest", "ReverseHourglassTest", "EarlyConservationTest") & is.null(modules))
    stop("Please specify the three modules: early, mid, and late using the argument 'module = list(early = ..., mid = ..., late = ...)'.", call. = FALSE)
  
  # if (!(is.element(transforms, c("log1p", "sqrt", "none", "log2", "log", "log10")))){
  #   stop("Please select the available transformations: 'log1p', 'sqrt', 'log2', 'log', 'log10' or 'none' using the argument 'transforms = c('log1p', 'sqrt', 'none')'.", call. = FALSE)
  # }
  
  # create placeholder vector
  vec_res <- NULL
  
  # output transforms in a vector
  for (i in transforms) {
    if(i == "none")
      tfExpressionSet <- ExpressionSet
    else if(i %in% c('log2', 'log', 'log10'))
      tfExpressionSet <- tf(ExpressionSet, FUN = i, pseudocount = pseudocount)
    else
      tfExpressionSet <- tf(ExpressionSet, FUN = i)
    test_function <- base::match.fun(TestStatistic)
    if(TestStatistic == "FlatLineTest")
      vec_res[i] <- test_function(tfExpressionSet,
                                  permutations       = permutations)[["p.value"]]
    else
      vec_res[i] <- test_function(tfExpressionSet,
                                  modules            = modules,
                                  permutations       = permutations)[["p.value"]]
  }
  return(vec_res)
}
