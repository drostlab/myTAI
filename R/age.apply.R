#' @title Age Category Specific apply Function
#' @description 
#' This function performs the split-apply-combine methodology on Phylostrata or Divergence Strata stored within the input ExpressionSet.
#' 
#' This function is very useful to perform any phylostratum or divergence-stratum specific analysis.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param FUN a function to be performed on the corresponding expression matrix of each phylostratum or divergence-stratum.
#' @param ... additional arguments of FUN.
#' @param as.list a boolean value specifying whether the output format shall be a matrix or a list object.
#' @details This function uses the \code{\link{split}} function to subset the expression matrix into
#' phylostratum specific sub-matrices. Internally using \code{\link{lapply}}, any function can
#' be performed to the sub-matrices. The return value of this function is a numeric matrix storing
#' the return values by \code{FUN} for each phylostratum and each developmental stage s. 
#' Note that the input \code{FUN} must be an function that can be applied to a matrix (e.g., \code{\link{colMeans}} or \code{\link{RE}}). 
#' In case you use an an anymous function you coud use \code{function(x) apply(x , 2 , var)} as an example to compute the variance of each phylostratum and each
#' developmental stage s.
#' @return Either a numeric matrix storing the return values of the applied function for each age class
#' or a numeric list storing the return values of the applied function for each age class in a list.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{split}}, \code{\link{tapply}}, \code{\link{lapply}}, \code{\link{RE}}, \code{\link{REMatrix}}
#' @examples
#'  
#'  # source the example dataset
#'  data(PhyloExpressionSetExample)
#'  
#' # Example 1
#' # get the relative expression profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, RE)
#'
#' # this is analogous to 
#' REMatrix(PhyloExpressionSetExample)

#' # Example 2
#' # compute the mean expression profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, colMeans)
#'
#' # Example 3
#' # compute the variance profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , var))
#'
#' # Example 4
#' # compute the range for each phylostratum
#' # Note: in this case, the range() function returns 2 values for each phylostratum
#' # and each developmental stage, hence one should use the argument 'as.list = TRUE'
#' # to make sure that the results are returned properly 
#' age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , range), as.list = TRUE)
#' 
#' @export
age.apply <- function(ExpressionSet,FUN, ... ,as.list = FALSE)
{
        
        is.ExpressionSet(ExpressionSet)
        
        f <- match.fun(FUN)
        ncols <- ncol(ExpressionSet)
        s <- split(ExpressionSet, ExpressionSet[ , 1])
        
        if(!as.list){
                res <- t(as.data.frame(lapply(s , function(x) f(as.matrix(x[ , 3:ncols]) , ...))))
                rownames(res) <- levels(as.factor(ExpressionSet[ , 1]))
        }
        
        if(as.list){
                res <- lapply(s , function(x) f(as.matrix(x[ , 3:ncols]) , ...))
                names(res) <- levels(as.factor(ExpressionSet[ , 1]))
        }
        
        return(res)
}

