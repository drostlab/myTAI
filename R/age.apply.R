#' @title Age Category Specific apply Function
#' @description 
#' This function performs the split-apply-combine methodology on Phylostrata or Divergence Strata stored within the input PhyloExpressionSet.
#' 
#' This function is very useful to perform any phylostratum or divergence-stratum specific analysis.
#' 
#' @param phyex_set a standard PhyloExpressionSet object.
#' @param FUN a function to be performed on the corresponding expression matrix of each phylostratum or divergence-stratum.
#' @param ... additional arguments of FUN.
#' @param as.list a boolean value specifying whether the output format shall be a matrix or a list object.
#' @details This function uses the \code{\link{split}} function to subset the expression matrix into
#' phylostratum specific sub-matrices. Internally using \code{\link{lapply}}, any function can
#' be performed to the sub-matrices. The return value of this function is a numeric matrix storing
#' the return values by \code{FUN} for each phylostratum and each developmental stage s. 
#' Note that the input \code{FUN} must be an function that can be applied to a matrix (e.g., \code{\link{colMeans}}). 
#' In case you use an anonymous function you could use \code{function(x) apply(x , 2 , var)} as an example to compute the variance of each phylostratum and each
#' developmental stage s.
#' @return Either a numeric matrix storing the return values of the applied function for each age class
#' or a numeric list storing the return values of the applied function for each age class in a list.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{split}}, \code{\link{tapply}}, \code{\link{lapply}}
#' @examples
#'  
#' # source the example dataset
#' data(example_phyex_set)
#'  
#' # Example 1
#' # get the relative expression profiles for each phylostratum
#' age.apply(example_phyex_set, relative_expression)
#'
#' # this is analogous to 
#' rel_exp_matrix(example_phyex_set)

#' # Example 2
#' # compute the mean expression profiles for each phylostratum
#' age.apply(example_phyex_set, colMeans)
#'
#' # Example 3
#' # compute the variance profiles for each phylostratum
#' age.apply(example_phyex_set, function(x) apply(x , 2 , var))
#'
#' # Example 4
#' # compute the range for each phylostratum
#' # Note: in this case, the range() function returns 2 values for each phylostratum
#' # and each developmental stage, hence one should use the argument 'as.list = TRUE'
#' # to make sure that the results are returned properly 
#' age.apply(example_phyex_set, function(x) apply(x , 2 , range), as.list = TRUE)
#' 
#' @export
age.apply <- function(phyex_set, FUN, ..., as.list = FALSE) {
    f <- match.fun(FUN)
    
    # Split the expression matrix by strata
    strata_list <- split(seq_len(nrow(phyex_set@expression_collapsed)), phyex_set@strata)
    
    # Apply function to each stratum's subset
    res <- lapply(strata_list, function(idx) {
        f(phyex_set@expression_collapsed[idx, , drop = FALSE], ...)
    })
    
    if (!as.list) {
        res <- t(as.data.frame(res))
        rownames(res) <- names(strata_list)
    }
    else {
        names(res) <- names(strata_list)
    }
        
    return(res)
}
