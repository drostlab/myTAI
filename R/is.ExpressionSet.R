#' @title Test ExpressionSet Standard
#' @description This function tests whether a given ExpressionSet follows the pre-defined PhyloExpressionSet or DivergenceExpressionSet standard.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet that shall be tested for format validity.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # read example PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#' 
#' is.ExpressionSet(PhyloExpressionSetExample)
#' 
#' @export
is.ExpressionSet <- function(ExpressionSet){
        
        ncols <- dim(ExpressionSet)[2]
        
        d.f_bool <- is.data.frame(ExpressionSet) | tibble::is.tibble(ExpressionSet)
        
        if (!d.f_bool)
                stop("ExpressionSet is not a data.frame or tibble.", call. = FALSE)
        
        age.vector_bool <- is.numeric(unlist(dplyr::select(ExpressionSet, 1)))
        
        if (!age.vector_bool)
                stop("The first column of the ExpressionSet needs to store numeric values.", call. = FALSE)
        
        gene.vector_bool <-
            ifelse(is.factor(unlist(dplyr::select(ExpressionSet, 2))),
                   is.character(levels(unlist(
                       dplyr::select(ExpressionSet, 2)
                   ))),
                   is.character(unlist(dplyr::select(ExpressionSet, 2))))
        
        if (!gene.vector_bool)
                stop("The second column of the ExpressionSet needs to store character values (e.g. gene ids).", call. = FALSE)
        
        expression.matrix_bool <- all(apply(dplyr::select(ExpressionSet, 3:ncol(ExpressionSet)), 2 , is.numeric))
        if (!expression.matrix_bool)
                stop("Columns 3 - ",ncol(ExpressionSet), " need to store numeric values (e.g. expression levels).", call. = FALSE)
        
        any.NA.values_bool <- !any(is.na(ExpressionSet))
        
        if (!any.NA.values_bool)
                stop("The ExpressionSet cannot contain 'NA' values. Please use na.omit() to remove 'NA' values from the input ExpressionSet.", call. = FALSE)
        
        if (all(c(d.f_bool,age.vector_bool,gene.vector_bool,expression.matrix_bool,any.NA.values_bool))) {
                return(TRUE)
        }
        
        else{
                stop("The present input object does not fulfill the ExpressionSet standard.", call. = FALSE)
        }
}
