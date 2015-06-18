#' @title Geometric Mean
#' @description 
#' This function computes the geometric mean of a numeric input vector \code{x}.
#' @param x a numeric vector for which geometric mean computations shall be performed.
#' @author Hajk-Georg Drost
#' @examples 
#' x <- 1:10
#' 
#' geom.mean(x)
#' 
#' @export

geom.mean <- function(x)
{
        if(is.numeric(x)){
                return(cpp_geom_mean(as.vector(x)))
        } else{
                stop("Please enter a numeric vector.")
        }
}