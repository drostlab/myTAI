#' @title Harmonic Mean
#' @description 
#' This function computes the harmonic mean of a numeric input vector \code{x}.
#' @param x a numeric vector for which harmonic mean computations shall be performed.
#' @author Hajk-Georg Drost
#' @examples 
#' x <- 1:10
#' 
#' harm.mean(x)
#' 
#' @export

harm.mean <- function(x)
{
        if(is.numeric(x)){
                return(cpp_harmonic_mean(as.vector(x)))
        } else{
                stop("Please enter a numeric vector.")
        }
}