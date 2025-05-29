#' @title Geometric Mean
#' @description 
#' This function computes the geometric mean of a numeric input vector \code{x}.
#' @param x a numeric vector for which geometric mean computations shall be performed.
#' @author Hajk-Georg Drost
#' @examples 
#' x <- 1:10
#' 
#' geom_mean(x)
#' 
#' @export
geom_mean <- function(x) {
    if (any(x < 0))
        NaN
    else if (any(x == 0))
        0
    else
        exp(mean(log(x)))
}



rowVars <- function(x, ...) {
    return(rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1))
}


