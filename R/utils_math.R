#' @title Geometric Mean
#' @description 
#' This function computes the geometric mean of a numeric input vector \code{x}.
#' @param x a numeric vector for which geometric mean computations shall be performed.
#' @author Hajk-Georg Drost
#' @examplesIf FALSE 
#' x <- 1:10
#' 
#' geom_mean(x)
#' 
#' @keywords internal
geom_mean <- function(x) {
    if (any(x < 0))
        NaN
    else if (any(x == 0))
        0
    else
        exp(mean(log(x)))
}



#' @title Row-wise Variance Calculation
#' @description Calculate the variance for each row of a matrix or data frame.
#' 
#' @param x A numeric matrix or data frame
#' @param ... Additional arguments passed to rowSums and rowMeans
#' 
#' @return A numeric vector containing the variance for each row
#' 
#' @details
#' This function computes the sample variance for each row using the formula:
#' var = sum((x - mean(x))^2) / (n - 1)
#' 
#' @examples
#' # Calculate row variances for a matrix
#' # mat <- matrix(1:12, nrow = 3)
#' # row_vars <- rowVars(mat)
#' 
#' @keywords internal
rowVars <- function(x, ...) {
    return(rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1))
}

#' @title Calculate Quantile Ranks
#' @description Calculate quantile ranks for a numeric vector, handling ties using average method.
#' @param x numeric vector for which to calculate quantile ranks
#' @return A numeric vector of quantile ranks between 0 and 1
#' @examples
#' # Calculate quantile ranks for a vector
#' ranks <- quantile_rank(c(1, 2, 3, 4, 5))
#' @export
quantile_rank <- function(x) {
    ranks <- base::rank(x, ties.method = "average")
    (ranks - 0.5) / length(x)
}

#' @title Format P-Value for Scientific Notation
#' @description Format p-values in scientific notation for plot annotations.
#' 
#' @param p Numeric p-value
#' @param sci_thresh Numeric threshold for using scientific notation (number of decimal places)
#' 
#' @return Expression object for use in plot annotations
#' 
#' @details
#' This function formats p-values in scientific notation using the format
#' "p = a × 10^b" which is suitable for ggplot2 annotations and maintains
#' proper mathematical formatting.
#' 
#' @examples
#' # Format p-value for plotting
#' expr <- exp_p(0.001)
#' # Use in ggplot annotation:
#' # annotate("text", x = 1, y = 1, label = expr, parse = TRUE)
#' 
#' @export
exp_p <- function(p, sci_thresh = 4) {
    if (!is.finite(p)) return("italic(p) == NA")
    if (p == 0) return("italic(p) == 0")
    parts <- strsplit(formatC(p, format = "e", digits = 3), "e")[[1]]
    base <- as.numeric(parts[1])
    expn <- as.integer(parts[2])
    if (abs(expn) >= sci_thresh) {
        paste0("italic(p) == ", base, " %*% 10^", expn)
    } else {
        paste0("italic(p) == ", signif(p, 3))
    }
}