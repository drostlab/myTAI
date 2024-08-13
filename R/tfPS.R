#' @title Transform age values
#' @description 
#' This function performs transformation of phylostratum values.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet object.
#' @param transform a character vector of any valid function that transforms PS values.
#' #' Possible values can be:
#' \itemize{
#' \item \code{transform} = \code{"qr"} (or \code{"quantilerank"}) : 
#' quantile rank transformation analogous to Julia function \code{StatsBase.quantilerank} 
#' using \code{method = :tied}.
#' }
#' @details This function transforms the phylostratum assignment. 
#' The return value of this function is a PhyloExpressionSet object with 
#' transformed phylostratum \code{tfPhylostratum} as the first column, satisfying
#' \code{\link{is.ExpressionSet}}. Note that the input \code{transform} must be an
#' available function, currently limited to only \code{"qr"} (or \code{"quantilerank"}).
#' @return A transformed PhyloExpressionSet object with \code{tfPhylostratum} as the first column.
#' @author Jaruwatana Sodai Lotharukpong
#' @seealso \code{\link{tf}}
#' @examples
#'  # source the example dataset
#'  data(PhyloExpressionSetExample)
#'  
#'  # get the relative expression profiles for each phylostratum
#'  tfPS(PhyloExpressionSetExample, transform = "qr")
#'
#' @export
tfPS <- function(ExpressionSet, transform){
  ExpressionSet <- as.data.frame(ExpressionSet)
  is.ExpressionSet(ExpressionSet)
  
  PS_vector <- ExpressionSet[,1]
  
  if(transform %in% c("qr", "quantilerank")){
    ranks <- base::rank(PS_vector, ties.method = "average")
    
    # -1 ensures that the ranks are scaled between 0 and 1
    # (when I checked using `data <- rnorm(100)`)
    
    message("Performing quantile rank transformation")
    
    tfPhylostratum <- (ranks - 1) / (length(PS_vector) - 1)
  }else{
    stop("Choose the following transformation functions: \"qr\"")
  }
  
  res <- base::cbind(
    base::as.data.frame(tfPhylostratum), 
    ExpressionSet[ , -1])
  
  return(res)
}
