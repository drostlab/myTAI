#' @title Transform Phylostratum Values
#' @description 
#' This function performs transformation of phylostratum values.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet object.
#' @param transform a character vector of any valid function that transforms PS values.
#' Possible values can be:
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
#' @return a standard PhloExpressionSet object storing transformed Phylostratum levels.
#' @author Jaruwatana Sodai Lotharukpong and Lukas Maischak
#' @seealso \code{\link{tf}}
#' @examples
#' # source the example dataset
#' data(PhyloExpressionSetExample)
#'  
#' # get the relative expression profiles for each phylostratum
#' tfPES <- tfPS(PhyloExpressionSetExample, transform = "qr")
#' head(tfPES)
#'
#' @export

tfPS <- function(ExpressionSet,transform){
  ExpressionSet <- as.data.frame(ExpressionSet)
  is.ExpressionSet(ExpressionSet)
  
  PS_vector <- ExpressionSet[,1]
  
  if(transform %in% c("qr", "quantilerank")){
    ranks <- base::rank(PS_vector, ties.method = "average")
    # using the tied method
    tfPhylostratum <- (ranks - 0.5) / length(PS_vector)
  }else{
    stop("Choose the following transformation functions: \"qr\"")
  }
  
  # We are only using "tied". 
  # Below are other methods for quantile rank transform,
  # from Julia's StatsBase.quantilerank
  
  # n <- length(PS_vector)
  # if (method == "inc") {tfPhylostratum <- (ranks - 1) / (n - 1)} 
  # else if (method == "exc") {tfPhylostratum <- ranks / (n + 1)}
  # else if (method == "compete") {tfPhylostratum <- (ranks - 1) / (n - 1)} 
  # else if (method == "tied") {tfPhylostratum <- (ranks - 0.5) / n}
  # else if (method == "strict") {tfPhylostratum <- (ranks - 1) / n}
  # else if (method == "weak") {tfPhylostratum <- ranks / n} 
  # else {
  #   stop("method=", method, " is not valid. Pass 'inc', 'exc', 'compete', 'tied', 'strict', or 'weak'.")
  # }
  
  res <- base::cbind(
    base::as.data.frame(tfPhylostratum), 
    ExpressionSet[ , -1])
  
  return(res)
}
