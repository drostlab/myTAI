
quant <- function(ExpressionMatrix,quantile = 0.9)
{
  nCols <- dim(ExpressionMatrix)[2]
  threshold <- vector(mode = "numeric",length = nCols)
  
  for(i in 1:nCols){
    threshold[i] <- as.numeric(quantile(ExpressionMatrix[ , i],probs = quantile))
  }
  return (threshold)
}

std.error <- function(x)
{
  if(is.numeric(x)){
    return(cpp_std_error(as.vector(x)))
  }
  else{
    stop("Please enter a numeric vector.")
  }
}

geom.mean <- function(x)
{
  
  if(is.numeric(x)){
    return(cpp_geom_mean(as.vector(x)))
  }
  
  else{
    stop("Please enter a numeric vector.")
  }
  
}

harmonic.mean <- function(x)
{
  
  N <- length(x)
  ratio <- 1/x
  harm <- 1 / ((1/N) * sum(ratio))
  return(harm)
  
}


re.colors <- function(n)
{
  
  colos <- c("black","red","green","brown","darkmagenta",
             "blue","darkred","darkblue","darkgreen", "orange",
             "azure4","gold4","greenyellow","hotpink4",
             "mediumorchid3","mediumorchid3","peachpuff4",
             "hotpink","lightgoldenrod", "peru", "slateblue3", "yellow4", "yellowgreen")
  
  return(colos[1:n])
  
}


#' @title A function to get a vector of length n storing a palette of colors for multiple bars in barplots.
#' @description A nice color palette for barplots with several bars.
#' @param n the number of colors to be in the palette. 
#' @return a character vector containing different color names that can be used for barplots.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{palette}}, \code{\link{re.colors}}
#' @examples \dontrun{
#' 
#' # get 5 different colors for 5 different bars
#' barplot_colors <- bar.colors(5)
#' 
#' 
#' }
bar.colors <- function(n)
{
  
  colos <- c("black","gray76","gray43","navy","lightskyblue3","palegreen4","seagreen1","lavenderblush1")
  return(colos[1:n])
  
}


# tests whether a given input data.frame 
# fulfills the ExpressionSet standard
is.ExpressionSet <- function(ExpressionSet){
  
  ncols <- dim(ExpressionSet)[2]
  
  d.f_bool <- is.data.frame(ExpressionSet)
  age.vector_bool <- is.numeric(ExpressionSet[ , 1])
  gene.vector_bool <- ifelse(is.factor(ExpressionSet[ , 2]),is.character(levels(ExpressionSet[ , 2])),is.character(ExpressionSet[ , 2]))
  expression.matrix_bool <- all(apply(ExpressionSet[ , 3:ncols] , 2 , is.numeric))
  any.NA.values_bool <- !any(apply(ExpressionSet , 2 , is.na))
  
  
  if(d.f_bool & age.vector_bool & gene.vector_bool & expression.matrix_bool & any.NA.values_bool){
    return(TRUE)
  }
  
  else{
    stop("The present input object does not fulfill the ExpressionSet standard.")
  }
}




