
pow <- function(x,power)
{
  return(x^power)
}


CollapseFromTo <- function(x,from,to,FUN, ...){  
        f <- match.fun(FUN); return(apply(x[ , from:to], 1 , f, ...))
        }


quant <- function(ExpressionMatrix,quantile = 0.9)
{
        nCols <- dim(ExpressionMatrix)[2]
        threshold <- vector(mode = "numeric",length = nCols)
        
        for(i in 1:nCols){
                threshold[i] <- as.numeric(stats::quantile(ExpressionMatrix[ , i],probs = quantile))
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


re.colors <- function(n)
{
        
#         colos <- c("black","red","green","brown","darkmagenta",
#                    "blue","darkred","darkblue","darkgreen", "orange",
#                    "azure4","gold4","greenyellow","hotpink4",
#                    "mediumorchid3","mediumorchid3","peachpuff4",
#                    "hotpink","lightgoldenrod", "peru", "slateblue3", "yellow4", "yellowgreen")
        
        return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(n) )
        
}


#' @title Color palette for barplots
#' @description A nice color palette for barplots with several bars.
#' @param n the number of colors to be in the palette. 
#' @return a character vector containing different color names that can be used for barplots.
#' @details This function can be used to select colors for bar plots. 
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # get 5 different colors for 5 different bars
#' barplot_colors <- bar.colors(5)
#' @export
bar.colors <- function(n)
{
        
        colos <- c("black","gray76","gray43","navy","lightskyblue3","palegreen4","seagreen1","lavenderblush1")
        return(colos[1:n])
        
}


get.contingency.tbl <- function(x, index){
        
        contig.tbl <- matrix(NA_real_,ncol = 2, nrow = 2)
        contig.tbl <- rbind(c(x[index, 1],x[index, 2]),
                 c(sum(x[-index, 1]),sum(x[-index, 2])))
        colnames(contig.tbl) <- c("BG","TestSet")
        rownames(contig.tbl) <- c(paste0("PS",index),paste0("-PS",index))
        return(contig.tbl)

}


# Function to determine the row indices
# GetColumnIndexFromTo(nrep = c(2,3,2))
# nrep is a variable number of columns (replicates per stage)
GetColumnIndexFromTo <- function(nrep){
        
        IndexOne <- vector("numeric")
        IndexTwo <- vector("numeric")
        
        IndexOne[1] <- 1
        IndexTwo[1] <- nrep[1]
        
        for (i in 2:length(nrep)){
                
                IndexOne[i] <-  IndexTwo[i - 1] + 1
                IndexTwo[i] <- IndexOne[i] - 1 + nrep[i]
        }
        
        res <- cbind(IndexOne,IndexTwo)
        colnames(res) <- c("From","To")
        rownames(res) <- paste0("CollapsedStage",1:length(nrep))
        return(res)
}



