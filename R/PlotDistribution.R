#' @title Plot the frequency distribution of Phylostrata or Divergence-Strata
#' @description This function plots the frequency distribution of genes within the 
#' corresponding \emph{phylostratigraphic map} or \emph{divergence map} and can be used to fastly visualize the PS or DS distribution of a given phylostratum vector or divergence-stratum vector.
#' @param PhyloExpressionSet a standard PhyloExpressionSet object.
#' @param plotText a boolean value specifying whether the total number of genes 
#' belonging to a corresponding phylostratum or divergence-stratum class
#' shall be plotted above each bar of the barplot.
#' @param as.ratio a boolean value specifying whether the relative frequencies
#' instead of absolute frequencies shall be plotted.
#' @param \dots default plot parameters.
#' @details 
#' The frequency distribution of all genes or a subset of genes might be of interest for subsequent analyses.
#'
#' For Example:
#'      
#' Filtering genes using gene cluster algorithms can result in different groups (classes) of genes.
#' For each gene group the phylostratum or divergence-stratum distribution can be visualized using this function
#' and can be compared between different groups.
#'
#' This analysis allows to compare different gene expression profiles (or gene groups in general) based
#' on their evolutionary origins or evolutionary relationships.
#' @return a barplot showing the phylostratum distribution or 
#' divergence-stratum distribution of a given numeric vector containing PS or DS values.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # load PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # plot the phylostratum distribution of a PhyloExpressionSet
#' PlotDistribution(PhyloExpressionSetExample, plotText = TRUE)
#'
#' # plot the relative frequency distribution of a PhyloExpressionSet
#' PlotDistribution(PhyloExpressionSetExample, plotText = TRUE, as.ratio = TRUE)
#'
#'
#' # a example for visualizing the PS distribution for a subset of genes
#' PlotDistribution(PhyloExpressionSetExample[sample(20000,5000) , ],
#'                  plotText = TRUE, as.ratio = TRUE)
#'
#' @export

PlotDistribution <- function(PhyloExpressionSet,plotText = TRUE,as.ratio = FALSE,...)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        
        AgeVector <- PhyloExpressionSet[ , 1]
        minPSVal <- min(AgeVector)
        maxPSVal <- max(AgeVector)
        
        nPS <- length(minPSVal:maxPSVal)
        
        # define arguments for different graphics functions
        barplot.args <- c("xlab","cex.lab","cex.axis","horiz","main","density","add","cex.names")
        text.args <- c("cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        if(as.ratio == FALSE){
                PStable <- table(factor(AgeVector,levels = 1:nPS))
                yLim <- ceiling((max(PStable) + (max(PStable)/2)))
                BarPlot <- do.call(graphics::barplot,c(list(PStable,ylab = "# Genes",ylim = c(0,yLim),border = "white"),
                                                       dots[!is.element(names(dots),c(text.args))]))
                
                if(plotText == TRUE){
                        do.call(graphics::text,c(list(x = BarPlot,y = (PStable + (PStable/3)),labels = ifelse(PStable > 0,PStable,"")),
                                                 dots[!is.element(names(dots),c(barplot.args))]))
                }
        }
        
        if(as.ratio == TRUE){
                
                PStable <- table(factor(AgeVector,levels = 1:nPS))
                relFreq <- as.numeric(format(PStable / sum(PStable),digits = 1))
                names(relFreq) <- names(PStable)
                yLim <- c(0,1)
                BarPlot <- do.call(graphics::barplot,c(list(relFreq,ylab = "rel. frequency of genes",ylim = yLim,border = "white"),
                                                       dots[!is.element(names(dots),c(text.args))]))
                
                if(plotText == TRUE){
                        do.call(graphics::text,c(list(x = BarPlot,y = (relFreq + .1),labels = ifelse(relFreq > 0,relFreq,"")),
                                                 dots[!is.element(names(dots),c(barplot.args))]))
                }
        }
}

