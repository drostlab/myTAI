### Function to plot the PS/DS distribution of a given PS/DS Vector
PlotDistribution <- function(PhyloExpressionSet,plotText=TRUE,as.ratio=FALSE,...)
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

