


PlotMeans <- function(ExpressionSet,Groups=list(NULL),legendName = NULL,...)
{
        
        is.ExpressionSet(ExpressionSet)
        
        if(any(sapply(Groups,is.null)))
                stop("Your Groups list does not store any items.")
        
        if(is.null(legendName))
                stop("Please specify the type of ExpressionSet you are working with: legendName = 'PS' or 'DS'.")
        
        ### getting the PS names available in the given expression set
        age_names <- as.character(names(table(ExpressionSet[ , 1])))
        
        # test whether all group elements are available in the age vector
        ra <- range(ExpressionSet[ , 1])
        if(!all(unlist(Groups) %in% as.numeric(age_names)))
                stop("There are items in your Group elements that are not available in the age column of your ExpressionSet.")
        
        ### getting the PS names available in the given expression set
        nPS <- length(age_names)
        nCols <- dim(ExpressionSet)[2]
        ### define and label the REmatrix that holds the rel. exp. profiles
        ### for the available PS
        MeanValsMatrix <- matrix(NA_real_,nPS,nCols-2)
        rownames(MeanValsMatrix) <- age_names
        colnames(MeanValsMatrix) <- names(ExpressionSet)[3:nCols]
        nGroups <- length(Groups)
        ### each PS class gets its corresponding color
        colos <- re.colors(nPS)
        nElements <- sapply(Groups,length)
        iterator <- 0
        ### for each phylostratum in the given dataset
        
        MeanValsMatrix <- age.apply(ExpressionSet, colMeans)
        
        ylim_min <- min(MeanValsMatrix) - (min(range(MeanValsMatrix))/6)
        ylim_max <- max(MeanValsMatrix) + (max(range(MeanValsMatrix))/5)
        
        
        # define arguments for different graphics functions
        plot.args <- c("lwd","col","lty","xlab","cex.lab","main")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        
        ### plot the rel. exp. levels in k different windows
        if(length(Groups) > 1)
                par(mfrow = rev(n2mfrow(nGroups)))
        
        for(j in 1:nGroups){
                
                if(j < 2){
                        
                        do.call(graphics::matplot,c(list(t(MeanValsMatrix[match(as.character(Groups[[j]]), rownames(MeanValsMatrix)) , ]),type = "l",
                                                         ylab = "Mean Expression Level",ylim = c(ylim_min,ylim_max), axes = FALSE, col = colos[match(as.character(Groups[[j]]), age_names)]),
                                                    dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(1,seq(1,nCols-2,1),names(ExpressionSet)[3:nCols]),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(2,seq(ylim_min,ylim_max,length.out = 5),format(seq(ylim_min,ylim_min,length.out = 5),digits = 2)),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::legend,c(list("top",legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
                                                        fill = colos[match(as.character(Groups[[j]]), age_names)],bty = "n",ncol = ceiling(nElements[j]/2)),
                                                   dots[!is.element(names(dots),c(axis.args,plot.args))]))
                        
                        iterator <- nElements[j]
                        
                }
                
                else{
                        
                        do.call(graphics::matplot,c(list(t(MeanValsMatrix[match(as.character(Groups[[j]]), rownames(MeanValsMatrix)) , ]),type = "l",
                                                         ylab = "Mean Expression Level",ylim = c(ylim_min,ylim_max), axes = FALSE, col = colos[match(as.character(Groups[[j]]), age_names)]),
                                                    dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(1,seq(1,nCols-2,1),names(ExpressionSet)[3:nCols]),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(2,seq(ylim_min,ylim_max,length.out = 5),format(seq(ylim_min,ylim_max,length.out = 5),digits = 2)),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::legend,c(list("top",legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
                                                        fill = colos[match(as.character(Groups[[j]]), age_names)],bty = "n",ncol = ceiling(nElements[j]/2)),
                                                   dots[!is.element(names(dots),c(axis.args,plot.args))]))
                        
                        iterator <- nElements[j]
                        
                }
        }
}


