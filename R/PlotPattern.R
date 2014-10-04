PlotPattern <- function(ExpressionSet,TestStatistic="FlatLineTest",
                        modules = NULL,permutations=1000,lillie.test = FALSE,
                        digits.ylab=3,p.value=TRUE,shaded.area=FALSE,y.ticks=4,...)
{
        
        is.ExpressionSet(ExpressionSet)
        
        if(!(is.element(TestStatistic, c("FlatLineTest","ReductiveHourglassTest")))){
                stop("Please enter a correct string for the test statistic: 'FlatLineTest' or 'ReductiveHourglassTest'.")
        }
        
        if((TestStatistic == "ReductiveHourglassTest") & is.null(modules))
                stop("Please specify the modules for the ReductiveHourglassTest: modules = list(early = ..., mid = ..., late = ...).")
        
        nCols <- dim(ExpressionSet)[2]
        resList <- vector("list", length = 2)
        age <- vector(mode = "numeric",length = (nCols-2))
        age <- TAI(ExpressionSet)
        
        # computing the standard error of the TAI/TDI pattern using bootstrap analyses
        x <- min(age);
        y <- max(age);
        
        if(TestStatistic == "FlatLineTest"){
                
                resList <- FlatLineTest(ExpressionSet,permutations = permutations)
                
        }
        
        if(TestStatistic == "ReductiveHourglassTest"){
                
                if(lillie.test == TRUE){
                        resList <- ReductiveHourglassTest(ExpressionSet,modules = modules, 
                                                          permutations = permutations, lillie.test = TRUE)
                }
                
                if(lillie.test == FALSE){
                        resList <- ReductiveHourglassTest(ExpressionSet,modules = modules, 
                                                          permutations = permutations, lillie.test = FALSE)
                }
                
        }
        
        # get p-value and standard deviation values from the test statistic 
        pval <- resList$p.value
        sd_vals <- resList$std.dev
        #random_mean <- resList$mean
        max_sd <- max(sd_vals)
        
        # plot the age pattern surrounded by the standard deviation
        
        # adjust the plotting margins 
        aspect_ratio <- ((max_sd + y) - (x - max_sd)) / 20
        ylim.range <- range((x - max_sd - aspect_ratio),(y + max_sd + aspect_ratio))  
        
        # define arguments for different graphics functions
        plot.args <- c("type","lwd","col","cex.lab","main","xlab","ylab")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        #
        #   plot phylotranscriptomic age
        # 
        
        if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) & (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                do.call(graphics::plot,c(list(x = age,ylim = c(ylim.range[1],ylim.range[2]),axes = FALSE), 
                                         dots[!is.element(names(dots),c(axis.args,legend.args))]))
        }
        
        # default: xlab = "Ontogeny" and ylab = "Age Index"
        else{
                do.call(graphics::plot,c(list(x = age,ylim = c(ylim.range[1],ylim.range[2]),axes = FALSE, 
                                              xlab = "Ontogeny",ylab = "Age Index" ), dots[!is.element(names(dots),c(axis.args,legend.args))]))  
        }
        
        do.call(graphics::axis,c(list(side = 1,at = 1:(nCols-2),
                                      labels = names(ExpressionSet)[3:nCols]),dots[!is.element(names(dots),c(plot.args,legend.args))]))
        
        do.call(graphics::axis,c(list(side = 2,at = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab),
                                      labels = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab)), 
                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
        
        # age + std.err
        lines(age + sd_vals,lwd = 2,col = "darkgrey")
        # age - std.err
        lines(age - sd_vals,lwd = 2,col = "darkgrey")
        
        if(p.value == TRUE){
                do.call(graphics::legend,c(x = "topleft",bty = "n",legend = paste("p = ",format(pval,digits = 3),sep = ""),
                                           dots[!is.element(names(dots),c(plot.args,axis.args))]))
                
        }
        
        # draw a dotted line passing the global minimum value of age
        #abline(h=min(age),col="black",lty="dotted");
        
        if(shaded.area == TRUE){
                #abline(v=c(mid[1],mid[length(mid)]),col="black",lty="dotted")
                col2alpha <- function (color, alpha = 0.2){
                        rgbCode <- col2rgb(color)[,1]
                        rgb(rgbCode[1], rgbCode[2], rgbCode[3], 255 * alpha, maxColorValue = 255)
                }
                usr <- par('usr')
                rect(modules[[2]][1], usr[3], modules[[2]][length(modules[[2]])], usr[4], col = col2alpha("midnightblue",alpha = 0.2)) 
        }
}