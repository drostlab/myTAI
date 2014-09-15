

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


 PlotRE <- function(ExpressionSet,Groups=list(NULL),legendName=NULL,...)
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
   
   nPS <- length(age_names)
   nCols <- dim(ExpressionSet)[2]
   ### define and label the REmatrix that holds the rel. exp. profiles
   ### for the available PS
   REmatrix <- matrix(NA_real_,nPS,nCols-2)
   rownames(REmatrix) <- age_names
   colnames(REmatrix) <- names(ExpressionSet)[3:nCols]
   nGroups <- length(Groups)
   ### each PS class gets its corresponding color
   colos <- re.colors(nPS)
   
   # number of items in each Groups element
   nElements <- sapply(Groups,length)
   
   # compute the relative expression matrix
   REmatrix <- REMatrix(ExpressionSet)
   
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
      do.call(graphics::matplot,c(list(x = t(REmatrix[match(as.character(Groups[[j]]), rownames(REmatrix)) , ]),
               type = "l",axes = FALSE, ylim = c(0,1.4),col = colos[match(as.character(Groups[[j]]), age_names)], ylab = "Relative Expression"), 
               dots[!is.element(names(dots),c(axis.args,legend.args))]))
       
      do.call(graphics::axis,c(list(side = 1,at = seq(1,nCols-2,1), labels = names(ExpressionSet)[3:nCols]), 
                                dots[!is.element(names(dots),c(plot.args,legend.args))]))
      
      do.call(graphics::axis,c(list(side = 2,at = seq(0,1.4,0.2), labels = seq(0,1.4,0.2)), 
                               dots[!is.element(names(dots),c(plot.args,legend.args))]))
      
      do.call(graphics::legend,c(list(x = "top",legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
               fill = colos[match(as.character(Groups[[j]]), age_names)],
               bty = "n",ncol = ceiling(nElements[j] / 2)), 
               dots[!is.element(names(dots),c(axis.args,plot.args))]))
      
     }

     else{
       do.call(graphics::matplot,c(list(x = t(REmatrix[match(as.character(Groups[[j]]), rownames(REmatrix)) , ]),type = "l",
               axes = FALSE,ylim = c(0,1.4),col = colos[match(as.character(Groups[[j]]), age_names)],
               ylab = "Relative Expression"), 
               dots[!is.element(names(dots),c(legend.args,axis.args))]))
       
       do.call(graphics::axis,c(list(side = 1, at = seq(1,nCols-2,1),labels = names(ExpressionSet)[3:nCols]), 
                                dots[!is.element(names(dots),c(plot.args,legend.args))]))
       
       do.call(graphics::axis,c(list(side = 2, at = seq(0,1.4,0.2),labels = seq(0,1.4,0.2)), 
                                dots[!is.element(names(dots),c(plot.args,legend.args))]))
       
       do.call(graphics::legend,c(list(x = "top",legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
              fill = colos[match(as.character(Groups[[j]]), age_names)],bty = "n",ncol = ceiling(nElements[j] / 2)), 
              dots[!is.element(names(dots),c(axis.args,plot.args))])) 
     } 
   }
 }



PlotBarRE <- function(ExpressionSet,Groups=list(NULL),wLength=0.1,ratio=FALSE,...)
{
   
   is.ExpressionSet(ExpressionSet)
   
   if(any(sapply(Groups,is.null)))
     stop("Your Groups list does not store any items.")
   
   if(any(sapply(Groups,function(x) length(x) < 2)))
     stop("Each Group class needs to store at least two items.")
   
   ### getting the PS names available in the given expression set
   age_names <- as.character(names(table(ExpressionSet[ , 1])))
   
   # test whether all group elements are available in the age vector
   ra <- range(ExpressionSet[ , 1])
   if(!all(unlist(Groups) %in% as.numeric(age_names)))
     stop("There are items in your Group elements that are not available in the age column of your ExpressionSet.")
   
   PS.Table <- age_names
   nPS <- length(PS.Table)
   PS.Names <- names(PS.Table)
   nCols <- dim(ExpressionSet)[2]
   nGroups <- length(Groups)
   MeanREClassValues <- matrix(NA_real_,nGroups,nCols-2)
   StdErr.RE.ClassValues <- matrix(NA_real_,nGroups,nCols-2)
   pValues <- 1
   AnovaListValues <- 1
   barColors <- bar.colors(nGroups)
   ### compute the relative expression profiles for all
   ### given phylostrata
   REmatrix <- matrix(NA_real_,ncol = nPS,nrow = nCols-2)
   REmatrix <- REMatrix(ExpressionSet)
  
   ### compute the mean relative expression levels for each PS-Group
   ### as well as the Std.Error of the relative expression levels
   ### included in each PS-Group
   for(i in 1:nGroups){
     MeanREClassValues[i , ] <- colMeans(REmatrix[match(as.character(Groups[[i]]), rownames(REmatrix)) , ])
     StdErr.RE.ClassValues[i , ] <- apply(REmatrix[match(as.character(Groups[[i]]), rownames(REmatrix)) , ],2,std.error)
   }   
   
   if(nGroups == 2){
     for(j in 1:(nCols-2)){
       
        testForConstantValues <- try(kruskal.test(list(REmatrix[match(as.character(Groups[[1]]), rownames(REmatrix)) , j], REmatrix[match(as.character(Groups[[2]]), rownames(REmatrix)) , j])), silent = FALSE)
        
        if(is(testForConstantValues, "try-error")){
          warning("Something went wrong with the Kruskal-Wallis Rank Sum Test... the p-value has been set to p = 1.")
          pValues[j] <- 1
        } 
          
        else
          pValues[j] <- testForConstantValues$p.value
       
     }

      FoldChangeOfMeanREValues <- apply(MeanREClassValues,2,function(x){return((x[1]) / (x[2]))})
      tmpFoldChanges <- FoldChangeOfMeanREValues
      REFoldChangeOfMeanREValues <- (tmpFoldChanges - min(tmpFoldChanges)) / (max(tmpFoldChanges) - min(tmpFoldChanges))
      
   }

   if(nGroups > 2){
     for(s in 1:(nCols-2)){
       for(k in 1:nGroups){
         AnovaListValues[k] <- list(REmatrix[match(as.character(Groups[[k]]), rownames(REmatrix)) , s])
       }
       
       pValues[s] <- kruskal.test(AnovaListValues)$p.value
       
     }
   }
   
   pValNames <- rep("",nCols-2)
   pValNames[which(pValues <= 0.05)] <- "*"
   pValNames[which(pValues <= 0.005)] <- "**"
   pValNames[which(pValues <= 0.0005)] <- "***"

   
   # define arguments for different graphics functions
   barplot.args <- c("xlab","cex.lab","cex.axis","horiz","main","density","add","cex.names")
   text.args <- c("cex")
   dots <- list(...)
   ellipsis.names <- names(dots)
   
   #print(pValues)
   if(ratio == TRUE){
     
     REBarPlot <- do.call(graphics::barplot,c(list(MeanREClassValues,beside = TRUE,ylim = c(0,1),
                          names.arg = names(ExpressionSet)[3:nCols],
                          col = barColors,border = "white"),
                          dots[!is.element(names(dots),c(text.args))]))
     
     do.call(graphics::text,c(list(apply(REBarPlot,2,mean),0.95,labels = pValNames),
                              dots[!is.element(names(dots),c(barplot.args))]))
     
     arrows(x0 = REBarPlot,y0 = ifelse(MeanREClassValues > 0,MeanREClassValues, (1/999)),x1 = REBarPlot,
            y1 = ifelse((StdErr.RE.ClassValues) == 0,MeanREClassValues + (1/999),
            MeanREClassValues + StdErr.RE.ClassValues),code = 2, angle = 90, length = wLength)
     
     par(xpd = TRUE)
     legend("topleft",inset = c(+0.2,0),legend = paste("Group ",1:length(Groups),sep = ""),
            fill = barColors,bty = "n",cex = 1.3,ncol = ceiling(nGroups/2))
     par(xpd = FALSE)
     if(nGroups == 2){
       lines(colMeans(REBarPlot), REFoldChangeOfMeanREValues,lty = 2,lwd = 5,col = "darkblue")
       do.call(graphics::axis,c(list(4,seq(0,1,0.2),format(seq(0,max(FoldChangeOfMeanREValues),length.out = 6),digits = 2)),
                                dots[!is.element(names(dots),c(text.args))]))
     }
   }

   if(ratio == FALSE){
     REBarPlot <- do.call(graphics::barplot,c(list(MeanREClassValues,beside = TRUE,ylim = c(0,1.6),
                          names.arg = names(ExpressionSet)[3:nCols],col = barColors,border = "white"),
                          dots[!is.element(names(dots),c(text.args))]))
     
     do.call(graphics::text,c(list(colMeans(REBarPlot),1.15,labels = pValNames),
                              dots[!is.element(names(dots),c(barplot.args))]))
     
     arrows(x0 = REBarPlot,y0 = ifelse(MeanREClassValues > 0,MeanREClassValues, (1/999)),x1 = REBarPlot,
            y1 = ifelse((StdErr.RE.ClassValues) == 0,MeanREClassValues + (1/999),
            MeanREClassValues + StdErr.RE.ClassValues),code = 2, angle = 90, length = wLength)
     
     legend("topleft",legend = paste("Group ",1:length(Groups),sep = ""),fill = barColors,bty = "n",ncol = ceiling(nGroups/2))
   }
   
}


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


### Function to plot the correlation between phylostrata and divergence-strata
### methods = "pearson" or "kendall" or "spearman"
PlotCorrelation <- function(PhyloExpressionSet,DivergenceExpressionSet,method="pearson",linearModel=F,main.text="",...)
{
  
  is.ExpressionSet(PhyloExpressionSet)
  is.ExpressionSet(DivergenceExpressionSet)
  
  colnames(PhyloExpressionSet)[2] <- "GeneID"
  colnames(DivergenceExpressionSet)[2] <- "GeneID"
  
  PS_DS.Subset <- merge(PhyloExpressionSet[ , 1:2], DivergenceExpressionSet[ , 1:2],by = "GeneID")
  
  CorrelationCoefficient <- cor(PS_DS.Subset[ , 2],PS_DS.Subset[ , 3],method = method)
  CorrCoeffasCharacter <- as.character(round(CorrelationCoefficient,3))
  
  nrows <- dim(PS_DS.Subset)[1]
  PS <- vector(mode = "numeric", length = nrows)
  DS <- vector(mode = "numeric", length = nrows)
  
  PS <- jitter(PS_DS.Subset[ , 2],1.5)
  DS <- jitter(PS_DS.Subset[ , 3],1.5)
  
  plot(PS,DS,main = paste(main.text,method," = ",CorrCoeffasCharacter,sep = ""),...)
  
  if(linearModel == TRUE)
    abline(lm(PS_DS.Subset[ , 3]~PS_DS.Subset[ , 2]),lwd = 5,col = "red")

}



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


