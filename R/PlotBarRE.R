#' @title Plot the Mean Relative Expression Levels as Barplot
#' @description This function takes a PhyloExpressionSet or DivergenceExpressionSet object as input and computes
#' for two or more defined phylostratum (divergence stratum) classes the statistical significance of
#' the differences of mean relative expression of these two (or more) corresponding phylostratum (divergence stratum) classes.
#' As test-statistic, the function performs a nonparametric \code{\link{kruskal.test}}
#'  based on relative expression values stored within each defined phylostratum class.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the phylostrata or divergence-strata that correspond to the same phylostratum class or divergence class.
#'  For ex. evolutionary old phylostrata: PS1-3 (Class 1) and evolutionary young phylostrata: PS4-12 (Class 2). 
#'  In this case, the \code{Groups} list could be assigned as, \code{Groups} = list(c(1:3), c(4:12)). 
#'  It is also possible to define more than 2 groups of evolutionary ages.
#'  For ex. \code{Groups} = list(c(1:3),c(4:8),c(9:12)) would perform a \code{\link{kruskal.test}} to determine the statistical significance
#'  of the evolutionary classes PS1-3, PS4-6, and PS9-12 based on their corresponding mean relative expression levels. 
#' @param wLength a numeric value defining the whiskers length above the bars. In case there are numerous different phylostratum classes
#'  a smaller \code{wLength} parameter should be used for better visualizations.
#' @param ratio a boolean value specifying whether the bars in the barplot represent the 
#'  mean relative expression level of phylostrata belonging to the same phylostratum class.
#'  In case \code{ratio} = TRUE, the ratio of the mean relative expression level of 
#'  the two phylostrata classes is plotted as lines within the barplot. This parameter can only be used for 2 class comparisons.
#' @param p.adjust.method correction method to adjust p-values for multiple comparisons (see \code{\link{p.adjust}} for possible methods).
#'   E.g., \code{p.adjust.method = "BH"} (Benjamini & Hochberg (1995)) or \code{p.adjust.method = "bonferroni"} (Bonferroni correction).
#' @param \dots default graphics parameters. 
#' @return A barplot comparing Phylostratum-Classes by its mean relative expression levels.
#' Significant stages are marked by '*' referring to statistically significant differences:
#'        
#' (1) '*'   = P-Value <= 0.05 
#'
#' (2) '**'  = P-Value <= 0.005  
#'
#' (3) '***' = P-Value <= 0.0005  
#' 
#' @details 
#' In case a large number of developmental stages is included in the input \code{ExpressionSet},
#'  p-values returned by \code{PlotBarRE} should be adjusted for multiple comparisons which can be done
#'  by specifying the \code{p.adjust.method} argument.
#'  
#' @references 
#' 
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' 
#' Domazet-Loso T. and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.
#' 
#' Myles Hollander and Douglas A. Wolfe (1973), Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 115-120.
#' 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{RE}}, \code{\link{REMatrix}},\code{\link{PlotRE}}, \code{\link{kruskal.test}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # example PhyloExpressionSet
#' PlotBarRE(ExpressionSet = PhyloExpressionSetExample,
#'           Groups        = list(c(1:3), c(4:12)))
#'
#'
#' # example DivergenceExpressionSet
#' PlotBarRE(ExpressionSet = DivergenceExpressionSetExample,
#'           Groups        = list(c(1:5), c(6:10)))
#'
#'
#' # Perform PlotBarRE() with p-value adjustment method Benjamini & Hochberg (1995)
#' PlotBarRE(ExpressionSet   = PhyloExpressionSetExample,
#'           Groups          = list(c(1:3), c(4:12)),
#'           p.adjust.method = "BH")
#'        
#'              
#' # Example: plot ratio
#' # the ratio curve visualizes the ratio between bar 1 / bar 2
#' # the z - axis shows the corresponding ratio value of bar 1 / bar 2
#' PlotBarRE(ExpressionSet = PhyloExpressionSetExample,
#'           Groups        = list(c(1:3), c(4:12)),
#'           ratio         = TRUE)
#' 
#' @export

PlotBarRE <- function(ExpressionSet,
                      Groups          = NULL,
                      wLength         = 0.1,
                      ratio           = FALSE,
                      p.adjust.method = NULL, ...)
{
        
        is.ExpressionSet(ExpressionSet)
        
        if(is.null(Groups))
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
        REmatrix <- age.apply(ExpressionSet = ExpressionSet, RE)
        
        ### compute the mean relative expression levels for each PS-Group
        ### as well as the Std.Error of the relative expression levels
        ### included in each PS-Group
        for(i in 1:nGroups){
                MeanREClassValues[i , ] <- colMeans(REmatrix[match(as.character(Groups[[i]]),
                                                                   rownames(REmatrix)) , ])
                
                StdErr.RE.ClassValues[i , ] <- apply(REmatrix[match(as.character(Groups[[i]]),
                                                                    rownames(REmatrix)) , ],2,std.error)
        }   
        
        if(nGroups == 2){
                for(j in 1:(nCols-2)){
                        
                        testForConstantValues <- try(stats::kruskal.test(list(REmatrix[match(as.character(Groups[[1]]),
                                                                                      rownames(REmatrix)) , j],
                                                                       REmatrix[match(as.character(Groups[[2]]),
                                                                                      rownames(REmatrix)) , j])),
                                                     silent = FALSE)
                        
                        if(methods::is(testForConstantValues, "try-error")){
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
                        
                        pValues[s] <- stats::kruskal.test(AnovaListValues)$p.value
                        
                }
        }
        
        if(!is.null(p.adjust.method)){
                
                pValues <- stats::p.adjust(pValues, method = p.adjust.method, n = nPS)
        }
        
        pValNames <- rep("",nCols-2)
        pValNames[which(pValues <= 0.05)] <- "*"
        pValNames[which(pValues <= 0.005)] <- "**"
        pValNames[which(pValues <= 0.0005)] <- "***"
        pValNames[which(is.na(pValNames))] <- ""
        
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
                
                suppressWarnings(graphics::arrows(x0 = REBarPlot,y0 = ifelse(MeanREClassValues > 0,MeanREClassValues, (1/999)),x1 = REBarPlot,
                       y1 = ifelse((StdErr.RE.ClassValues) == 0,MeanREClassValues + (1/999),
                                   MeanREClassValues + StdErr.RE.ClassValues),code = 2, angle = 90, length = wLength))
                
                graphics::par(xpd = TRUE)
                graphics::legend("topleft",
                                 inset  = c(+0.2,0),
                                 legend = paste("Group ",1:length(Groups),sep = ""),
                                 fill   = barColors,
                                 bty    = "n",
                                 cex    = 1.3,
                                 ncol   = ceiling(nGroups/2))
                
                graphics::par(xpd = FALSE)
                if(nGroups == 2){
                        graphics::lines(colMeans(REBarPlot), REFoldChangeOfMeanREValues,lty = 2,lwd = 5,col = "darkblue")
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
                
                suppressWarnings(graphics::arrows(x0 = REBarPlot,y0 = ifelse(MeanREClassValues > 0,MeanREClassValues, (1/999)),x1 = REBarPlot,
                       y1 = ifelse((StdErr.RE.ClassValues) == 0,MeanREClassValues + (1/999),
                                   MeanREClassValues + StdErr.RE.ClassValues),code = 2, angle = 90, length = wLength))
                
                graphics::legend("topleft",legend = paste("Group ",1:length(Groups),sep = ""),fill = barColors,bty = "n",ncol = ceiling(nGroups/2))
        }
        
}

