#' @title Plot the Mean Expression Profiles of Phylostrata or Divergence Strata
#' @description This function computes for each phylostratum or divergence-stratum the corresponding mean expression profile
#' and plots the profiles in N different windows corresponding to the given Phylostratum-Classes
#' that shall be compared.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the phylostrata or divergence-strata that correspond to the same phylostratum class or divergence class.
#' For ex. evolutionary old phylostrata: PS1-3 (Class 1) and evolutionary young phylostrata: PS4-12 (Class 2). 
#' In this case, the list could be assigned as, \code{Groups} = list(c(1:3), c(4:12)). 
#' It is also possible to define more than 2 groups of evolutionary ages.
#' For ex. \code{Groups} = list(c(1:3),c(4:8),c(9:12)).
#' @param legendName a character string specifying whether "PS" or "DS" are used.
#' @param colors colors for mean expression profiles. Default: \code{colors = NULL}, hence default colours are used.
#' @param \dots default graphics parameters.
#' @details 
#' 
#' This plot may be useful to compare the absolute mean expression        
#' levels of each phylostratum or divergence-stratum class.
#'
#'In different developmental processes different phylostratum or divergence-stratum
#' classes might be more expressed than others, hence contributing more to the overall
#' phylotranscriptomics pattern (\code{\link{TAI}} or \code{\link{TDI}}).
#' This plot can help to identify the phylostratum or divergence-stratum classes 
#' that contributes most to the overall transcriptome of the given developmental process.
#' @return a plot showing mean expression profiles of 
#' phylostrata or divergence-strata corresponding to the same phylostratum
#' class or divergence class.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{PlotBarRE}}, \code{\link{RE}}, \code{\link{REMatrix}}, \code{\link{PlotRE}}
#' @examples
#' 
#' # load PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # load PhyloExpressionSet
#' data(DivergenceExpressionSetExample)
#'
#' # plot evolutionary old PS (PS1-3) vs evolutionary young PS (PS4-12)
#' PlotMeans(PhyloExpressionSetExample,Groups = list(c(1:3), c(4:12)), 
#'           legendName = "PS", lty = 1, lwd = 5)
#'
#' # plot conserved DS (DS1-5) vs divergent DS (PS6-10)
#' # NOTE: DS are always defined in the range 1, 2, ... , 10.
#' # Hence, make sure that your groups are within this range!
#' PlotMeans(DivergenceExpressionSetExample,Groups = list(c(1:5), c(6:10)), 
#'           legendName = "DS", lty = 1, lwd = 5)
#'
#'
#'# adding custom colors for relative expression levels:
#' # -> colors should be ordered by PS/DS starting with PS1,2,3...
#' PlotMeans(PhyloExpressionSetExample,
#'        Groups     = list(c(1:3), c(4:12)), 
#'        legendName = "PS",
#'        colors     = c("black","red","green","brown","darkmagenta",
#'        "blue","darkred","darkblue","darkgreen", "orange",
#'        "azure4","gold4"), 
#'        lty        = 1, 
#'        lwd        = 5)
#' @export 

PlotMeans <- function(ExpressionSet,
                      Groups     = NULL,
                      legendName = NULL,
                      colors     = NULL, ...)
{
        
        is.ExpressionSet(ExpressionSet)
        
        if(is.null(Groups))
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
        if(!is.null(colors)){
                colos <- colors
        } else {
                colos <- re.colors(nPS)
        }
        
        nElements <- sapply(Groups,length)
        iterator <- 0
        ### for each phylostratum in the given dataset
        
        MeanValsMatrix <- age.apply(ExpressionSet, colMeans)
        
        ylim_min <- min(MeanValsMatrix) - (min(range(MeanValsMatrix))/6)
        ylim_max <- max(MeanValsMatrix) + (max(range(MeanValsMatrix))/5)
        
        
        # define arguments for different graphics functions
        plot.args <- c("lwd","col","lty","xlab","cex.lab","main","cex.main")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        
        ### plot the rel. exp. levels in k different windows
        if(length(Groups) > 1)
                graphics::par(mfrow = rev(grDevices::n2mfrow(nGroups)))
        
        for(j in 1:nGroups){
                
                if(j < 2){
                        
                        do.call(graphics::matplot,c(list(t(MeanValsMatrix[match(as.character(Groups[[j]]), rownames(MeanValsMatrix)) , ]),type = "l",
                                                         ylab = "Mean Expression Level",ylim = c(ylim_min,ylim_max), axes = FALSE, col = colos[match(as.character(Groups[[j]]), age_names)]),
                                                    dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(1,seq(1,nCols-2,1),names(ExpressionSet)[3:nCols]),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::axis,c(list(2,seq(ylim_min,ylim_max,length.out = 5),format(seq(ylim_min,ylim_max,length.out = 5),digits = 6)),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::legend,c(list("top", legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
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
                        
                        do.call(graphics::axis,c(list(2,seq(ylim_min,ylim_max,length.out = 5),format(seq(ylim_min,ylim_max,length.out = 5),digits = 6)),
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))
                        
                        do.call(graphics::legend,c(list("top",legend = paste(legendName,age_names[match(as.character(Groups[[j]]), age_names)],sep = ""),
                                                        fill = colos[match(as.character(Groups[[j]]), age_names)],bty = "n",ncol = ceiling(nElements[j]/2)),
                                                   dots[!is.element(names(dots),c(axis.args,plot.args))]))
                        
                        iterator <- nElements[j]
                        
                }
        }
}


