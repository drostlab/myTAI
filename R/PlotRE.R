#' @title Plot Relative Expression Levels
#' @description 
#' This function computes for each phylostratum or divergence stratum the corresponding relative expression profile
#' and plots the profiles in N different windows corresponding to the given phylostratum classes 
#' or divergence-stratum classes that shall be compared.
#' 
#' For each phylostratum or divergence-stratum the corresponding relative expression profile is being computed as follows:
#' 
#' \deqn{f_js = ( e_js - e_j min ) / ( e_j max - e_j min )}
#'
#' where \eqn{e_j min} and \eqn{e_j max} denote the minimum/maximum \code{\link{mean}} expression level 
#' of phylostratum j over  developmental stages s. This linear transformation corresponds to 
#' a shift by \eqn{e_j min} and a subsequent shrinkage by \eqn{e_j max - e_j min}. 
#' As a result, the relative expression level \eqn{f_js} of developmental stage s 
#' with minimum \eqn{e_js} is 0, the relative expression level \eqn{f_js} of the developmental 
#' stage s with maximum \eqn{e_js} is 1, and the relative expression levels \eqn{f_js} of 
#' all other stages s range between 0 and 1, accordingly.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the phylostrata or divergence strata that correspond 
#' to the same phylostratum class or divergence class.
#' For ex. evolutionary old phylostrata: PS1-3 (Class 1) 
#' and evolutionary young phylostrata: PS4-12 (Class 2). In this case, 
#' the list could be assigned as, \code{Groups} = list(c(1:3), c(4:12)). 
#' It is also possible to define more than 2 groups of evolutionary ages.
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param colors colors for relative expression profiles. Default: \code{colors = NULL}, hence default colours are used.
#' @param \dots default graphics parameters.
#' @details Studying the relative expression profiles of each phylostratum or divergence-stratum enables the detection
#' of common gene expression patterns shared by several phylostrata or divergence-strata.
#'
#' Finding similar relative expression profiles among phylostrata or divergence-strata suggests that 
#' phylostrata or divergence-strata sharing a similar relative expression profile are regulated by similar
#' gene regulatory elements. Hence, these common phylostrata or divergence-strata might govern similar processes in the given developmental time course. 
#' @return a plot showing the relative expression profiles of phylostrata or divergence-strata belonging to the same group.
#' @references 
#' Domazet-Loso T and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.
#'
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{PlotBarRE}}, \code{\link{RE}}, \code{\link{REMatrix}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet
#' PlotRE(PhyloExpressionSetExample,Groups = list(c(1:3), c(4:12)), 
#'        legendName = "PS", lty = 1, lwd = 5)
#'
#'
#' # or you can choose any combination of groups
#' PlotRE(PhyloExpressionSetExample,Groups = list(c(1,7,9), c(2:6,8,10:12)),
#'        legendName = "PS", lty = 1, lwd = 5)
#' 
#' # or multiple groups
#' PlotRE(PhyloExpressionSetExample,Groups = list(c(1,7,9), c(3:6,8),c(2,10:12)),
#'        legendName = "PS", lty = 1, lwd = 5)
#'        
#'        
#'        
#' # example DivergenceExpressionSet
#' PlotRE(DivergenceExpressionSetExample,Groups = list(c(1:5), c(6:10)), 
#'        legendName = "DS", lty = 1, lwd = 5)
#'
#'
#'
#' # adding custom colors for relative expression levels:
#' # -> colors should be ordered by PS/DS starting with PS1,2,3...
#' PlotRE(PhyloExpressionSetExample,
#'        Groups     = list(c(1:3), c(4:12)), 
#'        legendName = "PS",
#'        colors     = c("black","red","green","brown","darkmagenta",
#'        "blue","darkred","darkblue","darkgreen", "orange",
#'        "azure4","gold4"), 
#'        lty        = 1, 
#'        lwd        = 5)
#'   
#' @export

PlotRE <- function(ExpressionSet,
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
        
        nPS <- length(age_names)
        nCols <- dim(ExpressionSet)[2]
        ### define and label the REmatrix that holds the rel. exp. profiles
        ### for the available PS
        REmatrix <- matrix(NA_real_,nPS,nCols-2)
        rownames(REmatrix) <- age_names
        colnames(REmatrix) <- names(ExpressionSet)[3:nCols]
        nGroups <- length(Groups)
        ### each PS class gets its corresponding color
        
        if(!is.null(colors)){
                colos <- colors
        } else {
                colos <- re.colors(nPS)
        }
        
        # number of items in each Groups element
        nElements <- sapply(Groups,length)
        
        # compute the relative expression matrix
        REmatrix <- age.apply(ExpressionSet = ExpressionSet, RE)
        
        # define arguments for different graphics functions
        plot.args <- c("lwd","col","lty","xlab","cex.lab","main")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
        
        ### plot the rel. exp. levels in k different windows
        if(length(Groups) > 1)
                graphics::par(mfrow = rev(grDevices::n2mfrow(nGroups)))
        
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