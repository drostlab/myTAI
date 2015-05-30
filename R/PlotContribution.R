#' @title Plot the Phylostratum or Divergence Stratum Contribution to the Global Pattern
#' @description This function computes the cumulative contribution of each Phylostratum or Divergence Stratum to the global \code{\link{TAI}} or \code{\link{TDI}} profile.
#' 
#' 
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param colors a character vector specifying the colors that should be used for each Phylostratum/Divergence Stratum. Default: \code{colors = NULL}, meaning that default colors are used.
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param digits.ylab a numeric value specifying the number of digits shown for the TAI or TDI values on the y-axis.
#' @param y.ticks a numeric value specifying the number of ticks to be drawn on the y-axis.
#' @param ... additional \code{\link{plot}} parameters.
#' @details
#' Introduced by Domazet-Loso and Tautz (2010), this function allows users to visualize the cumulative contribution of each Phylostratum or Divergence Stratum to the global Transcriptome Age Index or Transcriptome Divergence Index profile to quantify how each Phylostratum or Divergence Stratum influences the profile of the global TAI or TDI pattern. 
#' 
#' @author Hajk-Georg Drost
#' @examples
#' 
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'  
#'  # visualize phylostratum contribution to global TAI
#'  PlotContribution(PhyloExpressionSetExample, legendName = "PS")
#'  
#'  # visualize divergence stratum contribution to global TDI
#'  PlotContribution(DivergenceExpressionSetExample, legendName = "DS")
#'  
#' @references
#' 
#' Domazet-Loso T. and Tautz D. (2010). A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature (468): 815-818.
#'    
#' @seealso \code{\link{pTAI}}, \code{\link{pTDI}}, \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{PlotPattern}}
#' @export         

PlotContribution <- function(ExpressionSet, 
                             colors      = NULL, 
                             legendName  = NULL, 
                             digits.ylab = 3,
                             y.ticks     = 5, ...){
        
        if(is.null(legendName))
                stop("Please specify whether your input ExpressionSet stores 'PS' or 'DS'.")
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- ncol(ExpressionSet)
        nPS <- length(table(ExpressionSet[ , 1]))
                
        if(is.null(colors)){
                
                colos <- re.colors(nPS)
                
        } 
        else if (!is.null(colors)) {
                if(length(colors) != nPS)
                        stop("The number of colors and number of PS/DS do not match...")
                
                colos <- colors
        }
        
        # define contribution matrix
        contrMatrix <- matrix(NA_real_,ncol = ncols,nrow = nPS)
        contrMatrix <- pTAI(ExpressionSet)
        
        ylim.range <- range(0,max(contrMatrix) + (max(contrMatrix)/5)) 
        
        # define arguments for different graphics functions
        plot.args <- c("type","lwd","col","cex.lab","main","xlab","ylab")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
                
        if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) || (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                
                do.call(graphics::matplot,c(list(x = t(contrMatrix), ylim = ylim.range, type = "l",lty = 1,col = colos, axes = FALSE), 
                                                 dots[!is.element(names(dots),c(axis.args,legend.args))]))
                }
                
                # default: xlab = "Ontogeny" and ylab = "Age Index"
                else {
                        if (legendName == "PS"){
                                do.call(graphics::matplot,c(list(x = t(contrMatrix),ylim = ylim.range, type = "l",lty = 1,col = colos,axes = FALSE, 
                                                              xlab = "Ontogeny",ylab = "Transcriptome Age Index" ), dots[!is.element(names(dots),c(axis.args,legend.args))]))  
                        }
                        
                        else if (legendName == "DS"){
                                do.call(graphics::matplot,c(list(x = t(contrMatrix),ylim = ylim.range, type = "l",lty = 1,col = colos, axes = FALSE, 
                                                              xlab = "Ontogeny",ylab = "Transcriptome Divergence Index" ), dots[!is.element(names(dots),c(axis.args,legend.args))])) 
                                
                        }
                        
                }
                
                do.call(graphics::axis,c(list(side = 1,at = 1:(ncols-2),
                                              labels = names(ExpressionSet)[3:ncols]),dots[!is.element(names(dots),c(plot.args,legend.args))]))
                
                do.call(graphics::axis,c(list(side = 2,at = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab),
                                              labels = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab)), 
                                         dots[!is.element(names(dots),c(plot.args,legend.args))]))
                
                do.call(graphics::legend,c(list(x = "top",bty = "n",legend = paste0(legendName,1:nPS),fill = colos, ncol = floor(nPS/2)),
                                           dots[!is.element(names(dots),c(plot.args,axis.args))]))
                

}







