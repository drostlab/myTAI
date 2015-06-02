#' @title Plot the Expression Profiles of a Gene Set
#' @description
#' This function simply visualizes the gene expression profiles of
#' a defined subset of genes stored in the input \code{ExpressionSet}.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param gene.set a character vector storing the gene ids for which gene expression profiles shall be visualized. 
#' @param get.subset a logical value indicating whether or not an \code{ExpressionSet} subset of the selected \code{gene.set} should be retuned. 
#' @param colors colors for gene expression profiles. Default: \code{colors = NULL}, hence default colours are used.
#' @param ... additional parameters passed to \code{\link{matplot}}.
#' @author Hajk-Georg Drost
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # png("test_png.png",700,400)
#' PlotGeneSet(ExpressionSet = PhyloExpressionSetExample, 
#'             gene.set      = PhyloExpressionSetExample[1:5, 2], 
#'             type          = "l", 
#'             lty           = 1, 
#'             lwd           = 4,
#'             xlab          = "Ontogeny",
#'             ylab          = "Expression Level")
#' 
#' # dev.off()
#' @export

PlotGeneSet <- function(ExpressionSet, 
                        gene.set, 
                        get.subset  = FALSE, 
                        colors      = NULL,
                        y.ticks     = 6,
                        digits.ylab = 4, ... ){
        
        is.ExpressionSet(ExpressionSet)
        
        GeneSubSet.indixes <- na.omit(match(tolower(gene.set), tolower(ExpressionSet[ , 2])))
        
#         if (is.na(GeneSubSet.indixes))
#                 stop ("None of your input gene ids could be found in the ExpressionSet.")
                                    
        if (length(GeneSubSet.indixes) != length(gene.set))
                warning ("Only ",length(GeneSubSet.indixes), " out of your ", length(gene.set), "gene ids could be found in the ExpressionSet.")
        
        GeneSubSet <- ExpressionSet[GeneSubSet.indixes , ]
        ncols <- ncol(GeneSubSet)
        
        if (is.null(colors))
                colors <- re.colors(length(GeneSubSet.indixes))
                                    
        if (length(colors) != length(gene.set))
                stop ("The number of colors and the number of genes do not match.")
                
        # define arguments for different graphics functions
        plot.args <- c("type","lwd","col","cex.lab","main","xlab","ylab")
        axis.args <- c("las", "cex.axis")
        legend.args <- c("border","angle","density","box.lwd","cex")
        dots <- list(...)
        ellipsis.names <- names(dots)
                
        ylim.range <- range(min(GeneSubSet[ , 3:ncols]),max(GeneSubSet[ , 3:ncols]))  
                
        if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) | (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                
                par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                        
                do.call(graphics::matplot,c(list(x = t(GeneSubSet[ , 3:ncols]), 
                                                 col  = colors[1:length(GeneSubSet.indixes)], 
                                                 axes = FALSE), 
                                            dots[!is.element(names(dots),c(axis.args,legend.args))]))
                                
        } else {      
                        do.call(graphics::matplot,c(list(x = t(GeneSubSet[ , 3:ncols]), 
                                                         col  = colors[1:length(GeneSubSet.indixes)], 
                                                         xaxt = "n", 
                                                         xlab = "Ontogeny",
                                                         ylab = "Expression Level",
                                                         axes = FALSE ), 
                                                    dots[!is.element(names(dots),c(axis.args,legend.args))]))
                }
                
                do.call(graphics::axis,c(list(side = 1,at = 1:(ncols-2),
                                              labels = names(ExpressionSet)[3:ncols]),dots[!is.element(names(dots),c(plot.args,legend.args))]))
                
                do.call(graphics::axis,c(list(side = 2,at = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab),
                                              labels = format(seq(ylim.range[1],ylim.range[2],length.out = y.ticks),digits = digits.ylab)), 
                                         dots[!is.element(names(dots),c(plot.args,legend.args))]))
                
                do.call(graphics::legend, c(list("topright", 
                       inset  = c(-0.2,0), 
                       legend = GeneSubSet[ , 2], 
                       title  = "Genes", 
                       col    = colors[1:length(GeneSubSet.indixes)],
                       lwd    = 4,
                       bty    = "n"), dots[!is.element(names(dots),c(plot.args,axis.args))]))

        if (get.subset)
                return(GeneSubSet)
}







