#' @title Plot the Expression Profiles of a Gene Set
#' @description
#' This function simply visualizes the gene expression profiles of
#' a defined subset of genes stored in the input \code{ExpressionSet}.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param gene.set a character vector storing the gene ids for which gene expression profiles shall be visualized. 
#' @param get.subset a logical value indicating whether or not an \code{ExpressionSet} subset of the selected \code{gene.set} should be retuned. 
#' @param use.only.map a logical value indicating whether instead of a standard \code{ExpressionSet} only a \code{Phylostratigraphic Map} or \code{Divergene Map} is passed to the function.
#' @param colors colors for gene expression profiles. Default: \code{colors = NULL}, hence default colours are used.
#' @param plot.legend a logical value indicating whether gene ids should be printed as legend next to the plot.
#' @param y.ticks a numeric value specifying the number of ticks to be drawn on the y-axis.
#' @param digits.ylab a numeric value specifying the number of digits shown for the expression levels on the y-axis.
#' @param ... additional parameters passed to \code{\link{matplot}}.
#' @author Hajk-Georg Drost
#' @details
#' 
#' This function simply visualizes or subsets the gene expression levels of a set of genes
#' that are stored in the input \code{ExpressionSet}.
#' 
#' @seealso \code{\link{SelectGeneSet}}, \code{\link{PlotEnrichment}}, \code{\link{DiffGenes}}  
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # the best parameter setting to visualize this plot:
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
#' 
#' # In case you would like to work with the expression levels
#' # of selected genes you can specify the 'get.subset' argument:
#' 
#' PlotGeneSet(ExpressionSet = PhyloExpressionSetExample, 
#'             gene.set      = PhyloExpressionSetExample[1:5, 2], 
#'             get.subset    = TRUE)
#' 
#' 
#' # get a gene subset using only a phylostratihraphic map
#' ExamplePSMap <- PhyloExpressionSetExample[ , 1:2]
#' 
#' PlotGeneSet(ExpressionSet = ExamplePSMap, 
#'             gene.set      = PhyloExpressionSetExample[1:5, 2], 
#'             get.subset    = TRUE,
#'             use.only.map  = TRUE)
#'             
#' @export

PlotGeneSet <- function(ExpressionSet, 
                        gene.set, 
                        get.subset   = FALSE,
                        use.only.map = FALSE,
                        colors       = NULL,
                        plot.legend  = TRUE,
                        y.ticks      = 6,
                        digits.ylab  = 4, ... ){
        
        if (!use.only.map)
                is.ExpressionSet(ExpressionSet)
        
        GeneSubSet.indixes <- stats::na.omit(match(tolower(gene.set), tolower(ExpressionSet[ , 2])))
        
        if (length(GeneSubSet.indixes) == 0)
                stop ("None of your input gene ids could be found in the ExpressionSet.")
                                    
        if (length(GeneSubSet.indixes) != length(gene.set))
                warning ("Only ",length(GeneSubSet.indixes), " out of your ", length(gene.set), " gene ids could be found in the ExpressionSet.")
        
        GeneSubSet <- ExpressionSet[GeneSubSet.indixes , ]
        ncols <- ncol(GeneSubSet)
        
        if(!is.null(colors)){
                if (length(colors) != length(gene.set))
                        stop ("The number of colors and the number of genes do not match.")
        }
        
        # http://www.compbiome.com/2010/12/r-using-rcolorbrewer-to-colour-your.html
        if (is.null(colors))
                 colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(GeneSubSet.indixes))                   
        if(!get.subset){
                # define arguments for different graphics functions
                plot.args <- c("type","lwd","col","cex.lab","main","xlab","ylab")
                axis.args <- c("las", "cex.axis")
                legend.args <- c("border","angle","density","box.lwd","cex")
                dots <- list(...)
                ellipsis.names <- names(dots)
                
                ylim.range <- range(min(GeneSubSet[ , 3:ncols]),max(GeneSubSet[ , 3:ncols]))  
                
                if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) | (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                        
                        if(plot.legend)
                                graphics::par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                        
                        do.call(graphics::matplot,c(list(x = t(GeneSubSet[ , 3:ncols]), 
                                                         col  = colors[1:length(GeneSubSet.indixes)], 
                                                         axes = FALSE), 
                                                    dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        
                } else {      
                        
                        if(plot.legend)
                                graphics::par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                        
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
                
                if(plot.legend){
                        do.call(graphics::legend, c(list("topright", 
                                                         inset  = c(-0.2,0), 
                                                         legend = GeneSubSet[ , 2], 
                                                         title  = "Genes", 
                                                         col    = colors[1:length(GeneSubSet.indixes)],
                                                         lwd    = 4,
                                                         bty    = "n"), dots[!is.element(names(dots),c(plot.args,axis.args))]))
                }
        }

        if (get.subset)
                return(GeneSubSet)
}








