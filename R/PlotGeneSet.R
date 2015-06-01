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
                        get.subset = FALSE, 
                        colors     = NULL, ... ){
        
        is.ExpressionSet(ExpressionSet)
        
        GeneSubSet.indixes <- na.omit(match(tolower(gene.set), tolower(ExpressionSet[ , 2])))
        
#         if (is.na(GeneSubSet.indixes))
#                 stop ("None of your input gene ids could be found in the ExpressionSet.")
                                    
        if (length(GeneSubSet.indixes) != length(gene.set))
                warning ("Only ",length(GeneSubSet.indixes), " out of your ", length(gene.set), "gene ids could be found in the ExpressionSet.")
        
        GeneSubSet <- ExpressionSet[GeneSubSet.indixes , ]
        ncols <- ncol(GeneSubSet)
        
        if (is.null(colors)){
                
                par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                matplot(t(GeneSubSet[ , 3:ncols]), 
                        col  = re.colors(length(GeneSubSet.indixes)),
                        xaxt = "n", ...)
                axis(1,at = seq(1,(ncols-2),1), labels = names(ExpressionSet)[3:ncols])
                
                legend("topright", 
                       inset  = c(-0.2,0), 
                       legend = GeneSubSet[ , 2], 
                       title  = "Genes", 
                       fill   = re.colors(length(GeneSubSet.indixes)),
                       bty    = "n")
        }
        
        if (!is.null(colors)){
                
                if (length(colors) != length(gene.set))
                        stop ("The number of colors and the number of genes do not match.")
                
                
                par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                matplot(t(GeneSubSet[ , 3:ncols]), 
                        col  = colors[1:length(GeneSubSet.indixes)], 
                        xaxt = "n", ...)
                axis(1,at = seq(1,(ncols-2),1), labels = names(ExpressionSet)[3:ncols])
                
                legend("topright", 
                       inset  = c(-0.2,0), 
                       legend = GeneSubSet[ , 2], 
                       title  = "Genes", 
                       fill   = colors[1:length(GeneSubSet.indixes)],
                       bty    = "n")
        }
        
        if (get.subset)
                return(GeneSubSet)
}








