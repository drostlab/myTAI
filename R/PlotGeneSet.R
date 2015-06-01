#' @title Plot the Expression Profiles of a Gene Set
#' @description
#' This function simply visualizes the gene expression profiles of
#' a defined subset of genes stored in the input \code{ExpressionSet}.
#' @param ExpressionSet
#' @param gene.set
#' @param get.subset
#' @param colors
#' @param ... additional parameters passed to \code{\link{matplot}}.
#' @author Hajk-Georg Drost
#' @examples
#' data(PhyloExpressioNSetExample)
#' 
#' @export

PlotGeneSet <- function(ExpressionSet, gene.set, get.subset = FALSE, colors = NULL, ... ){
        
        is.ExpressionSet(ExpressionSet)
        
        GeneSubSet.indixes <- match(tolower(gene.set), tolower(ExpressionSet[ , 2]))
        
#         if (is.na(GeneSubSet.indixes))
#                 stop ("None of your input gene ids could be found in the ExpressionSet.")
                                    
        if (length(GeneSubSet.indixes) != length(gene.set))
                warning ("Only ",length(GeneSubSet.indixes), " out of your ", length(gene.set), "gene ids could be found in the ExpressionSet.")
        
        GeneSubSet <- ExpressionSet[GeneSubSet.indixes , ]
        ncols <- ncol(GeneSubSet)
        
        if (is.null(colors)){
                
                par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
                matplot(t(GeneSubSet[ , 3:ncols]), col = re.colors(length(GeneSubSet.indixes)),...)
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
                        col = colors[1:length(GeneSubSet.indixes)], ...)
                
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








