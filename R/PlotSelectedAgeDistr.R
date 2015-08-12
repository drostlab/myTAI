#' @title Plot the PS or DS distribution of a selected set of genes
#' @description This function visualizes the PS or DS distribution of a selected set of genes
#' as histogram.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param gene.set a character vector storing the gene ids for which gene expression profiles shall be visualized.
#' @param legendName a character string specifying whether "PS" or "DS" are are visualized.
#' @param as.ratio logical value indicating whether or not relative frequencies shall be visualized.
#' @param use.only.map logical value indicating whether or not a Phylostratigraphic Map or Divergence Map should be passed to the \code{ExpressionSet} argument instead of a standard \code{ExpressionSet} object.
#' @param col colour of the bars.
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @author Hajk-Georg Drost
#' @examples 
#' data(PhyloExpressionSetExample)
#' 
#' # generate an example gene set
#' set.seed(123)
#' ExGeneSet <- sample(PhyloExpressionSetExample[ , 2], 5000)
#' 
#' # gene count example
#' PlotSelectedAgeDistr(ExpressionSet = PhyloExpressionSetExample,
#'                      gene.set      = ExGeneSet,
#'                      legendName    = "PS",
#'                      as.ratio      = TRUE)
#' 
#' # relative gene count example
#' PlotSelectedAgeDistr(ExpressionSet = PhyloExpressionSetExample,
#'                      gene.set      = ExGeneSet,
#'                      legendName    = "PS",
#'                      as.ratio      = FALSE)
#' 
#' @seealso \code{\link{PlotDistribution}}
#' @export


PlotSelectedAgeDistr <- function(ExpressionSet, gene.set,
                                 legendName   = NULL,
                                 as.ratio     = FALSE,
                                 use.only.map = FALSE,
                                 col          = "turquoise4",
                                 xlab         = NULL,
                                 ylab         = NULL){
        
        if (is.null(legendName))
                stop ("Please specify the type of ExpressionSet you are working with: legendName = 'PS' or 'DS'.")
        
        if (!use.only.map)
                is.ExpressionSet(ExpressionSet)
        
        AgeCategory <- GeneCount <- NULL
        
        GeneSubSet <- SelectGeneSet(ExpressionSet = ExpressionSet,
                                    gene.set      = gene.set,
                                    use.only.map  = use.only.map)
        
        names(GeneSubSet)[1] <- legendName
        
        nPS <- length(names(table(ExpressionSet[ , 1])))
        
        if (!as.ratio){
                
                if (is.null(xlab)){
                        xlab <- paste0("\n",legendName)
                }
                
                if (is.null(ylab)){
                        ylab <- "# Genes\n"
                }
                
                AgeDistr <- table(factor(GeneSubSet[ , 1], levels = 1:nPS))
        }
                
        if (as.ratio){
                
                if (is.null(xlab)){
                        xlab <- paste0("\n",legendName)
                }
                
                if (is.null(ylab)){
                        ylab <- "Rel. Number Of  Genes\n"
                }
                AgeDistr <- table(factor(GeneSubSet[ , 1], levels = 1:nPS))
                AgeDistr <- format(AgeDistr / sum(AgeDistr), digits = 2)
        }
                
        GeneSubSet.LongFormat <- data.frame(AgeCategory = paste0(legendName,names(AgeDistr)), GeneCount = as.numeric(AgeDistr))
        # get the correct order of PS or DS by specifying levels = paste0(legendName,names(AgeDistr))
        GeneSubSet.LongFormat[ , 1] <- factor(GeneSubSet.LongFormat[ , 1], levels = paste0(legendName,names(AgeDistr)))
        plot.res <- ggplot2::ggplot(GeneSubSet.LongFormat, ggplot2::aes(x = AgeCategory, y = GeneCount), order = FALSE) + ggplot2::geom_bar(stat="identity", fill = col) + ggplot2::labs(x = xlab, y = ylab)  + ggplot2::geom_text(ggplot2::aes( label = GeneCount, y =  GeneCount), stat= "identity", vjust = -.3) + ggplot2::theme_minimal()
        
        return (plot.res)
}



