#' @title Plot the Frequency Distribution of Phylostrata or Divergence Strata
#' @description This function plots the frequency distribution of genes within the 
#' corresponding \emph{phylostratigraphic map} or \emph{divergence map} and can be used to fastly visualize the PS or DS distribution of a given phylostratum vector or divergence-stratum vector.
#' @param PhyloExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param legendName a character string specifying whether "PS" or "DS" are visualized.
#' @param as.ratio a boolean value specifying whether the relative frequencies
#' instead of absolute frequencies shall be plotted.
#' @param use.only.map logical value indicating whether or not a Phylostratigraphic Map or Divergence Map should be passed to the \code{ExpressionSet} argument instead of a standard \code{ExpressionSet} object.
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @details 
#' The frequency distribution of all genes or a subset of genes might be of interest for subsequent analyses.
#'
#' For Example:
#'      
#' Filtering genes using gene cluster algorithms can result in different groups (classes) of genes.
#' For each gene group the phylostratum or divergence-stratum distribution can be visualized using this function
#' and can be compared between different groups.
#'
#' This analysis allows to compare different gene expression profiles (or gene groups in general) based
#' on their evolutionary origins or evolutionary relationships.
#' @return a barplot showing the phylostratum distribution or 
#' divergence-stratum distribution of a given numeric vector containing PS or DS values.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # load PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # plot the phylostratum distribution of a PhyloExpressionSet
#' PlotDistribution(PhyloExpressionSetExample)
#'
#' # plot the relative frequency distribution of a PhyloExpressionSet
#' PlotDistribution(PhyloExpressionSetExample, as.ratio = TRUE)
#'
#'
#' # a example for visualizing the PS distribution for a subset of genes
#' PlotDistribution(PhyloExpressionSetExample[sample(20000,5000) , ], as.ratio = TRUE)
#'
#' @seealso \code{\link{PlotSelectedAgeDistr}}
#' @export

PlotDistribution <- function(PhyloExpressionSet,
                             legendName   = "PS",
                             as.ratio     = FALSE,
                             use.only.map = FALSE,
                             xlab         = NULL,
                             ylab         = NULL)
{
        
        PlotSelectedAgeDistr(ExpressionSet = PhyloExpressionSet,
                             gene.set      = PhyloExpressionSet[ , 2],
                             legendName    = legendName,
                             as.ratio      = as.ratio,
                             use.only.map  = use.only.map,
                             xlab          = xlab,
                             ylab          = ylab)
        
        
}

