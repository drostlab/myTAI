#' @title Quantify the significant differences between gene expression distributions of PS or DS groups
#' @description This function performs a test to quantify the statistical significance between
#' the global expression level distributions of groups of PS or DS. It therefore, allows users to investigate
#' significant groups of PS or DS that significantly differ in their gene expression level distibution
#' within specific developmental stages or experiments.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the phylostrata or divergence strata that correspond 
#' to the same phylostratum class or divergence class.
#' For ex. evolutionary old phylostrata: PS1-3 (Class 1) 
#' and evolutionary young phylostrata: PS4-12 (Class 2). In this case, 
#' the list could be assigned as, \code{Groups} = list(c(1:3), c(4:12)).
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles. 
#' @param stat.test the statistical test to quantify PS or DS group differences.
#' @param gene.set a character vector storing the gene ids for which group specific differences shall be statistically quantified.
#' @param ... additional plot parameters.
#' @author Hajk-Georg Drost
#' @details 
#' The purpose of this function is to detect groups of PS or DS that significantly differ in their gene expression
#' level distributions on a global (transcriptome) level. Since relative expression levels (\code{\link{PlotRE}}) or
#' PS or DS specific mean expression levels (\code{\link{PlotMeans}}) are biased by highly expressed genes,
#' this function allows users to objectively test the significant difference of transcriptome expression between
#' groups of PS or DS in a specific developmental stage or experiment.
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' # perform a Wilcoxon Rank Sum test to statistically quantify the
#' # difference between PS-Group 1 expression levels versus PS-Group 2
#' # expression levels
#' GroupDiffs(ExpressionSet = PhyloExpressionSetExample,
#'            Groups       = list(group_1 = 1:3,group_2 = 4:12),
#'            legendName   = "PS")
#' 
#' # quantify the significant difference of a selected set of genes
#' set.seed(123)
#' ExampleGeneSet <- sample(PhyloExpressionSetExample[ , 2],5000)  
#'              
#' GroupDiffs(ExpressionSet = PhyloExpressionSetExample,
#'            Groups       = list(group_1 = 1:3,group_2 = 4:12),
#'            legendName   = "PS",
#'            gene.set     = ExampleGeneSet)               
#' 
#' 
#' @seealso \code{\link{PlotGroupDiffs}}, \code{\link{PlotMeans}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}, \code{\link{PlotCategoryExpr}}
#' @export

GroupDiffs <- function(ExpressionSet,
                       Groups     = NULL, 
                       legendName = NULL,
                       stat.test  = "wilcox.test",
                       gene.set   = NULL, ...){
        
        
        PlotGroupDiffs(ExpressionSet = ExpressionSet,
                       Groups        = Groups,
                       legendName    = legendName,
                       stat.test     = stat.test,
                       plot.p.vals   = FALSE,
                       gene.set      = gene.set, ...)
        
}



