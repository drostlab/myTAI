#' @title Select a Subset of Genes in an ExpressionSet
#' @description
#' Select a subset of genes stored in the input \code{ExpressionSet}.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param gene.set a character vector storing the gene ids for which gene expression profiles shall be visualized.
#' @param use.only.map a logical value indicating whether instead of a standard \code{ExpressionSet} only a \code{Phylostratigraphic Map} or \code{Divergene Map} is passed to the function. 
#' @author Hajk-Georg Drost
#' @details
#' 
#' This function selects a subset of genes specified in \code{gene.set} stored in the input \code{ExpressionSet} and returns a subset \code{ExpressionSet}.
#' 
#' This function is useful for studying the evolutionary \emph{properties} of a subset of genes
#' stored in the \code{ExpressionSet}. 
#' 
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # receive a subset ExpressionSet for the fist 5 genes stored in
#' # the PhyloExpressionSetExample
#' SelectGeneSet(ExpressionSet = PhyloExpressionSetExample,
#'             gene.set      = PhyloExpressionSetExample[1:5, 2])
#'             
#'             
#' # get a gene subset using only a phylostratihraphic map
#' ExamplePSMap <- PhyloExpressionSetExample[ , 1:2]
#' 
#' SelectGeneSet(ExpressionSet = ExamplePSMap,
#'               gene.set      = PhyloExpressionSetExample[1:5, 2],
#'               use.only.map  = TRUE)          
#'             
#' @seealso \code{\link{PlotGeneSet}}, \code{\link{PlotEnrichment}}, \code{\link{DiffGenes}}  
#' @export

SelectGeneSet <- function(ExpressionSet, gene.set, use.only.map = FALSE) {
        
        return( PlotGeneSet(ExpressionSet = ExpressionSet,
                            gene.set      = gene.set,
                            get.subset    = TRUE,
                            use.only.map  = use.only.map) )
        
        
}




