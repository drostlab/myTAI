#' @title Phylostratum or Divergence Stratum Enrichment of a given Gene Set based on Fisher's Test
#' @description 
#' This function computes the significance of enriched (over or underrepresented) Phylostrata or Divergence Strata within an input \code{test.set} based on the \code{\link{fisher.test}}. Please concult \code{\link{PlotEnrichment}} for details.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param test.set a character vector storing the gene ids for which PS/DS enrichment analyses should be performed.
#' @param use.only.map a logical value indicating whether instead of a standard \code{ExpressionSet} only a \code{Phylostratigraphic Map} or \code{Divergene Map} is passed to this function.
#' @param measure a character string specifying the measure that should be used to quantify over and under representation of PS/DS. Measures can either be \code{measure = "foldchange"} (odds) or \code{measure = "log-foldchange"} (log-odds). 
#' @param complete.bg a logical value indicating whether the entire background set of the input ExpressionSet should be considered when performing Fisher's exact test (\code{complete.bg = TRUE}) or whether genes that are stored in test.set should be excluded from the background set before performing Fisher's exact test (\code{complete.bg = FALSE}).
#' @param epsilon a small value to shift values by epsilon to avoid log(0) computations.
#' @author Hajk-Georg Drost
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' 
#' set.seed(123)
#' test_set <- sample(PhyloExpressionSetExample[ , 2],1000)
#' 
#' E.Result <- EnrichmentTest(ExpressionSet = PhyloExpressionSetExample,
#'                            test.set      = test_set ,
#'                            measure       = "log-foldchange")
#'                            
#' # get the log-fold change table
#' E.Result$enrichment.matrix
#' 
#' # get P-values for the enrichment significance for each Phylostratum
#' E.Result$p.values
#' 
#' @seealso \code{\link{PlotEnrichment}}, \code{\link{fisher.test}}
#'                            
#' @export 

EnrichmentTest <- function(ExpressionSet, 
                             test.set,
                             use.only.map = FALSE,
                             measure      = "log-foldchange",
                             complete.bg  = TRUE,
                             epsilon      = 1e-05){
        
        
        return( PlotEnrichment(ExpressionSet = ExpressionSet,
                               test.set      = test.set,
                               use.only.map  = use.only.map,
                               measure       = measure,
                               complete.bg   = complete.bg,
                               plot.bars     = FALSE) )
}