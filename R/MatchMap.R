#' @title Match a Phylostratigraphic Map or Divergence Map with a ExpressionMatrix
#' @description This function matches a \emph{Phylostratigraphic Map} or \emph{Divergence Map} only storing unique gene ids with a ExpressionMatrix
#' also storing only unique gene ids.
#' @param Map a standard \emph{Phylostratigraphic Map} or \emph{Divergence Map} object.
#' @param ExpressionMatrix  a standard ExpressionMatrix object.
#' @param remove.duplicates a logical value indicating whether duplicate gene ids should be removed from the data set.
#' @param accumulate an accumulation function such as \code{mean()}, \code{median()}, or \code{min()}
#' to accumulate multiple expression levels that map to the same unique gene id present in the \code{ExpressionMatrix}.
#' @details
#' 
#' In phylotranscriptomics analyses two major techniques are performed to
#' obtain standard \emph{Phylostratigraphic map} or \emph{Divergence map} objects.
#' 
#' To obtain a \emph{Phylostratigraphic Map}, \emph{Phylostratigraphy} (Domazet-Loso et al., 2007) has to be performed. To obtain a \emph{Divergence Map}, 
#' orthologous gene detection, Ka/Ks computations, and decilation (Quint et al., 2012; Drost et al., 2015) have to be performed.
#' 
#' The resulting standard \emph{Phylostratigraphic Map} or \emph{Divergence Map} objects consist of 2 colums storing the phylostratum assignment 
#' or divergence stratum assignment of a given gene in column one, and the corresponding gene id of that gene on columns two.
#' 
#' A standard ExpressionMatrix is a common gene expression matrix storing the gene ids or probe ids in the first column, and all experiments/stages/replicates in the following columns.
#' 
#' The \emph{MatchMap} function takes both standard datasets: \emph{Map} and \emph{ExpressionMatrix} to merge both data sets to obtain a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' 
#' This procedure is analogous to \code{\link{merge}}, but is customized to the \emph{Phylostratigraphic Map}, \emph{Divergence Map}, and \emph{ExpressionMatrix} standards to allow a faster and more intuitive usage. 
#' 
#' In case you work with an ExpressionMatrix that stores multiple expression levels for a unique gene id, you
#' can specify the \code{accumulation} argument to accumulate these multiple expression levels to obtain
#' one expression level for one unique gene.
#' 
#' 
#' @return a standard PhyloExpressionSet or DivergenceExpressionSet object. 
#' @references 
#' 
#'   Domazet-Loso T, Brajkovic J, Tautz D (2007) A phylostratigraphy approach to uncover the genomic history of major adaptations in metazoan lineages. Trends Genet. 23: 533-9.
#' 
#'   Domazet-Loso T, Tautz D (2010) A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature 468: 815-8.
#' 
#'   Quint M., Drost H.G., Gabel A., Ullrich K.K., Boenn M., Grosse I. (2012) A transcriptomic hourglass in plant embryogenesis. Nature 490: 98-101.
#' 
#'   Drost HG et al. (2015) \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' 
#' @author Hajk-Georg Drost
#' @examples
#'         
#' # load a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'         
#' # in a standard PhyloExpressionSet, 
#' # column one and column two denote a standard 
#' # phylostratigraphic map
#' PhyloMap <- PhyloExpressionSetExample[ , 1:2]
#'         
#' # look at the phylostratigraphic map standard
#' head(PhyloMap)
#'         
#' # in a standard PhyloExpressionSet, column two combined 
#' # with column 3 - N denote a standard ExpressionMatrix
#' ExpressionMatrixExample <- PhyloExpressionSetExample[ , c(2,3:9)]
#'         
#' # these two data sets shall illustrate an example 
#' # phylostratigraphic map that is returned
#' # by a standard phylostratigraphy run, and a expression set 
#' # that is the result of expression data analysis 
#' # (background correction, normalization, ...)
#'         
#' # now we can use the MatchMap function to merge both data sets
#' # to obtain a standard PhyloExpressionSet
#'         
#' PES <- MatchMap(PhyloMap, ExpressionMatrixExample)
#'         
#' # note that the function returns a head() 
#' # of the matched gene ids to enable
#' # the user to find potential mis-matches
#'         
#' # the entire procedure is analogous to merge() 
#' # with two data sets sharing the same gene ids 
#' # as column (primary key)
#' PES_merge <- merge(PhyloMap, ExpressionMatrixExample)
#'         
#'         
#'         
#'@export
MatchMap <- function(Map,ExpressionMatrix, remove.duplicates = FALSE, accumulate = NULL)
{
        
        names(ExpressionMatrix)[1] <- "GeneID"
        names(Map)[2] <- "GeneID"
        ExpressionMatrix[ , "GeneID"] <- tolower(ExpressionMatrix[ , "GeneID"])
        Map[ , "GeneID"] <- tolower(Map[ , "GeneID"])
        GeneID <- NULL
        
        if(remove.duplicates)
                Map <- Map[-which(duplicated(Map[ , "GeneID"])) , ]
        
        if(any(duplicated(Map[ , "GeneID"])))
                stop("You have duplicate Gene IDs in your Map. Please enter only unique Gene IDs.")
        
        if(!is.null(accumulate)){
                
                acc_fun <- match.fun(accumulate)
                ExpressionMatrix <- dplyr::summarise_each(dplyr::group_by(ExpressionMatrix, GeneID), dplyr::funs(acc_fun))
                
        }
        
        if(any(duplicated(ExpressionMatrix[ , "GeneID"])))
                stop("You have duplicate Gene IDs in your ExpressionMatrix. Please enter only unique Gene IDs, or specify the 'accumulate' argument.")
        
        joined_ExpressionMatrix <- dplyr::semi_join(ExpressionMatrix, Map, by = "GeneID")
        
        res_tbl <- merge(joined_ExpressionMatrix,Map, by = "GeneID")
        
        if(!any(duplicated(res_tbl[ , "GeneID"]))){
                
                return(res_tbl[ , c(ncol(res_tbl),1:(ncol(res_tbl)-1))])
                
        } else {
                stop("Something went wrong with matching Map and ExpressionMatrix! Please check for duplicate entries!")
        }
        
}
