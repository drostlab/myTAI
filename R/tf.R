#' @title Transform Gene Expression Levels
#' @description This function transforms the gene expression set stored in an input 
#' PhloExpressionSet or DivergenceExpressionSet object and returns a 
#' PhloExpressionSet or DivergenceExpressionSet object with transformed expression levels. 
#' The resulting transformed PhloExpressionSet or DivergenceExpressionSet 
#' object can then be used for subsequent analyses based on transformed expression levels.
#' @param ExpressionSet a standard PhloExpressionSet or DivergenceExpressionSet object.
#' @param FUN any valid function that transformes gene expression levels.
#' @details Motivated by the dicussion raised by Piasecka et al., 2013, the influence of
#' gene expression transformation on the global phylotranscriptomics pattern does not seem negligible.
#' Hence, different transformations can result in qualitatively different \code{\link{TAI}} or \code{\link{TDI}}
#' patterns.
#'
#' Initially, the \code{\link{TAI}} and \code{\link{TDI}} formulas were defined for absolute expression levels.
#' So using the initial \code{\link{TAI}} and \code{\link{TDI}} formulas with transformed expression levels
#' might turn out in qualitatively different patterns when compared with non-transformed expression levels,
#' but might also belong to a different class of models, since different valid expression level transformation functions result in different patterns.
#'
#' The purpose of this function is to allow the user to study the qualitative impact of different transformation functions on 
#' the global \code{\link{TAI}} and \code{\link{TDI}} pattern, or on any subsequent phylotranscriptomics analysis.
#'
#' The examples using the \emph{PhyloExpressionSetExample} data set show that using common gene expression 
#' transformation functions: \code{\link{log2}} (Quackenbush, 2001 and 2002), \code{\link{sqrt}} (Yeung et al., 2001), 
#' \code{\link[MASS]{boxcox}}, or \emph{inverse hyperbolic sine transformation}, each transformation results 
#' in qualitatively different patterns. 
#' Nevertheless, for each resulting pattern the statistical significance can be tested 
#' using either the \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}} (Drost et al., 2014) 
#' to quantify the significance of interest.
#' @return a standard PhloExpressionSet or DivergenceExpressionSet object storing transformed gene expression levels.
#' @references
#' Piasecka B, Lichocki P, Moretti S, et al. (2013) The hourglass and the early conservation models--co-existing patterns of developmental constraints in vertebrates. PLoS Genet. 9(4): e1003476.
#'
#' Quint M., Drost H.G., Gabel A., Ullrich K.K., Boenn M., Grosse I. (2012) A transcriptomic hourglass in plant embryogenesis. Nature 490: 98-101.
#'
#' Domazet-Loso T., Tautz D. (2010) A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature 468: 815-8.
#'
#' Drost H.G., Gabel A., Grosse I., Quint M. (2015) Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. Mol Biol Evol. 32 (5): 1221-1231. doi: 10.1093/molbev/msv012
#'
#' K.Y. Yeung et al.: Model-based clustering and data transformations for gene expression data. Bioinformatics 2001, 17:977-987
#'
#' K.Y. Yeung et al.:  Supplement to Model-based clustering and data transformations for gene expression data - Data Transformations and the Gaussian mixture assumption. Bioinformatics 2001, 17:977-987
#'
#' P.A.C. Hoen et al.: Deep sequencing-based expression analysis shows major advances in robustness, resolution and inter-lab portability over five microarray platforms. Nucleic Acids Research 2008, Vol. 36, No. 21
#'
#' H.H. Thygesen et al.: Comparing transformation methods for DNA microarray data. BMC Bioinformatics 2004, 5:77
#'
#' John Quackenbush: Microarray data normalization and transformation. Nature Genetics 2002, 32:496-501
#'
#' John Quackenbush: Computational Analysis of Microarray Data. Nature Reviews 2001, 2:418-427
#'
#' R. Nadon and J. Shoemaker: Statistical issues with microarrays: processing and analysis. TRENDS in Genetics 2002, Vol. 18 No. 5:265-271
#'
#' B.P. Durbin et al.: A variance-stabilizing transformation for gene-expression microarray data. Bioinformatics 2002, 18:S105-S110
#'
#' J. M. Bland et al.: Transforming data. BMJ 1996, 312:770
#'
#' \url{http://stats.stackexchange.com/questions/1444/how-should-i-transform-non-negative-data-including-zeros}
#'
#' \url{http://stats.stackexchange.com/questions/78929/how-can-i-estimate-theta-for-the-inverse-hyperbolic-sine-transformation}
#'
#' John B. Burbidge, Lonnie Magee and A. Leslie Robb (1988) Alternative Transformations to Handle Extreme Values of the Dependent Variable. Journal of the American Statistical Association, 83(401): 123-127.
#'
#' G. E. P. Box and D. R. Cox (1964) An Analysis of Transformations. Journal of the Royal Statistical Society. Series B (Methodological), 26(2): 211-252.
#' 
#' @author Hajk-Georg Drost
#' @seealso  \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#' 
#' # a simple example is to transform the gene expression levels
#' # of a given PhyloExpressionSet using a sqrt or log2 transformation
#' 
#' PES.sqrt <- tf(PhyloExpressionSetExample, sqrt)
#' 
#' PES.log2 <- tf(PhyloExpressionSetExample, log2)
#' 
#' # in case a given PhyloExpressionSet already stores gene expression levels
#' # that are log2 transformed and need to be re-transformed to absolute
#' # expression levels, to perform subsequent phylotranscriptomics analyses 
#' # (that are defined for absolute expression levels), one can re-transform
#' # a PhyloExpressionSet like this:
#' 
#' PES.absolute <- tf(PES.log2 , function(x) 2^x)
#' 
#' # which should be the same as  PhyloExpressionSetExample :
#' head(PhyloExpressionSetExample)
#' head(PES.absolute)
#' 
#' 
#' # plotting the TAI using log2 transformed expression levels
#' # and performing the Flat Line Test to obtain the p-value
#' PlotPattern(ExpressionSet = tf(PhyloExpressionSetExample, log2), 
#'             type          = "l", 
#'             lwd           = 5, 
#'             TestStatistic = "FlatLineTest")
#' 
#' @export

tf <- function(ExpressionSet, FUN){
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- dim(ExpressionSet)[2]
        f <- match.fun(FUN)
        res <- data.frame(ExpressionSet[ , 1:2] , apply(ExpressionSet[ , 3:ncols] , 2 , f))
        names(res) <- names(ExpressionSet)
        return(res)
}

