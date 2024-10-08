#' @title Transform Gene Expression Levels
#' @description This function transforms the gene expression set stored in an input 
#' PhloExpressionSet or DivergenceExpressionSet object and returns a 
#' PhloExpressionSet or DivergenceExpressionSet object with transformed expression levels. 
#' The resulting transformed PhloExpressionSet or DivergenceExpressionSet 
#' object can then be used for subsequent analyses based on transformed expression levels.
#' @param ExpressionSet a standard PhloExpressionSet or DivergenceExpressionSet object.
#' @param FUN any valid function that transforms gene expression levels.
#' @param pseudocount a numeric value to be added to the expression matrix prior to transformation.
#' @param integerise a boolean specifying whether the expression data should be rounded to the nearest integer.
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
#' Drost HG et al. (2015) Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
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
#' John B. Burbidge, Lonnie Magee and A. Leslie Robb (1988) Alternative Transformations to Handle Extreme Values of the Dependent Variable. Journal of the American Statistical Association, 83(401): 123-127.
#'
#' G. E. P. Box and D. R. Cox (1964) An Analysis of Transformations. Journal of the Royal Statistical Society. Series B (Methodological), 26(2): 211-252.
#' 
#' @author Hajk-Georg Drost
#' @seealso  \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}},
#' \code{\link{tfStability}}
#' @examples
#' \dontrun{
#' data(PhyloExpressionSetExample)
#' 
#' # a simple example is to transform the gene expression levels
#' # of a given PhyloExpressionSet using a sqrt or log2 transformation
#' 
#' PES.sqrt <- tf(PhyloExpressionSetExample, sqrt)
#' 
#' PES.log2 <- tf(PhyloExpressionSetExample, log2)
#' 
#' # plotting the TAI using log2 transformed expression levels
#' # and performing the Flat Line Test to obtain the p-value
#' PlotSignature(ExpressionSet = tf(PhyloExpressionSetExample, log2),
#'             permutations = 1000)
#' 
#' 
#' 
#' # in case the expression matrix contains 0s, a pseudocount can be added prior
#' # to certain transformations, e.g. log2(x+1) where 1 is the pseudocount.
#' 
#' PhyloExpressionSetExample[4,3] = 0
#' PES.log2 <- tf(PhyloExpressionSetExample, log2, pseudocount = 0)
#' 
#' # this should return -Inf at PES.log2[4,3] the issue here is that 
#' # -Inf cannot be used to compute the phylotranscriptomic profile.
#' 
#' PES.log2 <- tf(PhyloExpressionSetExample, log2, pseudocount = 1)
#' # log2 transformed expression levels can now be used in downstream analyses.
#' 
#' 
#' # to perform rank transformation
#' 
#' PES.rank <- tf(PhyloExpressionSetExample, FUN = function(x) apply(x, 2, base::rank))
#' 
#' 
#' # rlog and vst transformations are now also possible by loading the DESeq2 package
#' # and transforming the data with the parameter integerise = TRUE.
#' library(DESeq2) # make sure the DESeq2 version >= 1.29.15 for rlog
#' PES.vst <- tf(PhyloExpressionSetExample, vst, integerise = TRUE)
#' }
#' 
#' 
#' @export

tf <- function(ExpressionSet, FUN, pseudocount = 0, integerise = FALSE){
  if (!is.numeric(pseudocount) | !length(pseudocount) == 1) {
    stop("pseudocount must be a single numeric value")
  }
        ExpressionSet <- as.data.frame(ExpressionSet)
        is.ExpressionSet(ExpressionSet)
        
        ExpressionMatrix <- as.matrix(ExpressionSet[ , -c(1,2)] + pseudocount)
        # rownames(ExpressionMatrix) <- ExpressionSet[,2]
        
        if(integerise){
          ExpressionMatrix <- round(ExpressionMatrix, digits = 0)
        }
        
        f <- match.fun(FUN)
        
        res_mat <- f(ExpressionMatrix)
        
        res <- base::cbind(ExpressionSet[ , c(1,2)], base::as.data.frame(res_mat))
        return(res)
}
