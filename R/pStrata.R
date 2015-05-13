#' @title Compute Partial Strata Values
#' @description This function computes the partial \code{\link{TAI}} or \code{\link{TDI}}
#' values for all Phylostrata or Divergence Strata.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#' 
#' # compute partial TAI values for each Phylostratum
#' partialStrata <- pStrata(PhyloExpressionSetExample)
#' 
#' # show that colSums of pStrata is equavalent to the TAI values
#' all.equal(colSums(partialStrata),TAI(PhyloExpressionSetExample))
#' 
#' # show that colSums of pStrata is equavalent to colSums of pMatrix(PhyloExpressionSetExample)
#' all.equal(colSums(partialStrata),colSums(pMatrix(PhyloExpressionSetExample)))
#' 
#' 
#' @export

pStrata <- function(ExpressionSet){
        
        is.ExpressionSet(ExpressionSet)
        AGE <- NULL
        
        df <- as.data.frame(cbind(ExpressionSet[ , 1], pMatrix(ExpressionSet)))
        colnames(df)[1] <- "AGE"
        pStrataDF <- dplyr::summarise_each(dplyr::group_by(df,AGE),dplyr::funs(sum))
        colnames(pStrataDF)[1] <- names(ExpressionSet)[1]
        res <- as.matrix(pStrataDF[ , 2:(ncol(ExpressionSet)-1)])
        rownames(res) <- names(table(ExpressionSet[ , 1]))
        
        return(res)
}



