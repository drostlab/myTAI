#' @title Compute the Hourglass Score for the EarlyConservationTest
#' @description This function computes the EarlyConservationTest score for a given \code{\link{TAI}}
#' or \code{\link{TDI}} pattern.
#' 
#' The reductive early conservation test is a permutation test based on the following test statistic. 
#'
#' - A set of developmental stages is partitioned into three modules - early, mid, and late - based on prior biological knowledge.
#'
#' - The mean \code{\link{TAI}} or \code{\link{TDI}} value for each of the three modules T_early, T_mid, and T_late are computed. 
#'
#' - The two differences D1 = T_mid - T_early and D2 = T_late - T_early are calculated.
#'
#' - The minimum D_min of D1 and D2 is computed as final test statistic of the reductive hourglass test.
#'
#' This function \emph{ecScore} computes the \emph{D_min} value for a given \code{\link{TAI}} or \code{\link{TDI}}
#' stored in the \code{age_vals} argument.
#'
#' @param age_vals a numeric vector containing \code{\link{TAI}} or \code{\link{TDI}} values for each developmental stage s.
#' @param early a numeric vector storing the numeric stage values that correspond to the early phase of development.
#' @param mid a numeric vector storing the numeric stage values that correspond to the middle phase of development.
#' @param late a numeric vector storing the numeric stage values that correspond to the late phase of development.
#' @return a numeric value representing the early conservation score.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{EarlyConservationTest}}, \code{\link{TAI}}, \code{\link{TDI}}
#' @examples
#' 
#'  # read standard phylotranscriptomics data
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'
#'  # Example PhyloExpressionSet:
#'
#'  # compute the TAI profile
#'  TAIs <- TAI(PhyloExpressionSetExample)
#'
#'  # compute the early conservation score for the TAI profile
#'  ec_score <- ecScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7)
#'
#'
#'  # Example DivergenceExpressionSet:
#'
#'  # compute the TDI profile
#'  TDIs <- TDI(DivergenceExpressionSetExample)
#'
#'  # compute the early conservation score for the TDI profile
#'  ec_score <- ecScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7)
#'  
#'  # compute ecScore() vector from bootMatrix()
#'  apply(bootMatrix(PhyloExpressionSetExample,10),1,ecScore,early = 1:2,mid = 3:5,late = 6:7)
#'  
#'  
#' @export

ecScore <- function(age_vals,early,mid,late){
        
        D1 <- vector(mode = "numeric", length = 1)
        D2 <- vector(mode = "numeric", length = 1)
        D_min <- vector(mode = "numeric", length = 1)
        
        D1 <- mean(age_vals[mid]) - mean(age_vals[early])
        D2 <- mean(age_vals[late]) - mean(age_vals[early])
        
        D_min <- min(D1,D2)
        
        return(D_min)
}