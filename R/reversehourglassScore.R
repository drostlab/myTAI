#' @title Compute the Reverse Hourglass Score for the Reverse Hourglass Test
#' @description
#' 
#'         This function reduces the destruction of an hourglass shaped pattern to a single score value.
#'         
#'         Based on a given \code{\link{TAI}} or \code{\link{TDI}} pattern the given vector is being
#'         divided into three developmental modules: early, mid, and late. The corrisponding \code{\link{TAI}} or \code{\link{TDI}}
#'         values in each developmental module are accumulated using the \emph{scoringMethod} argument ("max-min" or "mean-mean").
#'         
#'         In more detail:
#'                 
#'         (1) for a given \code{\link{TAI}} or \code{\link{TDI}} vector \emph{tai_profile} or \emph{tdi_profile}, we classify each value of \emph{tai_profile} or \emph{tdi_profile} into its corresponding developmental module early, mid, or late.
#'         
#'         (2) accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using the arithmetic mean (\code{\link{mean}}) in case scoringMethod = "mean-mean", or accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using \code{\link{max}} for the early and late module and \code{\link{min}} for the mid module in case scoringMethod = "max-min".
#'         
#'         (3) then reduce the three values for each developmental module by computing the difference between: mid - early, and mid - late.
#'         
#'         (4) the two difference values are referred to as a_early and a_late. 
#v         
#'         
#'         Each developmental module now has an accumulated representation value which is being reduced to one value using the
#'         \emph{method} argument ("max", "min", or "mean"). 
#'         
#v         The "max", "min", or "mean" parameters refer to the following reduction procedure:
#'                 
#'                 Given the two accumulated values for each hourglass module: a_early and a_late,
#'         we reduce the two given values by:
#'                 
#' \itemize{
#' \item \emph{"max"}: \eqn{S = max{a_early,a_late}}
#' \item \emph{"min"}: \eqn{S = min{a_early,a_late}}
#' \item \emph{"mean"}: \eqn{S = mean{a_early,a_late}}
#' }
#'All together this results in a global score \emph{S}.
#'This global score \emph{S} is being returned by this function.
#' @param age_vals a numeric vector containing \code{\link{TAI}} or \code{\link{TDI}} values for each developmental stage s.
#' @param early a numeric vector including the numeric stage values that correspond to the early phase of development.
#' @param mid a numeric vector including the numeric stage values that correspond to the middle phase of development.
#' @param late a numeric vector including the numeric stage values that correspond to the late phase of development.
#' @param method to determine the two value reduction value, resulting in the global score S: "max", or "min", or "mean".
#' @param scoringMethod method to determine the module accumulation value: "max-min" or "mean-mean".
#' @param profile.warn a boolean value indicating whether a warning is printed when the low-high-low pattern isn't followed.
#' @details
#' 
#' The gpScore is a heuristic score enabling to  construct a test statistic to determine  
#' the significance of a present (phylotranscriptomic) hourglass pattern.
#' 
#' @return a numeric value representing the hourglass destruction score.
#' @author Hajk-Georg Drost
#' @references
#' Drost et al. (2015), Evidence for active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis. Mol Bio Evol.
#' 
#' @seealso \code{\link{ReverseHourglassTest}}, \code{\link{TAI}}, \code{\link{TDI}}
#' @examples
#' 
#'  # read standard phylotranscriptomics data
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'
#'  # example PhyloExpressionSet:
#'
#'  # compute the TAI profile
#'  TAIs <- TAI(PhyloExpressionSetExample)
#'
#'  # compute the global reverse hourglass destruction score 
#'  # for the TAIs profile using reduction method: mean(mean-mean)
#'  reversehourglass_score <- reversehourglassScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
#'                      method = "mean",scoringMethod = "mean-mean")
#'
#'
#'  # example DivergenceExpressionSet:
#'
#'  # compute the TDI profile
#'  TDIs <- TDI(DivergenceExpressionSetExample)
#'
#'  # compute the global reverse hourglass destruction score for the TDIs profile 
#'  # using reduction method: mean(mean-mean)
#'  reversehourglass_score <- reversehourglassScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7,
#'                      method = "mean",scoringMethod = "mean-mean")
#'                      
#'  # get warning if the expected pattern isn't followed
#'  reversehourglass_score <- reversehourglassScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
#'                      method = "mean",scoringMethod = "mean-mean",profile.warn=TRUE)
#'  
#' @export
reversehourglassScore <- function(age_vals,early,mid,late,method,scoringMethod,profile.warn=FALSE)
{
        
        Score.Early <- vector(mode = "numeric", length = 1)
        Score.Late <- vector(mode = "numeric", length = 1)
        age_valsEarly <- vector(mode = "numeric",length = length(early))
        age_valsMid <- vector(mode = "numeric",length = length(mid))
        age_valsLate <- vector(mode = "numeric",length = length(late))
        
        
        age_valsEarly <- age_vals[early]
        age_valsMid <- age_vals[mid]
        age_valsLate <- age_vals[late]
        
        
        
        if(scoringMethod == "max-min"){
                Score.Early <- min(age_valsMid) - max(age_valsEarly)
                Score.Late <- min(age_valsMid) - max(age_valsLate)
        }
        
        if(scoringMethod == "mean-mean"){
                Score.Early <- mean(age_valsMid) - mean(age_valsEarly)
                Score.Late <- mean(age_valsMid) - mean(age_valsLate)    
        }
        
        if(profile.warn){
          if(sign(Score.Early) == -1 | sign(Score.Late) == -1){
            message("The phylotranscriptomic pattern may not follow a reverse hourglass pattern (low-high-low).")
          } 
        }
        
        if(method == "max")
                return(max(Score.Early, Score.Late))
        if(method == "min")
                return(min(Score.Early, Score.Late))
        if(method == "mean")
                return(mean(Score.Early, Score.Late))
        
        
}


