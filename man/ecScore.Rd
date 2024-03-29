% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecScore.R
\name{ecScore}
\alias{ecScore}
\title{Compute the Hourglass Score for the EarlyConservationTest}
\usage{
ecScore(age_vals, early, mid, late, profile.warn = FALSE)
}
\arguments{
\item{age_vals}{a numeric vector containing \code{\link{TAI}} or \code{\link{TDI}} values for each developmental stage s.}

\item{early}{a numeric vector storing the numeric stage values that correspond to the early phase of development.}

\item{mid}{a numeric vector storing the numeric stage values that correspond to the middle phase of development.}

\item{late}{a numeric vector storing the numeric stage values that correspond to the late phase of development.}

\item{profile.warn}{a boolean value indicating whether a warning is printed when a low-mid-high pattern isn't followed.}
}
\value{
a numeric value representing the early conservation score.
}
\description{
This function computes the EarlyConservationTest score for a given \code{\link{TAI}}
or \code{\link{TDI}} pattern.

The reductive early conservation test is a permutation test based on the following test statistic. 

- A set of developmental stages is partitioned into three modules - early, mid, and late - based on prior biological knowledge.

- The mean \code{\link{TAI}} or \code{\link{TDI}} value for each of the three modules T_early, T_mid, and T_late are computed. 

- The two differences D1 = T_mid - T_early and D2 = T_late - T_early are calculated.

- The minimum D_min of D1 and D2 is computed as final test statistic of the reductive early conservation test.

This function \emph{ecScore} computes the \emph{D_min} value for a given \code{\link{TAI}} or \code{\link{TDI}}
stored in the \code{age_vals} argument.
}
\examples{

 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # Example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the early conservation score for the TAI profile
 ec_score <- ecScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7)


 # Example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the early conservation score for the TDI profile
 ec_score <- ecScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7)
 
 # compute ecScore() vector from bootMatrix()
 apply(bootMatrix(PhyloExpressionSetExample,10),1,ecScore,early = 1:2,mid = 3:5,late = 6:7)
 
 # get warning if the expected pattern isn't followed
 ec_score <- ecScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,profile.warn=TRUE)
 
}
\seealso{
\code{\link{EarlyConservationTest}}, \code{\link{TAI}}, \code{\link{TDI}}
}
\author{
Hajk-Georg Drost
}
