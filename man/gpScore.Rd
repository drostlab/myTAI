\name{gpScore}
\alias{gpScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the TAI or TDI destruction score of the global hourglass pattern based
on a pre-computed TAI or TDI vector.
}
\description{
This function reduces the destruction of an hourglass shaped pattern to a single score value.

Based on a given \code{\link{TAI}} or \code{\link{TDI}} pattern the given vector is being
divided into three developmental modules: early, mid, and late. The corrisponding \code{\link{TAI}} or \code{\link{TDI}}
values in each developmental module are accumulated using the \emph{scoringMethod} argument ("max-min" or "mean-mean").

In more detail:

(1) for a given \code{\link{TAI}} or \code{\link{TDI}} vector \emph{tai_profile} or \emph{tdi_profile}, we classify each value of \emph{tai_profile} or \emph{tdi_profile} into its corresponding developmental module early, mid, or late.

(2) accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using the arithmetic mean (\code{\link{mean}}) in case scoringMethod = "mean-mean", or accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using \code{\link{max}} for the early and late module and \code{\link{min}} for the mid module in case scoringMethod = "max-min".

(3) then reduce the three values for each developmental module by computing the difference between: early - mid, and late - mid.

(4) the two difference values are referred to as a_early and a_late. 


Each developmental module now has an accumulated representation value which is being reduced to one value using the
\emph{method} argument ("max", "min", or "mean"). 

The "max", "min", or "mean" parameters refer to the following reduction procedure:

Given the two accumulated values for each hourglass module: a_early and a_late,
we reduce the two given values by:

\emph{"max"}

\eqn{S = max{a_early,a_late}}

\emph{"min"}

\eqn{S = min{a_early,a_late}}

\emph{"mean"}

\eqn{S = mean{a_early,a_late}}

All together this results in a global score \emph{S}.
This global score \emph{S} is being returned by this function \code{\link{gpScore}}.

}
\usage{
gpScore(age_vals, early, mid, late, method, scoringMethod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age_vals}{
a numeric vector containing TAI or TDI values for each developmental stage s.
}
\item{early}{
a numeric vector including the numeric stage values that correspond to the early phase of development.
 }
  \item{mid}{
a numeric vector including the numeric stage values that correspond to the middle phase of development.
}
  \item{late}{
a numeric vector including the numeric stage values that correspond to the late phase of development.
}
  \item{method}{
method to determine the two value reduction value, resulting in the global score S: "max", or "min", or "mean".
}
  \item{scoringMethod}{
method to determine the module accumulation value: "max-min" or "mean-mean".
}
}
\details{
The gpScore is a heuristic score enabling to  construct a test statistic to determine 
the significance of a present (phylotranscriptomic) hourglass pattern.
}
\value{
a numeric value representing the hourglass destruction score.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ReductiveHourglassTest}}, \code{\link{TAI}}, \code{\link{TDI}}
}
\examples{
# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet:

# compute the TAI profile
TAIs <- TAI(PhyloExpressionSetExample)

# compute the global hourglass destruction score for the TAIs profile using reduction method: mean(mean-mean)
gScore <- gpScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,method = "mean",scoringMethod = "mean-mean")


# example DivergenceExpressionSet:

# compute the TDI profile
TDIs <- TDI(DivergenceExpressionSetExample)

# compute the global hourglass destruction score for the TDIs profile using reduction method: mean(mean-mean)
gScore <- gpScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7,method = "mean",scoringMethod = "mean-mean")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
