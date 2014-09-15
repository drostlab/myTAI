\name{REMatrix}
\alias{REMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the relative expression profiles of all phylostrata or divergence-strata
available in a PhyloExpressionSet or DivergenceExpressionSet.}
\description{
This function computes the relative expression profiles of all given phylostrata or divergence-strata within a given PhyloExpressionSet or DivergenceExpressionSet.
}
\usage{
REMatrix(ExpressionSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
}
\details{
For each phylostratum or divergence-stratum the corresponding relative expression profile is being computed as follows:

\deqn{f_js = ( e_js - e_j min ) / ( e_j max - e_j min )}

where \eqn{e_j min} and \eqn{e_j max} is the minimum/maximum mean expression level of phylostratum j over the  developmental stages s. This linear transformation corresponds to a shift by \eqn{e_j min} and a subsequent shrinkage by \eqn{e_j max - e_j min}. As a result, the relative expression level \eqn{f_js} of developmental stage s with minimum \eqn{e_js} is 0, the relative expression level \eqn{f_js} of the developmental stage s with maximum \eqn{e_js} is 1, and the relative expression levels \eqn{f_js} of all other stages s range between 0 and 1, accordingly.

}
\value{

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
\code{\link{RE}}
}
\examples{
# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
REMatrix(PhyloExpressionSetExample)

# example DivergenceExpressionSet
REMatrix(DivergenceExpressionSetExample)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
