\name{bootMatrix}
\alias{bootMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute TAI or TDI profiles based on N permuted PhyloExpressionSets or DivergenceExpressionSets.
}
\description{
This function computes the TAI for a row permutated PhyloExpressionSet or DivergenceExpressionSet.

One can determine the number of permutations. The function then returns a \code{\link{TAI}} or \code{\link{TDI}} matrix holding
the \code{\link{TAI}} or \code{\link{TDI}} profiles of the permutated PhyloExpressionSets or DivergenceExpressionSets. This procedure
can be used for Test-Statistics based on the TAI or TDI profiles.
}
\usage{
bootMatrix(ExpressionSet, permutations = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
  \item{permutations}{
a numeric value specifying the number of permutations to be performed.
}
}
\details{
The sampled \code{\link{TAI}} or \code{\link{TDI}} matrix samples the phylostratum or divergence-stratum vector of
a given PhyloExpressionSet or DivergenceExpressionSet and computes the corresponding TAI or TDI profiles
of the randomly assigned phylostrata or divergence-strata. This sampling is then performed N times, yielding N randomly sampled TAI or TDI profiles.
This random TAI or TDI profile matrix can furthermore be used to perform statistical tests (such as the \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}}) based on the significance of TAI or TDI patterns.

}
\value{
a numeric matrix representing N randomly permuted TAI or TDI profiles.
}
\references{
Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
}
\examples{

# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet using 1000 permutations
randomTAI.Matrix <- bootMatrix(PhyloExpressionSetExample, permutations = 1000)

# example DivergenceExpressionSet using 1000 permutations
randomTDI.Matrix <- bootMatrix(DivergenceExpressionSetExample, permutations = 1000)


}
