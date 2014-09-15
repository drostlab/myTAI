\name{MeanExpression}
\alias{MeanExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the absolute mean expression profiles for each phylostratum or divergence-stratum.}
\description{
This function computes the mean expression levels of all given phylostrata or divergence-strata of a given standard PhyloExpressionSet or DivergenceExpressionSet.
}
\usage{
MeanExpression(ExpressionSet, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
  \item{plot}{
a boolean value specifying whether a plot of the mean expression profiles shall be drawn.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

Either a plot showing the mean expression of a given PhyloExpressionSet or DivergenceExpressionSet, or a list cotaining:

\item{mean_expression}{a numeric vector storing the corresponding mean expression levels.}
\item{sd}{the standard deviation.}

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
\code{\link{PlotMeans}}
}
\examples{
# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# plot mean absolute expression levels for each phylostratum in a standard PhyloExpressionSet
MeanExpression(PhyloExpressionSetExample, plot = TRUE)

# plot mean absolute expression levels for each divergence-stratum in a standard DivergenceExpressionSet
MeanExpression(DivergenceExpressionSetExample, plot = TRUE)

MeanExpression(PhyloExpressionSetExample, plot = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
