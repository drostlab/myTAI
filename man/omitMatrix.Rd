\name{omitMatrix}
\alias{omitMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the TAI or TDI profiles without gene i for each single gene in a given PhyloExpressionSet or DivergenceExpressionSet object.}
\description{
For each gene i, exclude the corresponding gene i from the global PhyloExpressionSet or DivergenceExpressionSet and compute the TAI or TDI profile for the corresponding global PhyloExpressionSet or DivergenceExpressionSet with excluded gene i. This results in a TAI or TDI profile Matrix storing the TAI or TDI profile for each omitted gene i.

}
\usage{
omitMatrix(ExpressionSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
 
}
\details{
}
\value{
a numeric matrix storing TAI or TDI profile for each omitted gene i.
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

}
\examples{

# read standard phylotranscriptomics data
#data(PhyloExpressionSetExample)
#data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
#omMatrix <- omitMatrix(PhyloExpressionSetExample)

# example DivergenceExpressionSet
#omMatrix <- omitMatrix(DivergenceExpressionSetExample)



}
