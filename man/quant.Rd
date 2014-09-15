\name{quant}
\alias{quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the corresponding quantile value within a expression matrix.}
\description{
Function to compute the corresponding quantile value for a given quantile-threshold for
a given numeric matrix.}
\usage{
quant(ExpressionMatrix, quantile = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionMatrix}{
a numeric matrix representing a gene expression matrix for which the corresponding quantile value for each developmental stage or experiment shall be computed.
}
  \item{quantile}{
a numeric value specifying the quantile threshold that shall be applied to the gene expression matrix.}
}
\details{
This function computes for each developmental stage s or experiment s of a given expression matrix the corresponding quantile value
for a given quantile-threshold and returns a quantile-value vector.}
\value{
a numeric vector storing the corresponding quantile values of a given gene expression matrix.
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
\code{\link{quantile}}
}
\examples{
# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)

# get the 0.5 quantiles for each developmental stage within this PhyloExpressionSet
quant(PhyloExpressionSetExample[ , 3:9], quantile = 0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

