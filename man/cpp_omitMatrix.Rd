\name{cpp_omitMatrix}
\alias{cpp_omitMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the TAIs or TDIs for each omitted gene.
}
\description{
This function computes the TAI profile for each omitted gene within a given PhyloExpressionSet or DivergenceExpressionSet. 
}
\usage{
cpp_omitMatrix(ExpressionSet, AgeVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a numeric matrix storing the gene expression levels of the corresponding PhyloExpressionSet or DivergenceExpressionSet.
}
  \item{AgeVector}{
a numeric vector storing the phylostratum vector or divergence-stratum vector of the corresponding PhyloExpressionSet or DivergenceExpressionSet.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric matrix storing the omitted TAI or TDI profiles for each gene within a given PhyloExpressionSet or DivergenceExpressionSet.
}
\references{
Eddelbuettel D and Francois R. 2011 "Rcpp: Seamless R and C++ Integration". Journal of
Statistical Software, 40(8), 1-18. URL http://www.jstatsoft.org/v40/i08/.

Eddelbuettel D. 2013 "Seamless R and C++ Integration with Rcpp". Springer, New York. ISBN 978-1-4614-6867-7.
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

