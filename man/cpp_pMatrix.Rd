\name{cpp_pMatrix}
\alias{cpp_pMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the partial TAIs or TDIs for each gene.
}
\description{
The partial TAI or TDI matrix can be used to assign each gene its
corresponding TAI or TDI contribution profile. The overall sum of
partial TAI or TDI values over all genes results in the \emph{global} TAI or TDI profile that
can also be obtained from the \code{\link{TAI}} or \code{\link{TDI}} function.
}
\usage{
cpp_pMatrix(ExpressionSet, AgeVector)
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
a numeric matrix storing the partial TAI or TDI profiles for each gene within a given PhyloExpressionSet or DivergenceExpressionSet.
}
\references{
Domazet-Loso T and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.

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
\code{\link{pMatrix}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
