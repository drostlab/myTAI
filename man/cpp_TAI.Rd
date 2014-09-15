\name{cpp_TAI}
\alias{cpp_TAI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the Transcriptome Age Index or Transcriptome Divergence Index.
}
\description{
A C++ function called by the R function \code{\link{TAI}}.
}
\usage{
cpp_TAI(ExpressionSet,Phylostratum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
a numeric matrix containing gene expression levels.
}
\item{Phylostratum}{
a numeric vector containing PS or PS values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric vector containing the TAI or TDI values for all given developmental stages.
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
\code{\link{TAI}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
