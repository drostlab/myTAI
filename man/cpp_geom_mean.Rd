\name{cpp_geom_mean}
\alias{cpp_geom_mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the geometric mean of a vector.}
\description{
This function computes the geometric mean of a given numeric vector x.}
\usage{
cpp_geom_mean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numerical vector for which the geometric mean shall be computed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric value representing the corresponding geometric mean of the input vector.
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
\code{\link{harmonic.mean}}
}
\examples{
# simple vector
x <- 1:20

# geometric mean of the simple vector x
geom.mean(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
