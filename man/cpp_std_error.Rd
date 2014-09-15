\name{cpp_std_error}
\alias{cpp_std_error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the standard error of a numeric vector.
}
\description{
This function computes the standard error of a given numeric vector. Internally this function is written in C++.}
\usage{
cpp_std_error(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector for which a standard error should be computed.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric value representing the standard error of a given numeric vector.
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

# a simple numeric vector
x <- 1:20

# compute the standard error of x
std.error(x)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
