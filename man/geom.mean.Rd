\name{geom.mean}
\alias{geom.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the geometric mean of a vector.}
\description{
This function computes the geometric mean of a given numeric vector x.}
\usage{
geom.mean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numerical vector for which the geometric mean shall be computed.
}
}
\details{
The geometric mean is being computed as follows:

\deqn{( \prod x_i)^1/n  = exp( 1/n * \sum log(x_i))}

}
\value{
a numeric value representing the corresponding geometric mean of the input vector.
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

